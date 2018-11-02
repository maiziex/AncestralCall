#!/usr/bin/env python
import pdb
#pdb.set_trace()
from collections import defaultdict
import pickle
import numpy as np
import os
import sys
from argparse import ArgumentParser

parser = ArgumentParser(description="SV (INS) calling:")
parser.add_argument('--out_dir','-o_dir', help="Directory to store outputs", default='results_/')
parser.add_argument('--sample_name','-s', help="sample name")
parser.add_argument('--SV_len','-l',type=int,help="SV length threshold", default=1)

args = parser.parse_args()
use_chr_list = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]

def Extract_SV_info(hap_SV_file,output_file,SV_len):
    f = open(hap_SV_file,"r")
    contig_dict = defaultdict(list)
    curr = 0
    for line in f:
        #print(curr)
        curr += 1
        data = line.rsplit()
        if data[0] == "V":
            #print(data)
            mapq = int(data[5])
            repeat = int(data[4])
            if mapq >= 20 and repeat == 1: 
                chr_num = data[1]
                if chr_num in use_chr_list:
                    ref_start = int(data[2])
                    ref_end = int(data[3])
                    ref = data[6]
                    alt = data[7]
                    contig_num = int(data[8])
                    contig_start = int(data[9])
                    contig_end = int(data[10])
                    strand_dir = data[-1]
                    contig_dict[(chr_num,ref_start)].append([chr_num,ref_start, ref_end, ref, alt, contig_num,contig_start,contig_end,strand_dir])

    contig_dict_filter = defaultdict(list)
    count = 0
    for key, val in contig_dict.items():
        if len(val) > 1:
            #print(">1:")
            #print(val)
            count += 1
        else:
            contig_dict_filter[key] = val[0]


    print(len(contig_dict))
    print(len(contig_dict_filter))

    
    SV_dict = defaultdict(list)
    for key, val in contig_dict_filter.items():
        ref = val[3] 
        alt = val[4]
        len_var = abs(len(ref)-len(alt))
        if len_var >= SV_len or (alt == "-" and len_var == SV_len - 1) or (ref == "-" and len_var == SV_len - 1):
            SV_dict[key] = val
        
    pickle.dump(SV_dict, open(output_file,"wb"))

    print(len(SV_dict))
    return SV_dict


def compare_two_haploid_SV(SV_dict_contig_1_file, SV_dict_contig_2_file,lib_prefix):
    ins_homo_dict = defaultdict(int)
    ins_hetero_dict = defaultdict(int)
    ins_homo_sv = defaultdict(list)
    ins_hetero_sv = defaultdict(list)
    SV_dict_contig_1 = pickle.load(open(SV_dict_contig_1_file,"rb"))
    SV_dict_contig_2 = pickle.load(open(SV_dict_contig_2_file,"rb"))
    count_total_share = 0
    homo_key_dict = defaultdict(int)
    for key, val in SV_dict_contig_1.items():
        flag = 0
        if key in SV_dict_contig_2:
            #print(key)
            homo_key_dict[key] = 1
            val_2 = SV_dict_contig_2[key]
            chr_num = val[0]
            start_1 = val[1]
            end_1 = val[2]
            start_2 = val_2[1]
            end_2 = val_2[2]
            ref_1 = val[3]
            ref_2 = val_2[3]
            alt_1 = val[4]
            alt_2 = val_2[4]

            if start_1 == start_2 and end_1 == end_2 and ref_1 == ref_2 and alt_1 == alt_2:
                if ref_1 == "-": # ins
                    sv_len_1 = len(alt_1)
                    contig_num_1 = val[-4]
                    contig_start_1 = val[-3]
                    contig_end_1 = val[-2]
                    contig_strand_dir_1 = val[-1]
                    contig_num_2 = val_2[-4]
                    contig_start_2 = val_2[-3]
                    contig_end_2 = val_2[-2]
                    contig_strand_dir_2 = val_2[-1]
                    ins_homo_dict[key] = sv_len_1
                    #ins_homo_sv[key] = [chr_num,start_1,end_1,val[5],val[6],val[7],val_2[5],val_2[6],val_2[7]]
                    ins_homo_sv[key] = [chr_num,start_1,end_1,ref_1,alt_1,contig_num_1,contig_start_1,contig_end_1,contig_strand_dir_1,contig_num_2,contig_start_2,contig_end_2,contig_strand_dir_2]
            else:
                if ref_1 == "-": # ins
                    sv_len_1 = len(alt_1)
                    sv_len_2 = len(alt_2)
                    contig_num_1 = val[-4]
                    contig_start_1 = val[-3]
                    contig_end_1 = val[-2]
                    contig_strand_dir_1 = val[-1]
                    contig_num_2 = val_2[-4]
                    contig_start_2 = val_2[-3]
                    contig_end_2 = val_2[-2]
                    contig_strand_dir_2 = val_2[-1]
                    ins_homo_dict[key] = np.min([sv_len_1,sv_len_2])
                    #ins_homo_sv[key] = [chr_num,start_1,end_1,val[5],val[6],val[7],val_2[5],val_2[6],val_2[7]]
                    ins_homo_sv[key] = [chr_num,start_1,end_1,ref_1,alt_1,contig_num_1,contig_start_1,contig_end_1,contig_strand_dir_1,contig_num_2,contig_start_2,contig_end_2,contig_strand_dir_2]

        else:
            chr_num = val[0]
            start_1 = val[1]
            end_1 = val[2]
            ref_1 = val[3]
            alt_1 = val[4]
            sv_len = abs(len(ref_1) - len(alt_1)) + 1
            for _step in range(start_1 - 10, start_1 + 10):
                key_shift = (chr_num,_step)
                if key_shift in SV_dict_contig_2:
                    val_2 = SV_dict_contig_2[key_shift]
                    start_2 = val_2[1]
                    end_2 = val_2[2]
                    ref_2 = val_2[3]
                    alt_2 = val_2[4]
                    if abs(end_1 - end_2) <= 10:
                        flag = 1
                        homo_key_dict[key_shift] = 1
                        if ref_1 == "-": # ins
                            sv_len_1 = len(alt_1)
                            sv_len_2 = len(alt_2)
                            contig_num_1 = val[-4]
                            contig_start_1 = val[-3]
                            contig_end_1 = val[-2]
                            contig_strand_dir_1 = val[-1]
                            contig_num_2 = val_2[-4]
                            contig_start_2 = val_2[-3]
                            contig_end_2 = val_2[-2]
                            contig_strand_dir_2 = val_2[-1]
                            ins_homo_dict[key_shift] = np.min([sv_len_1,sv_len_2])
                            #ins_homo_sv[key_shift] = [chr_num,start_1,end_1,val[5],val[6],val[7],val_2[5],val_2[6],val_2[7]]
                            ins_homo_sv[key_shift] = [chr_num,start_1,end_1,ref_1,alt_1,contig_num_1,contig_start_1,contig_end_1,contig_strand_dir_1,contig_num_2,contig_start_2,contig_end_2,contig_strand_dir_2]

                        break
            if flag == 0:
                if ref_1 == "-": # ins
                    sv_len_1 = len(alt_1)
                    ins_hetero_dict[key] = sv_len_1
                    contig_start_1 = val[-3]
                    contig_end_1 = val[-2]
                    contig_num_1 = val[-4]
                    contig_strand_dir_1 = val[-1]
                    #ins_hetero_sv[key] = [chr_num,start_1,end_1,val[5],val[6],val[7],0,1,2]
                    ins_hetero_sv[key] = [chr_num,start_1,end_1,ref_1,alt_1,contig_num_1,contig_start_1,contig_end_1,contig_strand_dir_1,0,0,0,0]



    # processing the uniq in contig_2
    for key_2, val_2 in SV_dict_contig_2.items():
        if homo_key_dict[key_2] == 1:
            count_total_share += 1
        else:
            start_2 = val_2[1]
            end_2 = val_2[2]
            ref_2 = val_2[3]
            alt_2 = val_2[4]
            if ref_2 == "-": # ins
                sv_len_2 = len(alt_2)
                ins_hetero_dict[key_2] = sv_len_2
                contig_num_2 = val_2[-4]
                contig_start_2 = val_2[-3]
                contig_end_2 = val_2[-2]
                contig_strand_dir_2 = val_2[-1]
                #ins_hetero_sv[key_2] = [val_2[0],start_2,end_2,0,1,2,val_2[5],val_2[6],val_2[7]]
                ins_hetero_sv[key_2] = [val_2[0],start_2,end_2,ref_2,alt_2,0,0,0,0,contig_num_2,contig_start_2,contig_end_2,contig_strand_dir_2]



    print("done")
    pickle.dump(ins_homo_sv, open(lib_prefix + "_ins_homo_sv.p","wb"))
    pickle.dump(ins_hetero_sv, open(lib_prefix + "_ins_hetero_sv.p","wb"))



def Write_vcf_for_SV(ins_homo_file,ins_hetero_file,vcf_output,sample_name):
    fw = open(vcf_output,"w")
    fw.writelines("##fileformat=VCFv4.2\n")
    fw.writelines("##contig=<ID=chr1,length=248956422>\n")
    fw.writelines("##contig=<ID=chr2,length=242193529>\n")
    fw.writelines("##contig=<ID=chr3,length=198295559>\n")
    fw.writelines("##contig=<ID=chr4,length=190214555>\n")
    fw.writelines("##contig=<ID=chr5,length=181538259>\n")
    fw.writelines("##contig=<ID=chr6,length=170805979>\n")
    fw.writelines("##contig=<ID=chr7,length=159345973>\n")
    fw.writelines("##contig=<ID=chr8,length=145138636>\n")
    fw.writelines("##contig=<ID=chr9,length=138394717>\n")
    fw.writelines("##contig=<ID=chr10,length=133797422>\n")
    fw.writelines("##contig=<ID=chr11,length=135086622>\n")
    fw.writelines("##contig=<ID=chr12,length=133275309>\n")
    fw.writelines("##contig=<ID=chr13,length=114364328>\n")
    fw.writelines("##contig=<ID=chr14,length=107043718>\n")
    fw.writelines("##contig=<ID=chr15,length=101991189>\n")
    fw.writelines("##contig=<ID=chr16,length=90338345>\n")
    fw.writelines("##contig=<ID=chr17,length=83257441>\n")
    fw.writelines("##contig=<ID=chr18,length=80373285>\n")
    fw.writelines("##contig=<ID=chr19,length=58617616>\n")
    fw.writelines("##contig=<ID=chr20,length=64444167>\n")
    fw.writelines("##contig=<ID=chr21,length=46709983>\n")
    fw.writelines("##contig=<ID=chr22,length=50818468>\n")
    fw.writelines("##contig=<ID=chrX,length=156040895>\n")
    fw.writelines("##contig=<ID=chrY,length=57227415>\n")
    fw.writelines("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of SV:DEL=Deletion, INS=Insertion, SNP=snps\">\n")
    fw.writelines("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    fw.writelines("##FORMAT=<ID=Contig,Number=1,Type=String,Description=\"Contig information: Contig name, contig start, contig end\">\n")
    fw.writelines("#CHROM" + "\t" + "POS" + "\t" + "ID" + "\t" + "REF" + "\t" + "ALT" + "\t" + "QUAL" + "\t" + "FILTER" + "\t" + "INFO" + "\t" + "FORMAT" + "\t"  + sample_name + "\n")

    ins_homo = pickle.load(open(ins_homo_file,"rb"))
    ins_hetero = pickle.load(open(ins_hetero_file,"rb"))
    count_id = 1

    for key, val in ins_homo.items():
        print(key)
        chr_num =  str(val[0])
        start_ = val[1]
        end_ = val[2]
        ref = val[3]
        alt = val[4]
        GT = "1/1"
        contig_num_1 = val[5]
        contig_start_1 = val[6]
        contig_end_1 = val[7]
        contig_strand_dir_1 = val[8]
        contig_num_2 = val[9]
        contig_start_2 = val[10]
        contig_end_2 = val[11]
        contig_strand_dir_2 = val[12]
        contig_info = str(contig_num_1) + "_" + str(contig_start_1) + "_" + str(contig_end_1) + "_" + contig_strand_dir_1 + "_" +  str(contig_num_2) + "_" + str(contig_start_2) + "_" + str(contig_end_2) + "_" + contig_strand_dir_2
        fw.writelines(chr_num + "\t" + str(start_) + "\t" + "event" + str(count_id) + "\t" + ref + "\t" + alt + "\t" + "." + "\t" + "PASS" + "\t" + "SVTYPE=INS" + "\t" + "GT:Contig" + "\t"  + "1/1" + ":"  + contig_info +  "\n")

        count_id += 1


    for key, val in ins_hetero.items():
        print(key)
        chr_num =  str(val[0])
        start_ = val[1]
        end_ = val[2]
        ref = val[3]
        alt = val[4]
        GT = "0/1"
        contig_num_1 = val[5]
        contig_start_1 = val[6]
        contig_end_1 = val[7]
        contig_strand_dir_1 = val[8]
        contig_num_2 = val[9]
        contig_start_2 = val[10]
        contig_end_2 = val[11]
        contig_strand_dir_2 = val[12]
        contig_info = str(contig_num_1) + "_" + str(contig_start_1) + "_" + str(contig_end_1) + "_" + str(contig_strand_dir_1) + "_" +  str(contig_num_2) + "_" + str(contig_start_2) + "_" + str(contig_end_2) + "_" + str(contig_strand_dir_2)
        fw.writelines(chr_num + "\t" + str(start_) + "\t" + "event" + str(count_id) + "\t" + ref + "\t" + alt + "\t" + "." + "\t" + "PASS" + "\t" + "SVTYPE=INS" + "\t" + "GT:Contig" + "\t"  + "0/1" + ":"  + contig_info +  "\n")

        count_id += 1
 
    fw.close() 






if __name__ == "__main__":
    out_dir = args.out_dir
    if os.path.exists(out_dir):
        print("using existing output folder: " + out_dir)
    else:
        os.makedirs(out_dir)
    sample_name = args.sample_name
    SV_len = args.SV_len
    lib_prefix = out_dir + sample_name
    SV_dict_contig_1 = Extract_SV_info(lib_prefix + "_hp1.var.txt",lib_prefix + "_INS_dict_contig_1.p",SV_len)
    SV_dict_contig_2 = Extract_SV_info(lib_prefix + "_hp2.var.txt",lib_prefix + "_INS_dict_contig_2.p",SV_len)
    compare_two_haploid_SV(lib_prefix + "_INS_dict_contig_1.p", lib_prefix + "_INS_dict_contig_2.p",lib_prefix)
    Write_vcf_for_SV(lib_prefix + "_ins_homo_sv.p", lib_prefix + "_ins_hetero_sv.p", lib_prefix + "_INS.vcf",sample_name)








