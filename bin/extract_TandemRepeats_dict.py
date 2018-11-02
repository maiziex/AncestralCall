#!/usr/bin/env python
import pdb
#pdb.set_trace()
import gzip
from collections import defaultdict
import pickle
import os
import sys
from argparse import ArgumentParser


parser = ArgumentParser(description="Run depth all:")
parser.add_argument('--in_dir','-i_dir', help="Directory for inputs")
parser.add_argument('--out_dir','-o_dir', help="Directory to store outputs", default='results_/')
parser.add_argument('--sample_name','-s', help="sample name")
parser.add_argument('--variant_type','-v', help="variant type")
parser.add_argument('--in_vcf','-i_vcf', help="in_vcf")
parser.add_argument('--out_vcf','-o_vcf', help="out_vcf")

args = parser.parse_args()
script_path = os.path.dirname(os.path.abspath( __file__ ))
code_path = script_path + "/" 
other_tools_path = (os.path.abspath(os.path.join(code_path, '..'))) + "/"

def convert_vcf_to_fasta(vcf_file,output_fasta):
    f = open(vcf_file,"r")
    fw_fasta = open(output_fasta,"w")
    count = 0
    for line in f:
        data = line.rsplit()
        if data[0][0] != "#":
            ref = data[3]
            alt = data[4]
            GT = data[9].split(":")[0]
            count += 1
            event_name = data[2]
            if data[7] == "SVTYPE=INS":
                seq_used = data[4]
            elif data[7] == "SVTYPE=DEL":
                seq_used = data[3]
            fw_fasta.writelines(">" + event_name + "\n")
            fw_fasta.writelines(seq_used + "\n")

    print(count)
 
    fw_fasta.close()
    f.close()


def save_pickle_file(dict1,fname):
    for value in dict1:
        dict1[value] = dict(dict1[value])
    my_dict = dict(dict1)
    with open(fname, "wb") as f:
        pickle.dump(my_dict,f) 


def check_repeats_len(save_loci):
    repeats_len_dict = defaultdict(int)
    for loci in save_loci:
        _start = loci[0]
        _end = loci[1]
        for _step in range(_start, _end+1):
            repeats_len_dict[_step] = 1

    repeats_len = len(repeats_len_dict)
    return repeats_len


def Get_reads_match_tandem_repeats(tandem_repeats_file,lib_prefix,L_event_del_dict):
    qname_match_repeats = defaultdict(int)
    f = open(tandem_repeats_file,"r")
    count = 0
    count_repeats = 0
    for line in f:
        #print(count)
        data = line.rsplit()
        if count >= 8:
            if data != []:
                if data[0] == "Sequence:":
                    if count_repeats > 1:
                        repeats_len = check_repeats_len(save_loci)
                        repeats_percent = float(repeats_len)/len(L_event_del_dict[prev_qname][1])
                        print("repeats percent: " + str(repeats_percent))
                        qname_match_repeats[prev_qname] = repeats_percent
                        count_repeats = 0
                        save_loci = []
                    else:
                        count_repeats = 0
                        save_loci = []
                    prev_qname = data[1]
                    
                elif data[0] == "Parameters:":
                    count_repeats += 1
                else:
                    count_repeats += 1
                    _start = int(data[0])
                    _end = int(data[1])
                    save_loci.append([_start,_end])
        count += 1

    if count_repeats > 1:
        qname_match_repeats[prev_qname] = 1


    print("done~")
    return qname_match_repeats


def extract_SV_del_seq(vcf_file):
    f = open(vcf_file,"r")
    event_dict = defaultdict()
    for line in f:
        data = line.rsplit()
        if data[0][0] != "#":
            #print(data)
            event_name = data[2]
            event_dict[event_name] = [int(data[1]),data[3]]
    return event_dict


def extract_SV_ins_seq(vcf_file):
    f = open(vcf_file,"r")
    event_dict = defaultdict()
    for line in f:
        data = line.rsplit()
        if data[0][0] != "#":
            #print(data)
            event_name = data[2]
            event_dict[event_name] = [int(data[1]),data[4]]
    return event_dict


def add_repeats_tag_into_vcf(repeats_dict_flag,del_vcf,del_vcf2):
    f = open(del_vcf,"r")
    fw = open(del_vcf2,"w")
    for line in f:
        data = line.rsplit()
        if data[0][:2] == "##":
            fw.writelines(line)
        elif data[0][:2] == "#C":
            fw.writelines("##FORMAT=<ID=TandemRepeats,Number=1,Type=String,Description=\"Repeats flag: >0 or 0\">\n")
            fw.writelines(line)
        else:
            data[8] = "GT:Contig:Alu:TandemRepeats"
            event_name = data[2]
            if event_name in repeats_dict_flag:
                repeats_flag = repeats_dict_flag[event_name]
            else:
                repeats_flag = 0
            data[9] = data[9] + ":" + str(repeats_flag)
            fw.writelines("\t".join(data) + "\n")
    f.close()
    fw.close()
    



if __name__ == "__main__":
    in_dir = args.in_dir
    in_vcf = args.in_vcf
    out_dir = args.out_dir
    out_vcf = args.out_vcf
    sample_name = args.sample_name
    variant_type = args.variant_type
    if os.path.exists(out_dir):
        print("using existing output folder: " + out_dir)
    else:
        os.makedirs(out_dir)
    if variant_type == "DEL":
        del_vcf = in_dir + in_vcf
        del_vcf2 = in_dir + out_vcf
        del_fasta = out_dir + sample_name + "_DEL_with_alu.fasta"
        convert_vcf_to_fasta(del_vcf,del_fasta)
        try:
            trf_cmd = other_tools_path + "trf_tools/trf " + del_fasta + " 2 5 7 80 10 50 2000 -d -h" 
        except:
            trf_cmd = "trf " + del_fasta + " 2 5 7 80 10 50 2000 -d -h"  
        os.system(trf_cmd)
        os.system("mv " + sample_name + "_DEL_with_alu.fasta.2.5.7.80.10.50.2000.dat" + " " + out_dir)
        event_del_dict = extract_SV_del_seq(del_vcf)
        lib_prefix = out_dir + sample_name
        repeats_dict_flag = Get_reads_match_tandem_repeats(del_fasta + ".2.5.7.80.10.50.2000.dat",lib_prefix,event_del_dict)
        add_repeats_tag_into_vcf(repeats_dict_flag,del_vcf,del_vcf2)
    elif variant_type == "INS":
        ins_vcf = in_dir + in_vcf
        ins_vcf2 = in_dir + out_vcf
        ins_fasta = out_dir + sample_name + "_INS_with_alu.fasta"
        convert_vcf_to_fasta(ins_vcf,ins_fasta)
        try:
            trf_cmd = other_tools_path + "trf_tools/trf " + ins_fasta + " 2 5 7 80 10 50 2000 -d -h" 
        except:
            trf_cmd = "trf " + ins_fasta + " 2 5 7 80 10 50 2000 -d -h"  
        os.system(trf_cmd)
        os.system("mv " + sample_name + "_INS_with_alu.fasta.2.5.7.80.10.50.2000.dat" + " " + out_dir)
        event_ins_dict = extract_SV_ins_seq(ins_vcf)
        lib_prefix = out_dir + sample_name
        repeats_dict_flag = Get_reads_match_tandem_repeats(ins_fasta + ".2.5.7.80.10.50.2000.dat",lib_prefix,event_ins_dict)
        add_repeats_tag_into_vcf(repeats_dict_flag,ins_vcf,ins_vcf2)



