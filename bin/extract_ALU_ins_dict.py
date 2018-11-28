#!/usr/bin/env python
import pdb
#pdb.set_trace()
import gzip
from argparse import ArgumentParser
import os
import sys
from collections import defaultdict
import pickle
import subprocess
parser = ArgumentParser(description="Run depth all:")
parser.add_argument('--in_dir','-i_dir', help="Directory for inputs")
parser.add_argument('--out_dir','-o_dir', help="Directory to store outputs", default='results_/')
parser.add_argument('--sample_name','-s', help="sample name")
parser.add_argument('--alu_sequence','-as', help="alu fasta")
parser.add_argument('--out_vcf','-o_vcf', help="out_vcf")

args = parser.parse_args()
script_path = os.path.dirname(os.path.abspath( __file__ ))
code_path = script_path + "/" 
other_tools_path = (os.path.abspath(os.path.join(code_path, '..'))) + "/"

def extract_SV(vcf_file,vcf_output,output_fasta):
    f = open(vcf_file,"r")
    fw = open(vcf_output,"w")
    fw_fasta = open(output_fasta,"w")
    count = 0
    for line in f:
        data = line.rsplit()
        if data[0][0] == "#":
            fw.writelines(line)
        else:
            ref = data[3]
            alt = data[4]
            GT = data[9].split(":")[0]
            if abs(len(ref)-len(alt)) >= 250 and abs(len(ref) - len(alt)) <= 350:
                count += 1
                fw.writelines(line)
                event_name = data[2]
                if data[7] == "SVTYPE=INS":
                    seq_used = data[4]
                elif data[7] == "SVTYPE=DEL":
                    seq_used = data[3]
                fw_fasta.writelines(">" + event_name + "\n")
                fw_fasta.writelines(seq_used + "\n")

    #print(count)
    fw.close()
    fw_fasta.close()
    f.close()


def get_match_num_revised(cigar):
    cigar_len = len(cigar)
    cigar_list = []
    num_string = ""
    for i in range(cigar_len):
        letter = cigar[i]
        if letter.isalpha():
            cigar_list.append(num_string)
            num_string = ""
            cigar_list.append(letter)
        else:
            num_string += letter

    indices = [i for i, x in enumerate(cigar_list) if x == "M"]
    cumu_num = 0
    cumu_start = 0
    cumu_num_list = []
    cumu_start_list = []
    match_num_list = []
    for idx in indices:
        match_num = int(cigar_list[idx - 1])
        match_num_list.append(match_num)
        for i in range(0,idx,1):
            parameter = cigar_list[i+1]
            if parameter in ["M", "S", "H", "I"]:
                cumu_num += int(cigar_list[i])
            if parameter in ["M", "D"]:
                cumu_start += int(cigar_list[i])

        cumu_num_list.append(cumu_num-match_num)
        cumu_start_list.append(cumu_start-match_num)
        cumu_num = 0
        cumu_start = 0

    return (cumu_num_list,cumu_start_list,match_num_list)


def get_alu_dict(alu_bed_file):
    f_alu = open(alu_bed_file,"r")
    alu_dict = defaultdict(list)
    for line in f_alu:
        data = line.rsplit()
        chr_num = data[0]
        start_ = int(data[1])
        end_ = int(data[2])
        alu_dict[chr_num].append([start_,end_])
    pickle.dump(alu_dict,open("alu_dict.p","wb"))
    return alu_dict


def evaluate_alu_ins(vcf_file,lib_prefix,alu_event_save): 
    alu_dict_flag = defaultdict()
    f = open(vcf_file,"r")
    #fw = open(output_file,"w")
    count_alu = 0
    count_not_alu = 0
    count = 0
    for line in f:
        data = line.rsplit()
        if data[0][0] == "#":
            #fw.writelines(line)
            count += 1
        else:
            #print(data)
            chr_num = data[0]
            start_ = int(data[1])
            ref_len = len(data[3])
            alt_len = len(data[4])
            end_ = start_ + ref_len
            event_name = data[2]
            if ref_len < alt_len :
                if event_name in alu_event_save:
                    alu_flag = 1
                else:
                    alu_flag = 0
                if alu_flag == 1:
                    count_alu += 1
                    alu_dict_flag[event_name] = "alu_ins_1"
                else:
                    count_not_alu += 1
                    
    #pickle.dump(alu_dict_flag, open(lib_prefix + "_alu_dict_flag.p","wb"))
    print(count_alu,count_not_alu)
    return alu_dict_flag


def evaluate_sam_file(sam_file):
    f = open(sam_file,"r")
    alu_event_save = defaultdict(int)
    for line in f:
        data = line.rsplit()
        if data[0][0] != "@":
            #print(data)
            cigar = data[5]
            if cigar != "*":
                cumu_num_list,cumu_start_list,match_num_list = get_match_num_revised(cigar)
                total_match_num = sum(match_num_list)
                if float(total_match_num)/300 >= 0.8:
                    #print(total_match_num)
                    event_name = data[0]
                    alu_event_save[event_name] = 1
    return alu_event_save


def add_alu_tag_into_vcf(alu_dict_flag,ins_vcf,ins_vcf2):
    f = open(ins_vcf,"r")
    fw = open(ins_vcf2,"w")
    for line in f:
        data = line.rsplit()
        if data[0][:2] == "##":
            fw.writelines(line)
        elif data[0][:2] == "#C":
            fw.writelines("##FORMAT=<ID=Alu,Number=1,Type=String,Description=\"Alu flag: 1 or 0\">\n")
            fw.writelines(line)
        else:
            data[8] = "GT:Contig:Alu"
            event_name = data[2]
            if event_name in alu_dict_flag:
                if alu_dict_flag[event_name] == "alu_ins_1":
                    alu_flag = 1
                    
                else:
                    alu_flag = 0
            else:
                alu_flag = 0
            data[9] = data[9] + ":" + str(alu_flag)
            fw.writelines("\t".join(data) + "\n")
    f.close()
    fw.close()
    


if __name__ == "__main__":
    in_dir = args.in_dir
    out_dir = args.out_dir
    out_vcf = args.out_vcf
    sample_name = args.sample_name
    ALU_one_seq_fa = args.alu_sequence
    if os.path.exists(out_dir):
        print("using existing output folder: " + out_dir)
    else:
        os.makedirs(out_dir)
    ins_vcf = in_dir + sample_name + "_INS.vcf"
    ins_vcf2 = in_dir + out_vcf
    alu_vcf = out_dir + sample_name + "_alu_ins.vcf"
    alu_fasta = out_dir + sample_name + "_alu_ins.fasta"
    alu_sam = out_dir + sample_name + "_alu_ins.sam"

    extract_SV(ins_vcf,alu_vcf,alu_fasta)
    try:
        align_alu_cmd = other_tools_path + "minimap2/minimap2 -a " + ALU_one_seq_fa + " " + alu_fasta + " > " +  alu_sam
    except:
        align_alu_cmd = "minimap2 -a " + ALU_one_seq_fa + " " + alu_fasta + " > " +  alu_sam

    p = subprocess.Popen(align_alu_cmd,shell=True).wait()
    alu_event_save = evaluate_sam_file(alu_sam)
    lib_prefix = out_dir + sample_name
    alu_dict_flag = evaluate_alu_ins(alu_vcf,lib_prefix,alu_event_save)
    
    add_alu_tag_into_vcf(alu_dict_flag,ins_vcf,ins_vcf2)
