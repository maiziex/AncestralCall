#!/usr/bin/env python
import pdb
#pdb.set_trace()
import os
import sys
import numpy as np
from argparse import ArgumentParser
import pickle
from subprocess import Popen
from multiprocessing import Pool,cpu_count,active_children,Manager
import time
parser = ArgumentParser(description="Run depth all:")
parser.add_argument('--in_dir','-i_dir', help="Directory for inputs")
parser.add_argument('--out_dir','-o_dir', help="Directory to store outputs")
parser.add_argument('--ref_file','-r', help="reference fasta file")
parser.add_argument('--sample_name','-s', help="sample name")
parser.add_argument('--contig_hp1_fasta','-c_hp1', help="contig haplotype 1 fasta file")
parser.add_argument('--contig_hp2_fasta','-c_hp2', help="contig haplotype 2 fasta file")
parser.add_argument('--SV_len','-l',type=int,help="SV len threshold", default=20)
parser.add_argument('--num_threads','-nt', help="number of threads")
parser.add_argument('--sample_list','-sl', help='delimited list input', type=str)
parser.add_argument('--species_name_list','-ls', help='delimited list input', type=str)
parser.add_argument('--species_ref_list','-lr', help='delimited list input', type=str)
args = parser.parse_args()
sample_list = [item for item in args.sample_list.split(',')]
script_path = os.path.dirname(os.path.abspath( __file__ ))
code_path = script_path + "/" 
other_tools_path = (os.path.abspath(os.path.join(code_path, '..'))) + "/"

def read_ref(fasta_file,chr_num,out_dir):
    f = open(fasta_file,"r")
    count = 0
    ref_seq = ""
    for line in f:
        if count > 0:
            data = line.rsplit()
            ref_seq += data[0]
        count += 1
    print("total_len for chr" + str(chr_num))
    pickle.dump(ref_seq, open(out_dir + "ref_seq_chr" + str(chr_num) +  ".p","wb"))


def extract_ref_chr(ref_file,chr_num,out_dir):
    fw = open(out_dir + "genome_ref_chr" + str(chr_num) + ".fasta","w")
    f = open(ref_file,"r")
    flag = 0
    total_len = 0
    for line in f:
        data = line.rsplit()
        if data[0] == ">chr" + str(chr_num):
            fw.writelines(">" + str(chr_num) + "\n")
            flag = 1
        elif data[0][0] == ">" and flag == 1:
            break
        else:
            if flag == 1:
                total_len += len(data[0])
                fw.writelines(data[0] + "\n")
    print("chr" + str(chr_num) + ":")
    print(total_len)


def Call_one_sample(hp1_file,hp2_file,out_dir,ref_file,sample_name,SV_len,ref_list,species_list,xin):
    use_cmd = code_path + "Assembly_based_variants_call.py " + " --contig_hp1_fasta " + hp1_file + " --contig_hp2_fasta " + hp2_file + " --out_dir " + out_dir + " --ref_file " + ref_file + " --sample_name " + sample_name  + " --SV_len " + str(SV_len)  + " --species_ref_list " + ref_list + " --species_name_list " + species_list 
    Popen(use_cmd,shell=True).wait()


def Run_all(sample_list,out_dir,in_dir,ref_file,SV_len,species_list,ref_list,num_of_threads):
    all_cmd = ""
    count = 1
    total_num = len(sample_list)
    pool = Pool(processes=num_of_threads)
    for sample_name in sample_list:
        hp1_file = in_dir + sample_name + "_contig.1.fasta"
        hp2_file = in_dir + sample_name + "_contig.2.fasta"
        pool.apply_async(Call_one_sample,(hp1_file,hp2_file,out_dir,ref_file,sample_name,SV_len,ref_list,species_list,"xin"))
        count += 1
        if (count - 1)%num_of_threads == 0 or (count - 1) == total_num:
            pool.close()
            while len(active_children()) > 1:
                time.sleep(0.5)
            pool.join()
            pool = Pool(processes=num_of_threads)
    print("all done~")


def Merge_vcf_for_all_samples_del(out_dir):
    merge_cmd = "bcftools merge " 
    for sample_name in sample_list:
        origin_file = out_dir + sample_name + "_final_del.vcf"  
        sorted_file = out_dir + sample_name + "_final_del_sorted.vcf"  
        zip_file = out_dir + sample_name + "_final_del_sorted.vcf.gz"  
        sort_cmd = "bcftools sort " + origin_file + " > " + sorted_file
        zip_cmd = "bgzip -c " + sorted_file + " > " + zip_file
        idx_cmd = "tabix -p vcf " + zip_file
        Popen(sort_cmd,shell=True).wait()
        Popen(zip_cmd,shell=True).wait()
        Popen(idx_cmd,shell=True).wait()
        merge_cmd +=  zip_file + " " 
    merge_cmd += " -o " + out_dir + "Merged_del.vcf"
    Popen(merge_cmd,shell=True).wait()
    print(merge_cmd)
    print("all done")   


def Merge_vcf_for_all_samples_ins(out_dir):
    merge_cmd = "bcftools merge " 
    for sample_name in sample_list:
        origin_file = out_dir + sample_name + "_final_ins.vcf"  
        sorted_file = out_dir + sample_name + "_final_ins_sorted.vcf"  
        zip_file = out_dir + sample_name + "_final_ins_sorted.vcf.gz"  
        sort_cmd = "bcftools sort " + origin_file + " > " + sorted_file
        zip_cmd = "bgzip -c " + sorted_file + " > " + zip_file
        idx_cmd = "tabix -p vcf " + zip_file
        Popen(sort_cmd,shell=True).wait()
        Popen(zip_cmd,shell=True).wait()
        Popen(idx_cmd,shell=True).wait()
        merge_cmd +=  zip_file + " " 
    merge_cmd += " -o " + out_dir + "Merged_ins.vcf"
    Popen(merge_cmd,shell=True).wait()
    print(merge_cmd)
    print("all done")   


if __name__ == "__main__":
    out_dir = args.out_dir + "/"
    if os.path.exists(out_dir):
        print("using existing output folder: " + out_dir)
    else:
        os.makedirs(out_dir)
    in_dir = args.in_dir + "/"
    ref_file = args.ref_file
    species_ref_list = args.species_ref_list
    species_name_list = args.species_name_list
    SV_len = args.SV_len
    num_of_threads = int(args.num_threads)
    for chr_num in range(1,23):
        extract_ref_chr(ref_file,chr_num,out_dir)
    extract_ref_chr(ref_file,"X",out_dir)
    for chr_num in range(1,23):
        read_ref(out_dir + "genome_ref_chr" + str(chr_num) + ".fasta",chr_num,out_dir)
    read_ref(out_dir + "genome_ref_chrX.fasta",23,out_dir)
    Run_all(sample_list,out_dir,in_dir,ref_file,SV_len,species_name_list,species_ref_list,num_of_threads)
    Merge_vcf_for_all_samples_del(out_dir)
    Merge_vcf_for_all_samples_ins(out_dir)

