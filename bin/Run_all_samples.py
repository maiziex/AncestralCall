#!/usr/bin/env python
import pdb
#pdb.set_trace()
import os
import sys
import numpy as np
from argparse import ArgumentParser

parser = ArgumentParser(description="Run depth all:")
parser.add_argument('--in_dir','-i_dir', help="Directory for inputs")
parser.add_argument('--out_dir','-o_dir', help="Directory to store outputs")
parser.add_argument('--ref_dir','-r_dir', help="Directory to reference")
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


#sample_list = ['HG00250','HG00353','HG00512','HG00513','HG00731','HG00732','HG00733','HG00851','HG01971','HG03115','HG03838','hgp','NA12878','NA18552','NA19068','NA19238','NA19239','NA19440','NA19789','NA20587','NA21125','NA24143','NA24149']
#sample_list = ['HG00250','HG00353','HG00512']


def Run_all(sample_list,out_dir,in_dir,ref_file,SV_len,species_list,ref_list,num_of_threads,ref_dir):
    all_cmd = ""
    count = 1
    total_num = len(sample_list)
    for sample_name in sample_list:
        hp1_file = in_dir + sample_name + "_contig.1.fasta"
        hp2_file = in_dir + sample_name + "_contig.2.fasta"
        use_cmd = "Assembly_based_variants_call.py " + " --contig_hp1_fasta " + hp1_file + " --contig_hp2_fasta " + hp2_file + " --out_dir " + out_dir + " --ref_file " + ref_file + " --sample_name " + sample_name  + " --SV_len " + str(SV_len)  + " --species_ref_list " + ref_list + " --species_name_list " + species_list + " --ref_dir " + ref_dir  + " & "
        all_cmd += use_cmd
        count += 1
        if (count - 1)%num_of_threads == 0 or (count - 1) == total_num:
            all_cmd += " wait "
            os.system(all_cmd)
            print(all_cmd)
            all_cmd = ""

    print("all done~")


def Merge_vcf_for_all_samples(out_dir):
    merge_cmd = "bcftools merge " 
    for sample_name in sample_list:
        origin_file = out_dir + sample_name + "_final_del.vcf"  
        sorted_file = out_dir + sample_name + "_final_del_sorted.vcf"  
        zip_file = out_dir + sample_name + "_final_del_sorted.vcf.gz"  
        sort_cmd = "bcftools sort " + origin_file + " > " + sorted_file
        zip_cmd = "bgzip -c " + sorted_file + " > " + zip_file
        idx_cmd = "tabix -p vcf " + zip_file
        os.system(sort_cmd)
        os.system(zip_cmd)
        os.system(idx_cmd)
        merge_cmd +=  zip_file + " " 
    merge_cmd += " -o " + out_dir + "Merged.vcf"
    os.system(merge_cmd)
    print(merge_cmd)
    print("all done")   


if __name__ == "__main__":
    out_dir = args.out_dir
    ref_dir = args.ref_dir
    if os.path.exists(out_dir):
        print("using existing output folder: " + out_dir)
    else:
        os.makedirs(out_dir)
    in_dir = args.in_dir
    ref_file = args.ref_file
    species_ref_list = args.species_ref_list
    species_name_list = args.species_name_list
    SV_len = args.SV_len
    num_of_threads = int(args.num_threads)
    Run_all(sample_list,out_dir,in_dir,ref_file,SV_len,species_name_list,species_ref_list,num_of_threads,ref_dir)
    Merge_vcf_for_all_samples(out_dir)

