#!/usr/bin/env python
import pdb
#pdb.set_trace()
import os
from argparse import ArgumentParser
from collections import defaultdict
import subprocess
from subprocess import Popen, PIPE
from multiprocessing import Pool,cpu_count,active_children,Manager
import time
parser = ArgumentParser(description="Run depth all:")
parser.add_argument('--out_dir','-o_dir', help="Directory to store outputs")
parser.add_argument('--ref_file','-r', help="reference fasta file")
parser.add_argument('--sample_name','-s', help="sample name")
parser.add_argument('--contig_hp1_fasta','-c_hp1', help="contig haplotype 1 fasta file")
parser.add_argument('--contig_hp2_fasta','-c_hp2', help="contig haplotype 2 fasta file")
parser.add_argument('--SV_len','-l',type=int,help="SV len threshold", default=20)
parser.add_argument('--species_name_list','-ls', help='monkey list', type=str)
parser.add_argument('--species_ref_list','-lr', help='reference fasta list', type=str)
args = parser.parse_args()
species_list = [item for item in args.species_name_list.split(',')]
ref_list = [item for item in args.species_ref_list.split(',')]

script_path = os.path.dirname(os.path.abspath( __file__ ))
code_path = script_path + "/" 
other_tools_path = (os.path.abspath(os.path.join(code_path, '..'))) + "/"


def source_python_code():
    Popen("chmod +x " + code_path + "extract_ALU_del_dict.py",shell=True).wait()
    Popen("chmod +x " + code_path + "extract_ALU_ins_dict.py",shell=True).wait()
    Popen("chmod +x " + code_path + "extract_TandemRepeats_dict.py",shell=True).wait()
    Popen("chmod +x " + code_path + "extract_Flanking_seq_for_del.py",shell=True).wait()
    Popen("chmod +x " + code_path + "extract_Flanking_seq_for_ins.py",shell=True).wait()
    Popen("chmod +x " + code_path + "Evaluate_derived_sv_by_flankingseq_align_for_del.py",shell=True).wait()
    Popen("chmod +x " + code_path + "Evaluate_derived_sv_by_flankingseq_align_for_ins.py",shell=True).wait()


def get_paf(ref_file,hap_file,hap_paf,xin):
    use_cmd = other_tools_path + "minimap2/minimap2 -cx asm5 -t8 --cs " + ref_file + " " + hap_file  + " > " + hap_paf  
    Popen(use_cmd,shell=True).wait()


def sort_paf(hap_paf,hap_paf_sorted,xin):
    use_cmd =  "sort -k6,6 -k8,8n " + hap_paf + " > "  + hap_paf_sorted  
    Popen(use_cmd,shell=True).wait()


def get_var(hap_paf_sorted,hap_var_txt,xin):
    use_cmd = other_tools_path + "/k8-0.2.4/k8-Linux " + other_tools_path + "minimap2/misc/paftools.js " + " call -l 1 -L 1 -q 20 " +  hap_paf_sorted + " > " +  hap_var_txt 
    Popen(use_cmd,shell=True).wait()
    

def assembly_based_variants_call_paf(hap1_file,hap2_file,ref_file,sample_name,out_dir):
    hap1_paf = out_dir + sample_name + "_hp1.paf"
    hap2_paf = out_dir + sample_name + "_hp2.paf"
    pool = Pool(processes=2)
    pool.apply_async(get_paf,(ref_file,hap1_file,hap1_paf,"xin"))
    pool.apply_async(get_paf,(ref_file,hap2_file,hap2_paf,"xin"))
    pool.close()
    while len(active_children()) > 1:
        time.sleep(0.5)
    pool.join()
    print("paf done~")


def assembly_based_variants_call_sort(sample_name,out_dir):
    hap1_paf = out_dir + sample_name + "_hp1.paf"
    hap2_paf = out_dir + sample_name + "_hp2.paf"
    hap1_paf_sorted = out_dir + sample_name + "_hp1.paf.sorted"
    hap2_paf_sorted = out_dir + sample_name + "_hp2.paf.sorted" 
    pool = Pool(processes=2)
    pool.apply_async(sort_paf,(hap1_paf,hap1_paf_sorted,"xin"))
    pool.apply_async(sort_paf,(hap2_paf,hap2_paf_sorted,"xin"))
    pool.close()
    while len(active_children()) > 1:
        time.sleep(0.5)
    pool.join()
    print("sort done~")
 

def assembly_based_variants_call(sample_name,out_dir):
    hap1_paf_sorted = out_dir + sample_name + "_hp1.paf.sorted"
    hap2_paf_sorted = out_dir + sample_name + "_hp2.paf.sorted"
    hap1_var_txt = out_dir + sample_name + "_hp1.var.txt"
    hap2_var_txt = out_dir + sample_name + "_hp2.var.txt"
    pool = Pool(processes=2)
    pool.apply_async(get_var,(hap1_paf_sorted,hap1_var_txt,"xin"))
    pool.apply_async(get_var,(hap2_paf_sorted,hap2_var_txt,"xin"))
    pool.close()
    while len(active_children()) > 1:
        time.sleep(0.5)
    pool.join()
    print("var.txt done~")


def Call_SV_del_from_contigs(sample_name,out_dir,SV_len):
    try:
        all_cmd = "python " + code_path + "Extract_SV_info_from_contigs_use_overlap_for_del.py --SV_len " + str(SV_len) + " --sample_name " + sample_name + " --out_dir " + out_dir 
    except:
        all_cmd = "Extract_SV_info_from_contigs_use_overlap_for_del.py --SV_len " + str(SV_len) + " --sample_name " + sample_name + " --out_dir " + out_dir 
    Popen(all_cmd,shell=True).wait()
    print("DEL done~")


def Call_SV_ins_from_contigs(sample_name,out_dir,SV_len):
    try:
        all_cmd = "python " + code_path + "Extract_SV_info_from_contigs_use_shift_for_ins.py --SV_len " + str(SV_len) + " --sample_name " + sample_name + " --out_dir " + out_dir 
    except:
        all_cmd = "Extract_SV_info_from_contigs_use_shift_for_ins.py --SV_len " + str(SV_len) + " --sample_name " + sample_name + " --out_dir " + out_dir 
    Popen(all_cmd,shell=True).wait()
    print("INS done~")


def Call_alu_for_del(in_dir,out_dir,sample_name,ALU_one_seq_fa,out_vcf):
    try:
        called_cmd = "python " + code_path + "extract_ALU_del_dict.py --in_dir " + in_dir + " --out_dir " + out_dir + " --sample_name " + sample_name + " --alu_sequence " + ALU_one_seq_fa + " --out_vcf " + out_vcf
    except:
        called_cmd = "extract_ALU_del_dict.py --in_dir " + in_dir + " --out_dir " + out_dir + " --sample_name " + sample_name + " --alu_sequence " + ALU_one_seq_fa + " --out_vcf " + out_vcf
    Popen(called_cmd,shell=True).wait()

    print("alu done~")


def Call_alu_for_ins(in_dir,out_dir,sample_name,ALU_one_seq_fa,out_vcf):
    try:
        called_cmd = "python " + code_path + "extract_ALU_ins_dict.py --in_dir " + in_dir + " --out_dir " + out_dir + " --sample_name " + sample_name + " --alu_sequence " + ALU_one_seq_fa + " --out_vcf " + out_vcf
    except:
        called_cmd = "extract_ALU_ins_dict.py --in_dir " + in_dir + " --out_dir " + out_dir + " --sample_name " + sample_name + " --alu_sequence " + ALU_one_seq_fa  + " --out_vcf " + out_vcf
    Popen(called_cmd,shell=True).wait()
    print("alu done~")


def Call_repeats(in_dir,out_dir,sample_name,variant_type,in_vcf,out_vcf):
    try:
        called_cmd = "python " + code_path + "extract_TandemRepeats_dict.py --in_dir " + in_dir + " --out_dir " + out_dir + " --sample_name " + sample_name + " --variant_type " + variant_type + " --in_vcf " + in_vcf +  " --out_vcf " + out_vcf
    except:
        called_cmd = "extract_TandemRepeats_dict.py --in_dir " + in_dir + " --out_dir " + out_dir + " --sample_name " + sample_name + " --variant_type " + variant_type + " --in_vcf " + in_vcf + " --out_vcf " + out_vcf
    Popen(called_cmd,shell=True).wait()
    print("tandem repeats done~")


def extract_flanking_seq_del(in_dir,out_dir,sample_name):
    try:
        called_cmd = "python " + code_path + "extract_Flanking_seq_for_del.py  --in_dir " + in_dir + " --out_dir " + out_dir + " --sample_name " + sample_name 
    except:
        called_cmd = "extract_Flanking_seq_for_del.py  --in_dir " + in_dir + " --out_dir " + out_dir + " --sample_name " + sample_name 
    Popen(called_cmd,shell=True).wait()
    print("extract flanking seq done~")


def extract_flanking_seq_ins(in_dir,out_dir,sample_name):
    try:
        called_cmd = "python " + code_path + "extract_Flanking_seq_for_ins.py  --in_dir " + in_dir + " --out_dir " + out_dir + " --sample_name " + sample_name 
    except:
        called_cmd = "extract_Flanking_seq_for_ins.py  --in_dir " + in_dir + " --out_dir " + out_dir + " --sample_name " + sample_name 
    Popen(called_cmd,shell=True).wait()
    print("extract flanking seq done~")


def align_to_monkeys(monkey_ref,flanking_fasta,out_dir,sample_name,monkey_name,variant_type):
    flanking_sam = out_dir + sample_name + "_" + monkey_name + "_" + variant_type + "_for_flankingseq.sam"
    flanking_bam = out_dir + sample_name + "_" + monkey_name + "_" + variant_type + "_for_flankingseq.bam"
    flanking_bam_sorted = out_dir + sample_name + "_" + monkey_name + "_" + variant_type + "_for_flankingseq_sorted.bam"
    use_cmd = other_tools_path + "minimap2/minimap2 -a " + monkey_ref + " " + flanking_fasta + " > " + flanking_sam
    Popen(use_cmd,shell=True).wait()
    use_cmd = "samtools view -Sb " + flanking_sam + " > " + flanking_bam
    Popen(use_cmd,shell=True).wait()
    use_cmd =  "samtools sort " + flanking_bam + " -o " + flanking_bam_sorted
    Popen(use_cmd,shell=True).wait()
    use_cmd = "samtools index " +  flanking_bam_sorted
    Popen(use_cmd,shell=True).wait()


def evaluate_ancestral_del(flanking_bam_sorted,sample_name,monkey_name,out_dir,in_vcf,out_vcf):
    try:
        called_cmd = "python " + code_path +  "Evaluate_derived_sv_by_flankingseq_align_for_del.py  --out_dir " + out_dir + " --sample_name " + sample_name + " --monkey_name " + monkey_name + " --flanking_bam " + flanking_bam_sorted  + " --in_vcf " + in_vcf + " --out_vcf " + out_vcf
  
    except:
        called_cmd = "Evaluate_derived_sv_by_flankingseq_align_for_del.py  --out_dir " + out_dir + " --sample_name " + sample_name + " --monkey_name " + monkey_name + " --flanking_bam " + flanking_bam_sorted + " --in_vcf " + in_vcf + " --out_vcf " + out_vcf
    Popen(called_cmd,shell=True).wait()
    print("ancestral analysis for del done~")


def evaluate_ancestral_ins(flanking_bam_sorted,sample_name,monkey_name,out_dir,in_vcf,out_vcf):
    try:
        called_cmd = "python " + code_path +  "Evaluate_derived_sv_by_flankingseq_align_for_ins.py  --out_dir " + out_dir + " --sample_name " + sample_name + " --monkey_name " + monkey_name + " --flanking_bam " + flanking_bam_sorted  + " --in_vcf " + in_vcf + " --out_vcf " + out_vcf
  
    except:
        called_cmd = "Evaluate_derived_sv_by_flankingseq_align_for_ins.py  --out_dir " + out_dir + " --sample_name " + sample_name + " --monkey_name " + monkey_name + " --flanking_bam " + flanking_bam_sorted + " --in_vcf " + in_vcf + " --out_vcf " + out_vcf
    Popen(called_cmd,shell=True).wait()
    print("ancestral analysis for ins done~")



def Merge_all_variants(sample_name,out_dir):
    del_vcf = out_dir + sample_name + "_DEL.vcf"
    ins_vcf = out_dir + sample_name + "_INS.vcf"
    final_vcf = out_dir + sample_name + "_final.vcf"
    header_vcf = out_dir + "header.vcf"
    fw = open(out_dir + "header.vcf","w")
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

    fw.close()
    use_cmd = "cat " + header_vcf  + " "  + del_vcf + " " + ins_vcf + " > " + final_vcf 
    Popen(use_cmd,shell=True).wait()
    print(use_cmd)
    print("all done~")


if __name__ == "__main__":
    temp_del = defaultdict(str)
    temp_ins = defaultdict(str)
    source_python_code()
    out_dir = args.out_dir
    if os.path.exists(out_dir):
        print("using existing output folder: " + out_dir)
    else:
        os.makedirs(out_dir)
    ref_file = args.ref_file
    SV_len = args.SV_len
    sample_name = args.sample_name
    hap1_file = args.contig_hp1_fasta
    hap2_file = args.contig_hp2_fasta
    ALU_one_seq_fa = other_tools_path + "Alu_seq/ALU_one_seq.fa"
    assembly_based_variants_call_paf(hap1_file,hap2_file,ref_file,sample_name,out_dir)
    assembly_based_variants_call_sort(sample_name,out_dir)
    assembly_based_variants_call(sample_name,out_dir)
    Call_SV_del_from_contigs(sample_name,out_dir,SV_len)
    Call_SV_ins_from_contigs(sample_name,out_dir,SV_len)
    temp_del[1] = sample_name + "_temp_del_1.vcf"
    temp_ins[1] = sample_name + "_temp_ins_1.vcf"
    Call_alu_for_del(out_dir,out_dir,sample_name,ALU_one_seq_fa,temp_del[1])
    Call_alu_for_ins(out_dir,out_dir,sample_name,ALU_one_seq_fa,temp_ins[1])

    temp_del[2] = sample_name + "_temp_del_2.vcf"
    temp_ins[2] = sample_name + "_temp_ins_2.vcf"
    Call_repeats(out_dir,out_dir,sample_name,"DEL",temp_del[1],temp_del[2])
    Call_repeats(out_dir,out_dir,sample_name,"INS",temp_ins[1],temp_ins[2])

    extract_flanking_seq_del(out_dir,out_dir,sample_name)
    flanking_fasta = out_dir + sample_name + "_del_for_flankingseq.fasta"
    count = 0
    for monkey_name in species_list:
        monkey_ref = ref_list[count]
        align_to_monkeys(monkey_ref,flanking_fasta,out_dir,sample_name,monkey_name,"del")
        flanking_bam_sorted = out_dir + sample_name + "_" + monkey_name + "_del_for_flankingseq_sorted.bam"
        temp_del[count + 3] = sample_name + "_temp_del_" + str(count + 3) + ".vcf"
        evaluate_ancestral_del(flanking_bam_sorted,sample_name,monkey_name,out_dir,temp_del[count + 2], temp_del[count + 3])
        count += 1
    Popen("mv " + out_dir + temp_del[count + 2] + " " + out_dir + sample_name + "_final_del.vcf",shell=True).wait()
    extract_flanking_seq_ins(out_dir,out_dir,sample_name)
    flanking_fasta = out_dir + sample_name + "_ins_for_flankingseq.fasta"
    count = 0
    for monkey_name in species_list:
        monkey_ref = ref_list[count]
        align_to_monkeys(monkey_ref,flanking_fasta,out_dir,sample_name,monkey_name,"ins")
        flanking_bam_sorted = out_dir + sample_name + "_" + monkey_name + "_ins_for_flankingseq_sorted.bam"
        temp_ins[count + 3] = sample_name + "_temp_ins_" + str(count + 3) + ".vcf"
        evaluate_ancestral_ins(flanking_bam_sorted,sample_name,monkey_name,out_dir,temp_ins[count + 2],temp_ins[count + 3])
        count += 1
    
    Popen("mv " + out_dir + temp_ins[count + 2] + " " + out_dir + sample_name + "_final_ins.vcf",shell=True).wait()
    #os.sytems("rm " + out_dir + "*temp*.vcf")
    #os.sytems("rm " + out_dir + "*tmp*.vcf")
    #os.sytems("rm " + out_dir + "*_with_alu.vcf")
    #os.sytems("rm " + out_dir + "*_with_alu_with_repeats.vcf")
    #os.sytems("rm " + out_dir + "*.bam*")
    #os.sytems("rm " + out_dir + "*.sam")
    #os.sytems("rm " + out_dir + "*.p")

    


