#!/usr/bin/env python
import pdb
#pdb.set_trace()
import pysam
from collections import defaultdict
import pickle
from argparse import ArgumentParser
import os
import sys
parser = ArgumentParser(description="Run depth all:")
parser.add_argument('--out_dir','-o_dir', help="Directory to store outputs", default='results_/')
parser.add_argument('--monkey_name','-mn', help="monkey name")
parser.add_argument('--sample_name','-sn', help="sample name")
parser.add_argument('--flanking_bam','-bam', help="flanking bam")
parser.add_argument('--in_vcf','-i_vcf', help="in_vcf")
parser.add_argument('--out_vcf','-o_vcf', help="out_vcf")

args = parser.parse_args()

def get_match_num_revised_3(cigar):
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

    indices_match = [i for i, x in enumerate(cigar_list) if x == "M"]
    indices_del = [i for i, x in enumerate(cigar_list) if x == "D"]
    match_list = []
    del_list = []
    for idx in indices_match:
        match_list.append(int(cigar_list[idx-1]))
    for idx in indices_del:
        del_list.append(int(cigar_list[idx-1]))

    return (match_list,del_list)


def check_clipping_not_in_the_end(cigar):
    flag = 0
    if "H" not in cigar and "S" not in cigar:
        flag = 1
    else:
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
        if cigar_list[-1] == "S" or cigar_list[-1] == "H":
            flag = 0
        else:
            flag = 1
    return flag

    
def check_clipping_not_at_the_beginning(cigar):
    flag = 0
    if "H" not in cigar and "S" not in cigar:
        flag = 1
    else:
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
        if cigar_list[1] == "S" or cigar_list[1] == "H":
            flag = 0
        else:
            flag = 1
    return flag


def Evaluate_derived_sv(bam_file,lib_prefix):
    sam_file = pysam.AlignmentFile(bam_file, "rb")
    flanking_seq_dict = defaultdict(lambda: defaultdict(list))
    align_score_list = []
    event_align_chr_dict = defaultdict()
    for read in sam_file.fetch():
        tags = read.get_tags()
        AS_field = [s for s in tags if "AS" in s]
        if AS_field != []:
            align_score = int(AS_field[0][1])
            align_score_list.append(align_score)
            qname = read.qname
            start_pos = read.pos
            flanking_idx = int(qname.split("_")[0])
            event_chr = qname.split("_")[1].split("chr")[-1]
            event_name = qname.split("_")[2]
            event_size = int(qname.split("_")[-1])
            read_cigar = read.cigarstring
            read_chr = read.reference_name
            if read.is_reverse:
                read_reverse_flag = 1
            else:
                read_reverse_flag = 0
            if flanking_idx in flanking_seq_dict[event_name]:
                read_chr_reform = ''.join([n for n in read_chr if n.isdigit()])
                event_chr_reform = ''.join([n for n in event_chr if n.isdigit()])
                if read_chr_reform == event_chr_reform:
                    max_align_score = flanking_seq_dict[event_name][flanking_idx][-1]
                    if align_score > max_align_score:
                        flanking_seq_dict[event_name][flanking_idx] = [read_chr, start_pos, event_size,read_reverse_flag, read_cigar,align_score]
            else:
                read_chr_reform = ''.join([n for n in read_chr if n.isdigit()])
                event_chr_reform = ''.join([n for n in event_chr if n.isdigit()])
                if read_chr_reform == event_chr_reform:
                    flanking_seq_dict[event_name][flanking_idx] = [read_chr, start_pos, event_size,read_reverse_flag, read_cigar,align_score]

    print(len(flanking_seq_dict))

    count = 0
    count_1 = 0
    count_2 = 0
    count_3 = 0
    derived_del = defaultdict(int)
    derived_ins = defaultdict(int)
    SV_ref_coord = defaultdict(list)
    event_cigar = defaultdict(list)
    for event_name, val in flanking_seq_dict.items():
        if len(val) == 2:
            flanking_seq_1_info = val[1]
            flanking_seq_2_info = val[2]
            read_chr = val[1][0]
            chr_num = flanking_seq_1_info[0]
            reverse_flag_1 = flanking_seq_1_info[3]
            reverse_flag_2 = flanking_seq_2_info[3]
            align_score_1 = flanking_seq_1_info[-1]
            align_score_2 = flanking_seq_2_info[-1]
            if reverse_flag_1 == 0 and reverse_flag_2 == 0 and align_score_1 > 0 and align_score_2 > 0:
                count_1 += 1
                sv_size = flanking_seq_1_info[2]
                cigar_1 = flanking_seq_1_info[4]
                cigar_2 = flanking_seq_2_info[4]
                flag_1 = check_clipping_not_in_the_end(cigar_1)
                flag_2 = check_clipping_not_at_the_beginning(cigar_2)
                if flag_1 == 1 and flag_2 == 1:
                    start_1 = flanking_seq_1_info[1]
                    match_list, del_list = get_match_num_revised_3(cigar_1)
                    total_num = sum(match_list) + sum(del_list)
                    end_1 = start_1 + total_num
                    
                    start_2 = flanking_seq_2_info[1]
                    match_list, del_list = get_match_num_revised_3(cigar_2)
                    end_2 = start_2 + sum(match_list) + sum(del_list)
                    event_cigar[event_name] = [cigar_1,cigar_2]
                    if abs(start_2 - end_1) <= 2:
                        derived_ins[event_name] =  1
                        SV_ref_coord[event_name] = [chr_num,start_1, end_2, reverse_flag_1] 
                        event_align_chr_dict[event_name] = read_chr
                    elif float(start_2 - end_1)/sv_size >= 0.9 and float(start_2 - end_1)/sv_size <= 1.1:
                        derived_del[event_name] = 1
                        SV_ref_coord[event_name] = [chr_num,start_1, end_2, reverse_flag_1] 
                        event_align_chr_dict[event_name] = read_chr
                            
            elif reverse_flag_1 == 1 and reverse_flag_2 == 1 and align_score_1 > 0 and align_score_2 > 0 :
                count_2 += 1
                sv_size = flanking_seq_1_info[2]
                cigar_1 = flanking_seq_1_info[4]
                cigar_2 = flanking_seq_2_info[4]
                flag_1 = check_clipping_not_in_the_end(cigar_2)
                flag_2 = check_clipping_not_at_the_beginning(cigar_1)
                if flag_1 == 1 and flag_2 == 1:
                    start_2 =  flanking_seq_2_info[1]  
                    match_list, del_list = get_match_num_revised_3(cigar_2)
                    total_num = sum(match_list) + sum(del_list)
                    end_2 = start_2 + total_num  # == end_1

                    start_1 = flanking_seq_1_info[1]   # == start_2
                    match_list, del_list = get_match_num_revised_3(cigar_1)
                    end_1 = start_1 + sum(match_list) + sum(del_list)
                    event_cigar[event_name] = [cigar_1,cigar_2]

                    if abs(start_1 - end_2) <= 2:
                        derived_ins[event_name] =  1
                        SV_ref_coord[event_name] = [chr_num,start_2, end_1,reverse_flag_1] 
                        event_align_chr_dict[event_name] = read_chr

                    elif float(start_1 - end_2)/sv_size >= 0.9 and float(start_1 - end_2)/sv_size <= 1.1:
                        derived_del[event_name] = 1
                        SV_ref_coord[event_name] = [chr_num,start_2, end_1, reverse_flag_1] 
                        event_align_chr_dict[event_name] = read_chr

            else:
                #print(reverse_flag_1,reverse_flag_2)
                count_3 += 1
        else:
            count += 1
    pickle.dump(derived_del, open(lib_prefix + "_derived_del_by_flanking_dict_for_ins.p","wb"))
    pickle.dump(derived_ins, open(lib_prefix+ "_derived_ins_by_flanking_dict_for_ins.p","wb"))
    #pickle.dump(SV_ref_coord, open(lib_prefix + "SV_ref_coord.p","wb"))
    #pickle.dump(event_cigar, open(lib_prefix+ "event_cigar_dict.p","wb"))
    #pickle.dump(event_align_chr_dict, open(lib_prefix + "event_align_chr_dict.p","wb"))



def add_dr_tag_into_vcf(dr_dict_flag,ins_vcf,ins_vcf2,monkey_name):
    f = open(ins_vcf,"r")
    fw = open(ins_vcf2,"w")
    for line in f:
        data = line.rsplit()
        if data[0][:2] == "##":
            fw.writelines(line)
        elif data[0][:2] == "#C":
            fw.writelines("##FORMAT=<ID=del_on_ref_" + monkey_name + ",Number=1,Type=String,Description=\"del_on_ref flag: 1 or 0\">\n")
            fw.writelines(line)
        else:
            data[8] = data[8] + ":del_on_ref_" + monkey_name
            event_name = data[2]
            if event_name in dr_dict_flag:
                dr_flag = 1
            else:
                dr_flag = 0
            data[9] = data[9] + ":" + str(dr_flag)
            fw.writelines("\t".join(data) + "\n")
    f.close()
    fw.close()


def add_it_tag_into_vcf(it_dict_flag,ins_vcf,ins_vcf2,monkey_name):
    f = open(ins_vcf,"r")
    fw = open(ins_vcf2,"w")
    for line in f:
        data = line.rsplit()
        if data[0][:2] == "##":
            fw.writelines(line)
        elif data[0][:2] == "#C":
            fw.writelines("##FORMAT=<ID=ins_on_tar_" + monkey_name + ",Number=1,Type=String,Description=\"ins_on_tar flag: 1 or 0\">\n")
            fw.writelines(line)
        else:
            data[8] = data[8] + ":ins_on_tar_" + monkey_name
            event_name = data[2]
            if event_name in it_dict_flag:
                it_flag = 1
            else:
                it_flag = 0
            data[9] = data[9] + ":" + str(it_flag)
            fw.writelines("\t".join(data) + "\n")
    f.close()
    fw.close()




if __name__ == "__main__":
    out_dir = args.out_dir
    in_vcf = args.in_vcf
    out_vcf = args.out_vcf
    if os.path.exists(out_dir):
        print("using existing output folder: " + out_dir)
    else:
        os.makedirs(out_dir)

    flanking_bam = args.flanking_bam
    monkey_name = args.monkey_name
    sample_name = args.sample_name
    ins_vcf = out_dir + in_vcf
    ins_vcf2 = out_dir + sample_name + "_tmp_ins.vcf"
    ins_vcf3 = out_dir + out_vcf
    lib_prefix = out_dir + sample_name + "_" + monkey_name
    Evaluate_derived_sv(flanking_bam,lib_prefix)
    dr_dict_flag = pickle.load(open(lib_prefix + "_derived_del_by_flanking_dict_for_ins.p","rb"))
    it_dict_flag = pickle.load(open(lib_prefix + "_derived_ins_by_flanking_dict_for_ins.p","rb"))
    add_dr_tag_into_vcf(dr_dict_flag,ins_vcf,ins_vcf2,monkey_name)
    add_it_tag_into_vcf(it_dict_flag,ins_vcf2,ins_vcf3,monkey_name)


