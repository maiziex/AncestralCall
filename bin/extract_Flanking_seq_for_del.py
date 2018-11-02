#!/usr/bin/env python
import pdb
#pdb.set_trace()
import pickle
from collections import defaultdict
import os
import sys
from argparse import ArgumentParser

parser = ArgumentParser(description="Run depth all:")
parser.add_argument('--in_dir','-i_dir', help="Directory for inputs")
parser.add_argument('--out_dir','-o_dir', help="Directory to store outputs")
parser.add_argument('--ref_dir','-r_dir', help="Directory to store reference")
parser.add_argument('--sample_name','-s', help="sample name")

args = parser.parse_args()



def Extract_validated_SV_old(input_file,output_file,chr_num):
    f = open(input_file,"r")
    val_sv = defaultdict(list)
    for line in f:
        data = line.rsplit()
        #if data[9] != "0/0":
        print(data[9])
        cur_chr_num = data[1]
        if cur_chr_num == "chr" + str(chr_num):
            cur_name = data[0]
            _start = int(data[2])
            _end = int(data[3])
            _start_2 = _start - 500
            _end_2 = _end + 500
            val_sv[cur_name] = [_start, _end, _start_2, _end_2]
        #else:
            #print(data)
    pickle.dump(val_sv, open(output_file,"wb"))
    #print("done~")


def Extract_validated_SV(input_file,output_file,chr_num):
    f = open(input_file,"r")
    val_sv = defaultdict(list)
    for line in f:
        data = line.rsplit()
        if data[0][0] != "#":
            cur_chr_num = data[0]
            if cur_chr_num == "chr" + str(chr_num):
                cur_name = data[2]
                _start = int(data[1])
                _end = int(data[1]) + len(data[3])
                _start_2 = _start - 500
                _end_2 = _end + 500
                val_sv[cur_name] = [_start, _end, _start_2, _end_2]
    pickle.dump(val_sv, open(output_file,"wb"))
    f.close()


def write_validated_SV_fasta(input_file,ref_chr_dict,fw,chr_num):
    val_sv = pickle.load(open(input_file,"rb"))
    for key, val in val_sv.items():
        #print(val)
        _start = int(val[0])
        _end = int(val[1])
        _start_2 = int(val[2])
        _end_2 = int(val[3])
        seq_1_flankingseq = ref_chr_dict[_start_2:_start]
        seq_2_flankingseq = ref_chr_dict[_end:_end_2]
        fw.writelines(">1_chr" + str(chr_num) + "_" + key + "_" + str(_end - _start) + "\n")  # ref
        fw.writelines(seq_1_flankingseq + "\n")
        fw.writelines(">2_chr" + str(chr_num) + "_" + key + "_" + str(_end - _start) + "\n")  # del
        fw.writelines(seq_2_flankingseq + "\n")

    return fw





if __name__ == "__main__":
    in_dir = args.in_dir
    out_dir = args.out_dir
    ref_dir = args.ref_dir
    if os.path.exists(out_dir):
        print("using existing output folder: " + out_dir)
    else:
        os.makedirs(out_dir)
    sample_name = args.sample_name
    flanking_fasta = out_dir + sample_name + "_del_for_flankingseq.fasta"
    del_vcf = in_dir + sample_name + "_DEL_with_alu_with_repeats.vcf"
    fw = open(flanking_fasta,"w")
    for ii in range(1,23):
        print("chr" + str(ii))
        ref_chr_file = ref_dir + "ref_seq_chr" + str(ii) + ".p"
        ref_chr_dict = pickle.load(open(ref_chr_file,"rb"))
        flanking_chr_file = out_dir + sample_name + "_del_for_flankingseq_chr" + str(ii) + ".p"
        Extract_validated_SV(del_vcf,flanking_chr_file,ii)
        fw = write_validated_SV_fasta(flanking_chr_file,ref_chr_dict,fw,ii)


    print("chrX")
    ref_chr_file = ref_dir + "ref_seq_chr23.p"
    ref_chr_dict = pickle.load(open(ref_chr_file,"rb"))
    flanking_chr_file = out_dir + sample_name + "_del_for_flankingseq_chr23.p"
    Extract_validated_SV(del_vcf,flanking_chr_file,"X")
    fw = write_validated_SV_fasta(flanking_chr_file,ref_chr_dict,fw,"X")
    fw.close()






