# AncestralCall
Assemblies based structral variants calling for multiple samples, including anotation for Alu, Tandem Repeats, Ancestral call. 
output: merged.vcf for multiple samples, and single vcf for each sample. 
## Download:
```
git clone https://github.com/maiziex/AncestralCall.git
```

## Install:
```
cd AncestralCall
chmod +x install.sh
./install.sh
```

## Running The Code:
Put the "Ancestral/bin" in the ".bashrc" file, and source the ".bashrc" file
Or use the fullpath of "Run_all_samples.py"


### Step 1: Generate "merged.vcf" for all samples, and single vcf for each sample
```
Run_all_samples.py --ref_dir ./GRCh38_reference/ --in_dir ./diploid_contig/ --out_dir ./test/ --ref_file ./GRCh38_reference/genome.fa  --SV_len 20  --species_name_list "Chimp","Orang" --species_ref_list "pan_troglodytes_ref.fasta","pongo_abelii_ref.fasta" --num_threads 10 --sample_list 'HG00250','HG00353','HG00512'
```
--ref_dir: "./GRCh38_reference/" is the folder to store all reference files which would be downloaded by running "./install.sh".  <br />
--in_dir: "./diploid_contig/" is the input folder where you store the diploid assembled contig files. <br />
--out_dir: "./results/" is the folder name you can define to store the final results. <br />
--ref_file: "./GRCh38_reference/genome.fa" is the human reference fasta file which would be download by running "./install.sh"
--SV_len: "20" is the SV size you can define.  <br />
--species_name_list: "Chimp","Orang" are the species's name you can define. Each name is seperately by comma (","). <br />
--species_ref_list: "pan_troglodytes_ref.fasta","pongo_abelii_ref.fasta" are the reference fasta files for each species which you defined in the "--species_name_list", respectivley. Each reference file is seperately by comma (",") <br />
--num_threads: "10" is the number of threads you can define, which corresponds to number of samples.  <br />
--sample_list: 'HG00250','HG00353','HG00512' are the sample names corresponding to your contig files, which is the prefix of the contig files. <br />

### Step 2: Generate all the multiple-alignment files for each sample. 

