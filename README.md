# AncestralCall
Assemblies based structral variants calling for multiple samples, including anotation for Alu, Tandem Repeats, Ancestral call. 
output: merged.vcf for multiple samples, and single vcf for each sample. 

## Download:
```
git clone https://github.com/maiziex/AncestralCall.git
```

## Dependencies:
Ancestral utilizes <a href="https://www.python.org/downloads/">Python3</a>, <a href="https://github.com/lh3/minimap2/tree/master/misc">paftools (Called haploid assemblies based variants)</a>, <a href="https://tandem.bu.edu/trf/trf.html">trf (Tandem Repeats Finder)</a>, <a href="http://samtools.sourceforge.net/">SAMtools</a>, and <a href="https://github.com/lh3/minimap2">minimap2</a>. To be able to execute the above programs by typing their name on the command line, the program executables must be in one of the directories listed in the PATH environment variable (".bashrc"). <br />
Or you could just run "./install.sh" to install them, but make sure you have installed "conda" and "wget" first. 

## Install:
```
cd AncestralCall
chmod +x install.sh
./install.sh
```

## source folder:
After running "./install.sh", a folder "source" would be download, it includes human GRCh38 reference fasta file, and reference fasta files for Orangutan and Chimpanzee which you can use for ancestral call. You could also just download them by yourself from the corresponding official websites. 

## Running The Code:
Put the "Ancestral/bin" in the ".bashrc" file, and source the ".bashrc" file <br />
Or use the fullpath of "Run_all_samples.py"


### Step 1: Generate "merged.vcf" for all samples, and single vcf for each sample
```
Run_all_samples.py --in_dir ./diploid_contig/ --out_dir ./results/ --ref_file ./source/genome.fa  --SV_len 20  --species_name_list "Chimp","Orang" --species_ref_list "./source/pan_troglodytes_ref.fasta","./source/pongo_abelii_ref.fasta" --num_threads 10 --sample_list 'HG00250','HG00353','HG00512'
```

--in_dir: Required parameter. <br />
"./diploid_contig/" is the input folder where you store the diploid assembled contig files. <br />
<br />

--out_dir: Optional parameter, default = ./Ancestral_results/  <br />
"./results/" is the folder name you can define to store the final results.  <br />
<br />

--ref_file: Required augument. <br />
"./GRCh38_reference/genome.fa" is the human reference fasta file which can be download by running "./install.sh". <br />
<br /> 

--SV_len: Optional parameter, default = 20 <br />
"20" is the SV size you can define.<br />
--species_name_list: Required parameter. <br />
"Chimp","Orang" are the species's name you can define. Each name is seperately by comma (","). <br />
--species_ref_list: Required parameter.<br />
"pan_troglodytes_ref.fasta","pongo_abelii_ref.fasta" are the reference fasta files for each species which you defined in the "--species_name_list", respectivley. Each reference file is seperately by comma (",") <br />
--num_threads: Optional parameter, default = 2 <br />
"10" is the number of threads you can define, which corresponds to number of samples. default = 2 <br />
--sample_list: Required parameter. <br />
'HG00250','HG00353','HG00512' are the sample names corresponding to your contig files, which is the prefix of the contig files. <br />

### Step 2: Generate all the multiple-alignments files for each sample. 

