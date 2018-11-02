# change permission
cd bin
chmod +x Run_all_samples.py
cd ..

# install trf
cd trf_tools
chmod +x trf
cd ..

# install minimap2
git clone https://github.com/lh3/minimap2
cd minimap2 && make
cd ..

# download k8
curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf -

# download the reference file (GRCh38)
wget http://xinzhouneuroscience.org/wp-content/uploads/2018/10/source.tar.gz
tar -xvf source.tar.gz

conda install -c bioconda samtools 
conda install -c bioconda bcftools 
conda install -c bioconda htslib 





echo 'You have installed AncetralCall and the dependencies successfully!'
 



