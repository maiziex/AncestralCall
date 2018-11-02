# install minimap2
git clone https://github.com/lh3/minimap2
cd minimap2 && make
cd ..

# download k8
curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf -

# download the reference file (GRCh38)
wget http://xinzhouneuroscience.org/wp-content/uploads/2018/10/GRCh38_reference.tar.gz
tar -xvf GRCh38_reference.tar.gz

conda install -c bioconda bcftools 




echo 'You have installed AncetralCall and the dependencies successfully!'
 



