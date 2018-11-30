# change permission
cd bin
chmod +x *.py
cd ..

# install trf
cd trf_tools
chmod +x trf
cd ..

# install minimap2
git clone https://github.com/lh3/minimap2
cd minimap2 && make
cd ..

# download the reference file (GRCh38)
wget http://xinzhouneuroscience.org/wp-content/uploads/2018/10/source.tar.gz
tar -xvf source.tar.gz
rm source.tar.gz

# download k8
wget https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 
tar xvjf k8-0.2.4.tar.bz2
rm k8-0.2.4.tar.bz2


if ! [ -x "$(command -v samtools)" ];
then
    echo 'Error: samtools is not installed...'
    echo 'using Conda to install...'
    conda install -c bioconda samtools
else
    echo 'using existing samtools...'
fi

if ! [ -x "$(command -v minimap2)" ];
then
    echo 'Error: minimap2 is not installed...'
    echo 'using Conda to install...'
    conda install -c bioconda minimap2
else
    echo 'using existing minimap2...'
fi

if ! [ -x "$(command -v bcftools)" ];
then
    echo 'Error: bcftools is not installed...'
    echo 'using Conda to install...'
    conda install -c bioconda bcftools
else
    echo 'using existing bcftools...'
fi

if ! [ -x "$(command -v tabix)" ];
then
    echo 'Error: htslib is not installed...'
    echo 'using Conda to install...'
    conda install -c bioconda htslib
else
    echo 'using existing htslib...'
fi


echo 'You have installed AncetralCall dependencies and downloaded the source files successfully!'
 



