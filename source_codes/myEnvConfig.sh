#!/bin/bash
PREFIX=$1
echo "Your working directory is $PREFIX"

curl https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh -o $PREFIX/miniconda.sh
chmod +x $PREFIX/miniconda.sh
bash $PREFIX/miniconda.sh -b -p $PREFIX/miniconda2
sudo ln -s $PREFIX/miniconda2/etc/profile.d/conda.sh /etc/profile.d/conda.sh
 
source $PREFIX/miniconda2/bin/activate
conda install -c bioconda bwa=0.7.15 ensembl-vep=96.0 fgbio=0.8.1 picard=2.20.3 pysam=0.15.0.1 samtools=1.9 vardict=2019.06.04 -y
conda install -c r r=3.6.0 -y
#By default, VEP searches for caches in $HOME/.vep; to use a different directory when running VEP, use --dir_cache
#mkdir $PREFIX/vep/
#curl ftp://ftp.ensembl.org/pub/release-97/variation/indexed_vep_cache/homo_sapiens_vep_97_GRCh38.tar.gz -o $PREFIX/vep/homo_sapiens_vep_97_GRCh38.tar.gz
#tar xzf $PREFIX/vep/homo_sapiens_vep_97_GRCh38.tar.gz
 