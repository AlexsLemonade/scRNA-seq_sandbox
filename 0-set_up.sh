#!/bin/bash

# Install all the necssary packages, modules, and software for this pipeline:

# Install necessary bash commands
sudo apt-get install unzip
sudo apt-get install curl

# Install pip 
apt install python-pip

# Install other R packages 
apt-get install libgsl0-dev
apt-get install libxml2-dev
apt-get install libgit2-dev
install.packages(“devtools”)

# Install java developer kit - needed for FASTQC
sudo apt-get update
sudo apt-get install default-jdk

# Change directory so that these install in "bin"
cd /bin

# Get FASTQC
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip
unzip fastqc_v0.11.8.zip

# Install cutadapt -needed for Trim_galore, will install to python packages and /root/.local/bin
pip install cutadapt

# Download and unzip TrimGalore
curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.4.5.tar.gz -o trim_galore.tar.gz
tar xvzf trim_galore.tar.gz

# Install g++ dependency so that STAR can be installed
sudo apt-get install g++
sudo apt-get install make

# Install STAR for genomic alignment
wget https://github.com/alexdobin/STAR/archive/2.6.0a.tar.gz
tar -xzf 2.6.0a.tar.gz
cd STAR-2.6.0a
cd STAR/source
make STAR

# Set PATH so everything can be accessed
export PATH=$PATH:/usr/lib/rstudio-server/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/lib/rstudio-server/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/root/.local/lib/python2.7/site-packages:/root/.local/bin
