# adapted from rocker/tidyverse
FROM rocker/rstudio:3.5.1

ENV PATH=/usr/lib/rstudio-server/bin:/root/.local/lib/python2.7/site-packages:/root/.local/bin:$PATH

RUN apt-get update && apt-get install -y --no-install-recommends apt-utils

RUN apt-get update -qq && apt-get -y --no-install-recommends install \
  build-essential \
  libxml2-dev \
  libcairo2-dev \
  libsqlite-dev \
  libmariadbd-dev \
  libmariadb-client-lgpl-dev \
  libpq-dev \
  libssh2-1-dev \
  libsqlite-dev \
  python-dev \
  parallel \
  && R -e "source('https://bioconductor.org/biocLite.R')" \
  && install2.r --error \
    --deps TRUE \
    tidyverse \
    dplyr \
    devtools \
    optparse \
    fastqcr \
    rbamtools \
    rjson \
    caret \
    Rtsne \
    NMI \
  && R -e "BiocInstaller::biocLite(c('SRAdb', 'DBI', 'DESeq2', 'limma', 'edgeR'), suppressUpdates = TRUE)" 

# Install Rsubread by itself
RUN R -e 'BiocInstaller::biocLite("Rsubread")'

# Install other things
RUN apt-get update && apt-get install -y \
	dialog \
	unzip \
	curl \
	g++ \
	make 

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get -qq install default-jdk


RUN apt-get install -y python-pip libpq-dev python-dev

# Install the actual packages we will be using
# Cutadapt is needed for TrimGalore
# RSeQC is used for mapping quality control
RUN pip install \
	linux-utils \
	numpy \
	cutadapt \
	RSeQC

# Get FASTQC
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip && \
    unzip fastqc_v0.11.8.zip && \
    chmod 755 FastQC/fastqc 

# Download and unzip TrimGalore
RUN curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.4.5.tar.gz -o trim_galore.tar.gz && \
    tar xvzf trim_galore.tar.gz && \
    chmod 755 TrimGalore-0.4.5/trim_galore 

# Install HISAT2 for genomic alignment
RUN wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip && \
    unzip hisat2-2.1.0-Linux_x86_64.zip

# Set up samtools
RUN apt-get update && apt-get install -y --no-install-recommends libncurses5-dev && \
    wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 -O samtools.tar.bz2 && \
    tar -xjvf samtools.tar.bz2 && \
    cd samtools-1.3.1 && \
    make && make prefix=/usr/local/bin install

ENV PACKAGES git gcc make g++ cmake libboost-all-dev liblzma-dev libbz2-dev \
    ca-certificates zlib1g-dev curl unzip autoconf
ENV SALMON_VERSION 0.9.1

# salmon binary will be installed in /home/salmon/bin/salmon
# don't modify things below here for version updates etc.

WORKDIR /home

# installing cellRanger 
RUN curl -o cellranger-3.0.2.tar.gz "http://cf.10xgenomics.com/releases/cell-exp/cellranger-3.0.2.tar.gz?Expires=1547783967&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cDovL2NmLjEweGdlbm9taWNzLmNvbS9yZWxlYXNlcy9jZWxsLWV4cC9jZWxscmFuZ2VyLTMuMC4yLnRhci5neiIsIkNvbmRpdGlvbiI6eyJEYXRlTGVzc1RoYW4iOnsiQVdTOkVwb2NoVGltZSI6MTU0Nzc4Mzk2N319fV19&Signature=AzHE9HfWkBgaB9USR9Xumx2Smc7dRwv0NlnTBu6zvnlNpnkET09vm5-SgHOQlPbb4Ga-iZfKzZAm6CnpROmWKhG-TxSOagFnApG10NQJ6HjSzWU6ODvrsO6PQD5ROpMKT4ZmFoYBpxuz4jjcbUa0qgeK220qrbwO3U9bU1kDrI-qhWO3RQSRULVHeiPA4Dc01xEvYKHh7Q2L3a2pdhAlwI5AErxr3H-AddBunZOkXog4FgyW0PaL37iHhhkYaDdLv2QmbRqM1QgKZXYRSZBf4esp8Y46uDhwzS-vYA-fSgoGIAgPUt8AZNLBTedhPsF2NJvGGsek6DznZ9Dv557eZA__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA" && \
curl -O http://cf.10xgenomics.com/supp/cell-exp/ && \
refdata-cellranger-GRCh38-3.0.0.tar.gz && \
export PATH=/opt/cellranger-3.0.2:$PATH

RUN apt update && \
    apt install -y --no-install-recommends ${PACKAGES} && \
    apt clean

RUN curl -k -L https://github.com/COMBINE-lab/salmon/archive/v${SALMON_VERSION}.tar.gz -o salmon-v${SALMON_VERSION}.tar.gz && \
    tar xzf salmon-v${SALMON_VERSION}.tar.gz && \
    cd salmon-${SALMON_VERSION} && \
    mkdir build && \
    cd build && \
    cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local && make && make install

ENV PATH /home/salmon-${SALMON_VERSION}/bin:${PATH}

RUN salmon --version
