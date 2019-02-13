FROM rocker/tidyverse:3.4.3
MAINTAINER ccdl@alexslemonade.org
WORKDIR /rocker-build/

RUN apt-get update && apt-get install -y --no-install-recommends apt-utils
RUN apt-get install dialog apt-utils -y

RUN apt-get update -qq && apt-get -y --no-install-recommends install \
  build-essential \
  libxml2-dev \
  libsqlite-dev \
  libmariadbd-dev \
  libmariadbclient-dev \
  libmariadb-client-lgpl-dev \
  libpq-dev \
  libssh2-1-dev \
  pandoc

# Need this so Seurat will install
RUN R -e "devtools::install_github('UCSF-TI/fake-hdf5r', ref = 'c23358f4dd8b8135b9c60792de441eca3d867eba')"

RUN apt update && apt install -y dirmngr curl bash
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
  && install2.r --error \
  --deps TRUE \
  rjson \
  ggpubr \
  Seurat \
  optparse \
  tidyverse \
  dplyr \
  devtools \
  optparse \
  fastqcr \
  rjson \
  caret \
  Rtsne \
  NMI \
  && R -e "source('https://bioconductor.org/biocLite.R')" \
  && R -e "BiocInstaller::biocLite(c('ensembldb', 'DESeq2', 'qvalue', 'org.Hs.eg.db', 'org.Dr.eg.db', 'ComplexHeatmap', 'ConsensusClusterPlus', 'SRAdb', 'DBI', 'limma', 'edgeR', 'scran', 'scater'), suppressUpdates = TRUE)"

# Install R packages from github and urls
# Need most updated version of tximport so AlevinQC will install later
RUN R -e "devtools::install_github('mikelove/tximport', ref = 'b5b5fe11b0093b4b2784f982277b2aa66d2607f7')"
RUN R -e "devtools::install_github('csoneson/alevinQC', ref = '1fdf1c14b59eead3e239d7f99b607d59753e9420', dependencies = TRUE)"

# FastQC
RUN apt update && apt install -y fastqc

ENV PACKAGES git gcc make g++ libboost-all-dev liblzma-dev libbz2-dev \
    ca-certificates zlib1g-dev curl unzip autoconf

ENV SALMON_VERSION 0.12.0

# salmon binary will be installed in /home/salmon/bin/salmon
# don't modify things below here for version updates etc.

WORKDIR /home

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

# Apt doesn't have the latest version of cmake, so install it using their script.
RUN apt remove cmake cmake-data -y

RUN wget --quiet https://github.com/Kitware/CMake/releases/download/v3.13.3/cmake-3.13.3-Linux-x86_64.sh && \
    mkdir /opt/cmake && \
    sh cmake-3.13.3-Linux-x86_64.sh --prefix=/opt/cmake --skip-license && \
    sudo ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake

RUN curl -k -L https://github.com/COMBINE-lab/salmon/archive/v${SALMON_VERSION}.tar.gz -o salmon-v${SALMON_VERSION}.tar.gz && \
    tar xzf salmon-v${SALMON_VERSION}.tar.gz && \
    cd salmon-${SALMON_VERSION} && \
    mkdir build && \
    cd build && \
    cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local && make && make install

RUN apt-get update -y && apt-get install -y wget && \
    apt-get clean && rm -rf /var/cache/apt/archives/* /var/lib/apt/lists/* && \
    wget -P /usr/bin "https://raw.githubusercontent.com/inutano/pfastq-dump/master/bin/pfastq-dump" && \
    chmod +x /usr/bin/pfastq-dump && \
    wget -P / "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.2/sratoolkit.2.9.2-centos_linux64.tar.gz" && \
    tar zxf sratoolkit.2.9.2-centos_linux64.tar.gz && \
    cp -r sratoolkit.2.9.2-centos_linux64/bin/* /usr/bin && \
    rm -fr sratoolkit.2.9.2-centos_linux64*