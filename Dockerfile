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

# Install cellranger
RUN cd /tmp/ && \
	wget -O cellranger-3.0.2.tar.gz "http://cf.10xgenomics.com/releases/cell-exp/cellranger-3.0.2.tar.gz?Expires=1550715688&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cDovL2NmLjEweGdlbm9taWNzLmNvbS9yZWxlYXNlcy9jZWxsLWV4cC9jZWxscmFuZ2VyLTMuMC4yLnRhci5neiIsIkNvbmRpdGlvbiI6eyJEYXRlTGVzc1RoYW4iOnsiQVdTOkVwb2NoVGltZSI6MTU1MDcxNTY4OH19fV19&Signature=QXtAQYaZ5IJdrzA1jp78ETZWOEr2wzZltEYD2J8J~-apqM4E2z1P1P1QNMKkXSGGvWO5vx8jCDJDpluNw7M1jk9KzNMtmn1blQ0fPd8w-FRxkB5CicER1Me6oFaiepiHr3TV2vhFZlIfAvAytIBAL9EvTFgqJa8IYuYTpLOvbWhZrhk4nsG8bAoAD-RqYAHQcB2yyzoklOoJaKW1ihmL1def80az~xTLm3pmswyYhNNYlvhlbyVy2ameb1r5XqfKwW0wjp28Uj2c5cTs2GVKa9iWOLdOMwo5V5s6Dxi7NgZoE3JsPr44HfkCMbFFAKLy6KqIseXP7cO38uOXbYondQ__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA" && \	
	mv cellranger-3.0.2.tar.gz /opt/ && \
	cd /opt/ && \
	tar -xzvf cellranger-3.0.2.tar.gz && \
	rm -f cellranger-3.0.2.tar.gz

ENV PATH /opt/cellranger-3.0.2:$PATH

# Install bamtofastq convert
RUN wget http://cf.10xgenomics.com/misc/bamtofastq && \
    chmod 700 bamtofastq

