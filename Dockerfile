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
    remotes \
  && R -e "BiocInstaller::biocLite(c('SRAdb', 'DBI', "slalom", "SCnorm"), suppressUpdates = TRUE)" 

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
	numpy 

# Install Seurat for normalization
RUN R -e "remotes::install_github("UCSF-TI/fake-hdf5r");install.packages('Seurat')"

# Get FASTQC
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip && \
    unzip fastqc_v0.11.8.zip && \
    chmod 755 FastQC/fastqc 

ENV PACKAGES git gcc make g++ cmake libboost-all-dev liblzma-dev libbz2-dev \
    ca-certificates zlib1g-dev curl unzip autoconf
ENV SALMON_VERSION 0.9.1

# salmon binary will be installed in /home/salmon/bin/salmon
# don't modify things below here for version updates etc.

WORKDIR /home

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
