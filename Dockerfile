# adapted from rocker/tidyverse
FROM rocker/rstudio:3.5.1

ENV PATH=/usr/lib/rstudio-server/bin:/root/.local/lib/python2.7/site-packages:/root/.local/bin:$PATH

RUN apt-get update && apt-get install -y --no-install-recommends apt-utils

RUN apt-get update -qq && apt-get -y --no-install-recommends install \
  libxml2-dev \
  libcairo2-dev \
  libsqlite-dev \
  libmariadbd-dev \
  libmariadb-client-lgpl-dev \
  libpq-dev \
  libssh2-1-dev \
  libsqlite-dev \
  python-dev \
  && R -e "source('https://bioconductor.org/biocLite.R')" \
  && install2.r --error \
    --deps TRUE \
    tidyverse \
    dplyr \
    devtools \
    optparse \
    fastqcr \
  && R -e "BiocInstaller::biocLite(c('SRAdb', 'DBI'), suppressUpdates = TRUE)"

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
RUN pip install \
	linux-utils \
	cutadapt

# Get FASTQC
RUN  wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip
RUN  unzip fastqc_v0.11.8.zip
RUN chmod 755 FastQC/fastqc

# Download and unzip TrimGalore
RUN curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.4.5.tar.gz -o trim_galore.tar.gz
RUN tar xvzf trim_galore.tar.gz
RUN chmod 755 TrimGalore-0.4.5/trim_galore

# Install HISAT2 for genomic alignment
RUN wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip
RUN unzip hisat2-2.1.0-Linux_x86_64.zip

# Install Salmon
RUN wget https://cmake.org/files/v3.8/cmake-3.8.2-Linux-x86_64.tar.gz
RUN tar -zxvf cmake-3.8.2-Linux-x86_64.tar.gz
RUN export PATH=<SALMON_INSTALLATION_DIRECTORY>/cmake-3.8.2-Linux-x86_64/bin:$PATH
RUN cd cmake-3.8.2
RUN ./configure --prefix=INSTALLATION_DIRECTORY
RUN make
RUN make install
# RUN sudo apt-get install cmake
