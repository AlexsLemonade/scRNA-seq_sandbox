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

RUN apt-get install unzip
RUN apt-get install curl
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get -qq install default-jdk
RUN apt-get install g++
RUN apt-get install make
RUN apt-get install -y python-pip libpq-dev python-dev

# Install the actual packages we will be using
RUN pip install linux-utils
RUN pip install cutadapt

# Get FASTQC
RUN  wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip
RUN  unzip fastqc_v0.11.8.zip
RUN chmod 755 /bin/FastQC/fastqc

# Download and unzip TrimGalore
RUN curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.4.5.tar.gz -o trim_galore.tar.gz
RUN tar xvzf trim_galore.tar.gz
RUN chmod 755 /bin/TrimGalore-0.4.5/trim_galore

# Install STAR for genomic alignment
RUN wget https://github.com/alexdobin/STAR/archive/2.6.0a.tar.gz
RUN tar -xzf 2.6.0a.tar.gz