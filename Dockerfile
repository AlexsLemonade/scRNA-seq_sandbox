# adapted from rocker/tidyverse
FROM rocker/rstudio:3.5.1

ENV PATH=/usr/lib/rstudio-server/bin:/root/.local/lib/python2.7/site-packages:/root/.local/bin:$PATH

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
RUN pip install cutadapt
RUN apt-get install unzip
RUN apt-get install curl
RUN apt-get install default-jdk
RUN apt-get install g++
RUN apt-get install make
