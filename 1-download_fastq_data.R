# CCDL ALSF 2018
# C. Savonen 
# 
# Get raw reads for single cell glioblastoma data from Darmanis et al, 2017
#
#-------------------------- Get necessary packages-----------------------------#
if(!("SRAdb" %in% installed.packages())){
  source('http://bioconductor.org/biocLite.R')
  biocLite('SRAdb', dependencies = TRUE)
}
if(!("DBI" %in% installed.packages())){
  source('http://bioconductor.org/biocLite.R')
  biocLite('DBI')
}
library(SRAdb)
library(DBI)

#------------------- Connect to NCBI's SRA SQL database------------------------#
srafile = getSRAdbFile()
con = dbConnect(RSQLite::SQLite(), srafile)

# Get a list of the samples associated with the project we are interested in
files <- listSRAfile('SRP079058', con)

#-------------------------Get all the FASTQ files------------------------------#
if (!dir.exists("raw_data")){
  dir.create("raw_data")
}
getFASTQfile('SRP079058', con, destDir = "raw_data")
