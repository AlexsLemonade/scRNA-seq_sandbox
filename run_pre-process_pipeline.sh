#!/bin/bash
# C. Savonen
# CCDL for ALSF 2018

# Purpose: running the pre-processing steps for single cell RNA-seq data.

# Change your directory name, GEO ID, and SRP here. Then run the script.
dir=darmanis_data
GSE=GSE84465
SRP=SRP079058
label=darmanis

#----------------------------Make directories---------------------------------#
mkdir ${dir}
mkdir ${dir}/raw_data
mkdir ${dir}/aligned_reads
mkdir ${dir}/fastqc_reports
mkdir ${dir}/salmon_quants
mkdir results
mkdir ref_files

#--------------Will run setup only if it hasn't been ran before----------------#
# Will check for genome index first before running
if [ ! -e ref_files/human_index ]; then

  # Get the human transcriptome
  curl ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz \
    -o ref_files/Homo_sapiens.GRCh38.cdna.all.fa.gz

  # Index the human transcriptome
  salmon --threads=16 --no-version-check index \
    -t ref_files/Homo_sapiens.GRCh38.cdna.all.fa.gz \
    -i ref_files/human_index \
    -k 23
fi

#--------------------------- Download fastq data-------------------------------#
# Note: because running all ~3800 samples takes quite a bit of time, use -n option
# to set what number of randomly selected samples you would like
Rscript scripts/pre-processing/0-get_sample_download_list.R \
  -i ${SRP} \
  -o results \
  -d ${dir}/salmon_quants \
  -q ref_files/SRAmetadb.sqlite

# Change directory to the data
cd ${dir}

# Download each sample and run Salmon on it. Then remove the sample to save room.
for line in `cat ../results/files.to.download.txt`
  do
  # Download forward and reverse fastq files
  Rscript ../scripts/pre-processing/1-download_sra.R \
    -s $line \
    -q ../ref_files/SRAmetadb.sqlite \
    -d raw_data

# Run sequence quality control with FASTQC
/FastQC/fastqc raw_data/* --outdir fastqc_reports

  # For each fastq file pair run salmon for quantfication
  for f in `ls raw_data/*_1.fastq.gz | sed 's/_1.fastq.gz//' `
    do
    echo "Processing sample ${f}"

    # Run Salmon
    salmon quant -i ../ref_files/human_index -l A \
      -1 ${f}_1.fastq.gz \
      -2 ${f}_2.fastq.gz \
      -p 8 -o salmon_quants/$line \
      --gcBias --seqBias --biasSpeedSamp 5
    done
  rm raw_data/*
done

#-------------------------Obtain summary report of fastqc:---------------------#
cd ..

Rscript scripts/pre-processing/2-get_fastqc_reports.R \
  -d ${dir}/fastqc_reports \
  -o results

#-------------Make a gene matrix out of the Salmon quantification data---------#
Rscript scripts/pre-processing/3-make_gene_matrix.R \
  -d ${dir}/salmon_quants \
  -o ${dir} \
  -m 0.5 \
  -l ${dir}

#---------------------------Prep the data with this Rmd------------------------#
Rscript -e "rmarkdown::render('darmanis_data_prep.Rmd')"

#-------------------------------Run normalization------------------------------#
Rscript scripts/5-run_normalization.R \
  -d ${dir}/normalized_${label}/counts_${label}.tsv \
  -a all \
  -o ${dir}/normalized_${label} \
  -l ${label}
