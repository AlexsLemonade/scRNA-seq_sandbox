#!/bin/bash
# Note, change the directory to where you want these things to appear before running this code.

# Make directories
mkdir data
mkdir data/raw_data
mkdir data/fastqc_reports
mkdir data/fastqc_trimmed
mkdir data/aligned_reads
mkdir data/sorted_bams
mkdir data/fpkms
mkdir results
mkdir results/map_qc
mkdir results/salmon_quants
#--------------------------- Download fastq data-------------------------------#
# Note: because running all ~3800 samples takes quite a bit of time, this 
# example run is set to only run 10 samples. When you want to run the full set, 
# delete the -n argument.  
Rscript 0-download_fastq_data.R \
-i SRP079058 \
-d data/raw_data \
-n 20 

#------------------------Run fastqc and get reports----------------------------#
# Run sequence quality control with FASTQC
/FastQC/fastqc data/raw_data/* --outdir data/fastqc_reports

# Obtain summary report:
Rscript 1-get_fastqc_reports.R -d data/fastqc_reports -o results

#--------------Paired trimming of the adapters from these data:----------------#
# Go to the raw_data directory 
cd data/raw_data

# Do adapter trimming for all sample file pairs
for f in `ls *_1.fastq.gz | sed 's/_1.fastq.gz//'`
do
/TrimGalore-0.4.5/trim_galore \
--nextera \
--paired ${f}_1.fastq.gz ${f}_2.fastq.gz \
-o ../fastqc_trimmed
done

#--------------- Align trimmed reads to genome using HISAT2--------------------#
# Get premade human genome index:
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz -P ../
tar data/grch38.tar.gz

# Go to trimmed reads directory
cd ../fastqc_trimmed

# Do genome alignment for all sample pairs
for f in `ls -1 *_1_val_1.fq.gz | sed 's/_1_val_1.fq.gz//' `
do
/hisat2-2.1.0/hisat2 -x ../grch38/genome \
-1 ${f}_1_val_1.fq.gz \
-2 ${f}_2_val_2.fq.gz \
-S ../aligned_reads/${f}.bam
done

#--------------------Further prep the bam files with samtools------------------#
# This is a temporary fix until I figure out how to make the docker file do this
cd /samtools-1.3.1
make prefix=/bin/bash

cd ../aligned_reads
# Sort and index bam files
for f in `ls *.bam | sed 's/.bam//' `
do
/samtools-1.3.1/samtools sort ${f}.bam -o ${f}.sorted.bam
/samtools-1.3.1/samtools index ${f}.sorted.bam -o ${f}
done

# Can sort files in parallel to make it a tad faster (but will still take a couple hours)
# ls aligned_reads/*.bam | sed 's/.bam//' | parallel "/samtools-1.3.1/samtools sort {.}.bam -o {.}.sorted.bam"

#-------------------- Quantify with Salmon ------------------------------------#
# Now with Salmon instead! And we will compare the performance
# Get the human transcriptome
curl ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -o human.fa.gz

# Index the human transcriptome
salmon index -t human.fa.gz -i human_index

# Go to our trimmed read data
cd data/trimmed_reads

# Quantify each sample with salmon
for f in `ls *_1_val_1.fq.gz | sed 's/_1_val_1.fq.gz//' `
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -i human_index -l A \
-1 ${fn}/${samp}_1_val_1.fq.gz \
-2 ${fn}/${samp}_2_val_2.fq.gz \
-p 8 -o salmon_quants/${samp}_quant
done
