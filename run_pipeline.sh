#!/bin/bash
# Note, change the directory to where you want these things to appear before running this code.

# Make directories
mkdir data
mkdir data/raw_data
mkdir data/fastqc_reports
mkdir data/fastqc_trimmed
mkdir data/aligned_reads
mkdir data/salmon_quants
mkdir results
mkdir results/map_qc

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

#-------------------- Quantify with Salmon ------------------------------------#
# Now with Salmon instead! And we will compare the performance
cd ../data

# Get the human transcriptome
curl ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.pc_transcripts.fa.gz \
-o gencode.v29.pc_transcripts.fa.gz

#curl ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.lncRNA_transcripts.fa.gz \
# -o gencode.v29.lncRNA_transcripts.fa.gz

# curl ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.abinitio.fa.gz \
# -o Homo_sapiens.GRCh38.cdna.abinitio.fa.gz

# Index the human transcriptome
# Using abinitio because .all. appears to be too much to index and the process ends 
# up being "Killed." The option --perfectHash can help with this, but only to a certain extent. 
salmon --threads=16 --no-version-check index \
-t gencode.v29.pc_transcripts.fa.gz \
-i human_index \
--perfectHash -k 31 \
 --type quasi
 
# Go to our trimmed read data
cd fastqc_trimmed

# Quantify each sample with salmon
for f in `ls *_1_val_1.fq.gz | sed 's/_1_val_1.fq.gz//' `
do
echo "Processing sample ${f}"
salmon quant -i ../human_index -l A \
-1 ${f}_1_val_1.fq.gz \
-2 ${f}_2_val_2.fq.gz \
-p 8 -o ../salmon_quants/${f}_quant
--gcBias --seqBias --biasSpeedSamp 5
done
