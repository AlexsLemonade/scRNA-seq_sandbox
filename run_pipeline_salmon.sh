#!/bin/bash
# Note, change the directory to where you want these things to appear before running this code.

# Make directories
mkdir data
mkdir data/raw_data
mkdir data/fastqc_reports
mkdir data/salmon_quants
mkdir results

#--------------------------- Download fastq data-------------------------------#
# Note: because running all ~3800 samples takes quite a bit of time, use -n option 
# to set what number of randomly selected samples you would like  
Rscript 0-download_fastq_data.R \
-i SRP079058 \
-d data/raw_data 

#------------------------Run fastqc and get reports----------------------------#
# Run sequence quality control with FASTQC
/FastQC/fastqc data/raw_data/* --outdir data/fastqc_reports

# Obtain summary report:
Rscript 1-get_fastqc_reports.R -d data/fastqc_reports -o results

#-------------------- Quantify with Salmon ------------------------------------#
# Now with Salmon instead! And we will compare the performance
cd data

# Get the human transcriptome
#curl ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.transcripts.fa.gz \
#-o gencode.v29.transcripts.fa.gz

# Index the human transcriptome
salmon --threads=16 --no-version-check index \
-t Homo_sapiens.GRCh38.cdna.all.fa \
-i human_index \
-k 23
 
# Go to our trimmed read data
cd raw_data

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
