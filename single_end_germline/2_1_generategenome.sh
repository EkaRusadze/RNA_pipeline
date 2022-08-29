#!/bin/bash

#WORKFLOW FROM FASTQ TO VCF

#Step 1. Process the raw unmapped reads (fastq)
#1.1 Generate FastQC and MultiQC reports
#Optional step: if there are samples to be excluded, their names can be provided in a simple text file (each sample name on newline), and they will be moved to a different folder
#1.2 UMI extraction
#1.3 Trim the reads with cutadapt
#1.4 Generate FastQC and MultiQC reports 2

## 2. Mapping using STAR aligner  
        ##2.1 BUILD STAR genome



######################################

######## analysis_dir --- full path to the folder where all analysis files will be kept. They will be divided into different subfolders for each step.
######## softwares_folder --- full path to the softwares (STAR, Picard) folder
######## threads --- number of cores/cpus the process will use. If not provided, the default value is 4.
######## reference_genome38 --- full path to the reference genome file.
######## reference_gtf --- full path to the reference gtf (General Transfer Format) file.
######## temp_folder --- STAR needs a lot of memory to run and in most cases the default temp folder is not enough. It is STRONGLY RECOMMENDED to provide a full path to a new temp folder as STAR doesn't give an error if there isn't enough memory but the results will be incorrect!

#Optional: suggested reference files if working on human data (reference genome hg38) 
#### genome - http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
#### gtf - http://genome.ucsc.edu/cgi-bin/hgTables?

set -e

analysis_dir=/group/bioinf_biomarkers_rna/analysis_folder
softwares_folder=/group/bioinf_biomarkers_rna/Software_eka
threads=4
reference_genome38=/group/bioinf_biomarkers_rna/hg38/hg38.fa
reference_gtf=/group/bioinf_biomarkers_rna/hg38/hg38.gtf
temp_folder=/ssd2/users/eka/star
## STAR will throw an error if temp folder already exists
if [ -d "$temp_folder" ]; then
  rm -r $temp_folder
fi



############################################ 2.1 BUILD STAR genome #####################################

mkdir $analysis_dir/2_mapping_STAR/ -p
mkdir $analysis_dir/2_mapping_STAR/2.1_STAR_genomeDir/ -p

genomeDir=$analysis_dir/2_mapping_STAR/2.1_STAR_genomeDir

echo "Starting step 2.1"

#STAR manual https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf

$softwares_folder/STAR-2.7.10a/source/STAR --runMode genomeGenerate \
    --genomeDir $genomeDir \
    --genomeFastaFiles  $reference_genome38 \
    --sjdbGTFfile $reference_gtf \
    --sjdbOverhang 75 \
    --runThreadN $threads \
    --outTmpDir $temp_folder

echo "Finished step 2.1"



