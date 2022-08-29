#!/bin/bash

#WORKFLOW FROM FASTQ TO VCF

#Step 1. Process the raw unmapped reads (fastq)
#1.1 Generate FastQC and MultiQC reports


######################################


######## analysis_dir --- full path to the folder where all analysis files will be kept. They will be divided into different subfolders for each step.
######## fastQ_files_folder --- full path to the folder where the raw FASTQ samples to be analysed are located.
######## softwares_fastqc --- full path to the fastqc folder
######## threads --- number of cores/cpus the process will use. If not provided, the default value is 4.


set -e

analysis_dir=/group/bioinf_biomarkers_rna/analysis_folder
fastQ_files_folder=/group/bioinf_biomarkers_rna/EGAD00001003448
softwares_fastqc=/group/bioinf_biomarkers_rna/Software_eka/FastQC
threads=4

                                    
##################################### STEP 1.1 Run FASTQC ##################################### 

mkdir $analysis_dir/0_fastQ_data/ -p

#move all fastq files from the original folder into the analysis folder (0_fastQ_data)
 
cd $fastQ_files_folder
find . -type f -name  \*.fastq.gz  -printf "/%P\n" | while read FILE ; do filename=$(basename $FILE) ; mv ."$FILE" "$analysis_dir/0_fastQ_data"/$filename; done 
  
###### Now all raw data are here $analysis_dir/0_fastQ_data/*.fastq.gz 

#make an array of all files from dataset
files=$(ls $analysis_dir/0_fastQ_data/*.fastq.gz)


mkdir $analysis_dir/1_fromFastQ_toBAM/ -p
mkdir $analysis_dir/1_fromFastQ_toBAM/1.1_fastqc_reports/ -p

#FASTQC help/docs https://raw.githubusercontent.com/s-andrews/FastQC/master/INSTALL.txt

#### Run FASTQC

echo "Running FastQC on" ${#files[@]} "samples"

$softwares_fastqc/fastqc \
    --outdir $analysis_dir/1_fromFastQ_toBAM/1.1_fastqc_reports \
    --extract \
    --threads $threads \
   $files

mkdir $analysis_dir/1_fromFastQ_toBAM/1.1_multiqc_reports/ -p

##  for making a summary of all reports:

#multiqc help/docs https://multiqc.info/docs/

echo "Running MultiQC"

multiqc  $analysis_dir/1_fromFastQ_toBAM/1.1_fastqc_reports -o  $analysis_dir/1_fromFastQ_toBAM/1.1_multiqc_reports   --interactive -f

echo "MultiQC report is in" $analysis_dir/1_fromFastQ_toBAM/1.1_multiqc_reports

