#!/bin/bash

#WORKFLOW FROM FASTQ TO VCF

#Step 1. Process the raw unmapped reads (fastq)
#1.1 Generate FastQC and MultiQC reports
#Optional step: if there are samples to be excluded, their names can be provided in a simple text file (each sample name on newline), and they will be moved to a different folder
#1.2 UMI extraction
#1.3 Trim the reads with cutadapt
#1.4 Generate FastQC and MultiQC reports 2

## 2. Mapping using STAR aligner  
#2.1 BUILD STAR genome
#2.2 Alignment
#2.3 second-pass STAR
#2.4 Alignment with 2-pass STAR


######################################

######## analysis_dir --- full path to the folder where all analysis files will be kept. They will be divided into different subfolders for each step.
######## softwares_folder --- full path to the softwares (STAR, Picard) folder
######## threads --- number of cores/cpus the process will use. If not provided, the default value is 4.
######## reference_genome38 --- full path to the reference genome file.
######## reference_gtf --- full path to the reference gtf (General Transfer Format) file.
######## temp_folder --- STAR needs a lot of memory to run and in most cases the default temp folder is not enough. It is STRONGLY RECOMMENDED to provide a full path to a new temp folder as STAR doesn't give an error if there isn't enough memory but the results will be incorrect!
######## fastqfolder --- full path to the folder where the FASTQ files analyzed/processed in the last step are located (this path is printed out when running file 1.2).
######## overhang --- From the STAR manual: "Overhang specifies the length of the genomic sequence around the annotated junction to be used in constructing the splice junctions database. Ideally, this length should be equal to the ReadLength-1, where ReadLength is the length of the reads. For instance, for Illumina 2x100b paired-end reads, the ideal value is 100-1=99. In case of reads of varying length, the ideal value is max(ReadLength)-1. In most cases, the default value of 100 will work as well as the ideal value." If not provided, overhang will be set to 100.


set -e

analysis_dir=/group/bioinf_biomarkers_rna/analysis_folder
softwares_folder=/group/bioinf_biomarkers_rna/Software_eka
threads=4
reference_genome38=/group/bioinf_biomarkers_rna/hg38/hg38.fa
reference_gtf=/group/bioinf_biomarkers_rna/hg38/hg38.gtf
temp_folder=/ssd2/users/eka/star
fastqfolder=$analysis_dir/1_fromFastQ_toBAM/1.3_fastqc_removeAdapters
overhang=100
## STAR will throw an error if temp folder already exists
if [ -d "$temp_folder" ]; then
  rm -r $temp_folder
fi





files=$(ls $fastqfolder/*.fastq.gz)

############################################ 2.2 Alignment #####################################,

genomeDir=$analysis_dir/2_mapping_STAR/2.1_STAR_genomeDir

#STAR manual https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf



for file in $files; do
    sample_name=$(echo $file  | xargs -n 1 basename | cut -f 1 -d '.')


    mkdir $analysis_dir/2_mapping_STAR/2.2_Mapping/star_1pass/ -p
    mkdir $analysis_dir/2_mapping_STAR/2.2_Mapping/star_1pass/$sample_name/ -p
    cd $analysis_dir/2_mapping_STAR/2.2_Mapping/star_1pass/$sample_name/

    echo "Starting step 2.2 for sample "$sample_name

    $softwares_folder/STAR-2.7.10a/source/STAR  \
         --genomeDir  $genomeDir \
         --readFilesIn $file \
         --readFilesCommand zcat \
         --runThreadN $threads \
         --outTmpDir $temp_folder

    echo "Finished step 2.2 for sample "$sample_name

    ############################################ 2.3 second-pass STAR #####################################
    # For the 2-pass STAR, a new index is then created using splice junction information contained in the file SJ.out.tab from the first pass:

    mkdir $analysis_dir/2_mapping_STAR/2.1_STAR_genomeDir_pass/ -p
    mkdir $analysis_dir/2_mapping_STAR/2.1_STAR_genomeDir_pass/$sample_name/ -p

    echo "Starting step 2.3 for sample "$sample_name

    $softwares_folder/STAR-2.7.10a/source/STAR  --runMode genomeGenerate \
         --genomeDir  $analysis_dir/2_mapping_STAR/2.1_STAR_genomeDir_pass/$sample_name \
         --genomeFastaFiles  $reference_genome38 \
         --sjdbFileChrStartEnd  $analysis_dir/2_mapping_STAR/2.2_Mapping/star_1pass/$sample_name/SJ.out.tab \
         --sjdbOverhang $overhang \
         --runThreadN $threads \
         --outTmpDir $temp_folder

    echo "Finished step 2.3 for sample "$sample_name

    ############################################ 2.4 Alignment with 2-pass STAR #####################################
    ## The resulting index is then used to produce the final alignments as follows:
    mkdir  $analysis_dir/2_mapping_STAR/2.2_Mapping/star_2pass/ -p
    mkdir  $analysis_dir/2_mapping_STAR/2.2_Mapping/star_2pass/$sample_name/ -p
    cd  $analysis_dir/2_mapping_STAR/2.2_Mapping/star_2pass/$sample_name

    echo "Starting step 2.4 for sample "$sample_name    

    $softwares_folder/STAR-2.7.10a/source/STAR \
              --genomeDir $analysis_dir/2_mapping_STAR/2.1_STAR_genomeDir_pass/$sample_name \
              --readFilesIn $file \
              --readFilesCommand zcat \
              --runThreadN $threads \
              --outTmpDir $temp_folder

    echo "Finished step 2.4 for sample "$sample_name

    rm -rf $analysis_dir/2_mapping_STAR/2.2_Mapping/star_1pass/$sample_name
    rm -rf  $analysis_dir/2_mapping_STAR/2.1_STAR_genomeDir_pass/$sample_name
done

      

