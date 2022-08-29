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

## 3. Add read groups, sort, mark duplicates, and create index
#3.1 Sort based on coordinate
#3.2 AddOrReplaceReadGroups
#3.3 MarkDuplicates


######################################

######## analysis_dir --- full path to the folder where all analysis files will be kept. They will be divided into different subfolders for each step.
######## softwares_folder --- full path to the softwares (STAR, Picard) folder
######## threads --- number of cores/cpus the process will use. If not provided, the default value is 4.
######## reference_genome38 --- full path to the reference genome file.
######## reference_gtf --- full path to the reference gtf (General Transfer Format) file.
######## temp_folder --- STAR needs a lot of memory to run and in most cases the default temp folder is not enough. It is STRONGLY RECOMMENDED to provide a full path to a new temp folder as STAR doesn't give an error if there isn't enough memory but the results will be incorrect!
######## SAMfolder --- Full path to the folder where the output from second pass STAR is located (step 2.4). The correct path is already provided if no changes were done to the base code. 
######## memory --- allocate memory to java.

set -e

analysis_dir=/group/bioinf_biomarkers_rna/analysis_folder
softwares_folder=/group/bioinf_biomarkers_rna/Software_eka
threads=4
reference_genome38=/group/bioinf_biomarkers_rna/hg38/hg38.fa
reference_gtf=/group/bioinf_biomarkers_rna/hg38/hg38.gtf
temp_folder=/ssd2/users/eka/star
SAMfolder=$analysis_dir/2_mapping_STAR/2.2_Mapping/star_2pass
memory=-Xmx30g
#Picard requires an already existing temp folder
mkdir $temp_folder -p



folders=$(ls $SAMfolder)
                                      #####################################
                                      #####################################
                                      ### 3.  Mark for Duplicate ###
                                      #####################################
                                      #####################################
#. Add read groups, sort, mark duplicates, and create index
#The above step produces a SAM file, which is then put through the usual Picard processing steps: adding read group information, sorting, marking duplicates and indexing.

############################################ 3.1 Sort based on coordinate #####################################


#Picard docs https://broadinstitute.github.io/picard/command-line-overview.html#Overview

cd $SAMfolder
for sample_name in $folders; do
    # sort the sam files based on coordinate
    mkdir $analysis_dir/3_BAM_files/ -p
    mkdir $analysis_dir/3_BAM_files/3.1_bamSortedFiles/ -p
    mkdir $analysis_dir/3_BAM_files/3.2_deduppedBamSortedFiles/ -p
    
    echo "Started step 3.1 for sample "$sample_name

    java $memory -Djava.io.tmpdir=$temp_folder -jar $softwares_folder/picard-2.26.10/picard.jar SortSam I=$SAMfolder/$sample_name/Aligned.out.sam O=$analysis_dir/3_BAM_files/3.1_bamSortedFiles/$sample_name.sorted.sam SORT_ORDER=coordinate TMP_DIR=$temp_folder MAX_RECORDS_IN_RAM=5000000 

    echo "Finished step 3.1 for sample "$sample_name

    ############################################ 3.2 AddOrReplaceReadGroups #####################################

    echo "Started step 3.2 for sample "$sample_name

    # AddOrReplaceReadGroups
    java $memory -Djava.io.tmpdir=$temp_folder -jar $softwares_folder/picard-2.26.10/picard.jar AddOrReplaceReadGroups I=$analysis_dir/3_BAM_files/3.1_bamSortedFiles/$sample_name.sorted.sam  O= $analysis_dir/3_BAM_files/3.1_bamSortedFiles/$sample_name.addReplace.bam RGID=$sample_name RGPU=illumina RGPL=illumina RGLB=$sample_name RGSM=$sample_name SO=coordinate  VALIDATION_STRINGENCY=SILENT TMP_DIR=$temp_folder MAX_RECORDS_IN_RAM=5000000 

    echo "Finished step 3.2 for sample "$sample_name

    ############################################ 3.3 MarkDuplicates #####################################

    echo "Started step 3.3 for sample "$sample_name    
    
    java $memory -Djava.io.tmpdir=$temp_folder -jar  $softwares_folder/picard-2.26.10/picard.jar MarkDuplicates I=$analysis_dir/3_BAM_files/3.1_bamSortedFiles/$sample_name.addReplace.bam  O=$analysis_dir/3_BAM_files/3.2_deduppedBamSortedFiles/$sample_name.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT ASSUME_SORT_ORDER=coordinate M=$analysis_dir/3_BAM_files/3.2_deduppedBamSortedFiles/$sample_name.output.metrics TMP_DIR=$temp_folder MAX_RECORDS_IN_RAM=5000000 

    echo "Finished step 3.3 for sample "$sample_name

    rm $analysis_dir/3_BAM_files/3.1_bamSortedFiles/$sample_name.sorted.sam
    rm $analysis_dir/3_BAM_files/3.1_bamSortedFiles/$sample_name.addReplace.bam
    rm $analysis_dir/3_BAM_files/3.2_deduppedBamSortedFiles/$sample_name.output.metrics

done

                                      #####################################
                                      #####################################
                                      #####################################
                                      ##### End of the third Workflow #####
                                      #####################################
                                      #####################################
                                      ############### End  ###############
## the final ouput is stored in  $analysis_dir/3_BAM_files/3.2_deduppedBamSortedFiles/*.bam




