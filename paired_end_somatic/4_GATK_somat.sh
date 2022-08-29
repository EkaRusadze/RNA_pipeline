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

## 4. GATK best practice part 1
#4.1 Split read that contain Ns in their cigar string
#4.2 BQSR

######################################

######## analysis_dir --- full path to the folder where all analysis files will be kept. They will be divided into different subfolders for each step.
######## softwares_folder --- full path to the softwares (GATK, samtools) folder
######## threads --- number of cores/cpus the process will use. If not provided, the default value is 4.
######## reference_genome38 --- full path to the reference genome file.
######## known --- full path to the file that contains a set of known variants.
######## filesfolder --- full path to the folder where the output from the last step (step 3) is located.
######## ploidy --- ploidy of the organism. It's by default set to 2 (for human).
######## temp_folder --- full path to a specific temp folder (the default temp folder might not have enough memory for GATK to run without issues).


    #--known 
    #On 08.02.2022: I downloaded the --known-sites (DbSNP for GATK) from 
    #ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/GATK/All_20180418.vcf.gz

######################################
set -e
analysis_dir=/group/bioinf_biomarkers_rna/forMutect

softwares_GATK=/group/bioinf_biomarkers_rna/Software_eka/gatk-4.2.5.0/gatk-package-4.2.5.0-local.jar
softwares_samtools=/group/bioinf_biomarkers_rna/Software_eka/samtools-1.14

threads=27
reference_genome38=/group/bioinf_biomarkers_rna/hg38/hg38.fa
known=/group/bioinf_biomarkers_rna/hg38/All_20180418.vcf.gz
files=$(ls $analysis_dir/3_BAM_files/3.2_deduppedBamSortedFiles/*.bam)
ploidy=2
memory=-Xmx220g 
temp_folder=/ssd2/users/eka/gatk5


                                      #####################################
                                      #####################################
                                      ## 4.1 Splits reads that contain Ns ##
                                      #####################################
                                      #####################################

#Generate all needed files if not available

if [ ! -s $reference_genome38.fai ]; then
    echo "Generating .fai file for reference"
    cd $softwares_samtools 
    samtools faidx $reference_genome38
fi

filename=$(echo $reference_genome38 | rev | cut -d"." -f2-  | rev)

if [ ! -s $filename.dict ]; then
    echo "Generating .dict file for reference"
    java -jar $softwares_GATK CreateSequenceDictionary -R $reference_genome38
fi


if [ ! -s $known.tbi ]; then
    echo "Generating index file for known sites"
    java -jar $softwares_GATK IndexFeatureFile -I $known
fi

#Create folders 

mkdir $analysis_dir/4_processed_BAM/ -p
mkdir $analysis_dir/4_processed_BAM/4.1_SplitsCigar/ -p
mkdir $analysis_dir/4_processed_BAM/4.2_recal_data_table/ -p
mkdir $analysis_dir/4_processed_BAM/4.2_recal_reads/ -p



for i in $files; do
    sample_name=$(echo $i  | xargs -n 1 basename | cut -f 1 -d '.')

    echo 'Step 4.1. SplitNCigarReads:' $sample_name

    java $memory -jar $softwares_GATK  SplitNCigarReads \
    -R $reference_genome38 \
    -I  $analysis_dir/3_BAM_files/3.2_deduppedBamSortedFiles/$sample_name.bam \
    -O $analysis_dir/4_processed_BAM/4.1_SplitsCigar/$sample_name.bam \
    --tmp-dir $temp_folder





    ############################################  4.2 Recalibrate base quality scores (BQSR) #####################################
    #https://gatkforums.broadinstitute.org/gatk/discussion/2801/howto-recalibrate-base-quality-scores-run-bqsr
    #Recalibrate base quality scores in order to correct sequencing errors and other experimental artifacts.

    #4.2.a Analyze patterns of covariation in the sequence dataset
    echo "Step 4.2.a. BQSR:" $sample_name

    java $memory -jar $softwares_GATK BaseRecalibrator \
    -R $reference_genome38 \
    -I $analysis_dir/4_processed_BAM/4.1_SplitsCigar/$sample_name.bam \
    --known-sites $known  \
    -O $analysis_dir/4_processed_BAM/4.2_recal_data_table/$sample_name.recal_data_table \
    --tmp-dir $temp_folder

    #--known-sites $resources_folder/resources2/gatk/All_20170710_chrModif_new.vcf  \
    #On 08.02.2022: I downloaded the --known-sites (DbSNP for GATK) from 
    #ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/GATK/All_20180418.vcf.gz

    #4.2.b Apply the recalibls ration to your sequence data
    
    echo "Step 4.2.b. BQSR:" $sample_name

    java $memory -jar $softwares_GATK ApplyBQSR \
    -R $reference_genome38 \
    -I $analysis_dir/4_processed_BAM/4.1_SplitsCigar/$sample_name.bam \
    -bqsr $analysis_dir/4_processed_BAM/4.2_recal_data_table/$sample_name.recal_data_table \
    -O $analysis_dir/4_processed_BAM/4.2_recal_reads/$sample_name.bam \
    --tmp-dir $temp_folder

    # Expected Result
    #This creates a file called recal_reads.bam containing all the original reads, but now with exquisitely accurate base substitution, insertion and deletion quality scores. By default, the original quality scores are discarded in order to keep the file size down. However, you have the option to retain them by adding the flag –emit_original_quals to the PrintReads command, in which case the original qualities will also be written in the file, tagged OQ.


done




                ########################################################################################################################

                                 


