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
#4.3 Variant Calling

## 5. GATK best practice part 2
#5.1 Consolidate GVCFs for joint calling with GenotypeGVCFs
#5.2 VariantRecalibrators 
#5.3 Genotype Refinement
#5.4 VariantAnnotator


######################################

######## analysis_dir --- full path to the folder where all analysis files will be kept. They will be divided into different subfolders for each step.
######## softwares_folder --- full path to the softwares (GATK, samtools) folder
######## softwares_GATK --- full path to the GATK package jar file
######## threads --- number of cores/cpus the process will use. If not provided, the default value is 4.
######## reference_genome38 --- full path to the reference genome file.
######## filesfolder --- full path to the folder where the output from the last step (step 3) is located.
######## resource_folder --- full path to the folder where the resources for training variant recalibrator is located
######## annotator_resources --- full path to the folder the the annotator resource for Funcotator is located
######## interval_list --- full path to the file that lists one or more genomic intervals over which to operate (for GenomicsDBImport)
######## temp_folder --- full path to a specific temp folder (the default temp folder might not have enough memory for GATK to run without issues).
######## memory --- allocate memory to java.

### (note @me, the code will probably have to be reworked in variant recalibrator for later version re: resources)

### some help for resources: https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle

######################################
set -e

analysis_dir=/group/bioinf_biomarkers_rna/analysis_folder
softwares_folder=/group/bioinf_biomarkers_rna/Software_eka
softwares_GATK=/group/bioinf_biomarkers_rna/Software_eka/gatk-4.2.5.0/gatk-package-4.2.5.0-local.jar
threads=4
reference_genome38=/group/bioinf_biomarkers_rna/hg38/hg38.fa
filesfolder=$analysis_dir/4_processed_BAM/4.3_BP_RNA_files
resource_folder=/group/bioinf_biomarkers_rna/VCF_resources
annotator_resources=/group/bioinf_biomarkers_rna/VCF_resources/funcotator_dataSources.v1.7.20200521g
interval_list=/group/bioinf_biomarkers_rna/VCF_resources/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list
temp_folder=/ssd2/users/eka
memory=-Xmx30g 
                                 #####################################
 



files=$(ls $filesfolder/*.g.vcf)

jobName=AllSamples

mkdir $analysis_dir/5_VCFs/ -p
mkdir $analysis_dir/5_VCFs/5.1_GenotypeGVCFs/ -p


############################################ 5.1 Consolidate GVCFs for joint calling with GenotypeGVCFs  #####################################
#https://gatk.broadinstitute.org/hc/en-us/articles/360035889971--How-to-Consolidate-GVCFs-for-joint-calling-with-GenotypeGVCFs

#Here we need a code to make two columns tab delimited SampleName     path

for i in ${files[@]}; do
    filename=$(echo $i | xargs -n 1 basename | cut -f 1 -d '.')
    printf '%s\t%s\n' $filename $i
    #echo $filename $i
done > $analysis_dir/5_VCFs/5.1_GenotypeGVCFs/file_paths.map 


#5.1A GenomicsDBImport

java $memory -jar $softwares_GATK GenomicsDBImport \
-R $reference_genome38 \
--genomicsdb-workspace-path $analysis_dir/5_VCFs/5.1_GenotypeGVCFs/genomicsdb --sample-name-map $analysis_dir/5_VCFs/5.1_GenotypeGVCFs/file_paths.map \
-reader-threads 5 --intervals  $interval_list
--tmp-dir $temp_folder

#5.1B GenomicsDBImport
java $memory -jar $softwares_GATK GenotypeGVCFs \
-R $reference_genome38 \
-V gendb://$analysis_dir/5_VCFs/5.1_GenotypeGVCFs/genomicsdb \
-O $analysis_dir/5_VCFs/5.1_GenotypeGVCFs/2_Merged_$jobName.g.vcf
--tmp-dir $temp_folder

 
############################################ 5.2 VariantRecalibrators  #####################################


###### Regarding SNP and INDEL resources https://sites.google.com/a/broadinstitute.org/legacy-gatk-forum-discussions/frequently-asked-questions/1259-Which-training-sets-arguments-should-I-use-for-running-VQSR



mkdir $analysis_dir/5_VCFs/5.2_VariantRecalibrator/ -p

java $memory -jar $softwares_GATK VariantRecalibrator \
    -R $reference_genome38 \
    -V $analysis_dir/5_VCFs/5.1_GenotypeGVCFs/2_Merged_$jobName.g.vcf \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0  $resource_folder/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz \
    --resource:omni,known=false,training=true,truth=true,prior=12.0  $resource_folder/resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf.gz \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 $resource_folder/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $resource_folder/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
    -an QD \
	-an ReadPosRankSum \
	-an FS \
	-an SOR \
	-an DP  \
  	-mode SNP \
    -O $analysis_dir/5_VCFs/5.2_VariantRecalibrator/SNP_output.recal \
    --tranches-file $analysis_dir/5_VCFs/5.2_VariantRecalibrator/SNP_output.tranches \
    --rscript-file $analysis_dir/5_VCFs/5.2_VariantRecalibrator/SNP_output.plots.R
    --tmp-dir $temp_folder


# VariantRecalibrator INDEL
java $memory -jar $softwares_GATK VariantRecalibrator \
    -R $reference_genome38 \
    -V  $analysis_dir/5_VCFs/5.1_GenotypeGVCFs/2_Merged_$jobName.g.vcf \
    --resource:mills,known=false,training=true,truth=true,prior=12.0  $resource_folder/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $resource_folder/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
    -an QD \
    -an DP \
    -an FS \
    -an SOR \
    -an ReadPosRankSum \
    -mode INDEL  \
    -O $analysis_dir/5_VCFs/5.2_VariantRecalibrator/INDEL_output.recal \
    --tranches-file $analysis_dir/5_VCFs/5.2_VariantRecalibrator/INDEL_output.tranches \
    --rscript-file $analysis_dir/5_VCFs/5.2_VariantRecalibrator/INDEL_output.plots.R
    --tmp-dir $temp_folder




# ApplyVQSR SNPs
java $memory -jar $softwares_GATK ApplyVQSR \
    -R $reference_genome38 \
    -V $analysis_dir/5_VCFs/5.1_GenotypeGVCFs/2_Merged_$jobName.g.vcf \
    --truth-sensitivity-filter-level 99.0 \
    --recal-file $analysis_dir/5_VCFs/5.2_VariantRecalibrator/SNP_output.recal \
    --tranches-file $analysis_dir/5_VCFs/5.2_VariantRecalibrator/SNP_output.tranches \
    -mode SNP \
    -O  $analysis_dir/5_VCFs/5.2_VariantRecalibrator/recalibrated_$jobName.SNPs.vcf
    --tmp-dir $temp_folder





# ApplyVQSR INDEL
java $memory -jar $softwares_GATK ApplyVQSR \
    -R $reference_genome38 \
    -V  $analysis_dir/5_VCFs/5.2_VariantRecalibrator/recalibrated_$jobName.SNPs.vcf \
    --truth-sensitivity-filter-level 99.0 \
    --recal-file $analysis_dir/5_VCFs/5.2_VariantRecalibrator/INDEL_output.recal \
    --tranches-file $analysis_dir/5_VCFs/5.2_VariantRecalibrator/INDEL_output.tranches \
    -mode INDEL \
    -O  $analysis_dir/5_VCFs/5.2_VariantRecalibrator/recalibrated_$jobName.final.vcf
    --tmp-dir $temp_folder



############################################ 5.3 Genotype Refinement  #####################################

mkdir $analysis_dir/5_VCFs/5.3_GenotypeRefinement/ -p
nullSites=0.2

####### 

# We will keep only the sites that has missing variants below ($nullSites) and discard the rest


# Step 1: Filter low quality genotypes
### https://gatk.broadinstitute.org/hc/en-us/articles/360036350452-VariantFiltration

java $memory -jar $softwares_GATK VariantFiltration \
    -R $reference_genome38 \
    -V  $analysis_dir/5_VCFs/5.2_VariantRecalibrator/recalibrated_$jobName.final.vcf \
    -G-filter "GQ < 20.0" -G-filter-name lowGQ  \
    -O $analysis_dir/5_VCFs/5.3_GenotypeRefinement/Gfiltered_$jobName.vcf
    --tmp-dir $temp_folder


# Step 2: SelectVariants and apply the filters
###  https://gatk.broadinstitute.org/hc/en-us/articles/360036362532-SelectVariants

java $memory -jar $softwares_GATK SelectVariants \
    -R $reference_genome38 \
    --exclude-filtered \
    --exclude-non-variants \
    --set-filtered-gt-to-nocall  \
    -V $analysis_dir/5_VCFs/5.3_GenotypeRefinement/Gfiltered_$jobName.vcf \
    -O $analysis_dir/5_VCFs/5.3_GenotypeRefinement/SelectVariants1_$jobName.vcf
    --tmp-dir $temp_folder

# Step 3: SelectVariants and apply the filters
java $memory -jar $softwares_GATK SelectVariants \
    -R $reference_genome38 \
    --max-nocall-fraction $nullSites \
    -V $analysis_dir/5_VCFs/5.3_GenotypeRefinement/SelectVariants1_$jobName.vcf \
    -O $analysis_dir/5_VCFs/5.3_GenotypeRefinement/SelectVariants2_$jobName.$nullSites.vcf
    --tmp-dir $temp_folder

############################################ 5.4 VariantAnnotator  #####################################

## https://gatk.broadinstitute.org/hc/en-us/articles/360035889931-Funcotator-Information-and-Tutorial

mkdir $analysis_dir/5_VCFs/5.4_VariantAnnotator/ -p

java $memory -jar $softwares_GATK Funcotator \
    -R $reference_genome38 \
    -V $analysis_dir/5_VCFs/5.3_GenotypeRefinement/SelectVariants2_$jobName.$nullSites.vcf \
    --ref-version hg38 \
    --data-sources-path $annotator_resources \
    -O $analysis_dir/5_VCFs/5.4_VariantAnnotator/Funcotated_$jobName.$nullSites.vcf \
    --output-file-format VCF
    --tmp-dir $temp_folder






                                      #####################################
                                      #####################################
                                      #####################################
                                      ##### End of the third Workflow #####
                                      #####################################
                                      #####################################

                                      ############### End  ###############
## the final ouput is stored in $analysis_dir/5_VCFs/5.4_VariantAnnotator/Funcotated_$jobName.$nullSites.vcf


