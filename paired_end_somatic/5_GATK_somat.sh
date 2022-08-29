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

## 5. GATK best practice part 2
#5.1 Create Panel of Normals
#5.2 GenomicsDB and somatic panel of normals
#5.3 Call variants with Mutect2
#5.4 Estimate cross sample contamination
#5.5 Filter variant calls with FilterMutectCalls
#5.6 Genotype Refinement
#5.7 VariantAnnotator


######################################


######## analysis_dir --- full path to the folder where all analysis files will be kept. They will be divided into different subfolders for each step.
######## softwares_GATK --- full path to the GATK package jar file
######## softwares_bcftools --- full path to the bcftools folder
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



######################################
set -e
analysis_dir=/group/bioinf_biomarkers_rna/forMutect

softwares_GATK=/group/bioinf_biomarkers_rna/Software_eka/gatk-4.2.5.0/gatk-package-4.2.5.0-local.jar
resource_folder=/group/bioinf_biomarkers_rna/VCF_resources
softwares_bcftools=/group/bioinf_biomarkers_rna/Software_eka/bcftools
annotator_resources=/group/bioinf_biomarkers_rna/VCF_resources/funcotator_dataSources.v1.7.20200521s

threads=27
reference_genome38=/group/bioinf_biomarkers_rna/hg38/hg38.fa
ploidy=2
memory=-Xmx75g 
nullSites=0.2
jobName=mergeAll
interval_list=$resource_folder/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list
germline_resource_path=$resource_folder/somatic-hg38_af-only-gnomad.hg38.vcf.gz ##tbi needed too! #can be run without this resource
germline_resource=$(basename $germline_resource_path)
temp_folder=/ssd2/users/eka
                                 #####################################

if [ ! -s $germline_resource_path.tbi ]; then
    echo "Generating index file for germline resource"
    java -jar $softwares_GATK IndexFeatureFile -I $germline_resource_path
fi 

########### select the samples that availble in file_names only. These are the file that must be merged.


files=$(ls $analysis_dir/4_processed_BAM/4.2_recal_reads/*.bam)

modnames=()

for i in $files; do
    sample_name=$(echo $i  | xargs -n 1 basename | cut -f 1 -d '.')
    modname=${sample_name::-1}
    #echo $sample_name $modname
    modnames+=($modname)
done

uniqnames=($(echo "${modnames[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

echo ${#uniqnames[@]} ${uniqnames[@]}



################### Step 5.1. Create Panel of Normals (PON)
## https://gatk.broadinstitute.org/hc/en-us/articles/360035531132

mkdir $analysis_dir/5_VCFs/ -p
mkdir $analysis_dir/5_VCFs/5.1_PoN_VCF/ -p

for sample_name in ${uniqnames[@]}; do

    normal_sample=${sample_name}N
    #tumor_sample=${sample_name}T


    if [ -f $analysis_dir/4_processed_BAM/4.2_recal_reads/$normal_sample.bam ]; then
        echo 'Step 5.1. creating PoN VCF for sample:' $normal_sample
        java $memory -jar $softwares_GATK Mutect2 \
        -R $reference_genome38 \
        -I $analysis_dir/4_processed_BAM/4.2_recal_reads/$normal_sample.bam \
        --max-mnp-distance 0 \
        -O $analysis_dir/5_VCFs/5.1_PoN_VCF/$normal_sample.vcf.gz \
        --tmp-dir $temp_folder
    fi

done


#Here we need a code to make two column tab delimtee include SampleName     path

files=$(ls $analysis_dir/5_VCFs/5.1_PoN_VCF/*.vcf.gz)

mkdir $analysis_dir/5_VCFs/5.2_PoN/ -p

for i in ${files[@]}; do
    filename=$(echo $i | xargs -n 1 basename | cut -f 1 -d '.')
    printf '%s\t%s\n' $filename $i
    #echo $filename $i
done > $analysis_dir/5_VCFs/5.2_PoN/file_paths.map 





######5.2A GenomicsDBImport



java $memory -jar $softwares_GATK GenomicsDBImport \
-R $reference_genome38 \
--genomicsdb-workspace-path $analysis_dir/5_VCFs/5.2_PoN/pondb --sample-name-map $analysis_dir/5_VCFs/5.2_PoN/file_paths.map \
--intervals  $interval_list \
--tmp-dir $temp_folder

#######5.2B Create somatic PoN
java $memory -jar $softwares_GATK CreateSomaticPanelOfNormals \
-R $reference_genome38 \
-V gendb://$analysis_dir/5_VCFs/5.2_PoN/pondb \
-O $analysis_dir/5_VCFs/5.2_PoN/pon.vcf.gz \
--germline-resource $germline_resource_path \
--tmp-dir $temp_folder

#Outside the loop for step 5.4

mkdir $analysis_dir/5_VCFs/5.4_Calculate_Contamination/ -p

java $memory -jar $softwares_GATK SelectVariants \
-R $reference_genome38 \
--restrict-alleles-to BIALLELIC \
-V $germline_resource_path \
-O $analysis_dir/5_VCFs/5.4_Calculate_Contamination/biallelic_$germline_resource \
--tmp-dir $temp_folder


 
############################################ 5.3 Mutect2  #####################################
#########5.3 Call variants and generate a bamout file

# https://gatk.broadinstitute.org/hc/en-us/articles/360035889791#4
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531132

mkdir $analysis_dir/5_VCFs/5.3_Mutect2/ -p

for sample_name in ${uniqnames[@]}; do

    normal_sample=${sample_name}N
    tumor_sample=${sample_name}T

    java $memory -jar $softwares_GATK GetSampleName \
    -I $analysis_dir/4_processed_BAM/4.2_recal_reads/$normal_sample.bam \
    -O $analysis_dir/5_VCFs/5.3_Mutect2/${normal_sample}_name.txt

    normal_name=$(cat $analysis_dir/5_VCFs/5.3_Mutect2/${normal_sample}_name.txt)

    rm $analysis_dir/5_VCFs/5.3_Mutect2/${normal_sample}_name.txt

    java $memory -jar $softwares_GATK Mutect2 \
    -R $reference_genome38 \
    -I $analysis_dir/4_processed_BAM/4.2_recal_reads/$tumor_sample.bam \
    -I $analysis_dir/4_processed_BAM/4.2_recal_reads/$normal_sample.bam \
    --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
    -normal $normal_name \
    -pon $analysis_dir/5_VCFs/5.2_PoN/pon.vcf.gz \
    --germline-resource $germline_resource_path \
    -O $analysis_dir/5_VCFs/5.3_Mutect2/$sample_name.vcf.gz \
    --tmp-dir $temp_folder \
    -bamout $analysis_dir/5_VCFs/5.3_Mutect2/${sample_name}_tumornormal.bam 


    
 ##check if gene expression is better on original bam or bamout

###########Step 5.4. Estimate cross sample contamination

#Run GetPileupSummaries on tumor samples

    java $memory -jar $softwares_GATK GetPileupSummaries \
    -I $analysis_dir/4_processed_BAM/4.2_recal_reads/$tumor_sample.bam \
    -V $analysis_dir/5_VCFs/5.4_Calculate_Contamination/biallelic_$germline_resource \
    -L $analysis_dir/5_VCFs/5.4_Calculate_Contamination/biallelic_$germline_resource \
    --tmp-dir $temp_folder \
    --minimum-population-allele-frequency 0.01 \
    -O $analysis_dir/5_VCFs/5.4_Calculate_Contamination/${tumor_sample}_getpileupsummaries.table 

#Run GetPileupSummaries on normal samples

    java $memory -jar $softwares_GATK GetPileupSummaries \
    -I $analysis_dir/4_processed_BAM/4.2_recal_reads/$normal_sample.bam \
    -V $analysis_dir/5_VCFs/5.4_Calculate_Contamination/biallelic_$germline_resource \
    -L $analysis_dir/5_VCFs/5.4_Calculate_Contamination/biallelic_$germline_resource \
    --tmp-dir $temp_folder \
    -O $analysis_dir/5_VCFs/5.4_Calculate_Contamination/${normal_sample}_getpileupsummaries.table 


#Calculate contamination

    java $memory -jar $softwares_GATK CalculateContamination \
    -I $analysis_dir/5_VCFs/5.4_Calculate_Contamination/${tumor_sample}_getpileupsummaries.table \
    -matched $analysis_dir/5_VCFs/5.4_Calculate_Contamination/${normal_sample}_getpileupsummaries.table \
    --tmp-dir $temp_folder \
    -O $analysis_dir/5_VCFs/5.4_Calculate_Contamination/${tumor_sample}_contamination.table 




###########Step 5.5 Filter variant calls with FilterMutectCalls
    mkdir $analysis_dir/5_VCFs/5.5_FilteredCalls/ -p


    java $memory -jar $softwares_GATK FilterMutectCalls \
    -R $reference_genome38 \
    -V $analysis_dir/5_VCFs/5.3_Mutect2/$sample_name.vcf.gz \
    -stats $analysis_dir/5_VCFs/5.3_Mutect2/$sample_name.vcf.gz.stats \
    --contamination-table $analysis_dir/5_VCFs/5.4_Calculate_Contamination/${tumor_sample}_contamination.table \
    -O $analysis_dir/5_VCFs/5.5_FilteredCalls/$sample_name.vcf.gz \
    --tmp-dir $temp_folder

done


toMerge=$(ls $analysis_dir/5_VCFs/5.5_FilteredCalls/*.vcf.gz)

$softwares_bcftools/bcftools merge $toMerge -o $analysis_dir/5_VCFs/MergedMutectFiltered.vcf.gz ##can mutect2 output gvcf?



############################################ 5.6 Genotype Refinement  #####################################

mkdir $analysis_dir/5_VCFs/5.6_GenotypeRefinement/ -p
####### 
# We will keep only the sites that has missing variants below ($nullSites) and discard the rest, 
# then We will apply two scenarios
# 1- keep the missing GT as it is ./., then later in the python script we will replace it with 0/0 (or another technique).


# Step 1: Filter low quality genotypes
### https://gatk.broadinstitute.org/hc/en-us/articles/360036350452-VariantFiltration
java $memory -jar $softwares_GATK VariantFiltration \
    -R $reference_genome38 \
    -V $analysis_dir/5_VCFs/MergedMutectFiltered.vcf.gz \
    -G-filter "GQ < 20.0" -G-filter-name lowGQ  \
    -O $analysis_dir/5_VCFs/5.6_GenotypeRefinement/Gfiltered_$jobName.vcf


# Step 2: SelectVariants and apply the filters
###  https://gatk.broadinstitute.org/hc/en-us/articles/360036362532-SelectVariants
java $memory -jar $softwares_GATK SelectVariants \
    -R $reference_genome38 \
    --exclude-filtered \
    --exclude-non-variants \
    --set-filtered-gt-to-nocall  \
    -V $analysis_dir/5_VCFs/5.6_GenotypeRefinement/Gfiltered_$jobName.vcf \
    -O $analysis_dir/5_VCFs/5.6_GenotypeRefinement/SelectVariants1_$jobName.vcf


# It hence looks like --setFilteredGtToNocall and --maxNOCALLfraction have to be used sequentially, in two different runs of SelectVariants. Am I right?
# Yes, that is correct. They need to be run sequentially. I will make a ticket to add it into the docs.  Sheila
# Step 3: SelectVariants and apply the filteres
java $memory -jar $softwares_GATK SelectVariants \
    -R $reference_genome38 \
    --max-nocall-fraction $nullSites \
    -V $analysis_dir/5_VCFs/5.6_GenotypeRefinement/SelectVariants1_$jobName.vcf \
    -O $analysis_dir/5_VCFs/5.6_GenotypeRefinement/SelectVariants2_$jobName.$nullSites.vcf


############################################ 5.7_VariantAnnotator  #####################################

## https://gatk.broadinstitute.org/hc/en-us/articles/360035889931-Funcotator-Information-and-Tutorial

mkdir $analysis_dir/5_VCFs/5.7_VariantAnnotator/ -p

java $memory -jar $softwares_GATK Funcotator \
-R $reference_genome38 \
-V $analysis_dir/5_VCFs/5.6_GenotypeRefinement/SelectVariants2_$jobName.$nullSites.vcf \
--ref-version hg38 \
--data-sources-path $annotator_resources \
-O $analysis_dir/5_VCFs/5.7_VariantAnnotator/Funcotated_$jobName.vcf \
--output-file-format VCF



##########################################





