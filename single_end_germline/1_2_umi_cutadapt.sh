#!/bin/bash

#WORKFLOW FROM FASTQ TO VCF

#Step 1. Process the raw unmapped reads (fastq)
#1.1 Generate FastQC and MultiQC reports
#Optional step: if there are samples to be excluded, their names can be provided in a simple text file (each sample name on newline), and they will be moved to a different folder
#1.2 UMI extraction
#1.3 Trim the reads with cutadapt
#1.4 Generate FastQC and MultiQC reports 2



######################################


######## analysis_dir --- full path to the folder where all analysis files will be kept. They will be divided into different subfolders for each step.
######## fastQ_files_folder --- full path to the folder where the raw FASTQ samples to be analysed are located.
######## softwares_fastqc --- full path to the fastqc folder
######## threads --- number of cores/cpus the process will use. If not provided, the default value is 4.
######## exclude_samples_file --- full path to the text file of sample names that need to be excluded (each sample name on newline). If file is not provided, this step will be skipped.
######## UMI_pattern --- UMI pattern to be extracted. If pattern is not provided, this step will be skipped.
######## adapters_file --- full path to the text file of adapters that need to be trimmed (each adapter on newline). If file is not provided, this step will be skipped.

### Note: if both UMI extraction and adapter trimming steps are skipped, no new FastQC and MultiQC reports will be generated (as no changes have been applied to the FASTQ files). 



######################################
set -e

analysis_dir=/group/bioinf_biomarkers_rna/analysis_folder
fastQ_files_folder=/group/bioinf_biomarkers_rna/EGAD00001003448
softwares_fastqc=/group/bioinf_biomarkers_rna/Software_eka/FastQC
threads=4
exclude_samples_file=
UMI_pattern=""
adapters_file=


#check if exclude_samples file exists and save sample names in an array
count=0
if [ -s $exclude_samples_file ]; then
    excluded_samples=()
    while IFS= read -r line; do
        excluded_samples+=($line)
    done < $exclude_samples_file

    echo ${#excluded_samples[@]}" samples to exclude: " ${excluded_samples[@]}

    mkdir $analysis_dir/0_excluded_samples/ -p

#check if sample from fastq folder is in the excluded_samples array and move it to excluded_samples folder

    for sample in $analysis_dir/0_fastQ_data/*.fastq.gz; do
        samplename="$(basename "$sample")"
        if [[ " ${excluded_samples[*]} " == *" $samplename "* ]]; then
            mv $sample $analysis_dir/0_excluded_samples/$samplename
            echo "Moved "$samplename" to excluded_samples folder."
            (( count++ ))
        fi
    done
    
    echo "Excluded samples located in folder:" $analysis_dir/0_excluded_samples/

else
    echo "Couldn't find exclude_samples file."
fi

echo "Moved" $count "samples into excluded_samples folder."

#now we only have the necessary samples in 0_fastQ_data
files=$(ls $analysis_dir/0_fastQ_data/*.fastq.gz)

echo ${#files[@]} "samples used in further analysis"


##################################### 1.2  UMI extraction (UMI-tools) ##################################### 

#check if UMI_pattern file exists and use the pattern in UMI extraction
if [ ! -z $UMI_pattern ]; then

    echo "UMI pattern is "$UMI_pattern

    mkdir $analysis_dir/1_fromFastQ_toBAM/1.2_fastQ_data_UMI_Extracted/ -p

    #UMI tools help/docs https://umi-tools.readthedocs.io/en/latest/QUICK_START.html
    #UMI extrace usage/pattern https://umi-tools.readthedocs.io/en/latest/reference/extract.html

    ls -1 $files | xargs -n1 basename | parallel umi_tools extract --bc-pattern $UMI_pattern -I $analysis_dir/0_fastQ_data/{} -S $analysis_dir/1_fromFastQ_toBAM/1.2_fastQ_data_UMI_Extracted/{} 
    folder_nextstep=$analysis_dir/1_fromFastQ_toBAM/1.2_fastQ_data_UMI_Extracted

    echo "UMI extracted FASTQ files saved in" $folder_nextstep
else
    folder_nextstep=$analysis_dir/0_fastQ_data
    echo "No UMI pattern given. Files in the next step taken from folder" $folder_nextstep
fi

##################################### 1.3  Step removing adapters (5 prime) ##################################### 

### removes all sequences that have the adapter 

files=$(ls  $folder_nextstep/*.fastq.gz)

#check if adapters file exists and save adapters in an array
if [ -s $adapters_file ]; then
    adapters=()
    while IFS= read -r line; do
        adapters+=($line)
    done < $adapters_file
    echo "Adapters to trim: "${adapters[@]}

#use adapters for trimming    
    adapters_string=""
    
    for adapter in ${adapters[@]}; do
        adapters_string+="-b "$adapter" "
    done
    #echo $adapters_string

    mkdir $analysis_dir/1_fromFastQ_toBAM/1.3_fastqc_removeAdapters/ -p

    #cutadapt help/docs https://cutadapt.readthedocs.io/en/stable/ 
    #cutadapt quick guides https://cutadapt.readthedocs.io/en/stable/recipes.html

    ls -1 $files | xargs -n1 basename | parallel cutadapt --cores=0 $adapters_string -e 0.1 --discard-trimmed -o $analysis_dir/1_fromFastQ_toBAM/1.3_fastqc_removeAdapters/{}   $folder_nextstep/{}

    folder_nextstep=$analysis_dir/1_fromFastQ_toBAM/1.3_fastqc_removeAdapters
    
    echo "Trimmed FASTQ files saved in" $folder_nextstep

else
    echo "No adapters file found. Files in the next step taken from folder" $folder_nextstep
fi


                                     
################################ 1.4  Step run FASTQ  after removing adapters #####################################

## analysise all the fastq files after removing adapters 

#if changes have been done (i.e. UMIs extracted or adapters trimmed), run FastQC again

if [ "$folder_nextstep" != "$analysis_dir/0_fastQ_data" ]; then

    files=$(ls $folder_nextstep/*.fastq.gz)


    mkdir $analysis_dir/1_fromFastQ_toBAM/1.4_fastqc_reports_last/ -p

    $softwares_fastqc/fastqc \
        --outdir $analysis_dir/1_fromFastQ_toBAM/1.4_fastqc_reports_last \
        --extract \
        --threads $threads \
          $files


    mkdir $analysis_dir/1_fromFastQ_toBAM/1.4_multiqc_reports_last/ -p

    ##  for making a summary of all reports:
    multiqc $analysis_dir/1_fromFastQ_toBAM/1.4_fastqc_reports_last -o $analysis_dir/1_fromFastQ_toBAM/1.4_multiqc_reports_last  --interactive -f
    #multiqc $analysis_dir/1_fromFastQ_toBAM/1.3_fastqc_reports_afterRemoveAdapter --flat --export -f --pdf -o $analysis_dir/1_fromFastQ_toBAM/1.3_multiqc_reports_afterRemoveAdapter  --interactive -f

    echo "New MultiQC reports are located in" $analysis_dir/1_fromFastQ_toBAM/1.4_multiqc_reports_last

else
    echo "No changes done to the original FASTQ samples."
    echo "MultiQC reports are located in" $analysis_dir/1_fromFastQ_toBAM/1.1_multiqc_reports 
    echo "FASTQ samples for the next step are located in" $folder_nextstep

fi















