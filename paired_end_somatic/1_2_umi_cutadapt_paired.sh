#WORKFLOW FROM FASTQ TO DEGS


#General steps:
#1. Process the raw unmapped reads (fastq)
#2. Align the processed reads to genome using STAR algorithm (paired or single) (BAM)
#3. Read quantification (raw counts)
#4. Normalization and analysis with DEseq2
#5. GO term enrichment analysis with GSEA
#6. Results visualization



#Step 1. Process the raw unmapped reads (fastq)
#1.1 Generate quality report 1 (fastqc)
#1.2 Check for Phred score and overrepresent (exclude samples if needed)
#1.3 UMI extraction
#1.4 Trim the reads with cut adapter
#1.5 Generate quality report 2 (fastqc)
#1.6 Exclude samples if needed again



######################################
set -e

analysis_dir=
fastQ_files_folder=
softwares_folder=
exclude_samples_file=
UMI_pattern=""
UMI_pattern2=""
adapters_file=
threads=


#check if exclude_samples file exists and save sample names in an array
count=0
if [ -s $exclude_samples_file ]; then
    excluded_samples=()
    while IFS= read -r line; do
    #    echo $line
        excluded_samples+=($line)
    done < $exclude_samples_file

    echo ${#excluded_samples[@]}" samples to exclude: " ${excluded_samples[@]}

    mkdir $analysis_dir/0_excluded_samples/ -p

#check if sample from fastq folder is in the excluded_samples array and move it to excluded_samples folder

    for sample in $analysis_dir/0_fastQ_data/*.fastq.gz; do
        samplename="$(basename "$sample")"
    #    echo $samplename
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
files=()

for file in $analysis_dir/0_fastQ_data/*.fastq.gz; do
    files+=($file)
done
echo ${#files[@]}" samples used in further analysis"


files=$(echo ${files[*]}) ## convert array to strings with spaces (like a list)


###make a list with unique names

names=()

for file in $files; do
    #echo $file  
    #newName=$(basename $file | cut -d'_' -f1) #must be a single character
    basefile=$(basename $file)
    searchstring="_R"
    newName=${basefile%$searchstring*}
    #echo $newName
    names+=($newName)
    
done

#echo ${#names[@]} ${names[@]}

uniqnames=($(echo "${names[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

echo ${#uniqnames[@]} ${uniqnames[@]}


##################################### 1.2  UMI extraction (UMI-tools) ##################################### R
### Extract UMI barcode from a read and add it to the read name, leaving any sample barcode in place


#check if UMI_pattern exists and use the pattern in UMI extraction

if [ ! -z $UMI_pattern ] && [ ! -z $UMI_pattern2 ]; then

    echo "UMI pattern is "$UMI_pattern
    echo "UMI pattern2 is "$UMI_pattern2 

    mkdir $analysis_dir/1_fromFastQ_toBAM/1.2_fastQ_data_UMI_Extracted/ -p

    parallel umi_tools extract --bc-pattern $UMI_pattern --bc-pattern2 $UMI_pattern2 --read2-in $analysis_dir/0_fastQ_data/{}_R2_001.fastq.gz --read2-out $analysis_dir/1_fromFastQ_toBAM/1.2_fastQ_data_UMI_Extracted/{}_R2_001.fastq.gz -I $analysis_dir/0_fastQ_data/{}_R1_001.fastq.gz -S $analysis_dir/1_fromFastQ_toBAM/1.2_fastQ_data_UMI_Extracted/{}_R1_001.fastq.gz ::: ${uniqnames[@]}
    
    folder_nextstep=$analysis_dir/1_fromFastQ_toBAM/1.2_fastQ_data_UMI_Extracted

    echo "UMI extracted FASTQ files saved in" $folder_nextstep
else
    folder_nextstep=$analysis_dir/0_fastQ_data
    echo "At least one UMI pattern is not given! Files in the next step taken from folder" $folder_nextstep
fi

##################################### 1.3  Step removing adapters (5 prime) ##################################### R
### removes all sequences that have an adapter at stating at '5, with error of 10%
##########folder_nextstep=$analysis_dir/0_fastQ_data
files=$(ls  $folder_nextstep/*.fastq.gz)

names=()

for file in $files; do
    #echo $file  
    #newName=$(basename $file | cut -d'_' -f1) #must be a single character
    basefile=$(basename $file)
    searchstring="_R"
    newName=${basefile%$searchstring*}
    #echo $newName
    names+=($newName)
    
done

#echo ${#names[@]} ${names[@]}

uniqnames=($(echo "${names[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

echo ${#uniqnames[@]} ${uniqnames[@]}

#files=$(echo ls $in_dir/*.fastq.gz | xargs -n 1 basename)
#for i in  $files
#do
#cutadapt -g ^ACACTCTTTCCCTACACGACGCTCTTCCGATCT  -g ^GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -e 0.1 --cores=$threads --discard-trimmed -o $out_dir/fastQC_analysis/1.2_fastqc_removeAdapters/$i  $in_dir/$i
#done

#check if adapters file exists and save adapters in an array
if [ -s $adapters_file ]; then
    adapters=()
    while IFS= read -r line; do
        #echo $line
        adapters+=($line)
    done < $adapters_file
    echo "Adapters to trim: "${adapters[@]}

#use adapters for trimming    
    adapters_string=""
    
    for adapter in ${adapters[@]}; do
        adapters_string+="-b "$adapter" -B "$adapter" "
    done
    #echo $adapters_string

    mkdir $analysis_dir/1_fromFastQ_toBAM/1.3_fastqc_removeAdapters/ -p

    #ls -1 $files | xargs -n1 basename | parallel cutadapt --cores=0 -g ^ACACTCTTTCCCTACACGACGCTCTTCCGATCT  -g ^GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -e 0.1 --discard-trimmed -o $analysis_dir/1_fromFastQ_toBAM/1.3_fastqc_removeAdapters/{}   $folder_nextstep{}

    #cutadapt help/docs https://cutadapt.readthedocs.io/en/stable/ 
    #cutadapt quick guides https://cutadapt.readthedocs.io/en/stable/recipes.html
    echo "Starting cutadapt"

    parallel -P $threads cutadapt --cores=0 $adapters_string -e 0.1 --discard-trimmed -o $analysis_dir/1_fromFastQ_toBAM/1.3_fastqc_removeAdapters/{}_R1_001.fastq.gz -p $analysis_dir/1_fromFastQ_toBAM/1.3_fastqc_removeAdapters/{}_R2_001.fastq.gz $folder_nextstep/{}_R1_001.fastq.gz $folder_nextstep/{}_R2_001.fastq.gz ::: ${uniqnames[@]}


    folder_nextstep=$analysis_dir/1_fromFastQ_toBAM/1.3_fastqc_removeAdapters
    
    echo "Trimmed FASTQ files saved in" $folder_nextstep

else
    echo "No adapters file found. Files in the next step taken from folder" $folder_nextstep
fi


                                      #####################################
################################ 1.4  Step run FASTQ  after removing adapters #####################################
## analysise all the fastq files after removing adapters (5 prime):

#if changes have been done (i.e. UMIs extracted or adapters trimmed), run fastqc again

if [ "$folder_nextstep" != "$analysis_dir/0_fastQ_data" ]; then

    files=$(echo ls $folder_nextstep/*.fastq.gz)


    mkdir $analysis_dir/1_fromFastQ_toBAM/1.4_fastqc_reports_last/ -p

    $softwares_folder/FastQC/fastqc \
        --outdir $analysis_dir/1_fromFastQ_toBAM/1.4_fastqc_reports_last \
        --extract \
        --threads $threads \
          $files


    mkdir $analysis_dir/1_fromFastQ_toBAM/1.4_multiqc_reports_last/ -p

    ##  for making a summary of all reports:
    multiqc $analysis_dir/1_fromFastQ_toBAM/1.4_fastqc_reports_last -o $analysis_dir/1_fromFastQ_toBAM/1.4_multiqc_reports_last  --interactive -f
    #multiqc $analysis_dir/1_fromFastQ_toBAM/1.3_fastqc_reports_afterRemoveAdapter --flat --export -f --pdf -o $analysis_dir/1_fromFastQ_toBAM/1.3_multiqc_reports_afterRemoveAdapter  --interactive -f

    echo "MultiQC reports are located in" $analysis_dir/1_fromFastQ_toBAM/1.4_multiqc_reports_last

else
    echo "No changes done to the original FASTQ samples."
    echo "MultiQC reports are located in" $analysis_dir/1_fromFastQ_toBAM/1.1_multiqc_reports 
    echo "FASTQ samples for the next step are located in" $folder_nextstep

fi


                  












