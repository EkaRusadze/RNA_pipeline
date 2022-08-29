

###########################################
## Data Analysis
#########################################

## 1. Step FASTQC
        ##1.1 run FASTQC Analysis
        ##1.2 UMI extraction (UMI-tools)
        ##1.3 removing adapters (5 prime)
        ##1.4 run FASTQC  after removing adapters
        ##1.5 Manualy cut Poly A, Poly T tiles or different adapters # if it is needed

## 2. Mapping using STAR aligner  
        ##2.1 BUILD STAR genome
        ##2.2 Alignment
        ##2.3 second-pass STAR
        ##2.4 Alignment with 2-pass STAR

## 3. Add read groups, sort, mark duplicates, and create index
        ##3.1 Sort based on coordinate
        ##3.2 AddOrReplaceReadGroups
        ##3.3 MarkDuplicates

## 4. Gene expression analysis
        ##4.1 Make TxDbFromGFF
        ##4.2 Generate counts
        ##4.3 Sample-wise correlation analysis
        ##4.4 DEGs using DESeq


######################################

######## analysis_dir --- full path to the folder where you would like to keep all the analysed files (FASTQ, SAM, BAM). They will be divided into different subfolders for each step.
######## fastQ_files_folder --- full path to the folder where the FASTQ samples that you would like to analyse are located.
######## softwares_folder --- full path to the folder where you installed all the necessary software such as fastqc, STAR, etc.
######## threads --- number of cores/cpus the process will use. If not provided, the default value is 4.
######## reference_genome38 --- full path to your reference genome file.
######## reference_gtf --- full path to your reference gtf (General Transfer Format) file.
######## batchID --- number of the batch
######## batchfolder --- path to the folder where files for this batch are located
######## temp_folder --- STAR needs a lot of memory to run and in most cases the default temp folder is not enough. It is STRONGLY RECOMMENDED to provide a full path to a new temp folder as STAR doesn't give an error if there isn't enough memory but the results will be incorrect!
######## fastqfolder --- full path to the folder where the FASTQ files analyzed/processed in the last step are located.

set -e

analysis_dir=/group/bioinf_biomarkers_rna/testbatch
fastQ_files_folder=/group/bioinf_biomarkers_rna/testfullwf
softwares_folder=~/Software
threads=25

#genome downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
#gtf downloaded from http://genome.ucsc.edu/cgi-bin/hgTables?

reference_genome38=/group/bioinf_biomarkers_rna/hg38/hg38.fa
reference_gtf=/group/bioinf_biomarkers_rna/hg38/hg38.gtf

batchID=0
#fastqfolder=/group/bioinf_biomarkers_rna/startestdata
batchfolder=/group/bioinf_biomarkers_rna/batchdata/subdir$batchID
temp_folder=/ssd2/users/eka/batch$batchID
if [ -d "$temp_folder" ]; then
  rm -r $temp_folder
fi


                                      #####################################
                                      #####################################
                                      ### 2. Mapping using STAR aligner ###
                                      #####################################
                                      #####################################
############################################   Initialization #####################################




files=$(ls $batchfolder/*.fastq.gz)
echo ${files[@]}



############################################ 2.2 Alignment #####################################,
# Alignment jobs were executed as follows:
genomeDir=$analysis_dir/2_mapping_STAR/2.1_STAR_genomeDir

#STAR manual https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf


names=()

for file in $files; do
    echo $file  
    #newName=$(basename $file | cut -d'_' -f1) #must be a single character
    basefile=$(basename $file)
    searchstring="_R"
    newName=${basefile%$searchstring*}
    echo $newName
    names+=($newName)
    
done

echo ${#names[@]} ${names[@]}

uniqnames=($(echo "${names[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

echo ${#uniqnames[@]} ${uniqnames[@]}

for sample_name in ${uniqnames[@]}; do
    readone=${batchfolder}/${sample_name}_R1_001.fastq.gz
    readtwo=${batchfolder}/${sample_name}_R2_001.fastq.gz
    echo $readone
    echo $readtwo

    #sample_name=$(echo $i  | xargs -n 1 basename | cut -f 1 -d '.')
    #echo $sample_name

    #rm -rf $analysis_dir/2_mapping_STAR/2.2_Mapping/star_1pass/$sample_name/*
    #rm -rf  $analysis_dir/2_mapping_STAR/2.1_STAR_genomeDir_pass/$sample_name/*

    mkdir $analysis_dir/2_mapping_STAR/2.2_Mapping/star_1pass/ -p
    mkdir $analysis_dir/2_mapping_STAR/2.2_Mapping/star_1pass/batch$batchID -p
    mkdir $analysis_dir/2_mapping_STAR/2.2_Mapping/star_1pass/batch$batchID/$sample_name/ -p
    cd $analysis_dir/2_mapping_STAR/2.2_Mapping/star_1pass/batch$batchID/$sample_name/

    echo "starting step 2.2 for sample "$sample_name

    $softwares_folder/STAR-2.7.10a/source/STAR  \
         --genomeDir  $genomeDir \
         --readFilesIn $readone $readtwo \
         --readFilesCommand zcat \
         --runThreadN $threads \
         --outTmpDir $temp_folder

    echo "finished step 2.2 for sample "$sample_name

    ############################################ 2.3 second-pass STAR #####################################
    # For the 2-pass STAR, a new index is then created using splice junction information contained in the file SJ.out.tab from the first pass:

    mkdir $analysis_dir/2_mapping_STAR/2.1_STAR_genomeDir_pass/ -p
    mkdir $analysis_dir/2_mapping_STAR/2.1_STAR_genomeDir_pass/batch$batchID -p
    mkdir $analysis_dir/2_mapping_STAR/2.1_STAR_genomeDir_pass/batch$batchID/$sample_name/ -p

    echo "starting step 2.3 for sample "$sample_name

    $softwares_folder/STAR-2.7.10a/source/STAR  --runMode genomeGenerate \
         --genomeDir  $analysis_dir/2_mapping_STAR/2.1_STAR_genomeDir_pass/batch$batchID/$sample_name \
         --genomeFastaFiles  $reference_genome38 \
         --sjdbFileChrStartEnd  $analysis_dir/2_mapping_STAR/2.2_Mapping/star_1pass/batch$batchID/$sample_name/SJ.out.tab \
         --sjdbOverhang 75 \
         --limitSjdbInsertNsj 2000000 \
         --runThreadN $threads \
         --outTmpDir $temp_folder

    echo "finished step 2.3 for sample "$sample_name

    ############################################ 2.4 Alignment with 2-pass STAR #####################################
    ## The resulting index is then used to produce the final alignments as follows:
    mkdir  $analysis_dir/2_mapping_STAR/2.2_Mapping/star_2pass/ -p
    mkdir  $analysis_dir/2_mapping_STAR/2.2_Mapping/star_2pass/batch$batchID -p
    mkdir  $analysis_dir/2_mapping_STAR/2.2_Mapping/star_2pass/batch$batchID/$sample_name/ -p
    cd  $analysis_dir/2_mapping_STAR/2.2_Mapping/star_2pass/batch$batchID/$sample_name

    echo "starting step 2.4 for sample "$sample_name    

    $softwares_folder/STAR-2.7.10a/source/STAR --genomeDir $analysis_dir/2_mapping_STAR/2.1_STAR_genomeDir_pass/batch$batchID/$sample_name \
              --outSAMstrandField intronMotif \
              --readFilesIn $readone $readtwo \
              --readFilesCommand zcat \
              --runThreadN $threads \
              --outTmpDir $temp_folder

    echo "finished step 2.4 for sample "$sample_name

    #rm -rf $analysis_dir/2_mapping_STAR/2.2_Mapping/star_1pass/*
    #rm -rf  $analysis_dir/2_mapping_STAR/2.1_STAR_genomeDir_pass/*
done

                                      #####################################

                                      #####################################
                                      #####################################
                                      ##### End of the second Workflow #####
                                      #####################################
                                      #####################################
## the final ouput is stored in  $analysis_dir/2_mapping_STAR/2.2_Mapping/star_2pass/*


