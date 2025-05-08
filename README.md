## BIOL-505-ATAC-seq
# ATAC seq pipeline based off code from ENCODES atac-seq-pipeline

# Scripts
#-The scripts folder contains .sh files used for installing the conda environent for the project and converting the genome file to a tsv

# src
#The src folder contains python files used in the workflow

## Steps
#-Take in fasta or fastq file and convert it to a tsv using final_genome_builder.sh
#-Configure the atac_input.json file to take in the desired sequence
#-Run the atac.wdl file using atac_input.json as the input

## Workflow 

#Read in file and convert it to a tsv
#Trim ATAC-seq data
#Align sequences and convert it to a BAM file
#Filter sequences to remove duplicates and extract important regions
#Convert filtered sequence to TAG-ALIGN 
#Identifying open regions of the chromatin with peak calling
#Visualize Results

# Usage:
#./final_genome_builder.sh
#java -D.config.file=cromwell.conf -jar cromwell-50.jar run atac.wdl --inputs atac_input.json
