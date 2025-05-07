#! /bin/bash 

input_file_1=$1
input_file_2=$2

trimmomatic PE -phred33 ${input_file_1} ${input_file_2} \
trimmed_output_1.fastq.gz trimmed_output_1_unpaired.fastq.gz \
trimmed_output_2.fastq.gz trimmed_output_2_unpaired.fastq.gz \
ILLUMINACLIP:/kuhpc/sw/conda/latest/envs/trimmomatic/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

