#! /bin/bash


module load conda
module load fastqc
module load java/23

input_file_1=$1
input_file_2=$2

mkdir -p report_directory
fastqc ${input_file_1} -o fastqc_output/ 
fastqc ${input_file_2} -o fastqc_output/ 
