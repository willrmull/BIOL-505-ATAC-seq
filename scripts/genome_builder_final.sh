#!/bin/bash
#!!!!!!!!!!!!!!!!!NEEDS TO BE RUN IN AN ENVIRONMENT WITH SAMTOOLS AND BOWTIE2!!!!!!!!!!!!!!!!!!!!!!
set -e

# give user directions if they enter the wrong arugments
if [[ "$#" -lt 2 ]]; then
  echo "Usage: ./build_genome_data.sh [GENOME] [DEST_DIR]"
  echo "Example: ./build_genome_data.sh hg38 /your/genome/data/path/hg38"
  exit 2
fi

# assign first argument as the genome and assign the second variable as an absolute path to the output directory (even if entered as a relative path)
GENOME=$1
DEST_DIR=$(cd $(dirname $2) && pwd -P)/$(basename $2)
TSV=${DEST_DIR}/${GENOME}.tsv

#if the output directory doesn't exist, create it
#then navigate to the output directory
mkdir -p ${DEST_DIR}
cd ${DEST_DIR}

# Set reference FASTA URL and chromosome filtering regex
case "$GENOME" in
  hg38)
    REF_FA_URL="https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz"
    REGEX_BFILT_PEAK_CHR_NAME="chr[\dXY]+"
    ;;
  mm10)
    REF_FA_URL="https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/@@download/mm10_no_alt_analysis_set_ENCODE.fasta.gz"
    REGEX_BFILT_PEAK_CHR_NAME="chr[\dXY]+"
    ;;
  hg19)
    REF_FA_URL="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/referenceSequences/male.hg19.fa.gz"
    REGEX_BFILT_PEAK_CHR_NAME="chr[\dXY]+"
    ;;
#THIS IS OURS! It is a monarch butterfly genome from ncbi, the regex selects the regular numbered and Z chromosomes from the FASTA headers
  dp)
    REF_FA_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/135/715/GCF_018135715.1_MEX_DaPlex/GCF_018135715.1_MEX_DaPlex_genomic.fna.gz"
    REGEX_BFILT_PEAK_CHR_NAME="chromosome\s(\d+|Z).*"
    ;;
#if the genome isn't defined here, print an error 
  *)
    echo "Unsupported genome: $GENOME"
    exit 1
    ;;
esac

#Download the FASTA from the URL defined above
echo "=== Downloading reference FASTA for ${GENOME}..."
wget -c -O ref.fa.gz "$REF_FA_URL"
gunzip -f ref.fa.gz

echo "=== Indexing FASTA with samtools..."
samtools faidx ref.fa

#creates a file listing each chromosome
echo "=== Generating chrom sizes..."
cut -f1,2 ref.fa.fai > ${GENOME}.chrom.sizes

#create bowtie2 indexes
echo "=== Building Bowtie2 index..."
mkdir -p bowtie2_index
cd bowtie2_index
ln -s ../ref.fa .
bowtie2-build ref.fa ref.fa
cd ..

# creates a TSV file with the reference FASTA file, chromosome sizes file, regex filtering for chromosomes file, and indexed bowtie 2 file, alongside their respective paths
echo "=== Creating TSV file..."
cat > ${TSV} <<EOF
genome_name\t${GENOME}
ref_fa\t${DEST_DIR}/ref.fa
chrsz\t${DEST_DIR}/${GENOME}.chrom.sizes
regex_bfilt_peak_chr_name\t${REGEX_BFILT_PEAK_CHR_NAME}
bowtie2_index\t${DEST_DIR}/bowtie2_index
EOF

echo "=== Done! TSV file written to ${TSV}"

