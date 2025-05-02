#!/bin/bash
set -e

#give user instructions for how to input arguments if they don't input the correct number
if [[ "$#" -lt 2 ]]; then
  echo "Usage: ./build_genome_data.sh [GENOME] [DEST_DIR]"
  echo "Example: ./build_genome_data.sh hg38 /your/genome/data/path/hg38"
  exit 2
fi

# assign the first argument as the variable GENOME (is the code referenced later to specify which genome to use)
#assign the second argument as the output directory (makees it an absolute path regardless of if path entered is relative or absolute)
GENOME=$1
DEST_DIR=$(cd $(dirname $2) && pwd -P)/$(basename $2)
TSV=${DEST_DIR}/${GENOME}.tsv

#makes the output directory if it doesn't exist and navigates to it
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
#THIS IS THE ONE WE USE! It is a monarch butterfly genome from ncbi.  The regex selects the header line from only regular chromosomes (not fragments)
  dp)
    REF_FA_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/135/715/GCF_018135715.1_MEX_DaPlex/GCF_018135715.1_MEX_DaPlex_genomic.fna.gz"
    REGEX_BFILT_PEAK_CHR_NAME="^>NC_[^\s]+.*?(chromosome \d+|chromosome [A-Z])"
    ;;
  *)
    echo "Unsupported genome: $GENOME"
    exit 1
    ;;
esac

#Get the .gz genome from the URL specified above (NCBI in our case) and decompress
echo "=== Downloading reference FASTA for ${GENOME}..."
wget -c -O ref.fa.gz "$REF_FA_URL"
gunzip -f ref.fa.gz

#uses samtools to reference the reference genome.  It outputs a file titled ref.fa.fai containing:
##### sequence name (of each chromosome)
##### length of sequence in bp
##### the offset of the sequence (aka the position in the file (in bytes) that the sequence starts at, like a location marker)
##### number of base pairs per line
##### number of characters per line
echo "=== Indexing FASTA with samtools..."
samtools faidx ref.fa

#create a file with just sequence name and number of base pairs from ref.fa.fai
echo "=== Generating chrom sizes..."
cut -f1,2 ref.fa.fai > ${GENOME}.chrom.sizes

#Uses Bowtie2 to make an index of the reference genome
#indexes are pre-processed, compressed versions of the reference genome that makes alignment go faster
echo "=== Building Bowtie2 index..."
mkdir -p bowtie2_index
cd bowtie2_index
ln -s ../ref.fa .
bowtie2-build ref.fa ref.fa
cd ..

#creates a TSV file with the genome name, reference FASTA, chromosome size, regex pattern, and bowtie2 indexes, alongside their respective file paths
echo "=== Creating TSV file..."
cat > ${TSV} <<EOF
genome_name\t${GENOME}
ref_fa\t${DEST_DIR}/ref.fa
chrsz\t${DEST_DIR}/${GENOME}.chrom.sizes
regex_bfilt_peak_chr_name\t${REGEX_BFILT_PEAK_CHR_NAME}
bowtie2_index\t${DEST_DIR}/bowtie2_index
EOF

echo "=== Done! TSV file written to ${TSV}"
