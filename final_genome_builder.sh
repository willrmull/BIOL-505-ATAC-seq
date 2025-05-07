#!/bin/bash

# Exit immediately if a command exits with a non-zero status, if a variable is unset, or if a pipeline fails.
set -euo pipefail

# === User configuration ===
GENOME="danaus_plexippus"  # Name of the genome, used for naming outputs
REF_FA_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/135/715/GCF_018135715.1_MEX_DaPlex/GCF_018135715.1_MEX_DaPlex_genomic.fna.gz"  # URL to download the reference FASTA file
DEST_DIR="$(pwd)/genome_resources"  # Directory to store all generated files
REGEX_BFILT_PEAK_CHR_NAME="chr.*"  # Regular expression for chromosome names (used by downstream tools)
TSV="${DEST_DIR}/final_${GENOME}.tsv"  # Path to the final metadata TSV file

# === Create destination directory ===
mkdir -p "${DEST_DIR}"  # Make the destination directory if it doesn't exist

# === Download reference FASTA file ===
echo "=== Downloading reference FASTA file..."
curl -L ${REF_FA_URL} -o "${DEST_DIR}/ref.fa.gz"  # Download the gzipped reference genome FASTA file

# === Unzip FASTA file ===
echo "=== Unzipping reference FASTA file..."
gunzip "${DEST_DIR}/ref.fa.gz"  # Uncompress the FASTA file

# === Index FASTA ===
echo "=== Indexing FASTA with samtools..."
samtools faidx "${DEST_DIR}/ref.fa"  # Create FASTA index using samtools

# === Generate chromosome sizes ===
echo "=== Generating chrom.sizes..."
cut -f1,2 "${DEST_DIR}/ref.fa.fai" > "${DEST_DIR}/${GENOME}.chrom.sizes"  # Extract chromosome names and lengths

# === Build Bowtie2 index ===
echo "=== Building Bowtie2 index..."
mkdir -p "${DEST_DIR}/bowtie2_index"  # Create a directory for the Bowtie2 index
cd "${DEST_DIR}/bowtie2_index"  # Change to that directory
ln -sf ../ref.fa .  # Create a symbolic link to the FASTA file in the index directory
bowtie2-build ref.fa ref.fa  # Build the Bowtie2 index
cd ..  # Return to the previous directory

# === Create genome TSV ===
echo "=== Creating TSV file..."
cat > "${TSV}" <<EOF
genome_name	${GENOME}
ref_fa	${DEST_DIR}/ref.fa
chrsz	${DEST_DIR}/${GENOME}.chrom.sizes
regex_bfilt_peak_chr_name	${REGEX_BFILT_PEAK_CHR_NAME}
bowtie2_idx_tar	${DEST_DIR}/bowtie2_index
EOF

echo " Done! TSV file written to ${TSV}"  # Inform the user that the process is complete

