#!/usr/bin/env python

# Import necessary modules
import sys
import os
import re
import argparse
from encode_lib_common import (
    log, ls_l, mkdir_p, rm_f, run_shell_cmd, strip_ext_fastq, strip_ext_tar,
    untar)
from encode_lib_genomic import samtools_sort, bam_is_empty


# Argument parsing
def parse_arguments():
    parser = argparse.ArgumentParser(
        prog='ENCODE DCC bowtie2 aligner.',
        description='Align FASTQ files to a reference genome using Bowtie2 and convert to sorted BAM.'
    )

    # Positional arguments
    parser.add_argument('bowtie2_index_prefix_or_tar', type=str,
                        help='Prefix for Bowtie2 index or path to a .tar archive containing the index files.')

    parser.add_argument('fastqs', nargs='+', type=str,
                        help='Input FASTQ files (.gz compressed). Use 1 file for single-end, 2 for paired-end.')

    # Optional arguments
    parser.add_argument('--paired-end', action="store_true",
                        help='Set this flag if FASTQs are paired-end.')

    parser.add_argument('--local', action="store_true",
                        help='Use Bowtie2 local alignment mode (soft clipping). Default is end-to-end.')

    parser.add_argument('--multimapping', default=0, type=int,
                        help='Allow multimapping reads. Value passed to Bowtie2 -k parameter (value + 1).')

    parser.add_argument('--nth', type=int, default=1,
                        help='Number of threads to use.')

    parser.add_argument('--mem-gb', type=float,
                        help='Maximum memory in GB for samtools sort.')

    parser.add_argument('--out-dir', default='', type=str,
                        help='Directory where output files will be saved.')

    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'CRITICAL', 'ERROR'],
                        help='Set the log verbosity.')

    args = parser.parse_args()

    # Validate number of FASTQ files
    if args.paired_end and len(args.fastqs) != 2:
        raise argparse.ArgumentTypeError('Two FASTQs required for paired-end data.')
    if not args.paired_end and len(args.fastqs) != 1:
        raise argparse.ArgumentTypeError('One FASTQ required for single-end data.')

    # Set logging level
    log.setLevel(args.log_level)
    log.info(sys.argv)

    return args


# Alignment for single-end reads
def bowtie2_se(fastq, ref_index_prefix, multimapping, local, nth, mem_gb, out_dir):
    basename = os.path.basename(strip_ext_fastq(fastq))
    prefix = os.path.join(out_dir, basename)
    tmp_bam = f'{prefix}.bam'

    # Run bowtie2 and convert SAM to BAM with samtools
    run_shell_cmd(
        f'bowtie2 {"-k " + str(multimapping + 1) if multimapping else ""} '
        f'{"--local" if local else ""} --mm --threads {nth} -x {ref_index_prefix} '
        f'-U {fastq} | samtools view -bS - > {tmp_bam}'
    )

    # Sort BAM and return path
    bam = samtools_sort(tmp_bam, nth, mem_gb, out_dir)
    rm_f(tmp_bam)
    return bam


# Alignment for paired-end reads
def bowtie2_pe(fastq1, fastq2, ref_index_prefix, multimapping, local, nth, mem_gb, out_dir):
    basename = os.path.basename(strip_ext_fastq(fastq1))
    prefix = os.path.join(out_dir, basename)
    tmp_bam = f'{prefix}.bam'

    # Run bowtie2 with paired-end reads
    run_shell_cmd(
        f'bowtie2 {"-k " + str(multimapping + 1) if multimapping else ""} -X2000 '
        f'{"--local" if local else ""} --mm --threads {nth} -x {ref_index_prefix} '
        f'-1 {fastq1} -2 {fastq2} | samtools view -bS - > {tmp_bam}'
    )

    # Sort and return BAM
    bam = samtools_sort(tmp_bam, nth, mem_gb, out_dir)
    rm_f(tmp_bam)
    return bam

# Helper function: Check if Bowtie2 index exists
def chk_bowtie2_index(prefix):
    index_1 = f'{prefix}.1.bt2'
    index_2 = f'{prefix}.1.bt2l'
    if not (os.path.exists(index_1) or os.path.exists(index_2)):
        raise Exception(f"Bowtie2 index not found at prefix: {prefix}")

# Helper function: Find index prefix in directory
def find_bowtie2_index_prefix(d):
    if d == '':
        d = '.'
    for f in os.listdir(d):
        if f.endswith('.rev.1.bt2') or f.endswith('.rev.1.bt2l'):
            return re.sub(r'\.rev\.1\.(bt2|bt2l)$', '', f)
        elif f.endswith('.1.bt2') or f.endswith('.1.bt2l'):
            return re.sub(r'\.1\.(bt2|bt2l)$', '', f)
    return None

# Main workflow
def main():
    #
