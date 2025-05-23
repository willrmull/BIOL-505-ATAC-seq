#!/usr/bin/env python

# Import necessary libraries and modules
import sys  # Provides access to system-specific parameters and functions
import os  # For interacting with the operating system (file manipulation, paths, etc.)
import re  # Regular expression handling
import argparse  # For handling command-line arguments
from encode_lib_common import (
    get_num_lines, log, ls_l, mkdir_p, rm_f, run_shell_cmd, strip_ext_fastq,
    strip_ext_tar, untar, get_gnu_sort_param,
)
from encode_lib_genomic import (
    get_read_length, samtools_sort, bam_is_empty, get_samtools_res_param)

# Function to parse command-line arguments
def parse_arguments():
    # Create the argument parser
    parser = argparse.ArgumentParser(prog='ENCODE DCC bwa aligner.', description='')

    # Define the arguments the script will accept
    parser.add_argument('bwa_index_prefix_or_tar', type=str,
                        help='Path for prefix (or a tarball .tar) for reference bwa index. \
                            Prefix must be like [PREFIX].sa. \
                            TAR ball can have any [PREFIX] but it should not \
                            have a directory structure in it.')
    parser.add_argument('fastqs', nargs='+', type=str,
                        help='List of FASTQs (R1 and R2). \
                            FASTQs must be compressed with gzip (with .gz).')
    parser.add_argument(
        '--use-bwa-mem-for-pe', action="store_true",
        help='Use "bwa mem" for PAIRED-ENDED dataset with R1 FASTQ\'s read length >= --bwa-mem-read-len-limit. '
             'For shorter reads, bwa aln will be used. ')
    parser.add_argument(
        '--rescue-reads-for-bwa-mem', action="store_true",
        help='Use -P for "bwa mem" to rescue missing hits only (by using SW algorithm) '
             'but do not try to find hits that fit a proper pair.'
    )
    parser.add_argument('--bwa-mem-read-len-limit', type=int, default=70,
                        help='Read length limit for bwa mem (for paired-ended FASTQs only). '
                             'bwa aln will be used instead of bwa mem if R1 reads are shorter than this.')
    parser.add_argument('--paired-end', action="store_true",
                        help='Paired-end FASTQs.')
    parser.add_argument('--nth', type=int, default=1,
                        help='Number of threads to parallelize.')
    parser.add_argument('--mem-gb', type=float,
                        help='Max. memory for samtools sort and GNU sort -S '
                        '(half of this value will be used for GNU sort) in GB. '
                        'It should be total memory for this task (not memory per thread).')
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO',
                                 'WARNING', 'CRITICAL', 'ERROR',
                                 'CRITICAL'],
                        help='Log level')
    args = parser.parse_args()

    # Check if the fastqs have the correct dimension
    if args.paired_end and len(args.fastqs) != 2:
        raise argparse.ArgumentTypeError('Need 2 fastqs for paired end.')
    if not args.paired_end and len(args.fastqs) != 1:
        raise argparse.ArgumentTypeError('Need 1 fastq for single end.')

    # Validation for specific options
    if args.use_bwa_mem_for_pe and not args.paired_end:
        raise ValueError('--use-bwa-mem-for-pe is for paired-ended FASTQs only.')
    if not args.use_bwa_mem_for_pe and args.rescue_reads_for_bwa_mem:
        raise ValueError('--rescue-reads-for-bwa-mem is available only when --use-bwa-mem-for-pe is activated.')

    # Set up logging level
    log.setLevel(args.log_level)
    log.info(sys.argv)  # Log the command-line arguments used to run the script
    return args  # Return the parsed arguments

# Function for running bwa aln for single-end
def bwa_aln(fastq, ref_index_prefix, nth, out_dir):
    basename = os.path.basename(strip_ext_fastq(fastq))  # Get the base name without the extension
    prefix = os.path.join(out_dir, basename)  # Path to output directory
    sai = '{}.sai'.format(prefix)  # Output file for the SAI (sequential alignment index)

    # Build the BWA aln command
    cmd = 'bwa aln -q 5 -l 32 -k 2 -t {nth} {ref} {fastq} > {sai}'.format(
        nth=nth,
        ref=ref_index_prefix,
        fastq=fastq,
        sai=sai)
    run_shell_cmd(cmd)  # Run the command in the shell
    return sai  # Return the path to the SAI file

# Function for single-end alignment using bwa
def bwa_se(fastq, ref_index_prefix, nth, mem_gb, out_dir):
    basename = os.path.basename(strip_ext_fastq(fastq))  # Get the base name of the fastq file
    prefix = os.path.join(out_dir, basename)  # Path to output directory
    tmp_bam = '{}.bam'.format(prefix)  # Temporary BAM file

    # Run bwa aln to generate the SAI file
    sai = bwa_aln(fastq, ref_index_prefix, nth, out_dir)

    # Run bwa samse to generate a SAM file, then convert it to BAM format using samtools
    run_shell_cmd(
        'bwa samse {ref} {sai} {fastq} | '
        'samtools view -bS /dev/stdin {res_param} > {tmp_bam}'.format(
            ref=ref_index_prefix,
            sai=sai,
            fastq=fastq,
            res_param=get_samtools_res_param('view', nth=nth),
            tmp_bam=tmp_bam,
        )
    )
    rm_f(sai)  # Remove the SAI file

    # Sort the BAM file using samtools
    bam = samtools_sort(tmp_bam, nth, mem_gb)
    rm_f(tmp_bam)  # Remove the temporary BAM file

    return bam  # Return the path to the sorted BAM file

# Function for paired-end alignment using bwa
def bwa_pe(fastq1, fastq2, ref_index_prefix, nth, mem_gb, use_bwa_mem_for_pe,
           bwa_mem_read_len_limit, rescue_reads_for_bwa_mem, out_dir):
    basename = os.path.basename(strip_ext_fastq(fastq1))  # Get the base name of the first fastq file
    prefix = os.path.join(out_dir, basename)  # Path to output directory
    sam = '{}.sam'.format(prefix)  # SAM file
    badcigar = '{}.badReads'.format(prefix)  # File for bad CIGAR reads
    bam = '{}.bam'.format(prefix)  # Output BAM file

    temp_files = []  # Temporary files to delete later
    read_len = get_read_length(fastq1)  # Get the read length of the first fastq file

    log.info('Guessed read length of R1 FASTQ: {read_len}'.format(read_len=read_len))

    # Choose bwa mem or bwa aln based on read length
    if use_bwa_mem_for_pe and read_len >= bwa_mem_read_len_limit:
        log.info('Using bwa mem.')

        # Build and run the bwa mem command
        cmd = 'bwa mem -M {extra_param} -t {nth} {ref_index_prefix} {fastq1} {fastq2} | gzip -nc > {sam}'.format(
            extra_param='-P' if rescue_reads_for_bwa_mem else '',
            nth=nth,
            ref_index_prefix=ref_index_prefix,
            fastq1=fastq1,
            fastq2=fastq2,
            sam=sam,
        )
        temp_files.append(sam)
    else:
        log.info('Using bwa aln for each (R1 and R2) and then bwa sampe.')

        # Generate SAI files for both fastq1 and fastq2
        sai1 = bwa_aln(fastq1, ref_index_prefix, nth, out_dir)
        sai2 = bwa_aln(fastq2, ref_index_prefix, nth, out_dir)

        # Build and run the bwa sampe command
        cmd = 'bwa sampe {ref_index_prefix} {sai1} {sai2} {fastq1} {fastq2} | gzip -nc > {sam}'.format(
            ref_index_prefix=ref_index_prefix,
            sai1=sai1,
            sai2=sai2,
            fastq1=fastq1,
            fastq2=fastq2,
            sam=sam,
        )
        temp_files.extend([sai1, sai2, sam])

    run_shell_cmd(cmd)  # Run the bwa sampe or mem command

    # Check for bad CIGAR reads and remove them
    run_shell_cmd(
        'zcat -f {sam} | '
        'awk \'BEGIN {{FS="\\t" ; OFS="\\t"}} ! /^@/ && $6!="*" '
        '{{ cigar=$6; gsub("[0-9]+D","",cigar); '
        'n = split(cigar,vals,"[A-Z]"); s = 0; '
        'for (i=1;i<=n;i++) s=s+vals[i]; seqlen=length($10); '
        'if (s!=seqlen) print $1"\\t"; }}\' | '
        'sort {sort_param} | uniq > {badcigar}'.format(
            sam=sam,
            sort_param=get_gnu_sort_param(mem_gb * 1024 ** 3, ratio=0.5),
            badcigar=badcigar,
        )
    )

    # If bad CIGAR reads exist, remove them
    if get_num_lines(badcigar) > 0:
        run_shell_cmd(
            'zcat -f {sam} | grep -v -F -f {badcigar} | '
            'samtools view -Su /dev/stdin | samtools sort /dev/stdin -o {bam} -T {prefix} {res_param}'.format(
                sam=sam,
                badcigar=badcigar,
                bam=bam,
                prefix=prefix,
                res_param=get_samtools_res_param('sort', nth=nth, mem_gb=mem_gb),
            )
        )
    else:
        # If no bad CIGAR reads, sort directly
        run_shell_cmd(
            'samtools view -Su {sam} | samtools sort /dev/stdin -o {bam} -T {prefix} {res_param}'.format(
                sam=sam,
                bam=bam,
                prefix=prefix,
                res_param=get_samtools_res_param('sort', nth=nth, mem_gb=mem_gb),
            )
        )

    rm_f(temp_files)  # Remove all temporary files
    return bam  # Return the final BAM file

# Function to check if the BWA index exists
def chk_bwa_index(prefix):
    index_sa = '{}.sa'.format(prefix)
    if not os.path.exists(index_sa):
        raise Exception("BWA index does not exist. Prefix = {}".format(prefix))

# Function to find the BWA index prefix in a directory
def find_bwa_index_prefix(d):
    """
    Returns the prefix of BWA index, e.g., returns PREFIX if PREFIX.sa exists
    Args:
        d: directory to search for .sa file
    """
    if d == '':
        d = '.'
    for f in os.listdir(d):
        if f.endswith('.sa'):
            return re.sub('\.sa$', '', f)
    return None

# Main function to orchestrate the bwa alignment
def main():
    # Read parameters from command-line arguments
    args = parse_arguments()

    # Log the initialization process
    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # Declare temp arrays for files to delete later
    temp_files = []

    # If the bwa index is a tarball, unpack it
    if args.bwa_index_prefix_or_tar.endswith('.tar') or \
            args.bwa_index_prefix_or_tar.endswith('.tar.gz'):
        log.info('Unpacking bwa index tar...')
        tar = args.bwa_index_prefix_or_tar
        untar(tar, args.out_dir)  # Extract the tarball
        bwa_index_prefix = find_bwa_index_prefix(args.out_dir)  # Find the prefix in the output directory
        temp_files.append('{}*'.format(bwa_index_prefix))
    else:
        bwa_index_prefix = args.bwa_index_prefix_or_tar

    # Check if the BWA index exists
    chk_bwa_index(bwa_index_prefix)

    # Run bwa alignment
    log.info('Running bwa...')
    if args.paired_end:
        bam = bwa_pe(
            args.fastqs[0], args.fastqs[1],
            bwa_index_prefix, args.nth, args.mem_gb, args.use_bwa_mem_for_pe,
            args.bwa_mem_read_len_limit, args.rescue_reads_for_bwa_mem,
            args.out_dir)
    else:
        bam = bwa_se(
            args.fastqs[0],
            bwa_index_prefix, args.nth, args.mem_gb,
            args.out_dir)

    # Log removal of temporary files
    log.info('Removing temporary files...')
    rm_f(temp_files)

    # Check if BAM file is empty
    log.info('Checking if BAM file is empty...')
    if bam_is_empty(bam, args.nth):
        raise ValueError('BAM file is empty, no reads found.')

    # List all files in the output directory
    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')

# Ensure the script runs when executed directly
if __name__ == '__main__':
    main()
