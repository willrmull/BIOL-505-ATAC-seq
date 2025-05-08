#!/usr/bin/env python

# Importing necessary libraries and functions
import sys
import os
import argparse
from encode_lib_common import (
    log, ls_l, mkdir_p, read_tsv, run_shell_cmd,
    strip_ext_fastq)


# Function to parse command line arguments
def parse_arguments(debug=False):
    # Initialize argument parser for the script
    parser = argparse.ArgumentParser(prog='ENCODE DCC fastq merger.',
                                     description='Merge multiple fastq files.')
    # Argument to specify the fastq files or TSV file containing them
    parser.add_argument(
        'fastqs', nargs='+', type=str,
        help='TSV file path or list of FASTQs. '
             'FASTQs must be compressed with gzip (with .gz). '
             'Use TSV for multiple fastqs to be merged later. '
             'row=merge_id, col=end_id).')
    
    # Argument to specify if paired-end fastqs are used
    parser.add_argument('--paired-end', action="store_true",
                        help='Paired-end FASTQs.')
    
    # Argument to specify the number of threads for parallelization
    parser.add_argument('--nth', type=int, default=1,
                        help='Number of threads to parallelize.')
    
    # Argument to specify the output directory
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')
    
    # Argument to specify the log level for logging
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO',
                                 'WARNING', 'CRITICAL', 'ERROR',
                                 'CRITICAL'],
                        help='Log level')

    # Parse the arguments
    args = parser.parse_args()

    # Parse fastqs command line input: checking if it's a TSV or FASTQ file
    if args.fastqs[0].endswith('.gz') or args.fastqs[0].endswith('.fastq') or \
            args.fastqs[0].endswith('.fq'):  # If it's fastq
        args.fastqs = [[f] for f in args.fastqs]  # Convert to list of lists
    else:  # If it's TSV
        args.fastqs = read_tsv(args.fastqs[0])  # Read TSV into a list of lists

    # Validate the number of FASTQs based on whether they are paired-end or not
    for i, fastqs in enumerate(args.fastqs):
        if args.paired_end and len(fastqs) != 2:
            raise argparse.ArgumentTypeError(
                'Need 2 fastqs per replicate for paired end.')
        if not args.paired_end and len(fastqs) != 1:
            raise argparse.ArgumentTypeError(
                'Need 1 fastq per replicate for single end.')

    # Set the logging level for this script
    log.setLevel(args.log_level)
    log.info(sys.argv)  # Log the command line arguments
    return args  # Return the parsed arguments


# Function to merge fastq files
def merge_fastqs(fastqs, end, out_dir):
    """Merge fastq files into a single gzipped fastq file and save in out_dir/R1, out_dir/R2"""
    # Create the output directory for the corresponding read end (R1/R2)
    out_dir = os.path.join(out_dir, end)
    mkdir_p(out_dir)  # Create the directory if it does not exist
    # Define the prefix for the merged file using the base name of the first fastq file
    prefix = os.path.join(out_dir,
                          os.path.basename(strip_ext_fastq(fastqs[0])))

    # If there are multiple fastq files for the same end, merge them into one
    if len(fastqs) > 1:
        merged = '{}.merged.fastq.gz'.format(prefix)
    else:
        merged = '{}.fastq.gz'.format(prefix)

    # Use zcat to concatenate the fastq files, and gzip to compress the output
    cmd = 'zcat -f {} | gzip -nc > {}'.format(
        ' '.join(fastqs),  # Join the fastqs into a space-separated string
        merged)  # Output the merged fastq
    run_shell_cmd(cmd)  # Run the command in the shell
    return merged  # Return the path to the merged fastq


# Main function
def main():
    # Read command-line arguments
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)  # Ensure the output directory exists

    # Arrays to hold fastqs for R1 and R2 (for paired-end)
    fastqs_R1 = []
    fastqs_R2 = []
    
    # Populate the R1 and R2 lists
    for fastqs in args.fastqs:
        fastqs_R1.append(fastqs[0])  # First fastq is always R1
        if args.paired_end:  # If paired-end, second fastq is R2
            fastqs_R2.append(fastqs[1])

    log.info('Merging fastqs...')
    log.info('R1 to be merged: {}'.format(fastqs_R1))  # Log the R1 fastqs to be merged
    
    # Merge R1 fastqs and log the result
    merged_R1 = merge_fastqs(fastqs_R1, 'R1', args.out_dir)
    
    # If paired-end, merge R2 fastqs as well
    if args.paired_end:
        log.info('R2 to be merged: {}'.format(fastqs_R2))  # Log the R2 fastqs to be merged
        merged_R2 = merge_fastqs(fastqs_R2, 'R2', args.out_dir)

    # List all the files in the output directory
    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')  # Final log message


# Check if this script is being run as the main module and call main() if so
if __name__ == '__main__':
    main()
