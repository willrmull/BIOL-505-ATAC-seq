#!/usr/bin/env python

# Import necessary libraries for system operations, argument parsing, and the custom functions
import sys
import os
import argparse
from encode_lib_common import (
    mkdir_p, log, ls_l, rm_f, strip_ext_fastq)
from encode_lib_genomic import (
    get_read_length, remove_chrs_from_bam, samstat, samtools_index)


def parse_arguments():
    # Function to parse the command-line arguments
    parser = argparse.ArgumentParser(prog='ENCODE post align', description='')
    
    # Required input arguments: FASTQ file for R1 and BAM file
    parser.add_argument('fastq', type=str, help='Path for FASTQ R1')
    parser.add_argument('bam', type=str, help='Path for BAM')
    
    # Optional arguments
    parser.add_argument(
        '--chrsz', type=str,
        help='2-column chromosome sizes file. If not given, '
             'SAMstats on mito-free BAM will not be calculated.'
    )
    parser.add_argument('--mito-chr-name', default='chrM', help='Mitochondrial chromosome name.')
    parser.add_argument('--nth', type=int, default=1, help='Number of threads to parallelize.')
    parser.add_argument('--mem-gb', type=float,
                        help='Max. memory for samtools sort in GB. This should be total memory for the task, not per thread.')
    parser.add_argument('--out-dir', default='', type=str, help='Output directory.')
    parser.add_argument('--log-level', default='INFO', choices=['NOTSET', 'DEBUG', 'INFO',
                                 'WARNING', 'CRITICAL', 'ERROR', 'CRITICAL'],
                        help='Log level')

    args = parser.parse_args()  # Parse command-line arguments

    log.setLevel(args.log_level)  # Set the logging level
    log.info(sys.argv)  # Log the command-line invocation for tracking
    return args


def make_read_length_file(fastq, out_dir):
    '''
    Generate a text file containing the read length for the provided FASTQ file.
    '''
    basename = os.path.basename(strip_ext_fastq(fastq))  # Get the basename of the FASTQ file without extension
    prefix = os.path.join(out_dir, basename)  # Prefix for output files
    txt = '{}.read_length.txt'.format(prefix)  # Name of the output text file containing the read length
    
    # Get the read length using a utility function from the encode_lib_genomic module
    read_length = get_read_length(fastq)
    
    # Write the read length to a text file
    with open(txt, 'w') as fp:
        fp.write(str(read_length))
    
    return txt  # Return the file path


def main():
    # Main function to orchestrate the post-alignment tasks
    
    # Step 1: Read parameters from the command-line arguments
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)  # Create the output directory if it doesn't exist

    # Step 2: Generate the read length file for the FASTQ file
    log.info('Generating read length file...')
    make_read_length_file(args.fastq, args.out_dir)

    # Step 3: Run `samtools index` to index the BAM file
    log.info('Running samtools index...')
    samtools_index(args.bam, args.nth, args.out_dir)

    # Step 4: Run `SAMstat` to get statistics for the raw BAM file
    log.info('SAMstats on raw BAM...')
    samstat(args.bam, args.nth, args.mem_gb, args.out_dir)

    # Step 5: If chromosome sizes are provided, process the non-mitochondrial BAM
    if args.chrsz:
        log.info('SAMstats on non-mito BAM...')
        
        # Create a subdirectory to store the non-mitochondrial BAM results
        non_mito_out_dir = os.path.join(args.out_dir, 'non_mito')
        mkdir_p(non_mito_out_dir)

        # Remove mitochondrial chromosomes from the BAM file and generate the non-mito BAM
        non_mito_bam = remove_chrs_from_bam(args.bam, [args.mito_chr_name],
                                            args.chrsz, args.nth, non_mito_out_dir)

        # Run SAMstat on the non-mito BAM file
        samstat(non_mito_bam, args.nth, args.mem_gb, non_mito_out_dir)

        # Clean up by removing the non-mito BAM file after processing
        rm_f(non_mito_bam)

    # Step 6: List all files in the output directory to confirm the results
    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    # Final log message to indicate the process is complete
    log.info('All done.')


if __name__ == '__main__':
    main()  # Run the main function if the script is executed directly
