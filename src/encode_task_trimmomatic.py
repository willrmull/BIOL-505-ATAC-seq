#!/usr/bin/env python

# Import necessary modules
import sys
import os
import argparse
from encode_lib_common import (
    assert_file_not_empty, log, ls_l, mkdir_p, rm_f,
    run_shell_cmd, strip_ext_fastq)  # Functions for file handling, logging, etc.
from encode_lib_genomic import (
    locate_trimmomatic)  # Function to locate the Trimmomatic executable


# Function to parse command line arguments
def parse_arguments(debug=False):
    # Initialize argument parser
    parser = argparse.ArgumentParser(
        prog='ENCODE DCC Trimmomatic wrapper.')
    
    # Add expected command-line arguments with their descriptions
    parser.add_argument('--fastq1',
                        help='FASTQ R1 to be trimmed.')  # FASTQ file for read 1
    parser.add_argument('--fastq2',
                        help='FASTQ R2 to be trimmed.')  # FASTQ file for read 2 (for paired-end)
    parser.add_argument('--paired-end', action="store_true",
                        help='Paired-end FASTQs.')  # Flag for paired-end reads
    parser.add_argument('--crop-length', type=int, required=True,
                        help='Number of basepair to crop.'
                             'Trimmomatic\'s parameter CROP.')  # Length to crop reads
    parser.add_argument('--crop-length-tol', type=int, default=2,
                        help='Crop length tolerance to keep shorter reads '
                             'around the crop length. '
                             'Trimmomatic\'s parameter MINLEN will be --crop-length '
                             '- abs(--crop-length-tol).')  # Tolerance to avoid removing short reads
    parser.add_argument('--phred-score-format',
                        default='auto',
                        choices=['auto', 'phred33', 'phred64'],
                        help='Base encoding for Phred scores in FASTQs. '
                             'If it is not auto then -phred33 or -phred64 to '
                             'Trimmomatic\'s command line.')  # Phred score format (auto or specific)
    parser.add_argument('--out-dir-R1', default='', type=str,
                        help='Output directory for cropped R1 fastq.')  # Output directory for R1
    parser.add_argument('--out-dir-R2', default='', type=str,
                        help='Output directory for cropped R2 fastq.')  # Output directory for R2
    parser.add_argument('--trimmomatic-java-heap',
                        help='Trimmomatic\'s Java max. heap: java -jar Trimmomatic.jar '
                             '-Xmx[MAX_HEAP]')  # Optional Java heap size for Trimmomatic
    parser.add_argument('--nth', type=int, default=1,
                        help='Number of threads to parallelize.')  # Number of threads for parallel processing
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO',
                                 'WARNING', 'CRITICAL', 'ERROR',
                                 'CRITICAL'],
                        help='Log level')  # Logging level (default INFO)
    
    # Parse the arguments and return them
    args = parser.parse_args()

    # Validate crop length
    if not args.crop_length:
        raise ValueError('Crop length must be > 0.')

    log.setLevel(args.log_level)  # Set the log level
    log.info(sys.argv)  # Log the arguments
    return args


# Function to run Trimmomatic for single-end FASTQ
def trimmomatic_se(fastq1, crop_length, crop_length_tol,
                   phred_score_format, out_dir,
                   nth=1, java_heap=None):
    # Prefix for output file names (based on the input FASTQ file)
    prefix = os.path.join(out_dir,
                          os.path.basename(strip_ext_fastq(fastq1)))
    crop_length_tol = abs(crop_length_tol)  # Ensure tolerance is positive
    min_length = crop_length - crop_length_tol  # Minimum allowed read length
    # Output file name
    cropped = '{p}.crop_{cl}-{tol}bp.fastq.gz'.format(
        p=prefix, cl=crop_length, tol=crop_length_tol)

    # Set the Java heap size for Trimmomatic
    if java_heap is None:
        java_heap_param = '-Xmx6G'  # Default heap size
    else:
        java_heap_param = '-Xmx{}'.format(java_heap)

    # Handle Phred score format
    phred_score_format = phred_score_format.lower()
    if phred_score_format == 'auto':
        phred_score_param = ''
    elif phred_score_format == 'phred33':
        phred_score_param = '-phred33'
    elif phred_score_format == 'phred64':
        phred_score_param = '-phred64'
    else:
        raise ValueError('Wrong phred_score_format!')

    # Construct the Trimmomatic command and run it
    cmd = 'java -XX:ParallelGCThreads=1 {param} -jar {jar} SE -threads {nth} {phred_score_param} ' \
          '{fq1} {cropped} MINLEN:{ml} CROP:{cl}'.format(
        param=java_heap_param,
        jar=locate_trimmomatic(),  # Locate Trimmomatic executable
        nth=nth,
        phred_score_param=phred_score_param,
        fq1=fastq1,  # Input file
        cropped=cropped,  # Output cropped file
        ml=min_length,  # Minimum read length
        cl=crop_length,  # Crop length
    )
    run_shell_cmd(cmd)  # Execute the shell command

    return cropped


# Function to run Trimmomatic for paired-end FASTQs
def trimmomatic_pe(fastq1, fastq2, crop_length, crop_length_tol,
                   phred_score_format, out_dir_R1, out_dir_R2,
                   nth=1, java_heap=None):
    # Prefix for output file names for both read 1 and read 2
    prefix_R1 = os.path.join(
        out_dir_R1, os.path.basename(strip_ext_fastq(fastq1)))
    prefix_R2 = os.path.join(
        out_dir_R2, os.path.basename(strip_ext_fastq(fastq2)))

    # Ensure crop length tolerance is positive
    crop_length_tol = abs(crop_length_tol)
    min_length = crop_length - crop_length_tol  # Minimum allowed read length

    # Output file names for cropped reads
    cropped_R1 = '{p}.crop_{cl}-{tol}bp.fastq.gz'.format(
        p=prefix_R1, cl=crop_length, tol=crop_length_tol)
    cropped_R2 = '{p}.crop_{cl}-{tol}bp.fastq.gz'.format(
        p=prefix_R2, cl=crop_length, tol=crop_length_tol)

    # Temporary file names for Trimmomatic's internal use
    tmp_cropped_R1 = '{}.tmp'.format(cropped_R1)
    tmp_cropped_R2 = '{}.tmp'.format(cropped_R2)

    # Set the Java heap size for Trimmomatic
    if java_heap is None:
        java_heap_param = '-Xmx6G'  # Default heap size
    else:
        java_heap_param = '-Xmx{}'.format(java_heap)

    # Handle Phred score format
    phred_score_format = phred_score_format.lower()
    if phred_score_format == 'auto':
        phred_score_param = ''
    elif phred_score_format == 'phred33':
        phred_score_param = '-phred33'
    elif phred_score_format == 'phred64':
        phred_score_param = '-phred64'
    else:
        raise ValueError('Wrong phred_score_format!')

    # Construct the Trimmomatic command for paired-end reads and run it
    cmd = 'java -XX:ParallelGCThreads=1 {param} -jar {jar} PE -threads {nth} {phred_score_param} ' \
          '{fq1} {fq2} {cropped1} {tmp_cropped1} {cropped2} {tmp_cropped2} ' \
          'MINLEN:{ml} CROP:{cl}'.format(
        param=java_heap_param,
        jar=locate_trimmomatic(),  # Locate Trimmomatic executable
        nth=nth,
        phred_score_param=phred_score_param,
        fq1=fastq1,  # Input read 1
        fq2=fastq2,  # Input read 2
        cropped1=cropped_R1,  # Output for read 1
        tmp_cropped1=tmp_cropped_R1,  # Temporary file for read 1
        cropped2=cropped_R2,  # Output for read 2
        tmp_cropped2=tmp_cropped_R2,  # Temporary file for read 2
        ml=min_length,  # Minimum read length
        cl=crop_length,  # Crop length
    )
    run_shell_cmd(cmd)  # Execute the shell command

    # Remove temporary files
    rm_f([tmp_cropped_R1, tmp_cropped_R2])

    return cropped_R1, cropped_R2


# Main function to handle the workflow
def main():
    # Parse command-line arguments
    args = parse_arguments()

    # Create output directories if they do not exist
    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir_R1)
    if args.paired_end:
        mkdir_p(args.out_dir_R2)

    # Log the cropping parameters
    log.info(
        'Cropping fastqs with Trimmomatic... '
        'crop_length={cl}, crop_length_tol={clt}'.format(
            cl=args.crop_length,
            clt=args.crop_length_tol))
    
    # Call the appropriate function for paired or single-end reads
    if args.paired_end:
        cropped_R1, cropped_R2 = trimmomatic_pe(
            args.fastq1, args.fastq2,
            args.crop_length, args.crop_length_tol,
            args.phred_score_format,
            args.out_dir_R1, args.out_dir_R2,
            args.nth,
            args.trimmomatic_java_heap)
    else:
        cropped_R1 = trimmomatic_se(
            args.fastq1,
            args.crop_length, args.crop_length_tol,
            args.phred_score_format,
            args.out_dir_R1,
            args.nth,
            args.trimmomatic_java_heap)

    # List all files in the output directories
    log.info('List all files in output directory...')
    ls_l(args.out_dir_R1)
    if args.paired_end:
        ls_l(args.out_dir_R2)

    # Check if the output files are empty
    log.info('Checking if output is empty...')
    assert_file_not_empty(cropped_R1, help=
        'No reads in FASTQ after cropping. crop_length might be too high? '
        'While cropping, Trimmomatic (with MINLEN=crop_length-abs(crop_length_tol)) '
        'excludes all reads SHORTER than crop_length.')

    log.info('All done.')


# Run the main function if this script is executed directly
if __name__ == '__main__':
    main()
