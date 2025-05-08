#!/usr/bin/env python

# Standard libraries
import sys
import os
import argparse

# Import utility functions from ENCODE pipeline's shared library
from encode_lib_common import (
    assert_file_not_empty,  # Check that output file is not empty
    copy_f_to_f,            # Copy a file from one path to another
    log,                    # Logging utility
    ls_l,                   # List directory contents
    mkdir_p,                # Make directory if it doesn’t exist
    run_shell_cmd,          # Run a shell command with error handling
    strip_ext_fastq         # Strip .fastq/.fq/.gz extension from filename
)

def parse_arguments(debug=False):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        prog='ENCODE DCC fastq merger.',
        description='Trims FASTQ files to a fixed number of base pairs.')

    # Input FASTQ file
    parser.add_argument('fastq', type=str,
                        help='FASTQ file to be trimmed.')

    # Number of base pairs to keep
    parser.add_argument('--trim-bp', type=int, default=50,
                        help='Number of base pairs to retain after trimming.')

    # Output directory for the trimmed FASTQ
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')

    # Logging level
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'CRITICAL', 'ERROR'],
                        help='Log verbosity level.')

    # Parse arguments and initialize logging
    args = parser.parse_args()
    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args


def trim_fastq(fastq, trim_bp, out_dir):
    """
    Trim reads in the FASTQ to a specified length using external script `trimfastq.py`.

    Parameters:
        fastq (str): Path to input FASTQ file.
        trim_bp (int): Number of base pairs to retain.
        out_dir (str): Output directory for trimmed FASTQ.

    Returns:
        trimmed (str): Path to the trimmed output FASTQ file.
    """
    # Construct output file path
    prefix = os.path.join(out_dir, os.path.basename(strip_ext_fastq(fastq)))
    trimmed = '{}.trim_{}bp.fastq.gz'.format(prefix, trim_bp)

    # Run external Python script (`trimfastq.py`) and compress output
    cmd = 'python $(which trimfastq.py) {} {} | gzip -nc > {}'.format(
        fastq, trim_bp, trimmed)
    run_shell_cmd(cmd)

    # Check if any sequences were too short after trimming
    cmd2 = (
        'zcat -f {} | '
        '(grep "sequences shorter than desired length" || true) | wc -l'
    ).format(trimmed)

    # If all reads were too short, fall back to untrimmed input file
    if int(run_shell_cmd(cmd2)) > 0:
        copy_f_to_f(fastq, trimmed)

    return trimmed


def main():
    """Main function to coordinate trimming process."""
    # Parse command-line inputs
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # Trim reads to specified length
    log.info('Trimming FASTQ to {} bp...'.format(args.trim_bp))
    trimmed = trim_fastq(args.fastq, args.trim_bp, args.out_dir)

    # Validate output
    assert_file_not_empty(trimmed)

    # List contents of output directory for verification
    log.info('Listing output directory contents...')
    ls_l(args.out_dir)

    log.info('All done.')


# Entry point of the script
if __name__ == '__main__':
    main()
