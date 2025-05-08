#!/usr/bin/env python

# Standard libraries
import sys
import os
import argparse

# ENCODE utility functions
from encode_lib_common import (
    assert_file_not_empty,  # Ensure output file is non-empty
    get_num_lines,          # Count number of lines in a file
    log,                    # Logging utility
    ls_l,                   # List files in directory
    mkdir_p,                # Create directory if it doesn't exist
    rm_f,                   # Remove files
    run_shell_cmd,          # Run a shell command
    strip_ext_ta            # Strip extension from TAG-ALIGN file
)

# ENCODE genomic utilities for subsampling
from encode_lib_genomic import (
    subsample_ta_pe,  # Subsample paired-end TAG-ALIGN
    subsample_ta_se   # Subsample single-end TAG-ALIGN
)


def parse_arguments():
    """
    Parse command-line arguments.
    Returns a namespace with parsed arguments.
    """
    parser = argparse.ArgumentParser(
        prog='ENCODE DCC control TAG-ALIGN subsampler.',
        description=(
            'This script randomly subsamples a control TAG-ALIGN file. '
            'It does not validate whether the number of requested reads exceeds the available reads. '
            'If it does, the file will simply be shuffled without downsampling.')
    )

    # Required input: TAG-ALIGN file path
    parser.add_argument('ta', type=str,
                        help='Path to input TAG-ALIGN file.')

    # Optional flag: input is paired-end
    parser.add_argument('--paired-end', action="store_true",
                        help='Flag for paired-end TAG-ALIGN format.')

    # Required: number of reads to subsample
    parser.add_argument('--subsample', default=0, type=int,
                        help='Number of reads to subsample (must be positive).')

    # Output directory for subsampled result
    parser.add_argument('--out-dir', default='', type=str,
                        help='Directory where output will be saved.')

    # Logging verbosity
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO',
                                 'WARNING', 'CRITICAL', 'ERROR'],
                        help='Logging level.')

    args = parser.parse_args()

    # Validate that subsample value is greater than 0
    if not args.subsample:
        raise ValueError('--subsample should be a positive integer.')

    # Set logging level and log command-line arguments
    log.setLevel(args.log_level)
    log.info(sys.argv)

    return args


def main():
    """Main function to run subsampling process."""
    args = parse_arguments()

    log.info('Creating output directory (if needed)...')
    mkdir_p(args.out_dir)

    # Run the appropriate subsampling function depending on paired-end flag
    if args.paired_end:
        log.info('Subsampling paired-end TAG-ALIGN...')
        subsampled_ta = subsample_ta_pe(
            args.ta, args.subsample,
            non_mito=False,       # Do not exclude mitochondrial reads
            mito_chr_name=None,   # No specific mito chromosome name
            r1_only=False,        # Use both R1 and R2 reads
            out_dir=args.out_dir)
    else:
        log.info('Subsampling single-end TAG-ALIGN...')
        subsampled_ta = subsample_ta_se(
            args.ta, args.subsample,
            non_mito=False,       # Do not exclude mitochondrial reads
            mito_chr_name=None,   # No specific mito chromosome name
            out_dir=args.out_dir)

    # Check that the output file is not empty
    log.info('Verifying that output file is non-empty...')
    assert_file_not_empty(subsampled_ta)

    # List contents of output directory
    log.info('Listing output files...')
    ls_l(args.out_dir)

    log.info('Subsampling complete.')


# Run main if script is executed directly
if __name__ == '__main__':
    main()
