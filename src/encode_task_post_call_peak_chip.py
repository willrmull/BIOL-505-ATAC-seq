#!/usr/bin/env python


import sys
import os
import argparse
from encode_lib_common import (
    assert_file_not_empty,  # Function to assert if the file is not empty.
    log,                    # Logger utility for logging messages.
    ls_l,                   # List directory contents with detailed information.
    mkdir_p,                # Function to create directories (if they don't already exist).
)
from encode_lib_genomic import (
    peak_to_bigbed,         # Function to convert peak data to bigBed format.
    peak_to_hammock,        # Function to convert peak data to hammock format.
    get_region_size_metrics, # Function to calculate region size metrics (QC).
    get_num_peaks,          # Function to calculate the number of peaks in the data.
    peak_to_starch,         # Function to convert peak data to starch format.
)
from encode_lib_blacklist_filter import blacklist_filter  # Function for filtering peaks using a blacklist.
from encode_lib_frip import frip_shifted  # Function for shifted FRiP calculation.

def parse_arguments():
    """
    Parses command-line arguments using argparse. Provides a detailed explanation of each argument.
    
    Returns:
        args (Namespace): Parsed arguments from the command line.
    """
    parser = argparse.ArgumentParser(prog='ENCODE post_call_peak (chip)', description='')
    
    # Add arguments to the parser, each with a description of the argument's purpose.
    parser.add_argument('peak', type=str,
                        help='Path for PEAK file. Peak filename should be "*.*Peak.gz". '
                             'e.g. rep1.narrowPeak.gz')
    parser.add_argument('--ta', type=str,
                        help='TAG-ALIGN file.')
    parser.add_argument('--peak-type', type=str, required=True,
                        choices=['narrowPeak', 'regionPeak', 'broadPeak', 'gappedPeak'],
                        help='Peak file type.')
    parser.add_argument('--fraglen', type=int, required=True,
                        help='Fragment length.')
    parser.add_argument('--chrsz', type=str,
                        help='2-col chromosome sizes file.')
    parser.add_argument('--blacklist', type=str,
                        help='Blacklist BED file.')
    parser.add_argument('--regex-bfilt-peak-chr-name',
                        help='Keep chromosomes matching this pattern only '
                             'in .bfilt. peak files.')
    parser.add_argument('--mem-gb', type=float, default=4.0,
                        help='Max. memory for this job in GB. '
                             'This will be used to determine GNU sort -S (defaulting to 0.5 of this value). '
                             'It should be total memory for this task (not memory per thread).')
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'CRITICAL', 'ERROR', 'CRITICAL'],
                        help='Log level')
    
    # Parse the command-line arguments and return them.
    args = parser.parse_args()
    
    # If the blacklist argument is not provided or contains 'null', set it to an empty string.
    if args.blacklist is None or args.blacklist.endswith('null'):
        args.blacklist = ''
    
    # Set the logging level based on the provided argument.
    log.setLevel(args.log_level)
    
    # Log the command-line arguments for debugging purposes.
    log.info(sys.argv)
    
    return args


def main():
    """
    Main function that orchestrates the sequence of tasks by calling helper functions.
    """
    # Read and parse the command-line arguments.
    args = parse_arguments()

    # Log initialization and creation of the output directory.
    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # Log the blacklist-filtering of the peaks.
    log.info('Blacklist-filtering peaks...')
    bfilt_peak = blacklist_filter(
        args.peak, args.blacklist, args.regex_bfilt_peak_chr_name, args.out_dir)

    # Ensure that the output file is not empty.
    log.info('Checking if output is empty...')
    assert_file_not_empty(bfilt_peak)

    # Convert the filtered peak data to bigBed format.
    log.info('Converting peak to bigbed...')
    peak_to_bigbed(bfilt_peak, args.peak_type, args.chrsz,
                   args.mem_gb, args.out_dir)

    # Convert the filtered peak data to starch format.
    log.info('Converting peak to starch...')
    peak_to_starch(bfilt_peak, args.out_dir)

    # Convert the filtered peak data to hammock format.
    log.info('Converting peak to hammock...')
    peak_to_hammock(bfilt_peak, args.mem_gb, args.out_dir)

    # Perform shifted FRiP calculation with fragment length.
    log.info('Shifted FRiP with fragment length...')
    frip_qc = frip_shifted(args.ta, bfilt_peak,
                           args.chrsz, args.fraglen, args.out_dir)

    # Calculate region size metrics for the filtered peaks.
    log.info('Calculating (blacklist-filtered) peak region size QC/plot...')
    region_size_qc, region_size_plot = get_region_size_metrics(bfilt_peak)

    # Calculate the number of peaks in the filtered data.
    log.info('Calculating number of peaks (blacklist-filtered)...')
    num_peak_qc = get_num_peaks(bfilt_peak)

    # List all files in the output directory to check results.
    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    # Final log message indicating that all tasks are completed.
    log.info('All done.')


if __name__ == '__main__':
    # Run the main function when the script is executed directly.
    main()
