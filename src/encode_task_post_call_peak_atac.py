#!/usr/bin/env python

# Import necessary modules for system operations, argument parsing, and specific functions
import sys
import argparse
from encode_lib_common import (
    assert_file_not_empty,  # Custom function to check if a file is empty
    log,                    # Custom logging module
    ls_l,                   # Custom function to list files in a directory
    mkdir_p,                # Custom function to create directories
)
from encode_lib_genomic import (
    peak_to_bigbed,         # Convert peak data to BigBed format
    peak_to_hammock,        # Convert peak data to Hammock format
    get_region_size_metrics,  # Calculate region size metrics for peaks
    get_num_peaks,          # Get the number of peaks in the dataset
    peak_to_starch,         # Convert peaks to Starch format
)
from encode_lib_blacklist_filter import blacklist_filter  # Custom function for filtering peaks based on a blacklist
from encode_lib_frip import frip  # Custom function to calculate FRiP (Fraction of Reads in Peaks)


def parse_arguments():
    """
    Parse command-line arguments using argparse, allowing the user to specify input files, parameters, and options.
    """
    parser = argparse.ArgumentParser(prog='ENCODE post_call_peak (atac)', description='Process and analyze peak data')

    # Required argument: Path to the peak file
    parser.add_argument(
        'peak', type=str,
        help='Path for PEAK file. Peak filename should be "*.*Peak.gz". '
             'e.g. rep1.narrowPeak.gz'
    )

    # Optional argument: Path to the TAG-ALIGN file
    parser.add_argument('--ta', type=str, help='TAG-ALIGN file.')

    # Required argument: Type of peak file (narrowPeak, regionPeak, broadPeak, or gappedPeak)
    parser.add_argument('--peak-type', type=str, required=True,
                        choices=['narrowPeak', 'regionPeak', 'broadPeak', 'gappedPeak'],
                        help='Peak file type.')

    # Optional argument: Path to chromosome size file (2-column format)
    parser.add_argument('--chrsz', type=str, help='2-col chromosome sizes file.')

    # Optional argument: Path to a blacklist BED file
    parser.add_argument('--blacklist', type=str, help='Blacklist BED file.')

    # Optional argument: Regex pattern to filter peak chromosomes
    parser.add_argument('--regex-bfilt-peak-chr-name',
                        help='Keep chromosomes matching this pattern only in .bfilt. peak files.')

    # Optional argument: Memory limit for the task, used to determine sort memory
    parser.add_argument('--mem-gb', type=float, default=4.0,
                        help='Max. memory for this job in GB. '
                             'This will be used to determine GNU sort -S (defaulting to 0.5 of this value). '
                             'It should be total memory for this task (not memory per thread).')

    # Optional argument: Directory to store output files
    parser.add_argument('--out-dir', default='', type=str, help='Output directory.')

    # Optional argument: Log level to control the verbosity of logging
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'CRITICAL', 'ERROR', 'CRITICAL'],
                        help='Log level')

    # Parse the arguments
    args = parser.parse_args()

    # If no blacklist is specified, default it to an empty string
    if args.blacklist is None or args.blacklist.endswith('null'):
        args.blacklist = ''

    # Set logging level according to user input
    log.setLevel(args.log_level)
    
    # Log the command used to run the script
    log.info(sys.argv)
    
    return args


def main():
    """
    Main function that orchestrates the processing pipeline for peak files, including filtering, conversion, and analysis.
    """
    # Step 1: Parse the command-line arguments
    args = parse_arguments()

    # Step 2: Initialize and create the output directory if it doesn't exist
    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # Step 3: Perform blacklist filtering on the peak file
    log.info('Blacklist-filtering peaks...')
    bfilt_peak = blacklist_filter(
        args.peak, args.blacklist, args.regex_bfilt_peak_chr_name, args.out_dir
    )

    # Step 4: Check if the blacklist-filtered peak file is empty
    log.info('Checking if output is empty...')
    assert_file_not_empty(bfilt_peak)

    # Step 5: Convert the blacklist-filtered peak file to BigBed format
    log.info('Converting peak to bigbed...')
    peak_to_bigbed(bfilt_peak, args.peak_type, args.chrsz, args.mem_gb, args.out_dir)

    # Step 6: Convert the blacklist-filtered peak file to Starch format
    log.info('Converting peak to starch...')
    peak_to_starch(bfilt_peak, args.out_dir)

    # Step 7: Convert the blacklist-filtered peak file to Hammock format
    log.info('Converting peak to hammock...')
    peak_to_hammock(bfilt_peak, args.mem_gb, args.out_dir)

    # Step 8: Calculate FRiP (Fraction of Reads in Peaks) without fragment length
    log.info('FRiP without fragment length...')
    frip(args.ta, bfilt_peak, args.out_dir)

    # Step 9: Calculate and generate QC/plot for region size metrics of the blacklist-filtered peaks
    log.info('Calculating (blacklist-filtered) peak region size QC/plot...')
    get_region_size_metrics(bfilt_peak)

    # Step 10: Calculate the total number of peaks in the blacklist-filtered peak file
    log.info('Calculating number of peaks (blacklist-filtered)...')
    get_num_peaks(bfilt_peak)

    # Step 11: List all files in the output directory to confirm results
    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    # Final step: Log completion message
    log.info('All done.')


if __name__ == '__main__':
    # Run the main function when the script is executed
    main()
