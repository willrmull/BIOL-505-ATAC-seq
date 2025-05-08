#!/usr/bin/env python

# Import necessary libraries and functions
import sys
import os
import argparse
from encode_lib_common import (
    log, ls_l, mkdir_p, strip_ext)
from encode_lib_log_parser import parse_flagstat_qc


def parse_arguments():
    """
    Parses the command-line arguments for the script. These arguments 
    include paths to SAMstats log files (both non-mito and mito),
    output directory, and logging level.
    """
    # Initialize argument parser for the command-line interface
    parser = argparse.ArgumentParser(
        prog='ENCODE frac mito',  # Program name
        description='Calculates fraction of mito reads')  # Description of what the script does

    # Define required arguments: paths to SAMstats log files for non-mito and mito reads
    parser.add_argument('non_mito_samstat', type=str,
                        help='Path for SAMstats log file (non-mito reads)')
    parser.add_argument('mito_samstat', type=str,
                        help='Path for SAMstats log file (mito only reads)')

    # Optional argument: output directory for results
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory for the result files.')

    # Optional argument: log level for controlling verbosity
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO',
                                 'WARNING', 'CRITICAL', 'ERROR',
                                 'CRITICAL'],
                        help='Log level for the script.')

    # Parse the arguments from the command line
    args = parser.parse_args()

    # Set up logging with the specified log level
    log.setLevel(args.log_level)
    
    # Log the command line that was executed (for reproducibility)
    log.info(sys.argv)
    
    # Return the parsed arguments to be used in the main function
    return args


def frac_mito(non_mito_samstat, mito_samstat, out_dir):
    """
    This function calculates the fraction of mito reads based on the
    SAMstats log files. It reads the log files, extracts the relevant
    information (mapped reads for non-mito and mito), and computes
    the fraction of mito reads.

    Arguments:
    - non_mito_samstat: Path to SAMstats log file for non-mito reads
    - mito_samstat: Path to SAMstats log file for mito reads
    - out_dir: Output directory where the results will be saved

    Returns:
    - frac_mito_qc: Path to the output QC file with fraction of mito reads
    """
    # Define the prefix for output filenames based on the input non-mito SAMstats log file
    prefix = os.path.join(
        out_dir,
        os.path.basename(strip_ext(non_mito_samstat, 'non_mito.samstats.qc'))
    )
    
    # Define the path for the output file that will store the fraction of mito reads
    frac_mito_qc = '{}.frac_mito.qc'.format(prefix)

    # Parse the SAMstats log files to extract the relevant statistics
    non_mito_samstat_dict = parse_flagstat_qc(non_mito_samstat)  # Parse non-mito log file
    mito_samstat_dict = parse_flagstat_qc(mito_samstat)  # Parse mito log file

    # Determine the correct key for the "mapped" reads in the log files
    # This ensures compatibility with older versions of the log files
    if 'mapped' in non_mito_samstat_dict:
        key_mapped = 'mapped'
    elif 'mapped_reads' in non_mito_samstat_dict:
        key_mapped = 'mapped_reads'

    # Extract the number of mapped non-mito reads
    Rn = non_mito_samstat_dict[key_mapped]

    # Similarly, determine the correct key for "mapped" reads in the mito log file
    if 'mapped' in mito_samstat_dict:
        key_mapped = 'mapped'
    elif 'mapped_reads' in mito_samstat_dict:
        key_mapped = 'mapped_reads'

    # Extract the number of mapped mito reads
    Rm = mito_samstat_dict[key_mapped]

    # Calculate the fraction of mito reads
    frac = float(Rm) / float(Rn + Rm)

    # Write the calculated statistics (non-mito reads, mito reads, and fraction of mito reads)
    with open(frac_mito_qc, 'w') as fp:
        fp.write('non_mito_reads\t{}\n'.format(Rn))  # Write non-mito reads
        fp.write('mito_reads\t{}\n'.format(Rm))  # Write mito reads
        fp.write('frac_mito_reads\t{}\n'.format(frac))  # Write fraction of mito reads

    # Return the path to the output QC file
    return frac_mito_qc


def main():
    """
    Main function that orchestrates the script. It reads the command-line
    arguments, creates the necessary directories, computes the fraction
    of mito reads, and outputs the results to the specified directory.
    """
    # Parse the command-line arguments
    args = parse_arguments()

    # Log initialization and output directory creation
    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # Log the step of calculating the fraction of mito reads
    frac_mito_qc = frac_mito(args.non_mito_samstat,
                              args.mito_samstat,
                              args.out_dir)

    # Log the list of files generated in the output directory
    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    # Log that the script has completed
    log.info('All done.')


if __name__ == '__main__':
    main()  # Call the main function when the script is run
