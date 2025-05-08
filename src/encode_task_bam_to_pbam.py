#!/usr/bin/env python

# Import necessary libraries and modules
import sys  # Provides access to system-specific parameters and functions
import os  # For interacting with the operating system (file manipulation, paths, etc.)
import argparse  # For handling command-line arguments
from encode_lib_common import (
    log,  # Logging utility
    ls_l,  # Function to list files in a directory
    mkdir_p,  # Function to create directories if they do not exist
    rm_f,  # Function to remove files
)
from encode_lib_genomic import (
    bam_to_pbam,  # Function to convert BAM to PBAM format
)

# Function to parse command-line arguments
def parse_arguments():
    # Create the argument parser
    parser = argparse.ArgumentParser(prog='ENCODE bam to pbam', description='')

    # Define the arguments the script will accept
    parser.add_argument('bam', type=str, help='Path for BAM.')  # Path to the BAM file
    parser.add_argument('--ref-fa', type=str, help='Path for reference fasta.')  # Path to the reference FASTA file
    parser.add_argument('--delete-original-bam', action='store_true', help='Delete original BAM after conversion.')  # Flag to delete the original BAM file
    parser.add_argument('--out-dir', default='', type=str, help='Output directory.')  # Directory to save the output files
    parser.add_argument('--log-level', default='INFO', choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'CRITICAL', 'ERROR', 'CRITICAL'], help='Log level')  # Log level option

    # Parse the arguments
    args = parser.parse_args()

    # Set up logging based on the log level argument
    log.setLevel(args.log_level)
    log.info(sys.argv)  # Log the command used to run the script
    return args  # Return the parsed arguments

# Main function that controls the script workflow
def main():
    # Parse the command-line arguments
    args = parse_arguments()

    # Log the initialization process and create the output directory
    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)  # Create the output directory if it doesn't already exist

    # Log the conversion process from BAM to PBAM and call the function to perform the conversion
    log.info('Converting BAM into pBAM...')
    bam_to_pbam(args.bam, args.ref_fa, args.out_dir)  # Convert the BAM file to PBAM format

    # If the flag to delete the original BAM file is set, delete it after conversion
    if args.delete_original_bam:
        log.info('Deleting original BAM...')
        rm_f(args.bam)  # Remove the original BAM file

    # Log the list of all files in the output directory after the conversion
    log.info('List all files in output directory...')
    ls_l(args.out_dir)  # List the contents of the output directory

    # Log that the process is complete
    log.info('All done.')

# If this script is being run directly, execute the main function
if __name__ == '__main__':
    main()
