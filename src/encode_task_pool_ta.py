#!/usr/bin/env python

# Import necessary libraries and modules
import sys  # Provides access to system-specific parameters and functions
import os  # Provides functions to interact with the operating system
import argparse  # For parsing command-line arguments
from encode_lib_common import (
    log, ls_l, make_hard_link, mkdir_p, run_shell_cmd, strip_ext_ta
)  # Import custom utility functions from the encode_lib_common module

# Function to parse command-line arguments
def parse_arguments():
    # Initialize an argument parser with a description for the program
    parser = argparse.ArgumentParser(prog='ENCODE DCC TAGALIGN pooler.',
                                     description='Pool multiple TAGALIGN files into one.')

    # Define arguments for the script
    parser.add_argument('tas', nargs='+', type=str,
                        help='List of TAGALIGNs to be pooled.')  # Takes one or more TAGALIGN files as input
    parser.add_argument('--prefix', type=str,
                        help='Basename prefix.')  # Optional prefix for the output file name
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')  # Directory where the pooled file will be saved
    parser.add_argument('--col',
                        help='Number of columns to keep in a pooled TAGALIGN. '
                             'Keep all columns if not defined.')  # Optional, defines how many columns to keep in the pooled file
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'CRITICAL', 'ERROR', 'CRITICAL'],
                        help='Log level')  # Defines the log level for output (e.g., INFO, DEBUG, etc.)

    # Parse the command-line arguments
    args = parser.parse_args()

    # Set logging level based on parsed argument
    log.setLevel(args.log_level)

    # Log the command line arguments for debugging
    log.info(sys.argv)

    # Return the parsed arguments for use in other functions
    return args


# Function to pool multiple TAGALIGN files into one
def pool_ta(tas, col, basename_prefix, out_dir):
    # If there is more than one TAGALIGN file, proceed with pooling
    if len(tas) > 1:
        # Define the output file prefix
        if basename_prefix is not None:
            # Use the provided prefix if it exists
            prefix = os.path.join(out_dir, basename_prefix)
        else:
            # Otherwise, use the first TAGALIGN file's name (without extension) as the prefix
            prefix = os.path.join(out_dir,
                                  os.path.basename(strip_ext_ta(tas[0])))
        
        # Define the name for the pooled TAGALIGN file
        pooled_ta = '{}.pooled.tagAlign.gz'.format(prefix)

        # Build the shell command to pool the files
        cmd = 'zcat -f {} | '  # Decompress the TAGALIGN files
        if col is not None:
            # If a column number is specified, restrict the output to those columns
            cmd += 'cut -f 1-{} | '.format(col)
        
        # Pipe the output to gzip for compression
        cmd += 'gzip -nc > {}'
        
        # Format the command with the list of TAGALIGN files and the output file
        cmd = cmd.format(
            ' '.join(tas),  # Join the input files into a single string for the zcat command
            pooled_ta  # The name of the output file
        )

        # Run the shell command
        run_shell_cmd(cmd)

        # Return the path to the pooled TAGALIGN file
        return pooled_ta
    else:
        # Raise an error if there is less than two TAGALIGN files
        raise ValueError('Needs at least two TAs (or BEDs) to be pooled.')


# Main function that coordinates the execution of the script
def main():
    # Parse command-line arguments
    args = parse_arguments()

    # Log and create the output directory if it does not exist
    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # Log the pooling process
    log.info('Pooling TAGALIGNs...')
    pooled_ta = pool_ta(args.tas, args.col, args.prefix, args.out_dir)

    # Log the listing of all files in the output directory
    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    # Log the completion message
    log.info('All done.')


# Check if the script is being run directly (not imported as a module)
if __name__ == '__main__':
    main()
