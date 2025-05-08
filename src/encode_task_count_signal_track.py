#!/usr/bin/env python

# Import necessary libraries
import sys
import os
import argparse
from encode_lib_common import (
    log, ls_l, mkdir_p, rm_f, run_shell_cmd, strip_ext_ta,
    get_gnu_sort_param,
)


def parse_arguments():
    """
    Parses command-line arguments for the script.
    This function defines the input parameters like the TAGALIGN file,
    chromosome size file, output directory, memory options, and log level.
    """
    # Initialize argument parser for the command-line interface
    parser = argparse.ArgumentParser(
        prog='ENCODE DCC Count signal track generation')  # Program description

    # Define arguments to be passed in the command-line
    parser.add_argument('ta', type=str,
                        help='Path for TAGALIGN file.')  # Path to the TAGALIGN file (input data)
    parser.add_argument('--chrsz', type=str,
                        help='2-col chromosome sizes file.')  # Chromosome size file for genome coverage
    parser.add_argument('--mem-gb', type=float, default=4.0,
                        help='Max. memory for this job in GB. '
                             'This will be used to determine GNU sort -S (defaulting to 0.5 of this value). '
                             'It should be total memory for this task (not memory per thread).')  # Maximum memory for the job
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')  # Directory where the output files will be saved
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO',
                                 'WARNING', 'CRITICAL', 'ERROR',
                                 'CRITICAL'],
                        help='Log level')  # Log level for script execution (controls verbosity)

    # Parse the command-line arguments
    args = parser.parse_args()
    
    # Set up logging with the specified log level
    log.setLevel(args.log_level)
    
    # Log the arguments for reproducibility and debugging
    log.info(sys.argv)
    
    # Return the parsed arguments
    return args


def count_signal_track(ta, chrsz, mem_gb, out_dir):
    """
    Counts the signal from the provided TAGALIGN file by generating
    positive and negative strand BigWig files. It also manages the
    temporary files used in the process.
    
    Arguments:
    - ta: Path to the TAGALIGN file
    - chrsz: Path to the chromosome size file
    - mem_gb: Maximum memory for the job in GB
    - out_dir: Directory to store the output files
    
    Returns:
    - pos_bw: Path to the output positive strand BigWig file
    - neg_bw: Path to the output negative strand BigWig file
    """
    # Generate output file paths by removing the extension from the TAGALIGN file and appending relevant suffixes
    prefix = os.path.join(out_dir, os.path.basename(strip_ext_ta(ta)))
    pos_bw = '{}.positive.bigwig'.format(prefix)  # Output positive strand BigWig file
    neg_bw = '{}.negative.bigwig'.format(prefix)  # Output negative strand BigWig file

    # Temporary files to store intermediate results (bedgraph format)
    pos_bedgraph = '{}.positive.bedgraph'.format(prefix)
    neg_bedgraph = '{}.negative.bedgraph'.format(prefix)

    # List to store temporary files for cleanup later
    temp_files = []

    # Generate positive strand bedGraph from TAGALIGN using bedtools
    run_shell_cmd(
        'zcat -f {ta} | sort -k1,1 -k2,2n {sort_param} | '
        'bedtools genomecov -5 -bg -strand + -g {chrsz} -i stdin > {pos_bedgraph}'.format(
            ta=ta,  # Path to TAGALIGN file
            sort_param=get_gnu_sort_param(mem_gb * 1024 ** 3, ratio=0.5),  # Get sort parameters based on memory
            chrsz=chrsz,  # Path to chromosome size file
            pos_bedgraph=pos_bedgraph,  # Path to output positive strand bedgraph
        )
    )

    # Generate negative strand bedGraph from TAGALIGN using bedtools
    run_shell_cmd(
        'zcat -f {ta} | sort -k1,1 -k2,2n {sort_param} | '
        'bedtools genomecov -5 -bg -strand - -g {chrsz} -i stdin > {neg_bedgraph}'.format(
            ta=ta,  # Path to TAGALIGN file
            sort_param=get_gnu_sort_param(mem_gb * 1024 ** 3, ratio=0.5),  # Get sort parameters based on memory
            chrsz=chrsz,  # Path to chromosome size file
            neg_bedgraph=neg_bedgraph,  # Path to output negative strand bedgraph
        )
    )

    # Convert the positive strand bedgraph to BigWig format
    run_shell_cmd(
        'bedGraphToBigWig {pos_bedgraph} {chrsz} {pos_bw}'.format(
            pos_bedgraph=pos_bedgraph,  # Path to positive strand bedgraph
            chrsz=chrsz,  # Path to chromosome size file
            pos_bw=pos_bw,  # Path to output positive strand BigWig file
        )
    )

    # Convert the negative strand bedgraph to BigWig format
    run_shell_cmd(
        'bedGraphToBigWig {neg_bedgraph} {chrsz} {neg_bw}'.format(
            neg_bedgraph=neg_bedgraph,  # Path to negative strand bedgraph
            chrsz=chrsz,  # Path to chromosome size file
            neg_bw=neg_bw,  # Path to output negative strand BigWig file
        )
    )

    # Add the generated bedgraph files to the list of temporary files for cleanup
    temp_files.append(pos_bedgraph)
    temp_files.append(neg_bedgraph)

    # Clean up temporary files
    rm_f(temp_files)

    # Return the paths of the output BigWig files
    return pos_bw, neg_bw


def main():
    """
    Main function that orchestrates the workflow of the script.
    It handles reading command-line arguments, creating output directories,
    generating signal tracks, and logging the results.
    """
    # Parse command-line arguments
    args = parse_arguments()

    # Log the initialization step and create the output directory
    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # Log the signal track generation step
    log.info('Generating count signal tracks...')
    
    # Call the function to generate the signal tracks and capture the output BigWig file paths
    pos_bw, neg_bw = count_signal_track(
        args.ta,  # Path to the TAGALIGN file
        args.chrsz,  # Path to the chromosome size file
        args.mem_gb,  # Maximum memory for the job in GB
        args.out_dir  # Directory for saving the output files
    )

    # Log the listing of all files in the output directory
    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    # Log the completion of the task
    log.info('All done.')


if __name__ == '__main__':
    main()  # Call the main function when the script is run

