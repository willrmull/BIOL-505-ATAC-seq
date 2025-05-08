#!/usr/bin/env python

# Import necessary libraries and functions
import sys
import os
import argparse
from encode_lib_common import (
    log, ls_l, mkdir_p, pdf2png, rm_f, run_shell_cmd,
    strip_ext_ta)
from encode_lib_genomic import (
    subsample_ta_pe, subsample_ta_se)
from encode_lib_log_parser import parse_xcor_score


def parse_arguments():
    """
    Parses the command-line arguments for the script, allowing users to
    specify the input TAGALIGN file, output directory, and various
    parameters like subsampling, cross-correlation peak strandshift,
    and exclusion range.
    """
    # Initialize the argument parser
    parser = argparse.ArgumentParser(
        prog='ENCODE DCC cross-correlation analysis.',
        description='Performs cross-correlation analysis on TAGALIGN files.'
    )

    # Define required arguments
    parser.add_argument('ta', type=str,
                        help='Path for TAGALIGN file.')

    # Optional arguments for various parameters and options
    parser.add_argument('--mito-chr-name', default='chrM',
                        help='Mito chromosome name.')
    parser.add_argument('--subsample', type=int, default=0,
                        help='Subsample TAGALIGN reads.')
    parser.add_argument('--speak', type=int, default=-1,
                        help='User-defined cross-correlation peak strandshift \
                        (-speak= in run_spp.R). Disabled if -1.')
    parser.add_argument('--exclusion-range-min', type=int,
                        help='User-defined exclusion range minimum used for '
                             '-x=${xcor_exclusion_range_min}:${xcor_exclusion_range_max}')
    parser.add_argument('--exclusion-range-max', type=int,
                        help='User-defined exclusion range maximum used for '
                             '-x=${xcor_exclusion_range_min}:${xcor_exclusion_range_max}')
    parser.add_argument('--chip-seq-type', choices=['tf', 'histone'],
                        help='Type of ChIP-seq pipeline (histone or tf).')
    parser.add_argument('--paired-end', action="store_true",
                        help='Indicates paired-end TAGALIGN.')
    parser.add_argument('--nth', type=int, default=1,
                        help='Number of threads to parallelize.')
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO',
                                 'WARNING', 'CRITICAL', 'ERROR',
                                 'CRITICAL'],
                        help='Log level.')

    # Parse the arguments
    args = parser.parse_args()

    # Set the logging level and log the executed command
    log.setLevel(args.log_level)
    log.info(sys.argv)
    
    # Return the parsed arguments for further use
    return args


def get_exclusion_range_max(ta, chip_seq_type):
    """
    Estimates the maximum exclusion range for cross-correlation analysis.
    The function calculates the average read length from the TAGALIGN file
    and adjusts based on the ChIP-seq type (tf or histone).
    """
    # Estimate read length from TAGALIGN file
    cmd0 = "zcat -f {} > {}.tmp".format(ta, ta)
    run_shell_cmd(cmd0)

    # Calculate average read length
    cmd = "head -n 100 {}.tmp | awk 'function abs(v) "
    cmd += "{{return v < 0 ? -v : v}} BEGIN{{sum=0}} "
    cmd += "{{sum+=abs($3-$2)}} END{{print int(sum/NR)}}'"
    cmd = cmd.format(ta)
    read_len = int(run_shell_cmd(cmd))

    # Clean up temporary file
    rm_f(ta+'.tmp')

    # Adjust exclusion range based on ChIP-seq type
    if chip_seq_type == 'tf':
        return max(read_len + 10, 50)
    elif chip_seq_type == 'histone':
        return max(read_len + 10, 100)
    else:
        raise NotImplementedError('chip_seq_type not supported')


def xcor(ta, speak, mito_chr_name,
         nth, out_dir, chip_seq_type=None,
         exclusion_range_min=None, exclusion_range_max=None):
    """
    Performs the cross-correlation analysis using an R script (run_spp.R).
    It also processes the results to generate cross-correlation plots and
    compute fragment length estimations.
    """
    # Generate the output file paths
    prefix = os.path.join(out_dir,
                          os.path.basename(strip_ext_ta(ta)))
    xcor_plot_pdf = '{}.cc.plot.pdf'.format(prefix)
    xcor_score = '{}.cc.qc'.format(prefix)
    fraglen_txt = '{}.cc.fraglen.txt'.format(prefix)

    # Prepare exclusion range parameter if provided
    if chip_seq_type is not None and exclusion_range_min is not None:
        if exclusion_range_max is None:
            exclusion_range_max = get_exclusion_range_max(ta, chip_seq_type)

        exclusion_range_param = ' -x={}:{}'.format(
            exclusion_range_min, exclusion_range_max)
    else:
        exclusion_range_param = ''

    # Construct the command to run the R script for cross-correlation
    cmd1 = 'Rscript --max-ppsize=500000 $(which run_spp.R) -rf -c={} -p={} '
    cmd1 += '-filtchr="{}" -savp={} -out={} {}'
    cmd1 += exclusion_range_param
    cmd1 = cmd1.format(
        ta,
        nth,
        mito_chr_name,
        xcor_plot_pdf,
        xcor_score,
        '-speak={}'.format(speak) if speak >= 0 else '')  # Include speak parameter if needed
    run_shell_cmd(cmd1)

    # Clean up the xcor_score file by removing unwanted columns
    cmd2 = 'sed -r \'s/,[^\\t]+//g\' -i {}'
    cmd2 = cmd2.format(xcor_score)
    run_shell_cmd(cmd2)

    # Parse the cross-correlation score file to extract fragment length
    cmd3 = 'echo {} > {}'.format(
        parse_xcor_score(xcor_score)['estimated_fragment_len'],
        fraglen_txt)
    run_shell_cmd(cmd3)

    # Convert the cross-correlation plot from PDF to PNG format
    xcor_plot_png = pdf2png(xcor_plot_pdf, out_dir)
    
    # Return the generated files: plot PDF, PNG, score, and fragment length file
    return xcor_plot_pdf, xcor_plot_png, xcor_score, fraglen_txt


def main():
    """
    Main function that coordinates the execution of the script. It handles
    argument parsing, temporary file creation (subsampling TAGALIGN), 
    cross-correlation analysis, and file cleanup.
    """
    # Parse the command-line arguments
    args = parse_arguments()

    # Initialize output directory and log info
    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # Temporary file list to keep track of files that need to be removed later
    temp_files = []  

    # Subsample the TAGALIGN file based on the user’s input (if required)
    log.info('Subsampling TAGALIGN for cross-correlation...')
    if args.paired_end:
        ta_subsampled = subsample_ta_pe(
            args.ta, args.subsample, True,
            args.mito_chr_name, True, args.out_dir)
    else:
        ta_subsampled = subsample_ta_se(
            args.ta, args.subsample, True,
            args.mito_chr_name, args.out_dir)
    temp_files.append(ta_subsampled)  # Add the subsampled TAGALIGN to the temp list

    # Perform the cross-correlation analysis
    log.info('Performing cross-correlation analysis...')
    xcor_plot_pdf, xcor_plot_png, xcor_score, fraglen_txt = xcor(
        ta_subsampled, args.speak, args.mito_chr_name, args.nth, args.out_dir,
        args.chip_seq_type, args.exclusion_range_min, args.exclusion_range_max)

    # Remove temporary files generated during the process
    log.info('Removing temporary files...')
    rm_f(temp_files)

    # List all generated files in the output directory for the user
    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    # Log that the process is complete
    log.info('All done.')


if __name__ == '__main__':
    main()  # Run the main function when the script is executed
