#!/usr/bin/env python

# Import necessary libraries and modules
import sys
import os
import argparse
from encode_lib_common import (
    log, ls_l, mkdir_p, rm_f, run_shell_cmd, strip_ext_bam)
from encode_lib_genomic import (
    samtools_index)
from encode_lib_blacklist_filter import blacklist_filter_bam


def parse_arguments():
    """Parse command-line arguments."""
    # Initialize argument parser
    parser = argparse.ArgumentParser(
        prog='ENCODE DCC Fingerprint/JSD plot.')
    
    # Add positional argument for the list of filtered experiment BAM files
    parser.add_argument(
        'bams', nargs='+', type=str,
        help='List of paths for filtered experiment BAM files.')
    
    # Add optional argument for the filtered control BAM file
    parser.add_argument('--ctl-bam', type=str, default='',
                        help='Path for filtered control BAM file.')
    
    # Add optional argument for the blacklist BED file
    parser.add_argument('--blacklist', type=str, default='',
                        help='Blacklist BED file.')
    
    # Add argument for MAPQ threshold (default 30) for low MAPQ reads removal
    parser.add_argument('--mapq-thresh', default=30, type=int,
                        help='Threshold for low MAPQ reads removal.')
    
    # Add argument for the number of threads to use (default 1)
    parser.add_argument('--nth', type=int, default=1,
                        help='Number of threads to parallelize.')
    
    # Add optional argument for output directory
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')
    
    # Add argument for logging level (e.g., INFO, DEBUG)
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO',
                                 'WARNING', 'CRITICAL', 'ERROR',
                                 'CRITICAL'],
                        help='Log level')
    
    # Parse the arguments and return the parsed arguments object
    args = parser.parse_args()

    # Set logging level based on user input
    log.setLevel(args.log_level)
    log.info(sys.argv)  # Log the command used to run the script
    return args


def fingerprint(bams, ctl_bam, blacklist, mapq_thresh, nth, out_dir):
    """Generate Fingerprint plots and calculate Jensen-Shannon Divergence (JSD) between samples."""
    
    # Make BAM index (.bai) files first and filter BAMs with blacklist
    filtered_bams = []
    for bam in bams:
        # Apply blacklist filtering to each BAM file
        filtered_bam = blacklist_filter_bam(bam, blacklist, out_dir)
        # Index the filtered BAM file using samtools
        samtools_index(filtered_bam, nth)
        filtered_bams.append(filtered_bam)

    # If a control BAM is provided, filter and index it as well
    filtered_ctl_bam = None
    if ctl_bam:
        filtered_ctl_bam = blacklist_filter_bam(ctl_bam, blacklist, out_dir)
        samtools_index(filtered_ctl_bam, nth)

    # Define prefix for output files (using the first BAM file's base name)
    prefix = os.path.join(out_dir,
                          os.path.basename(strip_ext_bam(bams[0])))
    
    # Define output file names for the plot and temporary log
    plot_png = '{}.jsd_plot.png'.format(prefix)
    tmp_log = '{}.jsd.tmp'.format(prefix)

    # Initialize lists for sample labels, BAM paths, and JSD QC result file names
    labels = []
    bam_paths = []
    jsd_qcs = []

    # For each filtered experiment BAM, prepare the QC information and add labels
    for i, bam in enumerate(filtered_bams):
        prefix_ = os.path.join(out_dir,
                               os.path.basename(strip_ext_bam(bam)))
        jsd_qcs.append('rep{}.{}.jsd.qc'.format(i+1, prefix_))
        labels.append('rep{}'.format(i+1))  # Label as rep1, rep2, etc.
        bam_paths.append(bam)

    # If control BAM was provided, add it to the labels and paths
    if filtered_ctl_bam:
        labels.append('ctl1')
        bam_paths.append(filtered_ctl_bam)

    # Construct the command to run `plotFingerprint` to generate the fingerprint plot
    cmd = 'LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8 plotFingerprint -b {} '
    if filtered_ctl_bam:
        cmd += '--JSDsample {} '.format(filtered_ctl_bam)
    cmd += '--labels {} '
    cmd += '--outQualityMetrics {} '
    cmd += '--minMappingQuality {} '
    cmd += '-T "Fingerprints of different samples" '
    cmd += '--numberOfProcessors {} '
    cmd += '--plotFile {}'
    
    # Format the command with the appropriate arguments (BAM files, labels, etc.)
    cmd = cmd.format(
        ' '.join(bam_paths),
        ' '.join(labels),
        tmp_log,
        mapq_thresh,
        nth,
        plot_png)
    
    # Run the command using the helper function `run_shell_cmd`
    run_shell_cmd(cmd)

    # Clean up intermediate files (blacklist-filtered BAM files)
    if filtered_ctl_bam:
        rm_f(filtered_ctl_bam)
    rm_f(filtered_bams)

    # Parse the temporary log file to extract the JSD QC results for each experiment replicate
    with open(tmp_log, 'r') as fp:
        for i, line in enumerate(fp.readlines()):
            if i == 0:
                continue  # Skip header line
            if i > len(jsd_qcs):
                break  # Stop if there are more lines than expected
            with open(jsd_qcs[i-1], 'w') as fp2:
                # Write the JSD QC data for each replicate (excluding the first column)
                fp2.write('\t'.join(line.strip().split('\t')[1:]))
    
    # Remove the temporary log file after parsing
    rm_f(tmp_log)
    
    return plot_png, jsd_qcs


def main():
    """Main function to execute the fingerprint plotting and JSD calculation."""
    
    # Parse command-line arguments
    args = parse_arguments()

    # Log the initialization process and create the output directory if it doesn't exist
    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # Log the plotting and calculation process
    log.info('Plotting Fingerprint on BAMs and calculating JSD...')
    
    # Call the fingerprint function to generate the plot and QC results
    plot_png, jsd_qcs = fingerprint(
        args.bams, args.ctl_bam, args.blacklist, args.mapq_thresh,
        args.nth, args.out_dir)

    # Log the list of output files in the output directory
    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    # Log the completion of the process
    log.info('All done.')


if __name__ == '__main__':
    # Run the main function
    main()
