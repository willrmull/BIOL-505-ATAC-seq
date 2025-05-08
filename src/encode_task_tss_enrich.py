#!/usr/bin/env python

# Import necessary libraries for genomic data analysis, plotting, and logging
import matplotlib as mpl
mpl.use('Agg')  # Use Agg backend for matplotlib (for script-based plotting without a GUI)
import pybedtools  # For working with BED files
import numpy as np  # Numerical computations
from matplotlib import mlab  # For computing percentiles
from matplotlib import pyplot as plt  # For plotting
import sys  # System-specific parameters
import os  # For file system operations
import argparse  # Argument parsing for command-line options
from encode_lib_common import (
    strip_ext_bam, ls_l, log, logging, mkdir_p, rm_f)  # Utility functions
from encode_lib_genomic import (
    remove_read_group, samtools_index)  # Genomic utilities (e.g., removing read groups)
import metaseq  # A library for metagenomic signal processing
import warnings
warnings.filterwarnings("ignore")  # Ignore warnings (e.g., in plotting or data processing)


def parse_arguments():
    """
    Parse the command-line arguments to set up the parameters for the script.
    The arguments control various aspects like input files, output directories,
    and specific analysis options.
    """
    parser = argparse.ArgumentParser(prog='ENCODE TSS enrichment.')
    
    # Define command-line arguments
    parser.add_argument('--read-len-log', type=str,
                        help='Read length log file (from aligner task).')
    parser.add_argument('--read-len', type=int,
                        help='Read length (integer). This is ignored if '
                             '--read-len-log is defined.')
    parser.add_argument('--nodup-bam', type=str,
                        help='Raw BAM file (from task filter).')
    parser.add_argument('--chrsz', type=str,
                        help='2-col chromosome sizes file.')
    parser.add_argument('--tss', type=str, help='TSS definition bed file.')
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')
    parser.add_argument('--log-level', default='INFO', help='Log level',
                        choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING',
                                 'CRITICAL', 'ERROR', 'CRITICAL'])

    # Parse arguments from the command line
    args = parser.parse_args()

    # Ensure that either read_len_log or read_len is defined
    if args.read_len_log is None and args.read_len is None:
        raise ValueError('Either --read-len-log or --read-len must be defined.')

    # Set the logging level and log the command used to run the script
    log.setLevel(args.log_level)
    log.info(sys.argv)

    return args


def make_tss_plot(bam_file, tss, prefix, chromsizes,
                  read_len, bins=400, bp_edge=2000,
                  processes=8, greenleaf_norm=True):
    """
    Generate TSS enrichment plots. This function computes the TSS enrichment
    around the transcription start sites (TSS) and creates two plots: 
    one showing the aggregation plot and another showing the signal strength at each TSS.
    It also computes and returns the mean and standard deviation of the plot.

    Parameters:
        bam_file (str): The BAM file containing aligned reads.
        tss (str): The BED file defining the TSS positions.
        prefix (str): Output file prefix for generated plots.
        chromsizes (str): The chromosome sizes file (for binning).
        read_len (int): The length of the reads in the BAM file.
        bins (int): Number of bins to use in plotting.
        bp_edge (int): The number of base pairs to extend beyond the TSS.
        processes (int): Number of processes to parallelize the computation.
        greenleaf_norm (bool): Whether to use Greenleaf-style normalization.

    Returns:
        tuple: Paths to the generated plot files (small plot, large plot, log file).
    """
    logging.info('Generating TSS plot...')
    
    # Set file names for the TSS enrichment plots and QC log
    tss_plot_file = '{0}.tss_enrich.png'.format(prefix)
    tss_plot_large_file = '{0}.large_tss_enrich.png'.format(prefix)
    tss_log_file = '{0}.tss_enrich.qc'.format(prefix)

    # Load the TSS file using pybedtools (BED file)
    tss = pybedtools.BedTool(tss)
    tss_ext = tss.slop(b=bp_edge, g=chromsizes)  # Extend TSS regions by bp_edge on both sides

    # Load BAM file using metaseq to generate genomic signal (aligned read data)
    bam = metaseq.genomic_signal(bam_file, 'bam')
    
    # Shift the reads so that the center is aligned with the cut site of the TSS
    bam_array = bam.array(tss_ext, bins=bins, shift_width=-read_len/2,
                          processes=processes, stranded=True)

    # Apply Greenleaf normalization (optional) to adjust for noise in the bins
    if greenleaf_norm:
        # Use enough bins to cover 100 bp at the ends of the region (signal noise)
        num_edge_bins = int(100/(2*bp_edge/bins))
        bin_means = bam_array.mean(axis=0)
        avg_noise = (sum(bin_means[:num_edge_bins]) +
                     sum(bin_means[-num_edge_bins:]))/(2*num_edge_bins)
        bam_array /= avg_noise  # Normalize by average noise
    else:
        # If not Greenleaf norm, normalize by mapped read count
        bam_array /= bam.mapped_read_count() / 1e6

    # Plot the TSS enrichment
    fig = plt.figure()
    ax = fig.add_subplot(111)
    x = np.linspace(-bp_edge, bp_edge, bins)

    ax.plot(x, bam_array.mean(axis=0), color='r', label='Mean')
    ax.axvline(0, linestyle=':', color='k')  # TSS marker at 0

    # Find the peak of the TSS enrichment (max value)
    tss_point_val = max(bam_array.mean(axis=0))

    # Write the peak value to a log file
    with open(tss_log_file, 'w') as fp:
        fp.write(str(tss_point_val))

    # Set labels for the plot
    ax.set_xlabel('Distance from TSS (bp)')
    if greenleaf_norm:
        ax.set_ylabel('TSS Enrichment')
    else:
        ax.set_ylabel('Average read coverage (per million mapped reads)')
    ax.legend(loc='best')

    # Save the smaller TSS plot
    fig.savefig(tss_plot_file)

    # Generate a more detailed plot with more information
    upper_prct = 99
    if mlab.prctile(bam_array.ravel(), upper_prct) == 0.0:
        upper_prct = 100.0

    # Use metaseq to generate a detailed plot
    plt.rcParams['font.size'] = 8
    fig = metaseq.plotutils.imshow(bam_array,
                                   x=x,
                                   figsize=(5, 10),
                                   vmin=5, vmax=upper_prct, percentile=True,
                                   line_kwargs=dict(color='k', label='All'),
                                   fill_kwargs=dict(color='k', alpha=0.3),
                                   sort_by=bam_array.mean(axis=1))

    # Save the detailed plot
    fig.savefig(tss_plot_large_file)

    return tss_plot_file, tss_plot_large_file, tss_log_file


def main():
    """
    Main function that coordinates the execution of the script. It parses arguments,
    processes the input BAM file, generates TSS enrichment plots, and saves the results.
    """
    # Parse the command-line arguments
    args = parse_arguments()

    # Assign input arguments to variables for easier reference
    CHROMSIZES = args.chrsz
    TSS = args.tss if args.tss and os.path.basename(args.tss) != 'null' else ''
    FINAL_BAM = args.nodup_bam
    OUTPUT_PREFIX = os.path.join(
        args.out_dir,
        os.path.basename(strip_ext_bam(FINAL_BAM)))
    
    # Index the BAM file (to prepare for random access during processing)
    samtools_index(FINAL_BAM)

    # Remove the read group from the BAM file (clean up for analysis)
    RG_FREE_FINAL_BAM = remove_read_group(FINAL_BAM)

    # Initialize the output directory
    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # Retrieve the read length (either from log file or argument)
    if args.read_len_log:
        with open(args.read_len_log, 'r') as fp:
            read_len = int(fp.read().strip())
    elif args.read_len:
        read_len = args.read_len
    else:
        read_len = None

    # Generate TSS enrichment plots and log the results
    tss_plot, tss_large_plot, tss_enrich_qc = \
        make_tss_plot(FINAL_BAM,
                      TSS,
                      OUTPUT_PREFIX,
                      CHROMSIZES,
                      read_len)

    # Remove the temporary BAM file without read group
    rm_f(RG_FREE_FINAL_BAM)

    # List all files in the output directory
    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    # Indicate completion of the process
    log.info('All done.')


if __name__ == '__main__':
    main()  # Run the main function when the script is executed
