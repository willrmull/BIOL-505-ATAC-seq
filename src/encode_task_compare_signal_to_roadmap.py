#!/usr/bin/env python

# Import necessary libraries
import warnings
from matplotlib import pyplot as plt
import sys
import os
import argparse
from encode_lib_common import (
    strip_ext_bigwig, ls_l, log, mkdir_p)  # Utility functions for file handling
import numpy as np
import pandas as pd
import scipy.stats  # For statistical computations (Spearman correlation)
import matplotlib as mpl
mpl.use('Agg')  # Set up matplotlib to use non-GUI Agg backend

# Ignore warnings for cleaner output
warnings.filterwarnings("ignore")


def parse_arguments():
    """
    Parse the command-line arguments.
    """
    parser = argparse.ArgumentParser(prog='ENCODE compare signal to roadmap')
    
    # Add arguments for input files and parameters
    parser.add_argument('--bigwig', type=str,
                        help='BIGWIG file (from task macs2).')
    parser.add_argument('--dnase', type=str, help='DNase file.')
    parser.add_argument('--reg2map', type=str, help='Reg2map file.')
    parser.add_argument('--reg2map-bed', type=str, help='Reg2map bed file.')
    parser.add_argument('--roadmap-meta', type=str,
                        help='Roadmap metadata file.')
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')
    parser.add_argument('--log-level', default='INFO', help='Log level',
                        choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING',
                                 'CRITICAL', 'ERROR', 'CRITICAL'])
    
    # Parse arguments
    args = parser.parse_args()
    
    # Set logging level based on the user's input
    log.setLevel(args.log_level)
    log.info(sys.argv)  # Log the arguments for reference
    
    return args


def compare_to_roadmap(bw_file, regions_file, reg2map_file,
                       metadata, output_prefix):
    '''
    Takes a bigwig file and signal file, gets the bwAverageOverBed,
    then compares that signal with the signal in the Roadmap
    regions
    '''
    # Output filenames for results and logs
    out_file = '{0}.signal'.format(output_prefix)
    log_file = '{0}.roadmap_compare.log'.format(output_prefix)

    # Get the signal values for the peak regions using bigWigAverageOverBed
    # This command extracts signal values from the bigwig file over the regions specified in the bed file
    bw_average_over_bed = 'bigWigAverageOverBed {0} {1} {2}'.format(
        bw_file, regions_file, out_file)
    log.info(bw_average_over_bed)
    os.system(bw_average_over_bed)  # Run the command to get signal values

    # Read the output signal data
    sample_data = pd.read_table(out_file, header=None)
    sample_mean0_col = np.array(sample_data.iloc[:, 5])  # Extract signal values (6th column)

    # Read the roadmap signal data from the reg2map file
    roadmap_signals = pd.read_table(reg2map_file, compression='gzip')
    (nrow, ncol) = roadmap_signals.shape  # Get the shape of the roadmap signals

    # DataFrame to store correlation results
    results = pd.DataFrame(columns=('eid', 'corr'))
    
    # Open log file to write the results of correlation computation
    with open(log_file, 'w') as fp:
        for i in range(ncol):
            # Calculate Spearman correlation between each roadmap signal and the sample signal
            roadmap_i = roadmap_signals.iloc[:, i]
            spearman_corr = scipy.stats.spearmanr(np.array(roadmap_i),
                                                  sample_mean0_col)
            results.loc[i] = [roadmap_i.name, spearman_corr[0]]  # Store the correlation
            s = '{0}\t{1}'.format(roadmap_i.name, spearman_corr)  # Format output
            log.info(s)  # Log the result
            fp.write(s + '\n')  # Write result to log file

    # Read the metadata file to add more context to the plot (e.g., sample mnemonics)
    metadata = pd.read_table(metadata)
    metadata.columns = ['eid', 'mnemonic']  # Rename columns for clarity

    # Merge the correlation results with the metadata
    merged = pd.merge(metadata, results, on='eid')

    # Sort the results by correlation value
    sorted_results = merged.sort_values('corr', ascending=True)

    # Plot the results as a horizontal bar chart
    pos = np.array(range(ncol)) + 0.5
    fig = plt.figure(figsize=(5, int(ncol / 4)))  # Create figure for the plot
    plt.barh(pos, sorted_results['corr'], align='center', height=1.0)
    plt.yticks(pos, sorted_results['mnemonic'].tolist(), fontsize=7)  # Set y-axis labels (mnemonics)
    plt.xlabel('Spearmans correlation')  # Label for x-axis
    plt.title('Signal correlation to Roadmap DNase')  # Plot title
    plt.axis('tight')  # Remove extra space around the plot
    ax = plt.axes()  # Access axes for styling
    ax.yaxis.set_ticks_position('none')  # Remove y-axis ticks
    ax.spines['right'].set_visible(False)  # Remove right spine
    ax.spines['top'].set_visible(False)  # Remove top spine

    # Save the plot as a PNG file
    plot_img = output_prefix + '.roadmap_compare_plot.png'
    fig.savefig(plot_img, format='png', bbox_inches='tight')

    return plot_img


def main():
    # Read command-line parameters
    args = parse_arguments()
    
    # Extract filenames and parameters
    BIGWIG = args.bigwig
    DNASE = args.dnase
    OUTPUT_PREFIX = os.path.join(
        args.out_dir,
        os.path.basename(strip_ext_bigwig(BIGWIG)))  # Create output prefix from BigWig filename

    # Determine the Reg2Map file and bed file
    REG2MAP_BED = args.reg2map_bed if args.reg2map_bed and os.path.basename(
        args.reg2map_bed) != 'null' else DNASE
    REG2MAP = args.reg2map if args.reg2map and os.path.basename(
        args.reg2map) != 'null' else ''
    ROADMAP_META = args.roadmap_meta if args.roadmap_meta and os.path.basename(
        args.roadmap_meta) != 'null' else ''

    # Create output directory if it doesn't exist
    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # Compare the BigWig file to the Roadmap data
    compare_to_roadmap(BIGWIG, REG2MAP_BED, REG2MAP,
                       ROADMAP_META, OUTPUT_PREFIX)

    # List all files in the output directory for reference
    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    # Log that the process is complete
    log.info('All done.')


if __name__ == '__main__':
    main()  # Execute the main function
