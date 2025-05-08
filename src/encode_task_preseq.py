#!/usr/bin/env python

# Import necessary libraries for warnings, numerical computation, and plotting
import warnings
import numpy as np
from matplotlib import pyplot as plt
import sys
import os
import argparse
from encode_lib_common import (
    strip_ext_bam, ls_l, log, logging, rm_f, run_shell_cmd)
from encode_lib_genomic import (
    remove_read_group, locate_picard, samtools_sort)
import matplotlib as mpl
mpl.use('Agg')  # Use non-interactive Agg backend for matplotlib (useful for running scripts on a server)

# Ignore all warnings
warnings.filterwarnings("ignore")


def parse_arguments():
    # Function to parse command line arguments
    parser = argparse.ArgumentParser(prog='ENCODE preseq')
    parser.add_argument('--paired-end', action="store_true", help='Paired-end BAM.')
    parser.add_argument('--bam', type=str, help='Raw BAM file.')
    parser.add_argument('--picard-java-heap',
                        help='Picard\'s Java max. heap: java -jar picard.jar -Xmx[MAX_HEAP]')
    parser.add_argument('--nth', type=int, default=1, help='Number of threads to parallelize.')
    parser.add_argument('--mem-gb', type=float,
                        help='Max. memory for samtools sort in GB. It should be total memory for this task (not memory per thread).')
    parser.add_argument('--out-dir', default='', type=str, help='Output directory.')
    parser.add_argument('--log-level', default='INFO', help='Log level',
                        choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'CRITICAL', 'ERROR', 'CRITICAL'])
    args = parser.parse_args()  # Parse the arguments passed from the command line
    log.setLevel(args.log_level)  # Set the log level
    log.info(sys.argv)  # Log the command line arguments for tracking
    return args


def get_picard_complexity_metrics(aligned_bam, prefix, java_heap=None):
    '''
    Run Picard's EstimateLibraryComplexity to calculate the library complexity
    '''
    # Remove redundant or malformed read group information from the BAM file
    out_file = '{0}.picardcomplexity.qc'.format(prefix)
    
    # Set Java heap size for Picard; default is 6G if not provided
    if java_heap is None:
        java_heap_param = '-Xmx6G'
    else:
        java_heap_param = '-Xmx{}'.format(java_heap)
    
    # Command to run Picard's EstimateLibraryComplexity tool
    get_gc_metrics = (
        'mkdir -p tmp_java && java -Djava.io.tmpdir=$PWD/tmp_java '
        '{3} -XX:ParallelGCThreads=1 -jar '
        '{2} '
        'EstimateLibraryComplexity INPUT={0} OUTPUT={1} '
        'USE_JDK_DEFLATER=TRUE USE_JDK_INFLATER=TRUE '
        'VERBOSITY=ERROR '
        'QUIET=TRUE && rm -rf tmp_java').format(
        aligned_bam, out_file, locate_picard(), java_heap_param)
    
    # Execute the Picard command
    os.system(get_gc_metrics)

    # Extract the actual estimated library size from Picard's output
    header_seen = False
    est_library_size = 0
    with open(out_file, 'r') as fp:
        for line in fp:
            if header_seen:
                est_library_size = int(float(line.strip().split()[-1]))
                break
            if 'ESTIMATED_LIBRARY_SIZE' in line:
                header_seen = True

    return est_library_size


def run_preseq(bam_w_dups, prefix, nth=1, mem_gb=None):
    '''
    Runs the preseq tool to estimate sequencing complexity.
    '''
    # First sort the BAM file, since the sorted BAM is required for preseq
    sort_bam = samtools_sort(bam_w_dups, nth, mem_gb)

    logging.info('Running preseq...')
    preseq_data = '{0}.preseq.dat'.format(prefix)  # Output file for preseq data
    preseq_log = '{0}.preseq.log'.format(prefix)  # Output log for preseq run

    # Run preseq using the sorted BAM file and capture the data and log output
    run_shell_cmd(
        'preseq lc_extrap -P -B -o {preseq_data} {sort_bam} '
        '-seed 1 -v 2> {preseq_log}'.format(
            preseq_data=preseq_data,
            sort_bam=sort_bam,
            preseq_log=preseq_log,
        )
    )
    
    rm_f(sort_bam)  # Remove the sorted BAM file after the preseq run

    return preseq_data, preseq_log


def get_preseq_plot(data_file, prefix):
    '''
    Generates a plot based on the preseq data
    '''
    try:
        # Load the preseq data from the file, skipping the first row (header)
        data = np.loadtxt(data_file, skiprows=1)
    except IOError:
        return ''  # If the file doesn't exist, return an empty string
    
    data /= 1e6  # Scale the data to millions of reads

    # Create a new figure for the plot
    fig = plt.figure()

    # Plot the average expected yield (red line)
    plt.plot(data[:, 0], data[:, 1], 'r-')

    # Plot the confidence intervals (blue dashed lines)
    ci_lower, = plt.plot(data[:, 0], data[:, 2], 'b--')
    ci_upper, = plt.plot(data[:, 0], data[:, 3], 'b--')
    plt.legend([ci_lower], ['95% confidence interval'], loc=4)

    # Set plot title and labels
    plt.title('Preseq estimated yield')
    plt.xlabel('Sequenced fragments [ millions ]')
    plt.ylabel('Expected distinct fragments [ millions ]')

    # Save the plot as a PNG file
    plot_img = prefix + '.preseq.png'
    fig.savefig(plot_img, format='png')

    return plot_img  # Return the filename of the plot


def main():
    # Main function to orchestrate the analysis workflow
    args = parse_arguments()  # Read command line arguments

    # Prepare file paths and variables
    ALIGNED_BAM = args.bam
    OUTPUT_PREFIX = os.path.join(
        args.out_dir,
        os.path.basename(strip_ext_bam(ALIGNED_BAM)))  # Strip the .bam extension for output file naming
    RG_FREE_ALIGNED_BAM = remove_read_group(ALIGNED_BAM)  # Remove read group information for processing
    JAVA_HEAP = args.picard_java_heap

    # If paired-end, calculate library complexity using Picard's tool
    if args.paired_end:
        picard_est_lib_size = get_picard_complexity_metrics(
            RG_FREE_ALIGNED_BAM, OUTPUT_PREFIX, JAVA_HEAP)
    else:
        picard_est_lib_size = None  # Not applicable for single-end

    # Run preseq to estimate sequencing complexity
    preseq_data, preseq_log = run_preseq(
        ALIGNED_BAM, OUTPUT_PREFIX, args.nth, args.mem_gb)  # Sorted BAM passed to preseq

    # Generate the preseq plot
    get_preseq_plot(preseq_data, OUTPUT_PREFIX)

    # Write the Picard estimated library size to a file if calculated
    if picard_est_lib_size is not None:
        picard_est_lib_size_file = OUTPUT_PREFIX + '.picard_est_lib_size.qc'
        with open(picard_est_lib_size_file, 'w') as fp:
            fp.write(str(picard_est_lib_size)+'\n')

    # Clean up by removing the temporary BAM file without read group info
    rm_f(RG_FREE_ALIGNED_BAM)

    # List all the files in the output directory for logging
    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')  # Final log message


if __name__ == '__main__':
    main()  # Run the main function if the script is executed directly
