#!/usr/bin/env python

# Import required libraries
import warnings
import numpy as np
from collections import namedtuple
from scipy.signal import find_peaks_cwt
from matplotlib import pyplot as plt
import sys
import os
import argparse
from encode_lib_common import (
    strip_ext_bam, ls_l, log, rm_f, pdf2png)
from encode_lib_genomic import (
    remove_read_group, locate_picard)
import matplotlib as mpl
mpl.use('Agg')  # Use non-GUI backend for matplotlib (useful for server environments)

# Suppress warnings for cleaner output
warnings.filterwarnings("ignore")

# Define a namedtuple for storing the result of a QC check
QCResult = namedtuple('QCResult', ['metric', 'qc_pass', 'message'])

# Define an infinite constant to represent no upper or lower limit
INF = float("inf")


class QCCheck(object):
    """Base class for all QC checks"""
    def __init__(self, metric):
        self.metric = metric  # Store the metric name (e.g., 'Fraction of reads in NFR')

    def check(self, value):
        """Base check function (always returns True by default)"""
        return True

    def message(self, value, qc_pass):
        """Generate a message based on the check result"""
        return ('{}\tOK'.format(value) if qc_pass
                else '{}\tFailed'.format(value))

    def __call__(self, value):
        """Execute the check and return a QCResult object"""
        qc_pass = self.check(value)
        return QCResult(self.metric, qc_pass, self.message(value, qc_pass))


class QCIntervalCheck(QCCheck):
    """QC Check where value must fall within a specific interval (inclusive)"""
    def __init__(self, metric, lower, upper):
        super(QCIntervalCheck, self).__init__(metric)
        self.lower = lower
        self.upper = upper

    def check(self, value):
        """Check if the value is within the given range"""
        return self.lower <= value <= self.upper

    def message(self, value, qc_pass):
        """Message indicating whether the value falls within the specified range"""
        return ('{}\tOK'.format(value) if qc_pass else
                '{}\tout of range [{}, {}]'.format(value, self.lower, self.upper))


class QCLessThanEqualCheck(QCIntervalCheck):
    """QC Check where value must be less than or equal to the upper bound"""
    def __init__(self, metric, upper):
        super(QCLessThanEqualCheck, self).__init__(metric, -INF, upper)


class QCGreaterThanEqualCheck(QCIntervalCheck):
    """QC Check where value must be greater than or equal to the lower bound"""
    def __init__(self, metric, lower):
        super(QCGreaterThanEqualCheck, self).__init__(metric, lower, INF)


class QCHasElementInRange(QCCheck):
    """QC Check for the presence of an element within a specified range"""
    def __init__(self, metric, lower, upper):
        super(QCHasElementInRange, self).__init__(metric)
        self.lower = lower
        self.upper = upper

    def check(self, elems):
        """Check if any elements fall within the specified range"""
        return (len([elem for elem in elems
                     if self.lower <= elem <= self.upper]) > 0)

    def message(self, elems, qc_pass):
        """Generate a message indicating if elements are within the range"""
        return ('OK' if qc_pass else
                'Cannot find element in range [{}, {}]'.format(self.lower, self.upper))


def parse_arguments():
    """Parse command-line arguments"""
    parser = argparse.ArgumentParser(prog='ENCODE fragment length stat')
    parser.add_argument('--nodup-bam', type=str,
                        help='Raw BAM file (from task filter).')
    parser.add_argument('--picard-java-heap',
                        help='Picard\'s Java max. heap: java -jar picard.jar '
                             '-Xmx[MAX_HEAP]')
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')
    parser.add_argument('--log-level', default='INFO', help='Log level',
                        choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING',
                                 'CRITICAL', 'ERROR', 'CRITICAL'])
    args = parser.parse_args()  # Parse the arguments
    log.setLevel(args.log_level)  # Set logging level
    log.info(sys.argv)  # Log the command used to run the script
    return args


def read_picard_histogram(data_file):
    """Read Picard insert size distribution data"""
    with open(data_file) as fp:
        for line in fp:
            if line.startswith('## HISTOGRAM'):
                break
        data = np.loadtxt(fp, skiprows=1)

    return data


def get_insert_distribution(final_bam, prefix, java_heap=None):
    """Run Picard's CollectInsertSizeMetrics to get insert size distribution"""
    log.info('insert size distribution...')
    insert_data = '{0}.inserts.hist_data.log'.format(prefix)
    insert_plot = '{0}.inserts.hist_graph.pdf'.format(prefix)
    if java_heap is None:
        java_heap_param = '-Xmx6G'  # Default Java heap size
    else:
        java_heap_param = '-Xmx{}'.format(java_heap)
    
    # Construct the Picard command to calculate insert size metrics
    graph_insert_dist = ('java {4} -XX:ParallelGCThreads=1 -jar '
                         '{3} '
                         'CollectInsertSizeMetrics '
                         'INPUT={0} OUTPUT={1} H={2} '
                         'VERBOSITY=ERROR QUIET=TRUE '
                         'USE_JDK_DEFLATER=TRUE USE_JDK_INFLATER=TRUE '
                         'W=1000 STOP_AFTER=5000000').format(final_bam,
                                                             insert_data,
                                                             insert_plot,
                                                             locate_picard(),
                                                             java_heap_param)
    log.info(graph_insert_dist)  # Log the command being run
    os.system(graph_insert_dist)  # Execute the command
    return insert_data, insert_plot


def fragment_length_qc(data, prefix):
    """Perform fragment length quality control checks"""
    results = []

    NFR_UPPER_LIMIT = 150
    MONO_NUC_LOWER_LIMIT = 150
    MONO_NUC_UPPER_LIMIT = 300

    # Percentage of NFR (nucleosome free region) reads vs. total reads
    nfr_reads = data[data[:, 0] < NFR_UPPER_LIMIT][:, 1]
    percent_nfr = nfr_reads.sum() / data[:, 1].sum()
    results.append(
        QCGreaterThanEqualCheck('Fraction of reads in NFR', 0.4)(percent_nfr))

    # Ratio of NFR vs. mononucleosome reads
    mono_nuc_reads = data[
        (data[:, 0] > MONO_NUC_LOWER_LIMIT) &
        (data[:, 0] <= MONO_NUC_UPPER_LIMIT)][:, 1]
    percent_nfr_vs_mono_nuc = (nfr_reads.sum() / mono_nuc_reads.sum())
    results.append(
        QCGreaterThanEqualCheck('NFR / mono-nuc reads', 2.5)(
            percent_nfr_vs_mono_nuc))

    # Check for the presence of peaks in the data (e.g., NFR, Mono-Nuc, Di-Nuc)
    pos_start_val = data[0, 0]
    peaks = find_peaks_cwt(data[:, 1], np.array([25]))
    nuc_range_metrics = [
        ('Presence of NFR peak', 20 - pos_start_val, 90 - pos_start_val),
        ('Presence of Mono-Nuc peak', 120 - pos_start_val, 250 - pos_start_val),
        ('Presence of Di-Nuc peak', 300 - pos_start_val, 500 - pos_start_val)]
    for range_metric in nuc_range_metrics:
        results.append(QCHasElementInRange(*range_metric)(peaks))

    # Write results to file
    out = prefix + '.nucleosomal.qc'
    with open(out, 'w') as fp:
        for elem in results:
            fp.write(
                '\t'.join(
                    [elem.metric, str(elem.qc_pass), elem.message]) + '\n')

    return out


def fragment_length_plot(data_file, prefix, peaks=None):
    """Generate a plot of the fragment length distribution"""
    try:
        data = read_picard_histogram(data_file)  # Read data from file
    except IOError:
        return ''
    except TypeError:
        return ''

    fig = plt.figure()
    plt.bar(data[:, 0], data[:, 1])  # Plot the histogram
    plt.xlim((0, 1000))  # Set x-axis limits

    if peaks:  # If peaks are provided, highlight them on the plot
        peak_vals = [data[peak_x, 1] for peak_x in peaks]
        plt.plot(peaks, peak_vals, 'ro')  # Mark peaks with red dots

    plot_pdf = prefix + '.fraglen_dist.pdf'
    plot_png = prefix + '.fraglen_dist.png'

    # Save the plot as PDF and PNG
    fig.savefig(plot_pdf, format='pdf')
    pdf2png(plot_pdf, os.path.dirname(plot_pdf))  # Convert PDF to PNG
    rm_f(plot_pdf)  # Remove the temporary PDF file

    return plot_png


def main():
    """Main function to execute the QC checks and generate plots"""
    # Parse command-line arguments
    args = parse_arguments()

    FINAL_BAM = args.nodup_bam
    OUTPUT_PREFIX = os.path.join(
        args.out_dir,
        os.path.basename(strip_ext_bam(FINAL_BAM)))  # Output file prefix
    RG_FREE_FINAL_BAM = remove_read_group(FINAL_BAM)  # Remove read groups
    JAVA_HEAP = args.picard_java_heap  # Java heap size for Picard

    # Run Picard's CollectInsertSizeMetrics to get insert size distribution
    insert_data, insert_plot = get_insert_distribution(RG_FREE_FINAL_BAM,
                                                       OUTPUT_PREFIX,
                                                       JAVA_HEAP)

    # Run fragment length QC
    fragment_length_qc(read_picard_histogram(insert_data),
                       OUTPUT_PREFIX)

    # Generate the plot of the fragment length distribution
    fragment_length_plot(insert_data, OUTPUT_PREFIX)

    # Clean up temporary files
    rm_f(RG_FREE_FINAL_BAM)

    # Log files in the output directory
    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')  # Final log message


if __name__ == '__main__':
    main()  # Run the main function
