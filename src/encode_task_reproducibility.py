#!/usr/bin/env python

# Import standard Python libraries
import sys
import os
import argparse

# Import ENCODE library utilities for file handling and logging
from encode_lib_common import (
    copy_f_to_f,
    get_num_lines,
    infer_n_from_nC2,
    infer_pair_label_from_idx,
    log,
    mkdir_p,
)

# Import genomic QC and format conversion tools
from encode_lib_genomic import (
    peak_to_bigbed,
    peak_to_hammock,
    get_region_size_metrics,
    get_num_peaks,
    peak_to_starch,
)


# Argument parsing function
def parse_arguments():
    parser = argparse.ArgumentParser(
        prog='ENCODE DCC reproducibility QC.',
        description='Handles reproducibility QC on IDR or overlap peaks.'
    )

    # List of peak files from all true replicate pairs (can be empty for single replicate)
    parser.add_argument('peaks', type=str, nargs='*',
                        help='Sorted list of replicate pair peak files. Example for 4 reps: 0,1 0,2 ... 2,3.')

    # Required: Peak files from pseudo-replicates (e.g., rep1-pr1 vs rep1-pr2, etc.)
    parser.add_argument('--peaks-pr', type=str, nargs='+', required=True,
                        help='List of peak files from pseudo-replicates.')

    # Peak file from pooled pseudo-replicates (optional if no true replicate peaks)
    parser.add_argument('--peak-ppr', type=str,
                        help='Peak file from pooled pseudo-replicates.')

    # Peak file format
    parser.add_argument('--peak-type', type=str, default='narrowPeak',
                        choices=['narrowPeak', 'regionPeak', 'broadPeak', 'gappedPeak'],
                        help='Peak file type.')

    # Chromosome sizes file
    parser.add_argument('--chrsz', type=str,
                        help='Chromosome sizes file (2-column).')

    # Basename for outputs
    parser.add_argument('--prefix', type=str,
                        help='Output file prefix.')

    # Memory for sorting (in GB)
    parser.add_argument('--mem-gb', type=float, default=4.0,
                        help='Memory for GNU sort (GB). Used at 50% allocation.')

    # Output directory
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')

    # Logging verbosity
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'CRITICAL', 'ERROR'],
                        help='Log level')

    args = parser.parse_args()

    # Validation: ensure pseudo-replicate count matches expected combinations
    if len(args.peaks_pr) != infer_n_from_nC2(len(args.peaks)):
        raise argparse.ArgumentTypeError(
            'Mismatch in number of true replicate peak files and pseudo-replicate peaks.')

    # Set up logging
    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args


def main():
    # Parse input arguments
    args = parse_arguments()

    # Ensure output directory exists
    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    log.info('Reproducibility QC...')

    # Get number of peaks from each pseudo-replicate comparison
    N = [get_num_lines(peak) for peak in args.peaks_pr]

    # If multiple replicates provided
    if len(args.peaks):
        num_rep = infer_n_from_nC2(len(args.peaks))  # infer number of replicates
        num_peaks_tr = [get_num_lines(peak) for peak in args.peaks]

        Nt = max(num_peaks_tr)                       # max #peaks in true replicate comparisons
        Np = get_num_lines(args.peak_ppr)            # #peaks in pooled pseudo-replicate

        rescue_ratio = float(max(Np, Nt)) / float(min(Np, Nt))  # Reproducibility ratio: pooled vs true replicate
        self_consistency_ratio = float(max(N)) / float(min(N))  # Consistency across pseudo-replicates

        Nt_idx = num_peaks_tr.index(Nt)              # index of the peak file with the most peaks
        label_tr = infer_pair_label_from_idx(num_rep, Nt_idx)  # generate label for that comparison

        # Conservative peak set: the true replicate pair with most peaks
        conservative_set = label_tr
        conservative_peak = args.peaks[Nt_idx]
        N_conservative = Nt

        # Optimal peak set is the one with more peaks (between conservative and pooled-pseudo)
        if Nt > Np:
            optimal_set = conservative_set
            optimal_peak = conservative_peak
            N_optimal = N_conservative
        else:
            optimal_set = "pooled-pr1_vs_pooled-pr2"
            optimal_peak = args.peak_ppr
            N_optimal = Np
    else:
        # Single replicate fallback mode
        num_rep = 1
        Nt = 0
        Np = 0
        rescue_ratio = 0.0
        self_consistency_ratio = 1.0

        # Use the only pseudo-replicate peak file
        conservative_set = 'rep1-pr1_vs_rep1-pr2'
        conservative_peak = args.peaks_pr[0]
        N_conservative = N[0]
        optimal_set = conservative_set
        optimal_peak = conservative_peak
        N_optimal = N_conservative

    # Determine reproducibility status based on IDR criteria
    reproducibility = 'pass'
    if rescue_ratio > 2.0 or self_consistency_ratio > 2.0:
        reproducibility = 'borderline'
    if rescue_ratio > 2.0 and self_consistency_ratio > 2.0:
        reproducibility = 'fail'

    # Copy peak files to output with appropriate naming
    log.info('Writing optimal/conservative peak files...')
    optimal_peak_file = os.path.join(
        args.out_dir, '{}optimal_peak.{}.gz'.format(
            (args.prefix + '.') if args.prefix else '',
            args.peak_type))
    conservative_peak_file = os.path.join(
        args.out_dir, '{}conservative_peak.{}.gz'.format(
            (args.prefix + '.') if args.prefix else '',
            args.peak_type))
    copy_f_to_f(optimal_peak, optimal_peak_file)
    copy_f_to_f(conservative_peak, conservative_peak_file)

    # Convert peaks to multiple formats (bigBed, starch, hammock) if chrsz is provided
    if args.chrsz:
        log.info('Converting peak to bigbed...')
        peak_to_bigbed(optimal_peak_file, args.peak_type, args.chrsz, args.mem_gb, args.out_dir)
        peak_to_bigbed(conservative_peak_file, args.peak_type, args.chrsz, args.mem_gb, args.out_dir)

        log.info('Converting peak to starch...')
        peak_to_starch(optimal_peak_file, args.out_dir)
        peak_to_starch(conservative_peak_file, args.out_dir)

        log.info('Converting peak to hammock...')
        peak_to_hammock(optimal_peak_file, args.mem_gb, args.out_dir)
        peak_to_hammock(conservative_peak_file, args.mem_gb, args.out_dir)

    # Write reproducibility metrics to file
    log.info('Writing reproducibility QC log...')
    reproducibility_qc = '{}.reproducibility.qc'.format(args.prefix) if args.prefix else 'reproducibility.qc'
    reproducibility_qc = os.path.join(args.out_dir, reproducibility_qc)

    with open(reproducibility_qc, 'w') as fp:
        # Write header
        header = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
            'Nt',
            '\t'.join(['N{}'.format(i+1) for i in range(num_rep)]),
            'Np',
            'N_opt',
            'N_consv',
            'opt_set',
            'consv_set',
            'rescue_ratio',
            'self_consistency_ratio',
            'reproducibility',
        )
        # Write values
        line = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
            Nt,
            '\t'.join([str(i) for i in N]),
            Np,
            N_optimal,
            N_conservative,
            optimal_set,
            conservative_set,
            rescue_ratio,
            self_consistency_ratio,
            reproducibility)
        fp.write(header)
        fp.write(line)

    # Additional region-level QC
    log.info('Calculating (optimal) peak region size QC/plot...')
    region_size_qc, region_size_plot = get_region_size_metrics(optimal_peak_file)

    log.info('Calculating number of peaks (optimal)...')
    get_num_peaks(optimal_peak_file)

    log.info('All done.')


# Run the script
if __name__ == '__main__':
    main()
