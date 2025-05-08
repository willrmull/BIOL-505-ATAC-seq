#!/usr/bin/env python

# Import required libraries
import sys
import os
import argparse
import math

# Import shared ENCODE utility functions
from encode_lib_common import (
    assert_file_not_empty,  # Ensure output file is not empty
    log,                    # Logger setup
    ls_l,                   # List files in a directory
    mkdir_p,                # Make output directory
    rm_f,                   # Remove temporary file(s)
    run_shell_cmd,          # Execute shell commands
    get_gnu_sort_param      # Compute memory setting for GNU sort
)

# Import genomic-specific ENCODE utilities
from encode_lib_genomic import (
    peak_to_bigbed,         # Convert peak file to bigBed format
    peak_to_hammock,        # Convert peak file to Hammock format
    bed_clip,               # Clip peak coordinates to chromosome bounds
    peak_to_starch          # Convert peak file to Starch format
)

# Import additional ENCODE functions
from encode_lib_blacklist_filter import blacklist_filter  # Blacklist filtering
from encode_lib_frip import frip, frip_shifted            # FRiP score calculations


def parse_arguments():
    """
    Parse command-line arguments using argparse.
    """
    parser = argparse.ArgumentParser(
        prog='ENCODE DCC IDR.',
        description='Run IDR on narrowPeak or regionPeak files.'
    )

    # Required peak input files
    parser.add_argument('peak1', type=str, help='Replicate 1 peak file.')
    parser.add_argument('peak2', type=str, help='Replicate 2 peak file.')
    parser.add_argument('peak_pooled', type=str, help='Pooled peak file.')

    # Output prefix and peak type
    parser.add_argument('--prefix', default='idr', type=str, help='Output file prefix.')
    parser.add_argument('--peak-type', type=str, required=True,
                        choices=['narrowPeak', 'regionPeak', 'broadPeak', 'gappedPeak'],
                        help='Type of input peak files.')

    # IDR configuration
    parser.add_argument('--idr-thresh', default=0.1, type=float, help='IDR significance threshold.')
    parser.add_argument('--idr-rank', default='p.value', type=str,
                        choices=['p.value', 'q.value', 'signal.value'],
                        help='Column to rank peaks by.')

    # Optional enhancements
    parser.add_argument('--blacklist', type=str, help='BED file with blacklisted regions.')
    parser.add_argument('--regex-bfilt-peak-chr-name', help='Regex to filter chromosomes in peak files.')

    # FRiP score inputs
    parser.add_argument('--ta', type=str, help='TAGALIGN file for FRiP calculation.')
    parser.add_argument('--chrsz', type=str, help='Chromosome sizes file.')
    parser.add_argument('--fraglen', type=int, default=0,
                        help='Fragment length for shifted FRiP (used in ChIP-seq).')

    # Memory configuration
    parser.add_argument('--mem-gb', type=float, default=4.0,
                        help='Max memory for GNU sort (used for peak sorting).')

    # Output and logging
    parser.add_argument('--out-dir', default='', type=str, help='Output directory.')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'CRITICAL', 'ERROR'],
                        help='Logging verbosity.')

    args = parser.parse_args()

    # If blacklist file is not provided or set as 'null', ignore it
    if args.blacklist is None or args.blacklist.endswith('null'):
        args.blacklist = ''

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args


def get_npeak_col_by_rank(rank):
    """
    Return column index for ranking based on IDR ranking metric.
    """
    if rank == 'signal.value':
        return 7
    elif rank == 'p.value':
        return 8
    elif rank == 'q.value':
        return 9
    else:
        raise Exception('Invalid score ranking method')


def idr(basename_prefix, peak1, peak2, peak_pooled, peak_type, chrsz,
        thresh, rank, mem_gb, out_dir):
    """
    Main IDR pipeline using IDR framework and peak refinement.
    """
    # Construct all output filenames
    prefix = os.path.join(out_dir, basename_prefix)
    prefix += f'.idr{thresh}'
    idr_peak = f'{prefix}.{peak_type}.gz'
    idr_plot = f'{prefix}.unthresholded-peaks.txt.png'
    idr_stdout = f'{prefix}.log'
    idr_12col_bed = f'{peak_type}.12-col.bed.gz'
    idr_out = f'{prefix}.unthresholded-peaks.txt'
    idr_tmp = f'{idr_out}.tmp'
    idr_out_gz = f'{idr_out}.gz'

    # Run the IDR tool
    run_shell_cmd(
        f'idr --samples {peak1} {peak2} --peak-list {peak_pooled} --input-file-type narrowPeak '
        f'--output-file {idr_out} --rank {rank} --soft-idr-threshold {thresh} '
        f'--plot --use-best-multisummit-IDR --log-output-file {idr_stdout}'
    )

    # Clip peak coordinates that exceed chromosome sizes
    bed_clip(idr_out, chrsz, idr_tmp, no_gz=True)

    # Filter peaks above IDR threshold and sort
    col = get_npeak_col_by_rank(rank)
    neg_log10_thresh = -math.log10(thresh)
    run_shell_cmd(
        f"awk 'BEGIN{{OFS=\"\\t\"}} $12>={neg_log10_thresh} "
        f"{{if ($2<0) $2=0; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}' {idr_tmp} | "
        f"sort {get_gnu_sort_param(mem_gb * 1024**3, ratio=0.5)} | uniq | "
        f"sort -grk{col},{col} {get_gnu_sort_param(mem_gb * 1024**3, ratio=0.5)} | "
        f"gzip -nc > {idr_12col_bed}"
    )

    # Keep only first 10 columns (standard BED format)
    run_shell_cmd(
        f'zcat {idr_12col_bed} | '
        f'awk \'BEGIN{{OFS="\\t"}} {{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}\' | '
        f'gzip -nc > {idr_peak}'
    )

    # Compress the unthresholded full IDR output
    run_shell_cmd(f'cat {idr_tmp} | gzip -nc > {idr_out_gz}')

    # Clean up temporary files
    rm_f([idr_out, idr_tmp, idr_12col_bed])
    rm_f(f'{prefix}.*.noalternatesummitpeaks.png')

    return idr_peak, idr_plot, idr_out_gz, idr_stdout


def main():
    # Parse arguments and initialize
    args = parse_arguments()
    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # Run IDR and generate filtered peaks
    log.info('Running IDR...')
    idr_peak, idr_plot, idr_out_gz, idr_stdout = idr(
        args.prefix,
        args.peak1, args.peak2, args.peak_pooled, args.peak_type,
        args.chrsz,
        args.idr_thresh, args.idr_rank, args.mem_gb, args.out_dir,
    )

    # Check that resulting IDR peak file is not empty
    log.info('Checking if output is empty...')
    assert_file_not_empty(idr_peak, help=(
        'No IDR peaks found. IDR threshold might be too stringent '
        'or replicates have very poor concordance.')
    )

    # Filter out blacklisted regions from the IDR peaks
    log.info('Blacklist-filtering peaks...')
    bfilt_idr_peak = blacklist_filter(
        idr_peak, args.blacklist, args.regex_bfilt_peak_chr_name, args.out_dir)

    # Generate different formats for downstream use
    log.info('Converting peak to bigBed...')
    peak_to_bigbed(bfilt_idr_peak, args.peak_type, args.chrsz,
                   args.mem_gb, args.out_dir)

    log.info('Converting peak to Starch...')
    peak_to_starch(bfilt_idr_peak, args.out_dir)

    log.info('Converting peak to Hammock...')
    peak_to_hammock(bfilt_idr_peak, args.mem_gb, args.out_dir)

    # Optional: Calculate FRiP scores if tagAlign file is available
    if args.ta:
        if args.fraglen:
            log.info('Running shifted FRiP (ChIP-seq)...')
            frip_shifted(args.ta, bfilt_idr_peak,
                         args.chrsz, args.fraglen, args.out_dir)
        else:
            log.info('Running unshifted FRiP (ATAC-seq)...')
            frip(args.ta, bfilt_idr_peak, args.out_dir)

    # Final report of output directory
    log.info('Listing all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')


# Entrypoint
if __name__ == '__main__':
    main()
