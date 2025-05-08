#!/usr/bin/env python

# Standard libraries
import sys
import os
import argparse

# ENCODE DCC helper libraries
from encode_lib_common import (
    assert_file_not_empty, human_readable_number,
    log, ls_l, mkdir_p, rm_f, run_shell_cmd, strip_ext_ta,
    get_gnu_sort_param,
)
from encode_lib_genomic import bed_clip


# Function to parse command-line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(
        prog='ENCODE MACS2 callpeak',
        description='Wrapper for calling peaks using MACS2 with ENCODE-style filtering and formatting.'
    )

    # Required input tagAlign file
    parser.add_argument('ta', type=str,
                        help='Path for TAGALIGN file.')

    # Optional input files and settings
    parser.add_argument('--chrsz', type=str,
                        help='2-column chromosome sizes file.')
    parser.add_argument('--gensz', type=str,
                        help='Genome size (e.g., hs for human, ms for mouse, or total base pairs).')
    parser.add_argument('--pval-thresh', default=0.01, type=float,
                        help='P-value threshold for peak calling.')
    parser.add_argument('--smooth-win', default=150, type=int,
                        help='Smoothing window size for shifting reads (extsize).')
    parser.add_argument('--cap-num-peak', default=300000, type=int,
                        help='Maximum number of peaks to retain (top N by score).')
    parser.add_argument('--mem-gb', type=float, default=4.0,
                        help='Max memory (GB) to allocate for GNU sort.')
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory for results.')

    # Logging level
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'CRITICAL', 'ERROR'],
                        help='Set the logging level.')

    args = parser.parse_args()
    log.setLevel(args.log_level)  # Set the logger level
    log.info(sys.argv)  # Log the command-line arguments
    return args


# Function to call peaks using MACS2
def macs2(ta, chrsz, gensz, pval_thresh, smooth_win, cap_num_peak,
          mem_gb, out_dir):

    # Generate output file name prefix based on tagAlign filename
    prefix = os.path.join(out_dir, os.path.basename(strip_ext_ta(ta)))

    # Form final output file name (gzipped narrowPeak)
    npeak = '{}.{}.{}.narrowPeak.gz'.format(
        prefix,
        'pval{}'.format(pval_thresh),
        human_readable_number(cap_num_peak)
    )

    # Define temporary filenames used during processing
    npeak_tmp = '{}.tmp'.format(npeak)
    npeak_tmp2 = '{}.tmp2'.format(npeak)
    shiftsize = -int(round(float(smooth_win) / 2.0))  # Compute MACS2 shift size
    temp_files = []

    # Run MACS2 callpeak with the provided arguments
    run_shell_cmd(
        'macs2 callpeak '
        '-t {ta} -f BED -n {prefix} -g {gensz} -p {pval_thresh} '
        '--shift {shiftsize} --extsize {extsize} '
        '--nomodel -B --SPMR --keep-dup all --call-summits'.format(
            ta=ta,
            prefix=prefix,
            gensz=gensz,
            pval_thresh=pval_thresh,
            shiftsize=shiftsize,
            extsize=smooth_win
        )
    )

    # Sort peaks by score (column 8), rename them (Peak_1, Peak_2, ...), and fix boundaries
    run_shell_cmd(
        'LC_COLLATE=C sort -k 8gr,8gr {sort_param} "{prefix}_peaks.narrowPeak" | '
        'awk \'BEGIN{{OFS="\\t"}}'
        '{{$4="Peak_"NR; if ($2<0) $2=0; if ($3<0) $3=0; if ($10==-1) '
        '$10=$2+int(($3-$2+1)/2.0); print $0}}\' > {npeak_tmp}'.format(
            sort_param=get_gnu_sort_param(mem_gb * 1024 ** 3, ratio=0.5),
            prefix=prefix,
            npeak_tmp=npeak_tmp,
        )
    )

    # Retain only top N peaks (based on cap_num_peak)
    run_shell_cmd(
        'head -n {cap_num_peak} {npeak_tmp} > {npeak_tmp2}'.format(
            cap_num_peak=cap_num_peak,
            npeak_tmp=npeak_tmp,
            npeak_tmp2=npeak_tmp2,
        )
    )

    # Ensure peaks do not exceed chromosome boundaries
    bed_clip(npeak_tmp2, chrsz, npeak)

    # Remove temporary files
    rm_f([npeak_tmp, npeak_tmp2])
    temp_files.append("{prefix}_*".format(prefix=prefix))
    rm_f(temp_files)

    return npeak


# Main workflow
def main():
    # Parse user inputs
    args = parse_arguments()

    # Create output directory if it doesn't exist
    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # Call peaks using MACS2
    log.info('Calling peaks with MACS2...')
    npeak = macs2(
        args.ta,
        args.chrsz,
        args.gensz,
        args.pval_thresh,
        args.smooth_win,
        args.cap_num_peak,
        args.mem_gb,
        args.out_dir,
    )

    # Ensure peak output file is not empty
    log.info('Checking if output is empty...')
    assert_file_not_empty(npeak)

    # List contents of the output directory
    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')


# Run main only if script is executed directly
if __name__ == '__main__':
    main()
