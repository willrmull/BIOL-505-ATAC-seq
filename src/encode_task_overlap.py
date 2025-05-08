#!/usr/bin/env python

# Standard Python libraries
import sys
import os
import argparse

# ENCODE DCC utility functions
from encode_lib_common import (
    assert_file_not_empty,
    gunzip,
    log,
    ls_l,
    mkdir_p,
    rm_f,
    run_shell_cmd,
    get_gnu_sort_param,
)

# Genomics-related functions for conversions and filtering
from encode_lib_genomic import (
    peak_to_bigbed,
    peak_to_hammock,
    peak_to_starch,
)

from encode_lib_blacklist_filter import blacklist_filter
from encode_lib_frip import frip, frip_shifted


# Argument parsing function
def parse_arguments():
    parser = argparse.ArgumentParser(
        prog='ENCODE DCC Naive overlap.',
        description='Computes naive overlap between two replicate peak files and a pooled peak file.'
    )

    # Required peak files
    parser.add_argument('peak1', type=str, help='Replicate 1 peak file')
    parser.add_argument('peak2', type=str, help='Replicate 2 peak file')
    parser.add_argument('peak_pooled', type=str, help='Pooled replicate peak file')

    # Optional parameters
    parser.add_argument('--prefix', default='overlap', type=str, help='Basename prefix for output file')
    parser.add_argument('--peak-type', type=str, required=True,
                        choices=['narrowPeak', 'regionPeak', 'broadPeak', 'gappedPeak'],
                        help='Peak file format type')
    parser.add_argument('--nonamecheck', action='store_true',
                        help='Use if bedtools complains about name mismatches')
    parser.add_argument('--blacklist', type=str, help='Blacklist BED file')
    parser.add_argument('--regex-bfilt-peak-chr-name',
                        help='Regex to retain matching chromosome names in blacklist-filtered peaks')
    parser.add_argument('--ta', type=str, help='TAGALIGN file for FRiP calculation')
    parser.add_argument('--chrsz', type=str, help='Chromosome sizes file')
    parser.add_argument('--fraglen', type=int, default=0, help='Fragment length (for shifted FRiP)')
    parser.add_argument('--mem-gb', type=float, default=4.0, help='Memory in GB for sorting')
    parser.add_argument('--out-dir', default='', type=str, help='Output directory')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'CRITICAL', 'ERROR'],
                        help='Log verbosity level')

    args = parser.parse_args()

    # Set empty blacklist to blank if 'null' is provided
    if args.blacklist is None or args.blacklist.endswith('null'):
        args.blacklist = ''

    # Set up logging
    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args


# Main naive overlap logic
def naive_overlap(basename_prefix, peak1, peak2, peak_pooled, peak_type,
                  nonamecheck, mem_gb, out_dir):
    # Prepare output filename
    prefix = os.path.join(out_dir, basename_prefix) + '.overlap'
    overlap_peak = '{}.{}.gz'.format(prefix, peak_type)

    # Add bedtools option if nonamecheck is enabled
    nonamecheck_param = '-nonamecheck' if nonamecheck else ''

    # Set awk logic and cut fields depending on peak format
    if peak_type.lower() in ('narrowpeak', 'regionpeak'):
        awk_param = '{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}'
        cut_param = '1-10'
    elif peak_type.lower() == 'broadpeak':
        awk_param = '{s1=$3-$2; s2=$12-$11; if (($19/s1 >= 0.5) || ($19/s2 >= 0.5)) {print $0}}'
        cut_param = '1-9'
    elif peak_type.lower() == 'gappedpeak':
        awk_param = '{s1=$3-$2; s2=$18-$17; if (($31/s1 >= 0.5) || ($31/s2 >= 0.5)) {print $0}}'
        cut_param = '1-15'
    else:
        raise ValueError('Unsupported peak_type.')

    # Uncompress gzipped peaks (bedtools doesn't support gzipped input here)
    tmp1 = gunzip(peak1, 'tmp1', out_dir)
    tmp2 = gunzip(peak2, 'tmp2', out_dir)
    tmp_pooled = gunzip(peak_pooled, 'tmp_pooled', out_dir)

    # Intersect pooled vs replicate1 -> filter by 50% overlap -> intersect with replicate2 -> filter -> save
    run_shell_cmd(
        'intersectBed {nonamecheck_param} -wo '
        '-a {tmp_pooled} -b {tmp1} | '
        'awk \'BEGIN{{FS="\\t";OFS="\\t"}} {awk_param}\' | '
        'cut -f {cut_param} | sort {sort_param} | uniq | '
        'intersectBed {nonamecheck_param} -wo '
        '-a stdin -b {tmp2} | '
        'awk \'BEGIN{{FS="\\t";OFS="\\t"}} {awk_param}\' | '
        'cut -f {cut_param} | sort {sort_param} | uniq | gzip -nc > {overlap_peak}'.format(
            nonamecheck_param=nonamecheck_param,
            tmp_pooled=tmp_pooled,
            tmp1=tmp1,
            awk_param=awk_param,
            cut_param=cut_param,
            sort_param=get_gnu_sort_param(mem_gb * 1024 ** 3, ratio=0.5),
            tmp2=tmp2,
            overlap_peak=overlap_peak,
        )
    )

    # Clean up temp files
    rm_f([tmp1, tmp2, tmp_pooled])

    return overlap_peak


def main():
    # Parse CLI args
    args = parse_arguments()

    # Set up output directory
    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # Compute naive overlap between peaks
    log.info('Do naive overlap...')
    overlap_peak = naive_overlap(
        args.prefix, args.peak1, args.peak2, args.peak_pooled,
        args.peak_type, args.nonamecheck, args.mem_gb, args.out_dir,
    )

    # Filter overlapping peaks against blacklist
    log.info('Blacklist-filtering peaks...')
    bfilt_overlap_peak = blacklist_filter(
        overlap_peak, args.blacklist, args.regex_bfilt_peak_chr_name, args.out_dir)

    # Confirm the filtered result is not empty
    log.info('Checking if output is empty...')
    assert_file_not_empty(bfilt_overlap_peak)

    # Convert filtered peak to various formats
    log.info('Converting peak to bigbed...')
    peak_to_bigbed(bfilt_overlap_peak, args.peak_type, args.chrsz, args.mem_gb, args.out_dir)

    log.info('Converting peak to starch...')
    peak_to_starch(bfilt_overlap_peak, args.out_dir)

    log.info('Converting peak to hammock...')
    peak_to_hammock(bfilt_overlap_peak, args.mem_gb, args.out_dir)

    # Optionally compute FRiP score
    if args.ta:
        if args.fraglen:
            log.info('Shifted FRiP with fragment length...')
            frip_shifted(args.ta, bfilt_overlap_peak, args.chrsz, args.fraglen, args.out_dir)
        else:
            log.info('FRiP without fragment length...')
            frip(args.ta, bfilt_overlap_peak, args.out_dir)

    # List output contents
    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')


# Run main if executed as script
if __name__ == '__main__':
    main()
