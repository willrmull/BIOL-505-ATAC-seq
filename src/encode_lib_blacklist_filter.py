#!/usr/bin/env python

# This script filters out genomic regions (peaks) that overlap with blacklisted regions.
# Blacklisted regions are genomic intervals known to produce unreliable or artifactual signal.

import sys
import os
import argparse
from encode_lib_common import (
    get_ext, get_num_lines, gunzip, log, mkdir_p,
    rm_f, run_shell_cmd, strip_ext, strip_ext_bam)

# Function to parse command-line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC Blacklist filter.')
    parser.add_argument('peak', type=str, help='Peak file (BED or narrowPeak, optionally gzipped).')
    parser.add_argument('--blacklist', type=str, help='Blacklist BED file (gzipped or plain).')
    parser.add_argument('--regex-bfilt-peak-chr-name',
                        help='Only retain chromosomes matching this regex in the output.')
    parser.add_argument('--out-dir', default='', type=str, help='Directory to write filtered output.')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'CRITICAL', 'ERROR', 'CRITICAL'],
                        help='Set the logging verbosity.')
    args = parser.parse_args()

    # If blacklist is explicitly set to "null", treat it as not provided
    if args.blacklist is None or args.blacklist.endswith('null'):
        args.blacklist = ''

    # Set logging level and print initial log
    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args

# Function to filter peak file against a blacklist
def blacklist_filter(peak, blacklist, regex_bfilt_peak_chr_name, out_dir):
    # Generate output file prefix and suffix
    prefix = os.path.join(out_dir, os.path.basename(strip_ext(peak)))
    peak_ext = get_ext(peak)
    filtered = '{}.bfilt.{}.gz'.format(prefix, peak_ext)

    # Default to empty regex if none provided
    if regex_bfilt_peak_chr_name is None:
        regex_bfilt_peak_chr_name = ''

    # If blacklist or peak is missing or empty, just filter with regex (no blacklist removal)
    if blacklist == '' or get_num_lines(peak) == 0 or get_num_lines(blacklist) == 0:
        cmd = 'zcat -f {} | '.format(peak)
        cmd += 'grep -P \'{}\\b\' | '.format(regex_bfilt_peak_chr_name)
        cmd += 'gzip -nc > {}'.format(filtered)
        run_shell_cmd(cmd)
    else:
        # Unzip both peak and blacklist files for compatibility with bedtools
        tmp1 = gunzip(peak, 'tmp1', out_dir)
        tmp2 = gunzip(blacklist, 'tmp2', out_dir)

        # Remove peaks overlapping blacklist, cap score at
