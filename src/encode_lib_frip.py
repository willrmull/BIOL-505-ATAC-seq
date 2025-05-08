#!/usr/bin/env python

# Import necessary libraries
import sys
import os
import argparse

# Import utility functions from a shared ENCODE library
from encode_lib_common import (
    get_num_lines, gunzip, log, ls_l, mkdir_p, rm_f,
    run_shell_cmd, strip_ext, write_txt)

# ----------- Argument Parser -----------

def parse_arguments():
    """
    Parse command-line arguments for the FRiP script.
    Returns:
        args (Namespace): Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        prog='ENCODE DCC FRiP.',
        description='Compute Fraction of Reads in Peaks (FRiP) metric.')

    # Required arguments
    parser.add_argument('peak', type=str,
                        help='Peak file (can be .gz compressed).')
    parser.add_argument('ta', type=str,
                        help='TAGALIGN file (can be .gz compressed).')

    # Optional arguments
    parser.add_argument('--chrsz', type=str,
                        help='2-column chromosome sizes file. '
                             'Required for shifted FRiP.')
    parser.add_argument('--fraglen', type=int, default=0,
                        help='Fragment length for shifting reads. '
                             'If provided, perform shifted FRiP.')
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory for result files.')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'CRITICAL', 'ERROR'],
                        help='Set the logging level.')

    args = parser.parse_args()
    log.setLevel(args.log_level)  # Set logging level from args
    log.info(sys.argv)            # Log the command-line invocation
    return args

# ----------- Basic FRiP Calculation -----------

def frip(ta, peak, out_dir):
    """
    Compute standard FRiP metric: fraction of aligned reads in peak regions.
    Args:
        ta (str): TAGALIGN file.
        peak (str): Peak file.
        out_dir (str): Output directory.
    Returns:
        str: Path to the FRiP QC output file.
    """
    # Output file prefix based on peak name
    prefix = os.path.join(out_dir, os.path.basename(strip_ext(peak)))
    frip_qc = f'{prefix}.frip.qc'

    # Handle empty peak file case
    if get_num_lines(peak) == 0:
        val1 = 0.0
        tmp_files = []
    else:
        # Uncompress inputs (bedtools bug with gzipped input on both sides)
        tmp1 = gunzip(ta, 'tmp1', out_dir)
        tmp2 = gunzip(peak, 'tmp2', out_dir)

        # Intersect reads with peaks using bedtools, then count overlaps
        cmd = f'bedtools intersect -nonamecheck -a {tmp1} -b {tmp2} -wa -u | wc -l'
        val1 = run_shell_cmd(cmd)  # Number of reads in peaks
        tmp_files = [tmp1, tmp2]

    # Total number of reads in TAGALIGN file
    val2 = get_num_lines(ta)

    # Write FRiP score to file
    write_txt(frip_qc, str(float(val1) / float(val2)))

    # Cleanup temporary files
    rm_f(tmp_files)

    return frip_qc

# ----------- Shifted FRiP Calculation -----------

def frip_shifted(ta, peak, chrsz, fraglen, out_dir):
    """
    Compute FRiP after shifting reads (for ChIP-Seq data).
    Args:
        ta (str): TAGALIGN file.
        peak (str): Peak file.
        chrsz (str): Chromosome sizes file.
        fraglen (int): Estimated fragment length.
        out_dir (str): Output directory.
    Returns:
        str: Path to the FRiP QC output file.
    """
    prefix = os.path.join(out_dir, os.path.basename(strip_ext(peak)))
    frip_qc = f'{prefix}.frip.qc'

    # Calculate half fragment length for centering reads
    half_fraglen = (fraglen + 1) / 2

    # Handle empty peak file case
    if get_num_lines(peak) == 0:
        val1 = 0.0
    else:
        tmp2 = gunzip(peak, 'tmp2', out_dir)

        # Shift reads by half fragment length, filter valid entries, and intersect with peaks
        cmd = (
            f'bedtools slop -i {ta} -g {chrsz} -s -l {-half_fraglen} -r {half_fraglen} | '
            'awk \'{if ($2>=0 && $3>=0 && $2<=$3) print $0}\' | '
            f'bedtools intersect -nonamecheck -a stdin -b {tmp2} -wa -u | wc -l'
        )
        val1 = run_shell_cmd(cmd)
        rm_f(tmp2)

    val2 = get_num_lines(ta)

    # Write FRiP score to file
    write_txt(frip_qc, str(float(val1) / float(val2)))

    return frip_qc

# ----------- Main Execution -----------

def main():
    """
    Main function: parse arguments, run FRiP computation, and log output.
    """
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # Choose FRiP method depending on whether fragment length is provided
    if args.fraglen:
        frip_shifted(args.ta, args.peak, args.chrsz, args.fraglen, args.out_dir)
    else:
        frip(args.ta, args.peak, args.out_dir)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')

# ----------- Entry Point -----------

if __name__ == '__main__':
    main()
