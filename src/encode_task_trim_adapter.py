#!/usr/bin/env python

# Standard libraries
import sys
import os
import argparse
import copy

# Local utility to detect adapter sequences
from detect_adapter import detect_most_likely_adapter

# Common utility functions
from encode_lib_common import (
    copy_f_to_dir, log, ls_l, mkdir_p, read_tsv, rm_f,
    run_shell_cmd, strip_ext_fastq
)

# Function to merge multiple FASTQ files into one
from encode_task_merge_fastq import merge_fastqs


def parse_arguments(debug=False):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(prog='ENCODE DCC adapter trimmer.',
                                     description='Trims sequencing adapters from FASTQ files.')

    # Input: either a list of fastq.gz files or a TSV indicating groups to be merged
    parser.add_argument('fastqs', nargs='+', type=str,
                        help='List of FASTQ files or TSV file specifying sets of FASTQs to be merged.')

    # Option to auto-detect adapters (Illumina/Nextera/smallRNA)
    parser.add_argument('--auto-detect-adapter', action='store_true',
                        help='Automatically detect and trim adapters.')

    # Parameters passed to `cutadapt` for trimming
    parser.add_argument('--cutadapt-param', type=str, default='-e 0.1 -m 5',
                        help='cutadapt parameters (default: error rate 0.1, min length 5).')

    # One global adapter sequence to override all individual adapter inputs
    parser.add_argument('--adapter', type=str,
                        help='Global adapter to apply to all files.')

    # Per-file or TSV-specified adapter sequences
    parser.add_argument('--adapters', nargs='+', type=str,
                        help='List or TSV of adapter strings to use (per FASTQ).')

    # Paired-end flag
    parser.add_argument('--paired-end', action="store_true",
                        help='Enable paired-end mode.')

    # Number of threads
    parser.add_argument('--nth', type=int, default=1,
                        help='Number of threads (currently unused).')

    # Output directory
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')

    # Logging verbosity
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'CRITICAL', 'ERROR'],
                        help='Log verbosity level.')

    args = parser.parse_args()

    # If input is a list of FASTQs, wrap each in a list; otherwise read TSV as a matrix
    if args.fastqs[0].endswith('.gz') or args.fastqs[0].endswith('.fastq') or args.fastqs[0].endswith('.fq'):
        args.fastqs = [[f] for f in args.fastqs]
    else:
        args.fastqs = read_tsv(args.fastqs[0])  # Read TSV as matrix [merge_id][R1/R2]

    # If adapters are specified, parse as list or TSV
    if args.adapters:
        if os.path.exists(args.adapters[0]):
            args.adapters = read_tsv(args.adapters[0])  # TSV: matrix form
        else:
            args.adapters = [[a] for a in args.adapters]  # Wrap strings into matrix

    # Fill missing adapter values if not provided
    if args.adapter or not args.adapters:
        args.adapters = copy.deepcopy(args.fastqs)
        for i, adapters in enumerate(args.adapters):
            for j, adapter in enumerate(adapters):
                args.adapters[i][j] = args.adapter if args.adapter else ''

    # Check matrix dimensions match (fastqs vs adapters)
    if len(args.adapters) != len(args.fastqs):
        raise argparse.ArgumentTypeError('Mismatch in fastqs vs adapters count.')
    for i, fastqs in enumerate(args.fastqs):
        if args.paired_end and len(fastqs) != 2:
            raise argparse.ArgumentTypeError('Expecting 2 fastqs per group in paired-end mode.')
        if not args.paired_end and len(fastqs) != 1:
            raise argparse.ArgumentTypeError('Expecting 1 fastq per group in single-end mode.')
        if len(fastqs) != len(args.adapters[i]):
            raise argparse.ArgumentTypeError('Mismatch in fastqs and adapters for group {}.'.format(i))

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args


def trim_adapter_se(fastq, adapter, adapter_for_all, cutadapt_param, out_dir):
    """Trim adapter from a single-end FASTQ."""
    if adapter:
        prefix = os.path.join(out_dir, os.path.basename(strip_ext_fastq(fastq)))
        trimmed = '{}.trim.fastq.gz'.format(prefix)
        cmd = 'cutadapt {} -a {} {} | gzip -nc > {}'.format(
            cutadapt_param,
            adapter_for_all if adapter_for_all else adapter,
            fastq,
            trimmed)
        run_shell_cmd(cmd)
        return trimmed
    else:
        return copy_f_to_dir(fastq, out_dir)


def trim_adapter_pe(fastq1, fastq2, adapter1, adapter2, adapter_for_all,
                    cutadapt_param, out_dir):
    """Trim adapters from paired-end FASTQs."""
    if adapter1 and adapter2:
        prefix1 = os.path.join(out_dir, os.path.basename(strip_ext_fastq(fastq1)))
        prefix2 = os.path.join(out_dir, os.path.basename(strip_ext_fastq(fastq2)))
        trimmed1 = '{}.trim.fastq.gz'.format(prefix1)
        trimmed2 = '{}.trim.fastq.gz'.format(prefix2)

        cmd = 'cutadapt {} -a {} -A {} {} {} -o {} -p {}'.format(
            cutadapt_param,
            adapter_for_all if adapter_for_all else adapter1,
            adapter_for_all if adapter_for_all else adapter2,
            fastq1, fastq2,
            trimmed1, trimmed2)
        run_shell_cmd(cmd)
        return [trimmed1, trimmed2]
    else:
        # If no adapter given, just copy files
        fq1 = copy_f_to_dir(fastq1, out_dir)
        fq2 = copy_f_to_dir(fastq2, out_dir)
        return [fq1, fq2]


def main():
    # Parse and validate arguments
    args = parse_arguments()

    log.info('Initializing output directory...')
    mkdir_p(args.out_dir)

    temp_files = []  # List of temporary files to clean up

    # Adapter auto-detection step
    log.info('Detecting adapters...')
    for i in range(len(args.fastqs)):
        log.info('Detecting adapters for merge_id={}...'.format(i+1))
        fastqs = args.fastqs[i]
        adapters = args.adapters[i]

        if args.paired_end:
            if not args.adapter and args.auto_detect_adapter and not (adapters[0] and adapters[1]):
                args.adapters[i][0] = detect_most_likely_adapter(fastqs[0])
                args.adapters[i][1] = detect_most_likely_adapter(fastqs[1])
                log.info('Detected adapters for R1: {}, R2: {}'.format(
                         args.adapters[i][0], args.adapters[i][1]))
        else:
            if not args.adapter and args.auto_detect_adapter and not adapters[0]:
                args.adapters[i][0] = detect_most_likely_adapter(fastqs[0])
                log.info('Detected adapter for R1: {}'.format(args.adapters[i][0]))

    # Trimming step
    log.info('Trimming adapters...')
    trimmed_fastqs_R1 = []
    trimmed_fastqs_R2 = []
    for i in range(len(args.fastqs)):
        fastqs = args.fastqs[i]
        adapters = args.adapters[i]

        if args.paired_end:
            trimmed = trim_adapter_pe(
                fastqs[0], fastqs[1],
                adapters[0], adapters[1],
                args.adapter,
                args.cutadapt_param,
                args.out_dir)
            trimmed_fastqs_R1.append(trimmed[0])
            trimmed_fastqs_R2.append(trimmed[1])
        else:
            trimmed = trim_adapter_se(
                fastqs[0],
                adapters[0],
                args.adapter,
                args.cutadapt_param,
                args.out_dir)
            trimmed_fastqs_R1.append(trimmed)

    # Merge trimmed reads
    log.info('Merging fastqs...')
    log.info('Merging R1: {}'.format(trimmed_fastqs_R1))
    merge_fastqs(trimmed_fastqs_R1, 'R1', args.out_dir)
    if args.paired_end:
        log.info('Merging R2: {}'.format(trimmed_fastqs_R2))
        merge_fastqs(trimmed_fastqs_R2, 'R2', args.out_dir)

    # Record trimmed FASTQs for cleanup
    temp_files.extend(trimmed_fastqs_R1)
    temp_files.extend(trimmed_fastqs_R2)

    # Remove temp files
    log.info('Removing temporary files...')
    rm_f(temp_files)

    # Final listing of output directory contents
    log.info('Listing output directory...')
    ls_l(args.out_dir)

    log.info('All done.')


# Entry point
if __name__ == '__main__':
    main()
