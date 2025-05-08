#!/usr/bin/env python3

# This script filters multimapped reads from a SAM file based on a user-defined alignment threshold.
# It assumes the input is QNAME-sorted (i.e., all alignments for the same read appear consecutively).
# The script reads from standard input and writes to standard output.

import sys
import random
import argparse


# Function to parse command-line arguments
def parse_args():
    '''
    Parses command-line arguments for alignment threshold (-k) and paired-end flag (--paired-end).
    '''
    parser = argparse.ArgumentParser(
        description='Saves reads below an alignment threshold and discards all others.')

    # -k option: sets the alignment threshold per read (or per pair)
    parser.add_argument('-k', help='Alignment number cutoff')

    # --paired-end: flag to indicate paired-end sequencing data
    parser.add_argument('--paired-end', dest='paired_ended',
                        action='store_true', help='Data is paired-end')

    args = parser.parse_args()

    # Convert threshold to integer
    alignment_cutoff = int(args.k)
    paired_ended = args.paired_ended

    return alignment_cutoff, paired_ended


# Main execution logic
if __name__ == "__main__":
    '''
    Filters multimapped reads by discarding those with too many alignments,
    based on a cutoff. Outputs acceptable reads to stdout.
    '''

    # Parse arguments from command line
    alignment_cutoff, paired_ended = parse_args()

    # Double the cutoff for paired-end reads (because each fragment has 2 lines)
    if paired_ended:
        alignment_cutoff *= 2

    # Initialize read storage
    current_reads = []       # Stores all lines for the current read group
    current_qname = ''       # Keeps track of the current read's name

    # Read the SAM file line by line from stdin
    for line in sys.stdin:

        read_elems = line.strip().split('\t')  # Split line by tabs (SAM format)

        # Pass through header lines (start with '@') directly to output
        if read_elems[0].startswith('@'):
            sys.stdout.write(line)
            continue

        # If current read name matches the previous one, add to current group
        if read_elems[0] == current_qname:
            current_reads.append(line)
        else:
            # New read name encountered

            # If the previous read group has too many alignments, discard it
            if len(current_reads) > alignment_cutoff:
                # Start new group with current line
                current_reads = [line]
                current_qname = read_elems[0]
            elif len(current_reads) > 0:
                # Write previous read group to stdout (it passed the cutoff)
                for read in current_reads:
                    sys.stdout.write(read)

                # Start new group
                current_reads = [line]
                current_qname = read_elems[0]
            else:
                # First data line in the file, initialize the first group
                current_reads.append(line)
                current_qname = read_elems[0]

