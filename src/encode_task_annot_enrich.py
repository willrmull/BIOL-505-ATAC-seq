#!/usr/bin/env python

# Import necessary libraries and modules
import sys  # Provides access to system-specific parameters and functions
import os  # For interacting with the operating system (file manipulation, paths, etc.)
import argparse  # For handling command-line arguments
from encode_lib_common import (
    run_shell_cmd, strip_ext_ta,
    ls_l, get_num_lines, log)  # Import utility functions for various tasks
import warnings  # To manage warnings (such as silencing certain types)
warnings.filterwarnings("ignore")  # Ignore all warnings

# Function to parse command-line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(
        prog='ENCODE annot_enrich (fraction of reads in annotated regions)')
    
    # Define all the command-line arguments that the script accepts
    parser.add_argument(
        '--ta', type=str, help='TAG-ALIGN file (from task bam2ta).')  # Input TAG-ALIGN file
    parser.add_argument('--dnase', type=str, help='DNase definition bed file.')  # DNase region file
    parser.add_argument('--blacklist', type=str, help='Blacklist bed file.')  # Blacklist region file
    parser.add_argument('--prom', type=str,
                        help='Promoter definition bed file.')  # Promoter region file
    parser.add_argument('--enh', type=str,
                        help='Enhancer definition bed file.')  # Enhancer region file
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')  # Directory to save output files
    parser.add_argument('--log-level', default='INFO', help='Log level',
                        choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING',
                                 'CRITICAL', 'ERROR', 'CRITICAL'])  # Logging level
    
    # Parse the arguments and set the log level
    args = parser.parse_args()
    log.setLevel(args.log_level)  # Set logging level based on input
    log.info(sys.argv)  # Log the command used to execute the script
    
    return args  # Return the parsed arguments

# Function to compute the fraction of reads in specified regions (from a BED file)
def get_fract_reads_in_regions(reads_bed, regions_bed):
    """
    Function that takes in a BED file of reads and a BED file of regions and
    calculates the fraction of reads that overlap with the regions.
    """
    # Construct shell command to find the intersection of reads with the regions
    cmd = "bedtools sort -i {}  | "
    cmd += "bedtools merge -i stdin | "
    cmd += "bedtools intersect -u -nonamecheck -a {} -b stdin | "
    cmd += "wc -l"
    cmd = cmd.format(regions_bed, reads_bed)  # Format the command with provided files
    
    # Execute the command to count the number of intersecting reads
    intersect_read_count = int(run_shell_cmd(cmd))
    # Get the total number of reads in the reads BED file
    total_read_count = get_num_lines(reads_bed)
    # Calculate the fraction of reads in the regions
    fract_reads = float(intersect_read_count) / total_read_count

    return intersect_read_count, fract_reads  # Return the counts and the fraction

# Main function to control the script's workflow
def main():
    # Parse command-line arguments
    args = parse_arguments()

    # Assign input files and output directory based on arguments
    FINAL_BED = args.ta  # TAG-ALIGN file
    OUTPUT_PREFIX = os.path.join(
        args.out_dir,
        os.path.basename(strip_ext_ta(FINAL_BED)))  # Output prefix (strip extension from TAG-ALIGN file)

    # Assign other input files (optional) based on command-line arguments
    DNASE = args.dnase if args.dnase and os.path.basename(
        args.dnase) != 'null' else ''  # DNase regions (optional)
    BLACKLIST = args.blacklist if args.blacklist and os.path.basename(
        args.blacklist) != 'null' else ''  # Blacklist regions (optional)
    PROM = args.prom if args.prom and os.path.basename(
        args.prom) != 'null' else ''  # Promoter regions (optional)
    ENH = args.enh if args.enh and os.path.basename(args.enh) != 'null' else ''  # Enhancer regions (optional)

    result = []  # List to store the results

    # Check for DNase regions and calculate fraction of reads in those regions
    if DNASE:
        reads_dnase, fract_dnase = get_fract_reads_in_regions(FINAL_BED, DNASE)
        result.append(('fraction_of_reads_in_universal_DHS_regions',
                       str(reads_dnase), str(fract_dnase)))

    # Check for Blacklist regions and calculate fraction of reads in those regions
    if BLACKLIST:
        reads_blacklist, \
            fract_black_
