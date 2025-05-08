#!/usr/bin/env python

# Import necessary libraries and modules
import sys  # Provides access to system-specific parameters and functions
import os  # For interacting with the operating system (file manipulation, paths, etc.)
import argparse  # For handling command-line arguments
from encode_lib_common import (
    assert_file_not_empty, log, ls_l, mkdir_p, rm_f, run_shell_cmd,
    strip_ext_bam, strip_ext_ta)  # Import utility functions for common tasks
from encode_lib_genomic import (
    samtools_name_sort, subsample_ta_pe, subsample_ta_se)  # Import genomic functions for processing BAM and TAGALIGN

# Function to parse command-line arguments
def parse_arguments():
    # Create the argument parser
    parser = argparse.ArgumentParser(prog='ENCODE DCC BAM 2 TAGALIGN.', description='')

    # Define the arguments the script will accept
    parser.add_argument('bam', type=str, help='Path for BAM file.')  # Path to the BAM file
    parser.add_argument('--disable-tn5-shift', action="store_true", help='Disable TN5 shifting for DNase-Seq.')  # Option to disable TN5 shifting
    parser.add_argument('--mito-chr-name', default='chrM', help='Mito chromosome name.')  # Default mitochondrial chromosome name
    parser.add_argument('--subsample', type=int, default=0, help='Subsample TAGALIGN. This affects all downstream analysis.')  # Option to subsample TAGALIGN
    parser.add_argument('--paired-end', action="store_true", help='Paired-end BAM')  # Flag for paired-end BAM file
    parser.add_argument('--out-dir', default='', type=str, help='Output directory.')  # Directory to save the output files
    parser.add_argument('--nth', type=int, default=1, help='Number of threads to parallelize.')  # Number of threads for parallelization
    parser.add_argument('--mem-gb', type=float, help='Max. memory for samtools sort in GB.')  # Max memory for samtools sort (in GB)
    parser.add_argument('--log-level', default='INFO', choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'CRITICAL', 'ERROR', 'CRITICAL'], help='Log level')  # Log level option

    # Parse the arguments
    args = parser.parse_args()

    # Set up logging based on the log level argument
    log.setLevel(args.log_level)
    log.info(sys.argv)  # Log the command used to run the script
    return args  # Return the parsed arguments

# Function to convert BAM to TAGALIGN for single-end data
def bam2ta_se(bam, out_dir):
    # Set the output prefix by stripping the extension from the BAM file and adding the output directory
    prefix = os.path.join(out_dir, os.path.basename(strip_ext_bam(bam)))
    # Define the TAGALIGN output file path
    ta = '{}.tagAlign.gz'.format(prefix)

    # Construct the command to convert BAM to TAGALIGN using bedtools and awk
    cmd = 'bedtools bamtobed -i {} | '
    cmd += 'awk \'BEGIN{{OFS="\\t"}}{{$4="N";$5="1000";print $0}}\' | '
    cmd += 'gzip -nc > {}'
    cmd = cmd.format(bam, ta)

    # Run the command
    run_shell_cmd(cmd)
    return ta  # Return the path to the resulting TAGALIGN file

# Function to convert BAM to TAGALIGN for paired-end data
def bam2ta_pe(bam, nth, out_dir):
    # Set the output prefix by stripping the extension from the BAM file and adding the output directory
    prefix = os.path.join(out_dir, os.path.basename(strip_ext_bam(bam)))
    # Define the TAGALIGN output file path
    ta = '{}.tagAlign.gz'.format(prefix)
    # Define an intermediate BEDPE file path
    bedpe = '{}.bedpe.gz'.format(prefix)
    # Create a name-sorted BAM file using samtools
    nmsrt_bam = samtools_name_sort(bam, nth, out_dir)

    # First command to convert BAM to BEDPE and compress it
    cmd1 = 'LC_COLLATE=C bedtools bamtobed -bedpe -mate1 -i {} | '
    cmd1 += 'gzip -nc > {}'
    cmd1 = cmd1.format(nmsrt_bam, bedpe)
    run_shell_cmd(cmd1)
    rm_f(nmsrt_bam)  # Remove the name-sorted BAM file after use

    # Second command to process the BEDPE file and convert it to TAGALIGN
    cmd2 = 'zcat -f {} | '
    cmd2 += 'awk \'BEGIN{{OFS="\\t"}}'
    cmd2 += '{{printf "%s\\t%s\\t%s\\tN\\t1000\\t%s\\n'
    cmd2 += '%s\\t%s\\t%s\\tN\\t1000\\t%s\\n",'
    cmd2 += '$1,$2,$3,$9,$4,$5,$6,$10}}\' | '
    cmd2 += 'gzip -nc > {}'
    cmd2 = cmd2.format(bedpe, ta)

    # Run the second command
    run_shell_cmd(cmd2)
    rm_f(bedpe)  # Remove the intermediate BEDPE file

    return ta  # Return the path to the resulting TAGALIGN file

# Function to perform TN5 shifting on the TAGALIGN file
def tn5_shift_ta(ta, out_dir):
    # Set the output prefix by stripping the extension from the TAGALIGN file and adding the output directory
    prefix = os.path.join(out_dir, os.path.basename(strip_ext_ta(ta)))
    # Define the shifted TAGALIGN output file path
    shifted_ta = '{}.tn5.tagAlign.gz'.format(prefix)

    # Command to perform TN5 shifting on the TAGALIGN file (based on strand)
    cmd = (
        'zcat -f {ta} | awk \'BEGIN {{OFS = "\\t"}} {{'
        'if ($6 == "+") {{$2 = $2 + 4}} else if ($6 == "-") {{$3 = $3 - 5}} '
        'if ($2 >= $3) {{ if ($6 == "+") {{$2 = $3 - 1}} else {{$3 = $2 + 1}} }} '
        'print $0}}\' | gzip -nc > {shifted_ta}'
    ).format(ta=ta, shifted_ta=shifted_ta)

    # Run the command
    run_shell_cmd(cmd)
    return shifted_ta  # Return the path to the shifted TAGALIGN file

# Main function that controls the script workflow
def main():
    # Parse the command-line arguments
    args = parse_arguments()

    # Log the initialization process and create the output directory
    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # Declare a list to store temporary files that will be deleted later
    temp_files = []  # List of files to delete at the end

    # Log the conversion of BAM to TAGALIGN
    log.info('Converting BAM to TAGALIGN...')
    # Call bam2ta_pe or bam2ta_se based on whether paired-end data is provided
    if args.paired_end:
        ta = bam2ta_pe(args.bam, args.nth, args.out_dir)
    else:
        ta = bam2ta_se(args.bam, args.out_dir)

    # If subsampling is requested, perform subsampling on the TAGALIGN file
    if args.subsample:
        log.info('Subsampling TAGALIGN...')
        if args.paired_end:
            subsampled_ta = subsample_ta_pe(
                ta, args.subsample, False, args.mito_chr_name, False, args.out_dir)
        else:
            subsampled_ta = subsample_ta_se(
                ta, args.subsample, False, args.mito_chr_name, args.out_dir)
        temp_files.append(ta)  # Add the original TAGALIGN to the temp list for deletion
    else:
        subsampled_ta = ta  # If no subsampling, use the original TAGALIGN

    # If TN5 shifting is not disabled, perform TN5 shifting on the TAGALIGN file
    if args.disable_tn5_shift:
        shifted_ta = subsampled_ta  # If shifting is disabled, no change
    else:
        log.info("TN5-shifting TAGALIGN...")
        shifted_ta = tn5_shift_ta(subsampled_ta, args.out_dir)
        temp_files.append(subsampled_ta)  # Add the subsampled TAGALIGN to the temp list for deletion

    # Check if the shifted TAGALIGN file is empty
    log.info('Checking if output is empty...')
    assert_file_not_empty(shifted_ta)

    # Remove temporary files that were created during processing
    log.info('Removing temporary files...')
    rm_f(temp_files)

    # List all files in the output directory
    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    # Log that the script has completed
    log.info('All done.')

# If this script is being run directly, execute the main function
if __name__ == '__main__':
    main()
