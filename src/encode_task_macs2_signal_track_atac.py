#!/usr/bin/env python

# Import necessary libraries and modules
import sys  # Provides access to system-specific parameters and functions
import os  # Provides functions to interact with the operating system
import argparse  # For parsing command-line arguments
from encode_lib_common import (
    get_num_lines, log, ls_l, mkdir_p, rm_f, run_shell_cmd, strip_ext_ta,
    get_gnu_sort_param,
)  # Import custom utility functions from the encode_lib_common module


# Function to parse command-line arguments
def parse_arguments():
    # Initialize an argument parser with a description for the program
    parser = argparse.ArgumentParser(prog='ENCODE DCC MACS2 signal track',
                                     description='Generate MACS2 signal tracks (BigWig files) from TAGALIGN.')

    # Define arguments for the script
    parser.add_argument('ta', type=str,
                        help='Path for TAGALIGN file.')  # Input TAGALIGN file
    parser.add_argument('--chrsz', type=str,
                        help='2-col chromosome sizes file.')  # Chromosome sizes file
    parser.add_argument('--gensz', type=str,
                        help='Genome size (sum of entries in 2nd column of chr. sizes file, or hs for human, ms for mouse).')  # Genome size (either sum of chromosome sizes or a standard genome abbreviation)
    parser.add_argument('--pval-thresh', default=0.01, type=float,
                        help='P-Value threshold.')  # P-value threshold for peak calling
    parser.add_argument('--smooth-win', default=150, type=int,
                        help='Smoothing window size.')  # Size of smoothing window for signal generation
    parser.add_argument('--mem-gb', type=float, default=4.0,
                        help='Max. memory for this job in GB. This will be used to determine GNU sort -S (defaulting to 0.5 of this value).')  # Memory allocation for the job (in GB)
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')  # Output directory for the results
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'CRITICAL', 'ERROR', 'CRITICAL'],
                        help='Log level')  # Set the log level for messages (e.g., INFO, DEBUG)

    # Parse the command-line arguments and return them
    args = parser.parse_args()

    # Set logging level based on parsed argument
    log.setLevel(args.log_level)

    # Log the full command that was run (useful for debugging)
    log.info(sys.argv)

    return args


# Function to generate MACS2 signal tracks
def macs2_signal_track(ta, chrsz, gensz, pval_thresh, smooth_win, mem_gb, out_dir):
    # Define output file prefixes based on the input TAGALIGN file
    prefix = os.path.join(out_dir,
                          os.path.basename(strip_ext_ta(ta)))  # Create a prefix for output based on the input file name
    fc_bigwig = '{}.fc.signal.bigwig'.format(prefix)  # Output file for fold-change signal (BigWig format)
    pval_bigwig = '{}.pval.signal.bigwig'.format(prefix)  # Output file for p-value signal (BigWig format)

    # Temporary files for intermediate steps
    fc_bedgraph = '{}.fc.signal.bedgraph'.format(prefix)  # BedGraph file for fold-change signal
    fc_bedgraph_srt = '{}.fc.signal.srt.bedgraph'.format(prefix)  # Sorted BedGraph file for fold-change signal
    pval_bedgraph = '{}.pval.signal.bedgraph'.format(prefix)  # BedGraph file for p-value signal
    pval_bedgraph_srt = '{}.pval.signal.srt.bedgraph'.format(prefix)  # Sorted BedGraph file for p-value signal

    # Calculate the shift size for the MACS2 callpeak step
    shiftsize = -int(round(float(smooth_win)/2.0))  # Shift size for peak calling (half of smoothing window size)
    temp_files = []  # List to keep track of temporary files for cleanup

    # Run MACS2 callpeak to call peaks from the TAGALIGN file
    run_shell_cmd(
        'macs2 callpeak '
        '-t {ta} -f BED -n {prefix} -g {gensz} -p {pval_thresh} '
        '--shift {shiftsize} --extsize {extsize} '
        '--nomodel -B --SPMR '
        '--keep-dup all --call-summits '.format(
            ta=ta,
            prefix=prefix,
            gensz=gensz,
            pval_thresh=pval_thresh,
            shiftsize=shiftsize,
            extsize=smooth_win,
        )
    )

    # Use MACS2 bdgcmp to compute fold change (FE) from treatment and control pileups
    run_shell_cmd(
        'macs2 bdgcmp -t "{prefix}_treat_pileup.bdg" '
        '-c "{prefix}_control_lambda.bdg" '
        '--o-prefix "{prefix}" -m FE '.format(
            prefix=prefix,
        )
    )

    # Convert the fold-change signal from BedGraph format to BigWig format using bedtools
    run_shell_cmd(
        'bedtools slop -i "{prefix}_FE.bdg" -g {chrsz} -b 0 | '
        'bedClip stdin {chrsz} {fc_bedgraph}'.format(
            prefix=prefix,
            chrsz=chrsz,
            fc_bedgraph=fc_bedgraph,
        )
    )

    # Sort the BedGraph file and remove overlapping regions
    run_shell_cmd(
        'LC_COLLATE=C sort -k1,1 -k2,2n {sort_param} {fc_bedgraph} | '
        'awk \'BEGIN{{OFS="\\t"}}{{if (NR==1 || NR>1 && (prev_chr!=$1 '
        '|| prev_chr==$1 && prev_chr_e<=$2)) '
        '{{print $0}}; prev_chr=$1; prev_chr_e=$3;}}\' > {fc_bedgraph_srt}'.format(
            sort_param=get_gnu_sort_param(mem_gb * 1024 ** 3, ratio=0.5),
            fc_bedgraph=fc_bedgraph,
            fc_bedgraph_srt=fc_bedgraph_srt
        )
    )
    rm_f(fc_bedgraph)  # Remove the unsorted BedGraph file

    # Convert the sorted BedGraph file to BigWig format
    run_shell_cmd(
        'bedGraphToBigWig {fc_bedgraph_srt} {chrsz} {fc_bigwig}'.format(
            fc_bedgraph_srt=fc_bedgraph_srt,
            chrsz=chrsz,
            fc_bigwig=fc_bigwig,
        )
    )
    rm_f(fc_bedgraph_srt)  # Remove the sorted BedGraph file

    # Calculate the scaling value based on the number of tags per million
    sval = float(get_num_lines(ta)) / 1000000.0  # Scaling factor for p-value signal

    # Use MACS2 bdgcmp to compute p-values (ppois method)
    run_shell_cmd(
        'macs2 bdgcmp -t "{prefix}_treat_pileup.bdg" '
        '-c "{prefix}_control_lambda.bdg" '
        '--o-prefix {prefix} -m ppois -S {sval}'.format(
            prefix=prefix,
            sval=sval,
        )
    )

    # Convert the p-value signal from BedGraph format to BigWig format using bedtools
    run_shell_cmd(
        'bedtools slop -i "{prefix}_ppois.bdg" -g {chrsz} -b 0 | '
        'bedClip stdin {chrsz} {pval_bedgraph}'.format(
            prefix=prefix,
            chrsz=chrsz,
            pval_bedgraph=pval_bedgraph,
        )
    )

    # Sort the p-value BedGraph file and remove overlapping regions
    run_shell_cmd(
        'LC_COLLATE=C sort -k1,1 -k2,2n {sort_param} {pval_bedgraph} | '
        'awk \'BEGIN{{OFS="\\t"}}{{if (NR==1 || NR>1 && (prev_chr!=$1 '
        '|| prev_chr==$1 && prev_chr_e<=$2)) '
        '{{print $0}}; prev_chr=$1; prev_chr_e=$3;}}\' > {pval_bedgraph_srt}'.format(
            sort_param=get_gnu_sort_param(mem_gb * 1024 ** 3, ratio=0.5),
            pval_bedgraph=pval_bedgraph,
            pval_bedgraph_srt=pval_bedgraph_srt,
        )
    )
    rm_f(pval_bedgraph)  # Remove the unsorted p-value BedGraph file

    # Convert the sorted p-value BedGraph file to BigWig format
    run_shell_cmd(
        'bedGraphToBigWig {pval_bedgraph_srt} {chrsz} {pval_bigwig}'.format(
            pval_bedgraph_srt=pval_bedgraph_srt,
            chrsz=chrsz,
            pval_bigwig=pval_bigwig,
        )
    )
    rm_f(pval_bedgraph_srt)  # Remove the sorted p-value BedGraph file

    # Remove temporary files generated during the process
    temp_files.append("{prefix}_*".format(prefix=prefix))
    rm_f(temp_files)  # Clean up temporary files

    return fc_bigwig, pval_bigwig  # Return the paths to the generated BigWig files


# Main function to run the script
def main():
    # Parse command-line arguments
    args = parse_arguments()

    # Log and create the output directory if it does not exist
    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # Log the process of calling peaks and generating signal tracks
    log.info('Calling peaks and generating signal tracks with MACS2...')
    fc_bigwig, pval_bigwig = macs2_signal_track(
        args.ta,
        args.chrsz,
        args.gensz,
        args.pval_thresh,
        args.smooth_win,
        args.mem_gb,
        args.out_dir,
    )

    # List all files in the output directory
    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    # Log that the script has finished
    log.info('All done.')


# Check if the script is being run directly (not imported as a module)
if __name__ == '__main__':
    main()
