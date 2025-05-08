#!/usr/bin/env python

# Import necessary libraries and modules
import sys
import os
import argparse
from encode_lib_common import (
    assert_file_not_empty, human_readable_number,
    log, ls_l, mkdir_p, rm_f, run_shell_cmd, strip_ext_ta,
    get_gnu_sort_param,
)
from encode_lib_genomic import (
    subsample_ta_se, subsample_ta_pe, bed_clip,
)


def parse_arguments():
    """Parse command-line arguments."""
    # Initialize argument parser
    parser = argparse.ArgumentParser(
        prog='ENCODE DCC MACS2 callpeak')
    
    # Add positional argument for the TAGALIGN file (first) and control TAGALIGN file (second, optional)
    parser.add_argument(
        'tas', type=str, nargs='+',
        help='Path for TAGALIGN file (first) and '
             'control TAGALIGN file (second; optional).')
    
    # Add argument for fragment length (required)
    parser.add_argument('--fraglen', type=int, required=True,
                        help='Fragment length.')
    
    # Add optional argument for shift size in macs2 callpeak
    parser.add_argument('--shift', type=int, default=0,
                        help='macs2 callpeak --shift.')
    
    # Add optional argument for chromosome size file (used for peak clipping)
    parser.add_argument('--chrsz', type=str,
                        help='2-col chromosome sizes file.')
    
    # Add optional argument for genome size (either total chromosome sizes or specific abbreviations)
    parser.add_argument('--gensz', type=str,
                        help='Genome size (sum of entries in 2nd column of chr. sizes file, or hs for human, ms for mouse).')
    
    # Add argument for p-value threshold (default is 0.01)
    parser.add_argument('--pval-thresh', default=0.01, type=float,
                        help='P-Value threshold.')
    
    # Add argument for capping the number of peaks (default is 500000)
    parser.add_argument('--cap-num-peak', default=500000, type=int,
                        help='Capping number of peaks by taking top N peaks.')
    
    # Add optional argument for subsampling control file to a specific read depth
    parser.add_argument('--ctl-subsample', default=0, type=int,
                        help='Subsample control to this read depth (0: no subsampling).')
    
    # Add argument to specify if the control TAGALIGN file is paired-end
    parser.add_argument('--ctl-paired-end', action="store_true",
                        help='Paired-end control TA.')
    
    # Add argument for specifying max memory for the job (used to set sort memory)
    parser.add_argument('--mem-gb', type=float, default=4.0,
                        help='Max. memory for this job in GB. '
                             'This will be used to determine GNU sort -S (defaulting to 0.5 of this value). '
                             'It should be total memory for this task (not memory per thread).')
    
    # Add optional argument for specifying the output directory
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')
    
    # Add argument for setting log level (e.g., INFO, DEBUG)
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO',
                                 'WARNING', 'CRITICAL', 'ERROR',
                                 'CRITICAL'],
                        help='Log level')
    
    # Parse the arguments and return the parsed arguments object
    args = parser.parse_args()
    
    # If only one TAGALIGN file is provided, append an empty string for the control file
    if len(args.tas) == 1:
        args.tas.append('')
    
    # Set logging level based on the user's choice
    log.setLevel(args.log_level)
    log.info(sys.argv)  # Log the command used to run the script
    
    return args


def macs2(ta, ctl_ta, chrsz, gensz, pval_thresh, shift, fraglen, cap_num_peak,
          ctl_subsample, ctl_paired_end, mem_gb, out_dir):
    """Run MACS2 peak calling, sort peaks, and clip peaks to genome size."""
    
    # Extract the base name of the TAGALIGN file (excluding extension)
    basename_ta = os.path.basename(strip_ext_ta(ta))
    
    # If control TAGALIGN file is provided, handle subsampling if specified
    if ctl_ta:
        if ctl_subsample:
            if ctl_paired_end:
                ctl_ta = subsample_ta_pe(
                    ctl_ta, ctl_subsample,
                    non_mito=False, mito_chr_name=None, r1_only=False,
                    out_dir=out_dir)
            else:
                ctl_ta = subsample_ta_se(
                    ctl_ta, ctl_subsample,
                    non_mito=False, mito_chr_name=None,
                    out_dir=out_dir)

        # Set the base name for control TAGALIGN file
        basename_ctl_ta = os.path.basename(strip_ext_ta(ctl_ta))
        # Combine both experiment and control base names for a unique prefix
        basename_prefix = '{}_x_{}'.format(basename_ta, basename_ctl_ta)
        
        # If prefix exceeds 200 characters, truncate to avoid filesystem limits
        if len(basename_prefix) > 200:
            basename_prefix = '{}_x_control'.format(basename_ta)
    else:
        # If no control file is provided, just use the experiment base name
        basename_prefix = basename_ta
    
    # Define output prefix for all resulting files
    prefix = os.path.join(out_dir, basename_prefix)
    
    # Set the file names for the narrowPeak result and temporary files
    npeak = '{}.{}.{}.narrowPeak.gz'.format(
        prefix,
        'pval{}'.format(pval_thresh),
        human_readable_number(cap_num_peak))
    npeak_tmp = '{}.tmp'.format(npeak)
    npeak_tmp2 = '{}.tmp2'.format(npeak_tmp)
    
    # Initialize a list to store temporary files for cleanup
    temp_files = []

    # Run the MACS2 peak calling command
    run_shell_cmd(
        ' macs2 callpeak '
        '-t {ta} {ctl_param} -f BED -n {prefix} -g {gensz} -p {pval_thresh} '
        '--nomodel --shift {shiftsize} --extsize {extsize} --keep-dup all -B --SPMR'.format(
            ta=ta,
            ctl_param='-c {ctl_ta}'.format(ctl_ta=ctl_ta) if ctl_ta else '',
            prefix=prefix,
            gensz=gensz,
            pval_thresh=pval_thresh,
            shiftsize=0,
            extsize=fraglen,
        )
    )

    # Sort the peaks by score in descending order and handle peak formatting
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

    # Limit the number of peaks to the specified cap
    run_shell_cmd(
        'head -n {cap_num_peak} {npeak_tmp} > {npeak_tmp2}'.format(
            cap_num_peak=cap_num_peak,
            npeak_tmp=npeak_tmp,
            npeak_tmp2=npeak_tmp2,
        )
    )

    # Clip the peaks to the chromosome sizes to avoid out-of-bounds regions
    bed_clip(npeak_tmp2, chrsz, npeak)

    # Remove temporary files created during processing
    rm_f([npeak_tmp, npeak_tmp2])

    # Collect temporary files (e.g., the prefix-based files) for cleanup
    temp_files.append("{prefix}_*".format(prefix=prefix))
    rm_f(temp_files)

    return npeak


def main():
    """Main function to execute the MACS2 peak calling and processing."""
    
    # Parse command-line arguments
    args = parse_arguments()

    # Log the initialization process and create the output directory if it doesn't exist
    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # Log the peak calling process
    log.info('Calling peaks with macs2...')
    
    # Call the MACS2 function to perform peak calling and peak processing
    npeak = macs2(
        args.tas[0],
        args.tas[1],
        args.chrsz,
        args.gensz,
        args.pval_thresh,
        args.shift,
        args.fraglen,
        args.cap_num_peak,
        args.ctl_subsample,
        args.ctl_paired_end,
        args.mem_gb,
        args.out_dir,
    )

    # Log the output verification step to ensure the result file is not empty
    log.info('Checking if output is empty...')
    assert_file_not_empty(npeak)

    # Log the listing of all files in the output directory
    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    # Log completion of the process
    log.info('All done.')


if __name__ == '__main__':
    # Execute the main function
    main()
