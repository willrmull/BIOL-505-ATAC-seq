#!/usr/bin/env python

# Import necessary libraries
import sys
import os
import argparse
from encode_lib_common import (
    copy_f_to_f, get_num_lines, log, ls_l, mkdir_p, write_txt)  # Utility functions for file handling and I/O operations


def parse_arguments():
    """
    Parse command-line arguments to configure the script's execution.
    This includes specifying the input data (experiment and control replicates),
    output options, and control depth parameters.
    """
    # Create an argument parser for the command-line interface
    parser = argparse.ArgumentParser(
        prog='ENCODE DCC Choose control.',
        description='Choose appropriate control for each IP replicate.'
                    ' ctl_for_repN.tagAlign.gz will be generated for each '
                    'IP replicate on --out-dir. '
                    'This outputs a file with integers '
                    '(chosen control index for each replicate per line).')

    # Arguments related to input experiment and control data
    parser.add_argument('--tas', type=str, nargs='+', required=True,
                        help='List of experiment TAG-ALIGN per IP replicate.')  # IP replicate data (experiment)
    parser.add_argument('--ctl-tas', type=str, nargs='+', required=True,
                        help='List of control TAG-ALIGN per IP replicate.')  # Control replicate data
    parser.add_argument('--ta-pooled', type=str, nargs='*',
                        help='Pooled experiment TAG-ALIGN.')  # Option for pooled experiment data
    parser.add_argument('--ctl-ta-pooled', type=str, nargs='*',
                        help='Pooled control TAG-ALIGN.')  # Option for pooled control data

    # Control depth ratio and limit configurations
    parser.add_argument('--ctl-depth-ratio', type=float, required=True,
                        help='Control depth ratio (between any two controls).')  # Maximum allowable ratio between control depths
    parser.add_argument('--ctl-depth-limit', type=int, default=200000000,
                        help='Control depth limit. If read depth of chosen control exceeds this, it should be subsampled.')  # Max allowed control depth
    parser.add_argument('--exp-ctl-depth-ratio-limit', type=float, default=5.0,
                        help='Exp vs. control depth ratio limit.')  # Limit for the ratio of experiment depth to control depth

    # Flag for always using pooled control
    parser.add_argument('--always-use-pooled-ctl', action="store_true",
                        help='Always use pooled control for all IP replicates.')  # Flag for pooled control usage

    # Output file configurations
    parser.add_argument('--out-tsv-basename', default='chosen_ctl.tsv', type=str,
                        help='Output TSV basename (chosen control index per line for each replicate).')  # Output TSV file with chosen control index
    parser.add_argument('--out-tsv-subsample-basename', default='chosen_ctl_subsample.tsv', type=str,
                        help='Output TSV subsample basename (number of reads to subsample control per replicate).')  # Output TSV file for subsampling control
    parser.add_argument('--out-txt-subsample-pooled-basename', default='chosen_ctl_subsample_pooled.txt', type=str,
                        help='Output TXT subsample basename for pooled control.')  # Output TXT file for subsampling pooled control

    # Directory for saving the output files
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')  # Directory for results

    # Log level for script execution (helps with debugging and verbosity)
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'CRITICAL', 'ERROR', 'CRITICAL'],
                        help='Log level')  # Log level choices

    args = parser.parse_args()  # Parse the command-line arguments

    log.setLevel(args.log_level)  # Set the logging level based on user input
    log.info(sys.argv)  # Log the arguments used in the execution
    return args  # Return the parsed arguments


def main():
    """
    Main function to execute the script logic.
    This function selects appropriate controls based on specified depth ratios and conditions.
    It writes output files including control indices and subsampling information.
    """
    # Read parameters from command-line arguments
    args = parse_arguments()
    log.info('Initializing and making output directory...')

    # Make the output directory if it doesn't exist
    mkdir_p(args.out_dir)

    # Log reproducibility QC message
    log.info('Choosing appropriate control for each IP replicate...')
    
    # Get the number of replicates (IP experiments) and controls
    num_rep = len(args.tas)
    num_ctl = len(args.ctl_tas)

    # Calculate the number of lines (reads) in each experiment (IP) and control TAG-ALIGN files
    depths = [get_num_lines(ta) for ta in args.tas]  # Experiment depth
    depths_ctl = [get_num_lines(ctl_ta) for ctl_ta in args.ctl_tas]  # Control depth

    # Total depth for all IP and control replicates (pooled data)
    depth_rep_pooled = sum(depths)
    depth_ctl_pooled = sum(depths_ctl)

    # Create dictionaries to store the depth of each replicate and control
    depths = dict(enumerate(depths))
    depths_ctl = dict(enumerate(depths_ctl))

    # Add pooled data to the dictionaries
    depths[-1] = depth_rep_pooled
    depths_ctl[-1] = depth_ctl_pooled

    # Initialize list for chosen control indices (starting with 0 for all replicates)
    ctl_ta_idx = [0]*num_rep

    # Handling different cases based on the number of control replicates
    if num_ctl == 1:
        # If there is only one control, use it for all IP replicates
        pass
    elif args.always_use_pooled_ctl:
        # If the flag --always-use-pooled-ctl is set, use pooled control for all IP replicates
        ctl_ta_idx = [-1]*num_rep
    else:
        # If there are multiple controls, choose based on the depth ratio
        use_pooled_ctl = False  # Flag to decide if pooled control should be used

        # Check pairwise control depths and compare them against the depth ratio
        for i in range(num_ctl):
            for j in range(i+1, num_ctl):
                if depths_ctl[i]/float(depths_ctl[j]) > args.ctl_depth_ratio or \
                        depths_ctl[j]/float(depths_ctl[i]) > args.ctl_depth_ratio:
                    use_pooled_ctl = True
                    log.info(
                        'Number of reads in controls differ by a factor of {}.'
                        ' Using pooled controls.'.format(args.ctl_depth_ratio))
                    break

        # If the pooled control should be used, assign it for all replicates
        if use_pooled_ctl:
            ctl_ta_idx = [-1]*num_rep
        else:
            # Otherwise, compare individual control depths to experiment depths
            for i in range(num_rep):
                if i > num_ctl-1:
                    ctl_ta_idx[i] = -1  # Use pooled control for out-of-bound replicates
                elif depths_ctl[i] < depths[i]:
                    log.info(
                        'Fewer reads in control {} than experiment replicate '
                        '{}. Using pooled control for replicate {}.'.format(
                            i+1, i+1, i+1))
                    ctl_ta_idx[i] = -1  # Use pooled control for replicates with smaller control depth
                else:
                    ctl_ta_idx[i] = i  # Use the corresponding control for each replicate

    # Initialize subsampling information for control reads
    ctl_ta_subsample = [0] * num_rep
    ctl_ta_subsampled_pooled = 0
    
    if args.exp_ctl_depth_ratio_limit or args.ctl_depth_limit:
        # If depth ratio limit or control depth limit is specified, perform subsampling
        for rep in range(num_rep):
            chosen_ctl = ctl_ta_idx[rep]
            depth = depths[rep]
            depth_ctl = depths_ctl[chosen_ctl]
            limit = int(max(depth * args.exp_ctl_depth_ratio_limit, args.ctl_depth_limit))
            if depth_ctl > limit:
                ctl_ta_subsample[rep] = limit  # Subsample to the limit

        # Subsample pooled control if necessary
        limit = int(max(depth_rep_pooled * args.exp_ctl_depth_ratio_limit, args.ctl_depth_limit))
        if depth_ctl_pooled > limit:
            ctl_ta_subsampled_pooled = limit  # Subsample pooled control

    # Write the control index (selected for each replicate) to a TSV file
    log.info('Writing idx.txt...')
    out_txt = os.path.join(args.out_dir, args.out_tsv_basename)
    write_txt(out_txt, ctl_ta_idx)

    # Write the subsampling information for each replicate to a TSV file
    log.info('Writing subsample txt...')
    out_subsample_txt = os.path.join(args.out_dir, args.out_tsv_subsample_basename)
    write_txt(out_subsample_txt, ctl_ta_subsample)

    # Write the subsampling information for pooled control to a TXT file
    log.info('Writing subsample_pooled txt...')
    out_subsample_pooled_txt = os.path.join(args.out_dir, args.out_txt_subsample_pooled_basename)
    write_txt(out_subsample_pooled_txt, ctl_ta_subsampled_pooled)

    # List all files generated in the output directory
    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    # Log that the process is finished
    log.info('All done.')


if __name__ == '__main__':
    main()  # Execute the main function when the script is run
