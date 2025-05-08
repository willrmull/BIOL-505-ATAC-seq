#!/usr/bin/env python

# Import necessary modules for system operations, argument parsing, and the custom functions
import sys
import os
import argparse
from encode_lib_common import (
    assert_file_not_empty, get_num_lines, log, ls_l, mkdir_p, rm_f,
    run_shell_cmd, strip_ext_ta)
    

def parse_arguments():
    """
    Parse command-line arguments using argparse to provide flexibility in how the script is run.
    """
    parser = argparse.ArgumentParser(prog='ENCODE DCC pseudo replicator.')

    # Required argument: path to the TAGALIGN file
    parser.add_argument('ta', type=str, help='Path for TAGALIGN file.')

    # Optional argument: whether the data is paired-end (default is single-end)
    parser.add_argument('--paired-end', action="store_true", help='Paired-end TAGALIGN.')

    # Optional argument: random seed for pseudoreplication
    parser.add_argument('--pseudoreplication-random-seed', type=int, default=0,
                        help='Set it to 0 to use file\'s size (in bytes) as random seed.'
                             'Otherwise, this seed will be used for GNU shuf --random-source=sha256(seed).'
                             'Useful when random seed based on input file size doesn\'t work.')

    # Optional argument: output directory to store results
    parser.add_argument('--out-dir', default='', type=str, help='Output directory.')

    # Optional argument: log level to control verbosity of logs
    parser.add_argument('--log-level', default='INFO', choices=['NOTSET', 'DEBUG', 'INFO', 
                                                                'WARNING', 'CRITICAL', 'ERROR', 'CRITICAL'],
                        help='Log level')

    # Parse the arguments
    args = parser.parse_args()

    # Set logging level based on the user input
    log.setLevel(args.log_level)

    # Log the full command used to invoke the script
    log.info(sys.argv)
    
    return args


def spr_se(ta, pseudoreplication_random_seed, out_dir):
    """
    Perform single-end pseudoreplication on the TAGALIGN file.

    This function splits the input TAGALIGN file into two pseudo-replicates by shuffling the lines.
    It uses a random seed (either based on file size or provided by the user) to ensure reproducibility.

    Args:
    - ta: Path to the TAGALIGN file.
    - pseudoreplication_random_seed: Random seed for pseudoreplication (0 uses file size as seed).
    - out_dir: Directory where output files will be saved.

    Returns:
    - ta_pr1, ta_pr2: Paths to the two pseudo-replicates created by the function.
    """
    # Generate a prefix for output files based on the input file's name
    prefix = os.path.join(out_dir, os.path.basename(strip_ext_ta(ta)))

    # Temporary files for storing shuffled lines
    tmp_pr1 = '{}.00'.format(prefix)
    tmp_pr2 = '{}.01'.format(prefix)

    # Final output pseudo-replicates (gzipped TAGALIGN files)
    ta_pr1 = '{}.pr1.tagAlign.gz'.format(prefix)
    ta_pr2 = '{}.pr2.tagAlign.gz'.format(prefix)

    # Calculate the number of lines in the input file and divide by 2 for splitting
    nlines = int((get_num_lines(ta)+1)/2)

    # If no specific random seed is provided, use the file size as the seed
    if pseudoreplication_random_seed == 0:
        random_seed = run_shell_cmd('zcat -f {ta} | wc -c'.format(ta=ta))
        log.info('Using input file\'s size {random_seed} as random seed for pseudoreplication.'.format(random_seed=random_seed))
    else:
        random_seed = pseudoreplication_random_seed
        log.info('Using a fixed integer {random_seed} as random seed for pseudoreplication.'.format(random_seed=random_seed))

    # Shell command to shuffle the lines and split the file into two parts
    run_shell_cmd(
        'zcat {ta} | shuf --random-source=<(openssl enc '
        '-aes-256-ctr -pass pass:{random_seed} '
        '-nosalt </dev/zero 2>/dev/null) | '
        'split -d -l {nlines} - {prefix}.'.format(
            ta=ta, random_seed=random_seed, nlines=nlines, prefix=prefix)
    )

    # Gzip the resulting shuffled files into two pseudo-replicates
    run_shell_cmd('gzip -nc {tmp_pr1} > {ta_pr1}'.format(tmp_pr1=tmp_pr1, ta_pr1=ta_pr1))
    run_shell_cmd('gzip -nc {tmp_pr2} > {ta_pr2}'.format(tmp_pr2=tmp_pr2, ta_pr2=ta_pr2))

    # Remove temporary files after processing
    rm_f([tmp_pr1, tmp_pr2])

    # Return the paths to the two pseudo-replicates
    return ta_pr1, ta_pr2


def spr_pe(ta, pseudoreplication_random_seed, out_dir):
    """
    Perform paired-end pseudoreplication on the TAGALIGN file.

    This function performs pseudoreplication by shuffling paired-end reads. It splits the input TAGALIGN
    file into two pseudo-replicates (pr1 and pr2), while ensuring that the pairs are maintained.

    Args:
    - ta: Path to the TAGALIGN file.
    - pseudoreplication_random_seed: Random seed for pseudoreplication (0 uses file size as seed).
    - out_dir: Directory where output files will be saved.

    Returns:
    - ta_pr1, ta_pr2: Paths to the two pseudo-replicates created by the function.
    """
    # Generate a prefix for output files based on the input file's name
    prefix = os.path.join(out_dir, os.path.basename(strip_ext_ta(ta)))

    # Temporary files for storing shuffled lines
    tmp_pr1 = '{}.00'.format(prefix)
    tmp_pr2 = '{}.01'.format(prefix)

    # Final output pseudo-replicates (gzipped TAGALIGN files)
    ta_pr1 = '{}.pr1.tagAlign.gz'.format(prefix)
    ta_pr2 = '{}.pr2.tagAlign.gz'.format(prefix)

    # Calculate the number of lines for each pseudo-replicate
    nlines = int((get_num_lines(ta)/2+1)/2)

    # If no specific random seed is provided, use the file size as the seed
    if pseudoreplication_random_seed == 0:
        random_seed = run_shell_cmd('zcat -f {ta} | wc -c'.format(ta=ta))
        log.info('Using input file\'s size {random_seed} as random seed for pseudoreplication.'.format(random_seed=random_seed))
    else:
        random_seed = pseudoreplication_random_seed
        log.info('Using a fixed integer {random_seed} as random seed for pseudoreplication.'.format(random_seed=random_seed))

    # Shell command to shuffle the lines and split the file into two parts, handling paired-end data
    run_shell_cmd(
        'zcat -f {ta} | sed \'N;s/\\n/\\t/\' | '
        'shuf --random-source=<(openssl enc -aes-256-ctr '
        '-pass pass:{random_seed} -nosalt </dev/zero 2>/dev/null) | '
        'split -d -l {nlines} - {prefix}.'.format(
            ta=ta, random_seed=random_seed, nlines=nlines, prefix=prefix)
    )

    # Process the shuffled pairs and save them as gzipped files
    run_shell_cmd(
        'zcat -f {tmp_pr1} | '
        'awk \'BEGIN{{OFS="\\t"}} '
        '{{printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n'
        '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n",'
        '$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}\' | '
        'gzip -nc > {ta_pr1}'.format(
            tmp_pr1=tmp_pr1, ta_pr1=ta_pr1)
    )

    run_shell_cmd(
        'zcat -f {tmp_pr2} | '
        'awk \'BEGIN{{OFS="\\t"}} '
        '{{printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n'
        '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n",'
        '$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}\' | '
        'gzip -nc > {ta_pr2}'.format(
            tmp_pr2=tmp_pr2, ta_pr2=ta_pr2)
    )

    # Remove temporary files after processing
    rm_f([tmp_pr1, tmp_pr2])

    # Return the paths to the two pseudo-replicates
    return ta_pr1, ta_pr2


def main():
    """
    Main function that coordinates the pseudoreplication process based on whether the data is paired-end or single-end.
    """
    # Step 1: Parse command-line arguments
    args = parse_arguments()

    # Step 2: Initialize and create the output directory if it doesn't exist
    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # Step 3: Perform pseudoreplication
    log.info('Making self-pseudo replicates...')
    if args.paired_end:
        # For paired-end data, call spr_pe
        ta_pr1, ta_pr2 = spr_pe(args.ta, args.pseudoreplication_random_seed, args.out_dir)
    else:
        # For single-end data, call spr_se
        ta_pr1, ta_pr2 = spr_se(args.ta, args.pseudoreplication_random_seed, args.out_dir)

    # Step 4: List all files in the output directory
    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    # Step 5: Check if the output files are empty
    log.info('Checking if output is empty...')
    assert_file_not_empty(ta_pr1)
    assert_file_not_empty(ta_pr2)

    # Step 6: Log completion message
    log.info('All done.')


if __name__ == '__main__':
    # Run the main function if the script is executed directly
    main()
