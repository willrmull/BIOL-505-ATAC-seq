#!/usr/bin/env python

# Import required libraries
import math
import os
import gzip
import re
import subprocess

# Import custom utility functions from encode_lib_common module
from encode_lib_common import (
    get_gnu_sort_param, get_num_lines, get_peak_type, human_readable_number,
    rm_f, run_shell_cmd, strip_ext, strip_ext_bam, strip_ext_peak,
    strip_ext_ta, strip_ext_gz
)

# Default max memory for samtools per thread in MB
DEFAULT_SAMTOOLS_MAX_MEM_MB_PER_THREAD = 768

# Removes specific chromosomes from a BAM file using samtools and a filtered chrsz file
def remove_chrs_from_bam(bam, chrs, chrsz, nth=1, out_dir=''):
    if len(chrs) == 0:
        raise ValueError('There must be at least one chromosome, zero found.')

    # Create file paths for output BAM and temporary chromosome size file
    prefix = os.path.join(out_dir, os.path.basename(strip_ext_bam(bam)))
    suffix = 'no_{}'.format('_'.join(chrs))
    final_bam = '{}.{}.bam'.format(prefix, suffix)
    tmp_chrsz = '{}.{}.tmp.chrsz'.format(prefix, suffix)

    # Filter chromosomes and prepare temp chrsz file
    cmd0 = 'zcat -f {chrsz} | grep -v -P \'^({chrs})\\s\' | awk \'BEGIN{{OFS="\\t"}} {{print $1,0,$2}}\' > {tmp_chrsz}'.format(
        chrsz=chrsz,
        chrs='|'.join(chrs),
        tmp_chrsz=tmp_chrsz
    )
    run_shell_cmd(cmd0)

    # Extract the desired reads with samtools view
    cmd1 = 'samtools view -b -L {tmp_chrsz} {bam} {res_param} > {final_bam}'.format(
        tmp_chrsz=tmp_chrsz,
        bam=bam,
        res_param=get_samtools_res_param('view', nth=nth),
        final_bam=final_bam
    )
    run_shell_cmd(cmd1)
    rm_f(tmp_chrsz)  # Clean up temporary file

    return final_bam

# Run SAMstats tool on a BAM file and generate QC metrics
def samstat(bam, nth=1, mem_gb=None, out_dir=''):
    prefix = os.path.join(out_dir, os.path.basename(strip_ext_bam(bam)))
    samstat_qc = '{}.samstats.qc'.format(prefix)

    # Sort BAM by name and pass to SAMstats
    run_shell_cmd(
        'samtools sort -n {bam} -T {prefix}.tmp {res_param} -O sam | SAMstats --sorted_sam_file - --outf {samstat_qc}'.format(
            bam=bam,
            prefix=prefix,
            res_param=get_samtools_res_param('sort', nth=nth, mem_gb=mem_gb),
            samstat_qc=samstat_qc,
        )
    )
    return samstat_qc

# Create BAM index using samtools
def samtools_index(bam, nth=1, out_dir=''):
    bai = '{}.bai'.format(bam)
    run_shell_cmd(
        'samtools index {bam} {res_param}'.format(
            bam=bam,
            res_param=get_samtools_res_param('index', nth=nth),
        )
    )
    if os.path.abspath(out_dir) != os.path.abspath(os.path.dirname(bam)):
        return os.path.join(out_dir, os.path.basename(bai))
    else:
        return bai

# Generate samtools resource parameter string (-@, -m)
def get_samtools_res_param(subcmd, nth=1, mem_gb=None):
    res_param = ''
    if subcmd == 'index':
        res_param += '-@ {num_total_threads} '.format(num_total_threads=nth)
    else:
        res_param += '-@ {num_additional_threads} '.format(num_additional_threads=nth - 1)

    if subcmd == 'sort' and nth and mem_gb:
        mem_mb_per_thread = min(
            math.floor(mem_gb * 1024.0 / nth),
            DEFAULT_SAMTOOLS_MAX_MEM_MB_PER_THREAD
        )
        res_param += '-m {mem}M '.format(mem=mem_mb_per_thread)

    return res_param

# Sort BAM by position
def samtools_sort(bam, nth=1, mem_gb=None, out_dir=''):
    prefix = os.path.join(out_dir, os.path.basename(strip_ext_bam(bam)))
    srt_bam = '{}.srt.bam'.format(prefix)

    run_shell_cmd(
        'samtools sort {bam} -o {srt_bam} -T {prefix} {res_param}'.format(
            bam=bam,
            srt_bam=srt_bam,
            prefix=prefix,
            res_param=get_samtools_res_param('sort', nth=nth, mem_gb=mem_gb),
        )
    )
    return srt_bam

# Sort BAM by read name
def samtools_name_sort(bam, nth=1, mem_gb=None, out_dir=''):
    prefix = os.path.join(out_dir, os.path.basename(strip_ext_bam(bam)))
    nmsrt_bam = '{}.nmsrt.bam'.format(prefix)

    run_shell_cmd(
        'samtools sort -n {bam} -o {nmsrt_bam} -T {prefix} {res_param}'.format(
            bam=bam,
            nmsrt_bam=nmsrt_bam,
            prefix=prefix,
            res_param=get_samtools_res_param('sort', nth=nth, mem_gb=mem_gb),
        )
    )
    return nmsrt_bam

# Check if BAM is empty
def bam_is_empty(bam, nth=1):
    cmd = 'samtools view -c {bam} {res_param}'.format(
        bam=bam,
        res_param=get_samtools_res_param('view', nth=nth),
    )
    return int(run_shell_cmd(cmd)) == 0

# Locate picard.jar or conda-installed picard
def locate_picard():
    try:
        cmd = 'which picard.jar'
        return run_shell_cmd(cmd)
    except:
        try:
            cmd = 'which picard'
            picard = run_shell_cmd(cmd)
            jar_path = os.path.realpath(picard) + '.jar'
            if os.path.isfile(jar_path) and os.access(jar_path, os.R_OK):
                return jar_path
            raise Exception('Picard found, but jar file missing.')
        except:
            raise Exception('Cannot find picard.jar or conda installation of Picard tools')

# Locate trimmomatic.jar or conda-installed trimmomatic
def locate_trimmomatic():
    try:
        cmd = 'which trimmomatic.jar'
        return run_shell_cmd(cmd)
    except:
        try:
            cmd = 'which trimmomatic'
            trimmomatic = run_shell_cmd(cmd)
            jar_path = os.path.realpath(trimmomatic) + '.jar'
            if os.path.isfile(jar_path) and os.access(jar_path, os.R_OK):
                return jar_path
            raise Exception('Trimmomatic found, but jar file missing.')
        except:
            raise Exception('Cannot find trimmomatic.jar or conda installation of trimmomatic')

