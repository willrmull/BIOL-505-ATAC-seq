#!/usr/bin/env python

# Imports
import sys
import os
import re
import csv
import logging
import subprocess
import math
import signal
import time
import argparse

# Configure logging format and level
logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout)
log = logging.getLogger(__name__)

BIG_INT = 99999999  # A very large integer for placeholder usage

# ---------- Utility Functions ----------

def get_ticks():
    """Return high-resolution current time in seconds for benchmarking."""
    return getattr(time, 'perf_counter', getattr(time, 'time'))()

# ----------- File Extension Stripping Utilities -----------

# Strip common extensions from file names
def strip_ext_fastq(fastq): return re.sub(r'\.(fastq|fq|Fastq|Fq)\.gz$', '', str(fastq))
def strip_ext_bam(bam): return re.sub(r'\.(bam|Bam)$', '', str(bam))
def strip_ext_tar(tar): return re.sub(r'\.(tar|tar\.gz)$', '', str(tar))
def strip_ext_ta(ta): return re.sub(r'\.(tagAlign|TagAlign|ta|Ta)\.gz$', '', str(ta))
def strip_ext_bed(bed): return re.sub(r'\.(bed|Bed)\.gz$', '', str(bed))
def strip_ext_npeak(npeak): return re.sub(r'\.(narrowPeak|NarrowPeak)\.gz$', '', str(npeak))
def strip_ext_rpeak(rpeak): return re.sub(r'\.(regionPeak|RegionPeak)\.gz$', '', str(rpeak))
def strip_ext_gpeak(gpeak): return re.sub(r'\.(gappedPeak|GappedPeak)\.gz$', '', str(gpeak))
def strip_ext_bpeak(bpeak): return re.sub(r'\.(broadPeak|BroadPeak)\.gz$', '', str(bpeak))
def strip_ext_bigwig(bw): return re.sub(r'\.(bigwig|bw)$', '', str(bw))
def strip_ext_gz(f): return re.sub(r'\.gz$', '', str(f))

def get_ext(f):
    """Return the extension of a file (ignoring .gz)."""
    f_wo_gz = re.sub(r'\.gz$', '', str(f))
    return f_wo_gz.split('.')[-1]

def strip_ext(f, ext=''):
    """Remove file extension and optionally a .gz suffix."""
    if ext == '':
        ext = get_ext(f)
    return re.sub(r'\.({}|{}\.gz)$'.format(ext, ext), '', str(f))

# ----------- Peak Type Utilities -----------

def get_peak_type(peak):
    """Determine peak type based on its extension."""
    if strip_ext_npeak(peak) != peak:
        return 'narrowPeak'
    elif strip_ext_rpeak(peak) != peak:
        return 'regionPeak'
    elif strip_ext_bpeak(peak) != peak:
        return 'broadPeak'
    elif strip_ext_gpeak(peak) != peak:
        return 'gappedPeak'
    else:
        raise Exception('Unsupported peak type for {}'.format(peak))

def strip_ext_peak(peak):
    """Return filename without extension for known peak types."""
    peak_type = get_peak_type(peak)
    return strip_ext(peak, peak_type)

# ----------- Read/Write and Formatting Utilities -----------

def human_readable_number(num):
    """Format a number with human-readable suffix (e.g., K, M)."""
    for unit in ['', 'K', 'M', 'G', 'T', 'P']:
        if abs(num) < 1000:
            return '{}{}'.format(num, unit)
        num = int(num / 1000.0)
    return '{}E'.format(num)

def human_readable_filesize(num):
    """Format a byte count into a human-readable string."""
    for unit in ['', 'KB', 'MB', 'GB', 'TB', 'PB']:
        if abs(num) < 1024.0:
            return '{}{}'.format(num, unit)
        num = int(num / 1024.0)
    return '{}EB'.format(num)

def read_tsv(tsv):
    """Read a tab-separated file into a list of lists."""
    result = []
    with open(tsv, 'r') as fp:
        for row in csv.reader(fp, delimiter='\t'):
            result.append(row if row else [''])
    return result

def write_tsv(tsv, arr):
    """Write a list of lists to a tab-separated file."""
    with open(tsv, 'w') as fp:
        for i, arr2 in enumerate(arr):
            line = '\t'.join([str(a) for a in arr2])
            line += '\n' if i < len(arr) - 1 else ''
            fp.write(line)

def write_txt(f, s):
    """Write a string or list of strings to a text file."""
    with open(f, 'w') as fp:
        arr = [s] if type(s) != list else s
        for a in arr:
            fp.write(str(a) + '\n')

# ----------- Filesystem Utilities -----------

def mkdir_p(dirname):
    """Create a directory if it does not exist."""
    if dirname and not os.path.exists(dirname):
        os.makedirs(dirname)

def rm_f(files):
    """Remove one or more files."""
    if files:
        cmd = 'rm -f {}'.format(' '.join(files)) if isinstance(files, list) else f'rm -f {files}'
        run_shell_cmd(cmd)

def ls_l(d):
    """Run ls -l on a directory."""
    cmd = 'ls -l {}'.format(d)
    run_shell_cmd(cmd)

def touch(f):
    """Create an empty file or update its timestamp."""
    run_shell_cmd('touch {}'.format(f))

# ----------- Compression Utilities -----------

def gunzip(f, suffix, out_dir):
    """Uncompress .gz file to out_dir with optional suffix."""
    if not f.endswith('.gz'):
        raise Exception('Cannot gunzip a file without .gz extension.')
    gunzipped = os.path.join(out_dir, os.path.basename(strip_ext_gz(f)))
    if suffix:
        gunzipped += f'.{suffix}'
    run_shell_cmd(f'zcat -f {f} > {gunzipped}')
    return gunzipped

def untar(tar, out_dir):
    """Extract tar or tar.gz archive to a given directory."""
    cmd = 'tar zxvf {} --no-same-owner -C {}'.format(tar, out_dir or '.') if tar.endswith('.gz') \
        else 'tar xvf {} --no-same-owner -C {}'.format(tar, out_dir or '.')
    run_shell_cmd(cmd)

# ----------- File Checks -----------

def get_num_lines(f):
    """Return number of lines in a (possibly gzipped) file."""
    return int(run_shell_cmd(f'zcat -f {f} | wc -l'))

def assert_file_not_empty(f, help=''):
    """Check if file exists and is not empty."""
    if not os.path.exists(f):
        raise Exception(f'File does not exist ({f}). Help: {help}')
    elif get_num_lines(f) == 0:
        raise Exception(f'File is empty ({f}). Help: {help}')

# ----------- Linking and Copying -----------

def hard_link(f, link):
    """Create a hard link to a file."""
    if os.path.abspath(f) == os.path.abspath(link):
        raise Exception(f'Trying to hard-link itself. {f}')
    os.link(f, link)
    return link

def make_hard_link(f, out_dir):
    """Create a hard link in a specified directory."""
    linked = os.path.join(out_dir, os.path.basename(f))
    rm_f(linked)
    return hard_link(f, linked)

def soft_link(f, link):
    """Create a symbolic (soft) link."""
    if os.path.abspath(f) == os.path.abspath(link):
        raise Exception(f'Trying to soft-link itself. {f}')
    os.symlink(f, link)
    return link

def make_soft_link(f, out_dir):
    """Create a soft link in the specified directory."""
    linked = os.path.join(out_dir, os.path.basename(f))
    rm_f(linked)
    return soft_link(f, linked)

def copy_f_to_f(f, dest):
    """Copy a file to a destination path."""
    if os.path.abspath(f) == os.path.abspath(dest):
        raise Exception(f'Trying to copy to itself. {f}')
    run_shell_cmd(f'cp -f {f} {dest}')
    return dest

def copy_f_to_dir(f, out_dir):
    """Copy a file into a directory."""
    dest = os.path.join(out_dir, os.path.basename(f))
    return copy_f_to_f(f, dest)

# ----------- Sorting Helpers -----------

def get_gnu_sort_param(max_mem_job, ratio=0.5):
    """Compute memory size parameter for GNU sort based on max memory."""
    mem_mb = int(math.ceil(max_mem_job * ratio / (1024 * 1024)))
    return f'-S {mem_mb}M'

# ----------- Time and Conversion -----------

def now():
    """Return current timestamp as string."""
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

def pdf2png(pdf, out_dir):
    """Convert first page of a PDF file to PNG image."""
    prefix = os.path.join(out_dir, os.path.basename(strip_ext(pdf)))
    png = f'{prefix}.png'
    cmd = ('gs -dFirst
