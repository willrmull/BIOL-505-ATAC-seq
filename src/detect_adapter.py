import sys
import gzip

# Global flag to control verbose output
VERBOSE = False

# Dictionary of common adapter sequences used in various library preparations
# Keys are names; values are the byte-encoded sequences
adapters = {
    'Illumina': b'AGATCGGAAGAGC',
    'Nextera ': b'CTGTCTCTTATA',
    'smallRNA': b'TGGAATTCTCGG'
}

# Function to open a file (gzip-compressed or not)
def open_gz(fname):
    # If file ends in '.gz', open with gzip, else open as binary
    return gzip.open(fname) if fname.endswith('.gz') else open(fname, 'rb')

# Detect which adapter sequences are present in reads
def detect_adapters_and_cnts(fname, max_n_lines=1000000):
    # Initialize counter for each adapter type
    adapter_cnts = {
        'Illumina': 0,
        'Nextera ': 0,
        'smallRNA': 0
    }

    # Open input FASTQ file (gzip or plain text)
    with open_gz(fname) as fp:
        # Iterate through lines up to max_n_lines
        for seq_index, line in enumerate(fp):
            if seq_index >= max_n_lines:
                break

            # In a FASTQ file, every 4th line starting at 1 is a sequence line
            if seq_index % 4 != 1:
                continue

            # Check if any known adapter appears in this


