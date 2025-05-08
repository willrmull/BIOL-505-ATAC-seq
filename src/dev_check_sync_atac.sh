#!/bin/bash

# This script compares the versions of various source code files
# between a local working directory and the official pipeline source directory.
# It's useful for identifying any local modifications made to the ATAC-seq pipeline scripts.

# 'diff' is used to show line-by-line differences between the original and local versions.
# Each command compares the source version from the official pipeline (../../atac-seq-pipeline/src/)
# against a potentially modified version in the current working directory.

# Compare the BAM to tagAlign task implementation
diff ../../atac-seq-pipeline/src/encode_task_bam2ta.py encode_task_bam2ta.py 

# Compare the blacklist filtering module
diff ../../atac-seq-pipeline/src/encode_lib_blacklist_filter.py encode_lib_blacklist_filter.py 

# Compare genomic helper functions (e.g., for working with BAMs, chromosomes)
diff ../../atac-seq-pipeline/src/encode_lib_genomic.py encode_lib_genomic.py 

# Compare the log parsing utility module
diff ../../atac-seq-pipeline/src/encode_lib_log_parser.py encode_lib_log_parser.py 

# Compare the common utility functions used throughout the pipeline
diff ../../atac-seq-pipeline/src/encode_lib_common.py encode_lib_common.py 

# Compare the Bowtie2 alignment task implementation
diff ../../atac-seq-pipeline/src/encode_task_bowtie2.py encode_task_bowtie2.py 

# Compare the filtering task (e.g., removing low-quality or duplicate reads)
diff ../../atac-seq-pipeline/src/encode_task_filter.py encode_task_filter.py 

# Compare the post-alignment processing step
diff ../../atac-seq-pipeline/src/encode_task_post_align.py encode_task_post_align.py 

# Compare FRiP score calculation module
diff ../../atac-seq-pipeline/src/encode_lib_frip.py encode_lib_frip.py 

# Compare the IDR (Irreproducible Discovery Rate) peak comparison task
diff ../../atac-seq-pipeline/src/encode_task_idr.py encode_task_idr.py 

# Compare the peak overlap and reproducibility module
diff ../../atac-seq-pipeline/src/encode_task_overlap.py encode_task_overlap.py 

# Compare the task for pooling tagAlign files
diff ../../atac-seq-pipeline/src/encode_task_pool_ta.py encode_task_pool_ta.py 

# Compare the QC report generator script
diff ../../atac-seq-pipeline/src/encode_task_qc_report.py encode_task_qc_report.py 

# Compare the reproducibility analysis script
diff ../../atac-seq-pipeline/src/encode_task_reproducibility.py encode_task_reproducibility.py 

# Compare the SPR (Self-consistency Peak Recall) task implementation
diff ../../atac-seq-pipeline/src/encode_task_spr.py encode_task_spr.py 

# Compare the cross-correlation analysis script
diff ../../atac-seq-pipeline/src/encode_task_xcor.py encode_task_xcor.py

# Compare the Jensen-Shannon divergence (JSD) task script
diff ../../atac-seq-pipeline/src/encode_task_jsd.py encode_task_jsd.py

# Compare the GC bias analysis task script
diff ../../atac-seq-pipeline/src/encode_task_gc_bias.py encode_task_gc_bias.py


