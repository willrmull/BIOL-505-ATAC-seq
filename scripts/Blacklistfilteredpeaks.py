#!/usr/bin/env bash
# -*- coding: utf-8 -*-

# Define filenames
PEAK="${PREFIX}.narrowPeak.gz"
FILTERED_PEAK="${PREFIX}.narrowPeak.filt.gz"

# Blacklist filtering and cleanup
bedtools intersect -v -a ${PEAK} -b ${BLACKLIST} \
    | awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' \
    | grep -P 'chr[\dXY]+[ \t]' \
    | gzip -nc > ${FILTERED_PEAK}


