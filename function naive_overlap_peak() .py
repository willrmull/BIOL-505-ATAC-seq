#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#Identify peaks in the pooled sample that are supported by both Rep1 and Rep2
#The result: Peaks in the pooled sample that have at least 50% reciprocal overlap with both Rep1 and Rep2.
intersectBed -wo -a Pooled.narrowPeak.gz -b Rep1.narrowPeak.gz | 
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort | uniq | 
intersectBed -wo -a stdin -b Rep2.narrowPeak.gz | 
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort | uniq > PooledInRep1AndRep2.narrowPeak.gz


#Filter out blacklisted regions and cap peak scores
#intersect -v: Removes peaks that overlap with blacklisted regions.
bedtools intersect -v -a PooledInRep1AndRep2.narrowPeak.gz -b ${BLACKLIST} | awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' | grep -P 'chr[\dXY]+[ \t]'  | gzip -nc > PooledInRep1AndRep2.filt.narrowPeak.gz

"""

