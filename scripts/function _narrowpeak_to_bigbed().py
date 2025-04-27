#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#This code takes a gzipped narrowPeak file and a chromosome size file, 
#filters and formats them properly, then converts the peak file into a BigBed file for genome browser visualization.
#Filter chromosome sizes
cat ${CHRSZ} | grep -P 'chr[\dXY]+[ \t]' > ${BIGBED}.chrsz.tmp
#Sort peak file
zcat ${PEAK} | sort -k1,1 -k2,2n > ${BIGBED}.tmp
# Clip peaks to chromosome sizes
#bedClip: Ensures no peak entries extend beyond chromosome boundaries.
#Input: Sorted BED file and chromosome size file.
#Output: Clipped BED file: ${BIGBED}.tmp2.
bedClip ${BIGBED}.tmp ${BIGBED}.chrsz.tmp ${BIGBED}.tmp2
#Convert BED to BigBed
bedToBigBed -type=bed6+4 -as=narrowPeak.as  ${BIGBED}.tmp2 ${BIGBED}.chrsz.tmp ${BIGBED}
#Clean up temporary files
rm -f ${BIGBED}.tmp ${BIGBED}.tmp2 ${BIGBED}.chrsz.tmp



