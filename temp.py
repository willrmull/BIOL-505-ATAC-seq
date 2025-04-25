# -*- coding: utf-8 -*-

##-t $tag: input ChIP-seq reads in BED format.

##-f BED: specifies BED input format.

##-n "$prefix": output prefix.

##-g "$gensz": effective genome size.

##-p $pval_thresh: p-value threshold for peak calling.

##--shift, --extsize: manual fragment size modeling.

##--nomodel: disables MACS2's internal model building.

##-B --SPMR: outputs signal tracks in bedGraph format (normalized per million reads).

##--keep-dup all: keeps all duplicate reads.

##--call-summits: refines peak regions to highest point.

macs2 callpeak \
    -t $tag -f BED -n "$prefix" -g "$gensz" -p $pval_thresh \
    --shift $shiftsize  --extsize $smooth_window --nomodel -B --SPMR --keep-dup all --call-summits
    
# Sort by Col8 in descending order and replace long peak names in Column 4 with Peak_<peakRank>

##Sorts peaks by column 8 (e.g. p-value or score) in descending order.

##Replaces the peak name in column 4 with Peak_1, Peak_2, etc.

##Outputs only the top ${NPEAKS} peaks to $peakfile, compresses with gzip.

sort -k 8gr,8gr "$prefix"_peaks.narrowPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' 
| head -n ${NPEAKS} | gzip -nc > $peakfile 

#Removes files not needed after sorting and filtering peaks.

rm -f "$prefix"_peaks.narrowPeak
rm -f "$prefix"_peaks.xls
rm -f "$prefix"_summits.bed

##Computes Poisson p-value-based enrichment signal track.

##Calculates fold enrichment (FE) by comparing treatment vs control.

##Produces a new .bdg file.

macs2 bdgcmp -t "$prefix"_treat_pileup.bdg -c "$prefix"_control_lambda.bdg\
    --o-prefix "$prefix" -m FE

##Ensures coordinates don’t exceed chromosome lengths (slopBed and bedClip).

slopBed -i "$prefix"_FE.bdg -g "$chrsz" -b 0 | bedClip stdin "$chrsz" $fc_bedgraph
rm -f "$prefix"_FE.bdg

##Sorts and converts the FE bedGraph to a BigWig file for efficient visualization.

sort -k1,1 -k2,2n $fc_bedgraph > $fc_bedgraph_srt
bedGraphToBigWig $fc_bedgraph_srt "$chrsz" "$fc_bigwig"

##Removes intermediate files.

rm -f $fc_bedgraph $fc_bedgraph_srt

# sval counts the number of tags per million in the (compressed) BED file

sval=$(wc -l <(zcat -f "$tag") | awk '{printf "%f", $1/1000000}')

##Computes Poisson p-value-based enrichment signal track.

macs2 bdgcmp\
    -t "$prefix"_treat_pileup.bdg -c "$prefix"_control_lambda.bdg\
    --o-prefix "$prefix" -m ppois -S "${sval}"
    
## Processes and converts the p-value signal track to BigWig.   
    
slopBed -i "$prefix"_ppois.bdg -g "$chrsz" -b 0 | bedClip stdin "$chrsz" $pval_bedgraph
rm -f "$prefix"_ppois.bdg
sort -k1,1 -k2,2n $pval_bedgraph > $pval_bedgraph_srt
bedGraphToBigWig $pval_bedgraph_srt "$chrsz" "$pval_bigwig"
rm -f $pval_bedgraph $pval_bedgraph_srt

## Deletes temporary signal tracks.

rm -f "$prefix"_treat_pileup.bdg "$prefix"_control_lambda.bdg


