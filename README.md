# BIOL-505-ATAC-seq
#ATAC-seq pipeline is an application that uses both bash and python scripts to fetch genome and ATAC-sequencing data, then using a workflow to process #scripts that process that data into peaks and visualization. 
#Features
#-gets genome and ATAC-sequence data
#-trims ATAC-seq data and aligns it to the genome
#-filters the aligned sequences
#-identifying open regions of the chromatin with peak calling
#-visualization

#Usage:
#./final_genome_builder.sh
#java -D.config.file=cromwell.conf -jar cromwell-50.jar run atac.wdl --inputs atac_input.json
