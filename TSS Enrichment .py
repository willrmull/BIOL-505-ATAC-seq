#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#bam_file: Path to a BAM file containing aligned reads.

#tss: A BED file with TSS positions.

#prefix: Prefix for output files.

#chromsizes: Genome sizes file (for proper slop operations).

#read_len: Length of sequencing reads.

#bins: Number of bins across the TSS region (default 400).

#bp_edge: Number of base pairs upstream and downstream of TSS to include (default 2000 bp).

#processes: Number of parallel processes to use.

#greenleaf_norm: Whether to apply Greenleaf normalization (TSS enrichment normalization

def make_tss_plot(bam_file, tss, prefix, chromsizes, read_len, bins=400, bp_edge=2000,
                  processes=8, greenleaf_norm=True):
    '''
    Take bootstraps, generate tss plots, and get a mean and
    standard deviation on the plot. Produces 2 plots. One is the
    aggregation plot alone, while the other also shows the signal
    at each TSS ordered by strength.
    '''
    logging.info('Generating tss plot...')
#    Sets the output file names for:

#    A basic line plot.

#    A detailed heatmap-style plot.

#    A QC file with the TSS enrichment value.

    tss_plot_file = '{0}.tss_enrich.png'.format(prefix)
    tss_plot_large_file = '{0}.large_tss_enrich.png'.format(prefix)
    tss_log_file = '{0}.tss_enrich.qc'.format(prefix)

    #  Load and Extend TSS Regions
    tss = pybedtools.BedTool(tss)
    tss_ext = tss.slop(b=bp_edge, g=chromsizes)

    # Load BAM File and Extract Signal
    bam = metaseq.genomic_signal(bam_file, 'bam') # Need to shift reads and just get ends, just load bed file?
    bam_array = bam.array(tss_ext, bins=bins, shift_width = -read_len/2, # Shift to center the read on the cut site
                          processes=processes, stranded=True)

    # Actually first build an "ends" file
    #get_ends = '''zcat {0} | awk -F '\t' 'BEGIN {{OFS="\t"}} {{if ($6 == "-") {{$2=$3-1; print}} else {{$3=$2+1; print}} }}' | gzip -c > {1}_ends.bed.gz'''.format(bed_file, prefix)
    #print(get_ends)
    #os.system(get_ends)

    #bed_reads = metaseq.genomic_signal('{0}_ends.bed.gz'.format(prefix), 'bed')
    #bam_array = bed_reads.array(tss_ext, bins=bins,
    #                      processes=processes, stranded=True)
 # Normalization (Greenleaf style): Find the avg height
    # at the end bins and take fold change over that
    #Normalization
    if greenleaf_norm:
        # Use enough bins to cover 100 bp on either end
        num_edge_bins = int(100/(2*bp_edge/bins))
        bin_means = bam_array.mean(axis=0)
        avg_noise = (sum(bin_means[:num_edge_bins]) +
                     sum(bin_means[-num_edge_bins:]))/(2*num_edge_bins)
        bam_array /= avg_noise
    else:
        bam_array /= bam.mapped_read_count() / 1e6

    # Generate a line plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    x = np.linspace(-bp_edge, bp_edge, bins)

    ax.plot(x, bam_array.mean(axis=0), color='r', label='Mean')
    ax.axvline(0, linestyle=':', color='k')

    # Note the middle high point (TSS)
    tss_point_val = max(bam_array.mean(axis=0))

    # write tss_point_val to file
    with open(tss_log_file, 'w') as fp:
        fp.write(str(tss_point_val))

    ax.set_xlabel('Distance from TSS (bp)')
    if greenleaf_norm:
        ax.set_ylabel('TSS Enrichment')
    else:
        ax.set_ylabel('Average read coverage (per million mapped reads)')
    ax.legend(loc='best')

    fig.savefig(tss_plot_file)

    # Print a more complicated plot with lots of info

    # Find a safe upper percentile - we can't use X if the Xth percentile is 0
    upper_prct = 99
    if mlab.prctile(bam_array.ravel(), upper_prct) == 0.0:
        upper_prct = 100.0

    plt.rcParams['font.size'] = 8
    #Detailed Heatmap Plot
    fig = metaseq.plotutils.imshow(bam_array,
                                   x=x,
                                   figsize=(5, 10),
                                   vmin=5, vmax=upper_prct, percentile=True,
                                   line_kwargs=dict(color='k', label='All'),
                                   fill_kwargs=dict(color='k', alpha=0.3),
                                   sort_by=bam_array.mean(axis=1))

    # And save the file
    fig.savefig(tss_plot_large_file)

    return tss_plot_file, tss_plot_large_file, tss_log_file



