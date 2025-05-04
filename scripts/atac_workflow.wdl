#This input block configures how the pipeline should run, including:
#What software environment to use (docker, conda, etc.).
#What reference genome and experimental data to use.
#Optional tuning parameters for trimming, alignment, peak calling, etc.
version 1.0
#Provide container or environment info
struct RuntimeEnvironment {
    String docker
    String singularity
    String conda
}
            pipeline_metadata: {
                title: 'Pipeline metadata',
                description: 'Metadata for a pipeline (e.g. title and description).'
            },
            reference_genome: {
                title: 'Reference genome',
                description: 'Genome specific files. e.g. reference FASTA, bowtie2 index, chromosome sizes file.',
                help: 'Choose one atac.genome_tsv file that defines all genome specific parameters in it or define each genome specific parameter in input JSON to override those defined in genome TSV file. If you use Caper then use https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v1/[GENOME]_caper.tsv. Caper will automatically download/install all files defined in such TSV. Otherwise download genome TSV file by using a shell script (scripts/download_genome_data.sh [GENOME] [DEST_DIR]). Supported genomes are hg38, hg19, mm10 and mm9. See pipeline documentation if you want to build genome database from your own FASTA file. If some genome data are missing then analyses using such data will be skipped.'
            },
            input_genomic_data: {
                title: 'Input genomic data',
                description: 'Genomic input files for experiment.',
                help: 'Pipeline can start with any types of experiment data (e.g. FASTQ, BAM, NODUP_BAM, TAG-ALIGN, PEAK). Choose one type and leave others empty. FASTQs have a variable for each biological replicate. e.g. atac.fastqs_rep1_R1 and atac.fastqs_rep2_R1. You can define up to 10 experiment replicates. For other types, there is an array to define file for each biological replicate. e.g. atac.bams: ["rep1.bam", "rep1.bam"]. Define sequential endedness with atac.paired_end, if you have mixed SE and PE replicates then define atac.paired_ends instead for each replicate. e.g. atac.paired_ends: [false, true].'
            },
            adapter_trimming: {
                title: 'Adapter trimming',
                description: 'Parameters for adapter trimming.',
                help: 'Use atac.auto_detect_adapter to automatically detect/trim 3 adapters (Illumina: AGATCGGAAGAGC, Nextera: CTGTCTCTTATA, smallRNA: TGGAATTCTCGG) or manually define adapter sequence to be trimmed (atac.adapter or atac.adapters_repX_RY). Leave all parameters undefined/empty if your FASTQs are already trimmed.'
            },
            pipeline_parameter: {
                title: 'Pipeline parameter',
                description: 'Pipeline type and flags to turn on/off analyses.',
                help: 'Pipeline can run as DNAse-seq mode. The only difference is TN5-shifting of read in ATAC-seq mode. Use atac.align_only to align FASTQs without peak calling.'
            },
            alignment: {
                title: 'Alignment',
                description: 'Parameters for alignment.',
                help: 'Pipeline calculates mitochondrial fraction of reads in raw BAM. But after that it filters out mitochondrial reads (e.g. chrM, MT) from NODUP_BAMs (filtered/deduped). It is controlled by atac.filter_chrs array. If you want to keep mitochondrial reads then make this array empty.'
            },
            peak_calling: {
                title: 'Peak calling',
                description: 'Parameters for peak calling.',
                help: 'This group includes statistical thresholds for peak-calling or post-peak-calling analyses: p-val, FDR, IDR.'
            },
            resource_parameter: {
                title: 'Resource parameter',
                description: 'Number of CPUs (threads), max. memory and walltime for tasks.',
                help: 'Resource settings are used for determining an instance type on cloud backends (e.g. GCP, AWS) and used for submitting tasks to a cluster engine (e.g. SLURM, SGE, ...). Walltime (atac.*_time_hr) is only used for cluster engines. Other tasks default to use 1 CPU and 4GB of memory.'
            }
        }
    }
    input {
        # group: runtime_environment
        String docker = 'encodedcc/atac-seq-pipeline:v2.2.2'
        String singularity = 'https://encode-pipeline-singularity-image.s3.us-west-2.amazonaws.com/atac-seq-pipeline_v2.2.2.sif'
        String conda = 'encd-atac'
        String conda_macs2 = 'encd-atac-macs2'
        String conda_spp = 'encd-atac-spp'
        String conda_python2 = 'encd-atac-py2'

        # group: pipeline_metadata
        String title = 'ATAC-seq on Danaus plexippus'
        String description = 'ATAC-seq on mid-fifth instar Danaus plexippus caterpillars'

        # group: reference_genome
        #Choose either a to reference genome TSV file or to define them all manually
        File? genome_tsv = final_danaus_plexippus.tsv
        String? genome_name = dp
        File? ref_fa
        File? ref_mito_fa
        File? bowtie2_idx_tar
        File? bowtie2_mito_idx_tar
        File? chrsz
        File? blacklist
        File? blacklist2
        String? regex_bfilt_peak_chr_name
        String? gensz
        File? tss
        File? dnase
        File? prom
        File? enh
        File? reg2map
        File? reg2map_bed
        File? roadmap_meta

        # group: input_genomic_data
        #Choose one data type to start with (FASTQ, BAM, etc.).
        Boolean? paired_end 
        Array[Boolean] paired_ends = [true]
        Array[File] fastqs_rep1_R1 = [s3://sra-pub-src-16/SRR13566303/ATAC028_DpM5thHD2_R1.fastq.gz.1]
        Array[File] fastqs_rep1_R2 = [s3://sra-pub-src-16/SRR13566303/ATAC028_DpM5thHD2_R2.fastq.gz.1]
        Array[File] bams = []
        Array[File] nodup_bams = []
        Array[File] tas = []
        Array[File] peaks = []
        Array[File] peaks_pr1 = []
        Array[File] peaks_pr2 = []
        File? peak_pooled
        File? peak_ppr1
        File? peak_ppr2

        # group: set pipeline_parameters
        String pipeline_type = 'atac'
        Boolean align_only = false
        Boolean true_rep_only = false
        Boolean enable_xcor = false
        Boolean enable_count_signal_track = false
        Boolean enable_idr = false
        Boolean enable_preseq = false
        Boolean enable_fraglen_stat = true
        Boolean enable_tss_enrich = true
        Boolean enable_gc_bias = true

        # group: adapter_trimming
        String cutadapt_param = '-e 0.1 -m 5'
        #Automatically detect adapters
        Boolean auto_detect_adapter = true
        String? adapter
        #Or define adapters manually
        Array[String] adapters_rep1_R1 = []
        Array[String] adapters_rep1_R2 = []
        # group: alignment
        Int multimapping = 4
        String dup_marker = 'picard'
        Boolean no_dup_removal = false
        Int mapq_thresh = 30
        filter_chrs = []
        Int subsample_reads = 0
        Int xcor_subsample_reads = 25000000
        Array[Int?] read_len = []
        Int pseudoreplication_random_seed = 0

        # group: peak_calling
        #limits the maximum number of peaks kept from peak calling.
        Int cap_num_peak = 300000
        #Sets the p-value threshold for MACS2 peak calling
        Float pval_thresh = 0.01
        #sets the smoothing window size (in base pairs) used when computing signal for peak calling.
        #Smoothing helps reduce noise by averaging signals over nearby bases.
        Int smooth_win = 150
        #Sets the IDR (Irreproducible Discovery Rate) threshold for filtering reproducible peaks between replicates.
        #there’s a ≤5% chance that a peak is due to noise.
        Float idr_thresh = 0.05

        # group: resource_parameter
        Int align_cpu = 6
        Float align_mem_factor = 0.15
        Int align_time_hr = 48
        Float align_disk_factor = 8.0

        Int filter_cpu = 4
        Float filter_mem_factor = 0.4
        Int filter_time_hr = 24
        Float filter_disk_factor = 8.0

        Int bam2ta_cpu = 2
        Float bam2ta_mem_factor = 0.3
        Int bam2ta_time_hr = 12
        Float bam2ta_disk_factor = 4.0

        Float spr_mem_factor = 20.0
        Float spr_disk_factor = 30.0

        Int jsd_cpu = 4
        Float jsd_mem_factor = 0.1
        Int jsd_time_hr = 12
        Float jsd_disk_factor = 2.0

        Int xcor_cpu = 2
        Float xcor_mem_factor = 1.0
        Int xcor_time_hr = 6
        Float xcor_disk_factor = 4.5

        Int call_peak_cpu = 2
        Float call_peak_mem_factor = 4.0
        Int call_peak_time_hr = 24
        Float call_peak_disk_factor = 30.0

        Float macs2_signal_track_mem_factor = 12.0
        Int macs2_signal_track_time_hr = 24
        Float macs2_signal_track_disk_factor = 80.0

        Float preseq_mem_factor = 0.5
        Float preseq_disk_factor = 5.0

        String? filter_picard_java_heap
        String? preseq_picard_java_heap
        String? fraglen_stat_picard_java_heap
        String? gc_bias_picard_java_heap
    }
    #This block of code configures the resource requirements for various tasks in the ENCODE ATAC-seq pipeline. 
    #It helps the pipeline schedule and run jobs on a compute backend (e.g. a cloud provider, SLURM, or a local workstation) by specifying:
    #How many CPU cores each task should use
    #How much memory to allocate (based on genome size)
    #How much disk space to reserve
    #Maximum time a task is allowed to run
    #Optional Java heap size overrides for Picard-based tools
    parameter_meta {
        docker: {
            description: 'Default Docker image URI to run WDL tasks.',
            group: 'runtime_environment',
            example: 'ubuntu:20.04'
        }
        singularity: {
            description: 'Default Singularity image URI to run WDL tasks. For Singularity users only.',
            group: 'runtime_environment',
            example: 'docker://ubuntu:20.04'
        }
        conda: {
            description: 'Default Conda environment name to run WDL tasks. For Conda users only.',
            group: 'runtime_environment',
            example: 'encd-atac'
        }
        conda_macs2: {
            description: 'Conda environment name for task macs2. For Conda users only.',
            group: 'runtime_environment',
            example: 'encd-atac-macs2'
        }
        conda_spp: {
            description: 'Conda environment name for tasks spp/xcor. For Conda users only.',
            group: 'runtime_environment',
            example: 'encd-atac-spp'
        }
        conda_python2: {
            description: 'Conda environment name for tasks with python2 wrappers (tss_enrich). For Conda users only.',
            group: 'runtime_environment',
            example: 'encd-atac-py2'
        }
        title: {
            description: 'Experiment title.',
            group: 'pipeline_metadata',
            example: 'ENCSR356KRQ (subsampled 1/400)'
        }
        description: {
            description: 'Experiment description.',
            group: 'pipeline_metadata',
            example: 'ATAC-seq on primary keratinocytes in day 0.0 of differentiation (subsampled 1/400)'
        }
        genome_tsv: {
            description: 'Reference genome database TSV.',
            group: 'reference_genome',
            help: 'This TSV files includes all genome specific parameters (e.g. reference FASTA, bowtie2 index). You can still invidiaully define any parameters in it. Parameters defined in input JSON will override those defined in genome TSV.',
            example: 'https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v1/hg38_caper.tsv'
        }
        genome_name: {
            description: 'Genome name.',
            group: 'reference_genome'
        }
        ref_fa: {
            description: 'Reference FASTA file.',
            group: 'reference_genome'
        }
        bowtie2_idx_tar: {
            description: 'BWA index TAR file.',
            group: 'reference_genome'
        }
        chrsz: {
            description: '2-col chromosome sizes file.',
            group: 'reference_genome'
        }
        blacklist: {
            description: 'Blacklist file in BED format.',
            group: 'reference_genome',
            help: 'Peaks will be filtered with this file.'
        }
        blacklist2: {
            description: 'Secondary blacklist file in BED format.',
            group: 'reference_genome',
            help: 'If it is defined, it will be merged with atac.blacklist. Peaks will be filtered with merged blacklist.'
        }
        regex_bfilt_peak_chr_name: {
            description: 'Reg-ex for chromosomes to keep while filtering peaks.',
            group: 'reference_genome',
            help: 'Chromosomes defined here will be kept. All other chromosomes will be filtered out in .bfilt. peak file. This is done along with blacklist filtering peak file.'
        }
        gensz: {
            description: 'Genome sizes. "hs" for human, "mm" for mouse or sum of 2nd columnin chromosome sizes file.',
            group: 'reference_genome'
        }
        tss: {
            description: 'TSS file in BED format.',
            group: 'reference_genome'
        }
        dnase: {
            description: 'Open chromatin regions file in BED format.',
            group: 'reference_genome'
        }
        prom: {
            description: 'Promoter regions file in BED format.',
            group: 'reference_genome'
        }
        enh: {
            description: 'Enhancer regions file in BED format.',
            group: 'reference_genome'
        }
        reg2map: {
            description: 'Cell type signals file.',
            group: 'reference_genome'
        }
        reg2map_bed: {
            description: 'File of regions used to generate reg2map signals.',
            group: 'reference_genome'
        }
        roadmap_meta: {
            description: 'Roadmap metadata.',
            group: 'reference_genome'
        }
        paired_end: {
            description: 'Sequencing endedness.',
            group: 'input_genomic_data',
            help: 'Setting this on means that all replicates are paired ended. For mixed samples, use atac.paired_ends array instead.',
            example: true
        }
        paired_ends: {
            description: 'Sequencing endedness array (for mixed SE/PE datasets).',
            group: 'input_genomic_data',
            help: 'Whether each biological replicate is paired ended or not.'
        }
        fastqs_rep1_R1: {
            description: 'Read1 FASTQs to be merged for a biological replicate 1.',
            group: 'input_genomic_data',
            help: 'Define if you want to start pipeline from FASTQ files. Pipeline can start from any type of inputs (e.g. FASTQs, BAMs, ...). Choose one type and fill paramters for that type and leave other undefined. Especially for FASTQs, we have individual variable for each biological replicate to allow FASTQs of technical replicates can be merged. Make sure that they are consistent with read2 FASTQs (atac.fastqs_rep1_R2). These FASTQs are usually technical replicates to be merged.',
            example: [
                "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep1/pair1/ENCFF341MYG.subsampled.400.fastq.gz",
                "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep1/pair1/ENCFF106QGY.subsampled.400.fastq.gz"
            ]
        }
        fastqs_rep1_R2: {
            description: 'Read2 FASTQs to be merged for a biological replicate 1.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read1 FASTQs (atac.fastqs_rep1_R1). These FASTQs are usually technical replicates to be merged.',
            example: [
                "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep1/pair2/ENCFF248EJF.subsampled.400.fastq.gz",
                "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep1/pair2/ENCFF368TYI.subsampled.400.fastq.gz"
            ]
        }
        bams: {
            description: 'List of unfiltered/raw BAM files for each biological replicate.',
            group: 'input_genomic_data',
            help: 'Define if you want to start pipeline from BAM files. Unfiltered/raw BAM file generated from aligner (e.g. bowtie2). Each entry for each biological replicate. e.g. [rep1.bam, rep2.bam, rep3.bam, ...].'
        }
        nodup_bams: {
            description: 'List of filtered/deduped BAM files for each biological replicate',
            group: 'input_genomic_data',
            help: 'Define if you want to start pipeline from filtered BAM files. Filtered/deduped BAM file. Each entry for each biological replicate. e.g. [rep1.nodup.bam, rep2.nodup.bam, rep3.nodup.bam, ...].'
        }
        tas: {
            description: 'List of TAG-ALIGN files for the biological replicate.',
            group: 'input_genomic_data',
            help: 'Define if you want to start pipeline from TAG-ALIGN files. TAG-ALIGN is in a 6-col BED format. It is a simplified version of BAM. Each entry for each biological replicate. e.g. [rep1.tagAlign.gz, rep2.tagAlign.gz, ...].'
        }
        peaks: {
            description: 'List of NARROWPEAK files (not blacklist filtered) for each biological replicate.',
            group: 'input_genomic_data',
            help: 'Define if you want to start pipeline from PEAK files. Each entry for each biological replicate. e.g. [rep1.narrowPeak.gz, rep2.narrowPeak.gz, ...]. Define other PEAK parameters (e.g. atac.peaks_pr1, atac.peak_pooled) according to your flag settings (e.g. atac.true_rep_only) and number of replicates. If you have more than one replicate then define atac.peak_pooled, atac.peak_ppr1 and atac.peak_ppr2. If atac.true_rep_only flag is on then do not define any parameters (atac.peaks_pr1, atac.peaks_pr2, atac.peak_ppr1 and atac.peak_ppr2) related to pseudo replicates.'
        }    
        pipeline_type: {
            description: 'Pipeline type. atac for ATAC-Seq or dnase for DNase-Seq.',
            group: 'pipeline_parameter',
            help: 'The only difference of two types is that TN5 shifting of TAG-ALIGN is done for atac. TAG-ALIGN is in 6-col BED format. It is a simplified version of BAM.',
            choices: ['atac', 'dnase'],
            example: 'atac'
        }
        align_only: {
            description: 'Align only mode.',
            group: 'pipeline_parameter',
            help: 'Reads will be aligned but there will be no peak-calling on them.'
        }
        true_rep_only: {
            description: 'Disables all analyses related to pseudo-replicates.',
            group: 'pipeline_parameter',
            help: 'Pipeline generates 2 pseudo-replicate from one biological replicate. This flag turns off all analyses related to pseudos (with prefix/suffix pr, ppr).'
        }
        enable_count_signal_track: {
            description: 'Enables generation of count signal tracks.',
            group: 'pipeline_parameter'
        }
        enable_fraglen_stat: {
            description: 'Enables calculation of fragment length distribution/statistics.',
            group: 'pipeline_parameter'
        }
        enable_tss_enrich: {
            description: 'Enables TSS enrichment plot generation.',
            group: 'pipeline_parameter'
        } 
        cutadapt_param: {
            description: 'Parameters for cutadapt.',
            group: 'adapter_trimming',
            help: 'It is -e 0.1 -m 5 by default (err_rate=0.1, min_trim_len=5). You can define any parameters that cutadapt supports.'
        }
        auto_detect_adapter: {
            description: 'Auto-detect/trim adapter sequences.',
            group: 'adapter_trimming',
            help: 'Can detect/trim three types of adapter sequences. Illumina: AGATCGGAAGAGC, Nextera: CTGTCTCTTATA, smallRNA: TGGAATTCTCGG.',
            example: true
        }
        adapter: {
            description: 'Adapter for all FASTQs.',
            group: 'adapter_trimming',
            help: 'Define if all FASTQs have the same adapter sequence. Otherwise define adapter sequence for individual FASTQ in atac.adapters_repX_R1 and atac.adapters_repX_R2 instead. Use atac.auto_detect_adapter if you want to detect adapters automatically. If all of your FASTQs are already trimmed then leave all adapter-related parameters undefined/empty.'
        }
        adapters_rep1_R1: {
            description: 'Adapter sequences for read1 FASTQs to be merged for a biological replicate 1.',
            group: 'adapter_trimming',
            help: 'Make sure that they are consistent with read2 FASTQs (atac.adapters_rep1_R2). You can combine this with atac.auto_detect_adapter. Pipeline will auto-detect/trim adapter sequences for null entry in this list. e.g. ["AAGGCCTT", null, "AAGGCCTT"].'
        }
        adapters_rep1_R2: {
            description: 'Adapter sequences for read2 FASTQs to be merged for a biological replicate 1.',
            group: 'adapter_trimming',
            help: 'Make sure that they are consistent with read1 FASTQs (atac.adapters_rep1_R1).'
        }
        multimapping: {
            description: 'Number of multimappers.',
            group: 'alignment',
            help: 'It is 4 by default. Set it to 0 if your sample does not have multimappers.'
        }
        dup_marker: {
            description: 'Marker for duplicate reads. picard or sambamba.',
            group: 'alignment',
            help: 'picard for Picard MarkDuplicates or sambamba for sambamba markdup.'
        }
        no_dup_removal: {
            description: 'Disable removal of duplicate reads during filtering BAM.',
            group: 'alignment',
            help: 'Duplicate reads are filtererd out during filtering BAMs to gerenate NODUP_BAM. This flag will keep all duplicate reads in NODUP_BAM. This flag does not affect naming of NODUP_BAM. NODUP_BAM will still have .nodup. suffix in its filename.'
        }
        mapq_thresh: {
            description: 'Threshold for low MAPQ reads removal.',
            group: 'alignment',
            help: 'Low MAPQ reads are filtered out while filtering BAM.'
        }
        filter_chrs: {
            description: 'List of chromosomes to be filtered out while filtering BAM.',
            group: 'alignment',
            help: 'It is ["chrM", "MT"] by default. Therefore, mitochondrial reads will be filtered out while filtering. Make it empty if you want to keep all reads.'
        }
        subsample_reads: {
            description: 'Subsample reads. Shuffle and subsample reads.',
            group: 'alignment',
            help: 'This affects all downstream analyses after filtering BAM. (e.g. all TAG-ALIGN files, peak-calling). Reads will be shuffled only if actual number of reads in BAM exceeds this number.  0 means disabled.'
        }
        read_len: {
            description: 'Read length per biological replicate.',
            group: 'alignment',
            help: 'Pipeline can estimate read length from FASTQs. If you start pipeline from other types (BAM, NODUP_BAM, TA, ...) than FASTQ. Then provide this for some analyses that require read length (e.g. TSS enrichment plot).'
        }
        cap_num_peak: {
            description: 'Upper limit on the number of peaks.',
            group: 'peak_calling',
            help: 'Called peaks will be sorted in descending order of score and the number of peaks will be capped at this number by taking first N peaks.'
        }
        pval_thresh: {
            description: 'p-value Threshold for MACS2 peak caller.',
            group: 'peak_calling',
            help: 'macs2 callpeak -p'
        }
        smooth_win: {
            description: 'Size of smoothing windows for MACS2 peak caller.',
            group: 'peak_calling',
            help: 'This will be used for both generating MACS2 peaks/signal tracks.'
        }
        align_cpu: {
            description: 'Number of cores for task align.',
            group: 'resource_parameter',
            help: 'Task align merges/crops/maps FASTQs.'
        }
        align_mem_factor: {
            description: 'Multiplication factor to determine memory required for task align.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of FASTQs to determine required memory of instance (GCP/AWS) or job (HPCs).'
        }
        align_time_hr: {
            description: 'Walltime (h) required for task align.',
            group: 'resource_parameter',
            help: 'This is for HPCs only. e.g. SLURM, SGE, ...'
        }
        align_disk_factor: {
            description: 'Multiplication factor to determine persistent disk size for task align.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of FASTQs to determine required disk size of instance on GCP/AWS.'
        }
        filter_cpu: {
            description: 'Number of cores for task filter.',
            group: 'resource_parameter',
            help: 'Task filter filters raw/unfilterd BAM to get filtered/deduped BAM.'
        }
        filter_mem_factor: {
            description: 'Multiplication factor to determine memory required for task filter.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of BAMs to determine required memory of instance (GCP/AWS) or job (HPCs).'
        }
        filter_time_hr: {
            description: 'Walltime (h) required for task filter.',
            group: 'resource_parameter',
            help: 'This is for HPCs only. e.g. SLURM, SGE, ...'
        }
        filter_disk_factor: {
            description: 'Multiplication factor to determine persistent disk size for task filter.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of BAMs to determine required disk size of instance on GCP/AWS.'
        }
        bam2ta_cpu: {
            description: 'Number of cores for task bam2ta.',
            group: 'resource_parameter',
            help: 'Task bam2ta converts filtered/deduped BAM in to TAG-ALIGN (6-col BED) format.'
        }
        bam2ta_mem_factor: {
            description: 'Multiplication factor to determine memory required for task bam2ta.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of filtered BAMs to determine required memory of instance (GCP/AWS) or job (HPCs).'
        }
        bam2ta_time_hr: {
            description: 'Walltime (h) required for task bam2ta.',
            group: 'resource_parameter',
            help: 'This is for HPCs only. e.g. SLURM, SGE, ...'
        }
        bam2ta_disk_factor: {
            description: 'Multiplication factor to determine persistent disk size for task bam2ta.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of filtered BAMs to determine required disk size of instance on GCP/AWS.'
        }
        spr_mem_factor: {
            description: 'Multiplication factor to determine memory required for task spr.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of filtered BAMs to determine required memory of instance (GCP/AWS) or job (HPCs).'
        }
        spr_disk_factor: {
            description: 'Multiplication factor to determine persistent disk size for task spr.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of filtered BAMs to determine required disk size of instance on GCP/AWS.'
        }
        jsd_cpu: {
            description: 'Number of cores for task jsd.',
            group: 'resource_parameter',
            help: 'Task jsd plots Jensen-Shannon distance and metrics related to it.'
        }
        jsd_mem_factor: {
            description: 'Multiplication factor to determine memory required for task jsd.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of filtered BAMs to determine required memory of instance (GCP/AWS) or job (HPCs).'
        }
        jsd_time_hr: {
            description: 'Walltime (h) required for task jsd.',
            group: 'resource_parameter',
            help: 'This is for HPCs only. e.g. SLURM, SGE, ...'
        }
        jsd_disk_factor: {
            description: 'Multiplication factor to determine persistent disk size for task jsd.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of filtered BAMs to determine required disk size of instance on GCP/AWS.'
        }
        call_peak_cpu: {
            description: 'Number of cores for task call_peak. MACS2 is single-thread. No more than 2 is required.',
            group: 'resource_parameter',
            help: 'Task call_peak call peaks on TAG-ALIGNs by using MACS2 peak caller.'
        }
        call_peak_mem_factor: {
            description: 'Multiplication factor to determine memory required for task call_peak.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of TAG-ALIGNs (BEDs) to determine required memory of instance (GCP/AWS) or job (HPCs).'
        }
        call_peak_time_hr: {
            description: 'Walltime (h) required for task call_peak.',
            group: 'resource_parameter',
            help: 'This is for HPCs only. e.g. SLURM, SGE, ...'
        }
        call_peak_disk_factor: {
            description: 'Multiplication factor to determine persistent disk size for task call_peak.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of TAG-ALIGNs (BEDs) to determine required disk size of instance on GCP/AWS.'
        }
        macs2_signal_track_mem_factor: {
            description: 'Multiplication factor to determine memory required for task macs2_signal_track.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of TAG-ALIGNs (BEDs) to determine required memory of instance (GCP/AWS) or job (HPCs).'
        }
        macs2_signal_track_time_hr: {
            description: 'Walltime (h) required for task macs2_signal_track.',
            group: 'resource_parameter',
            help: 'This is for HPCs only. e.g. SLURM, SGE, ...'
        }
        macs2_signal_track_disk_factor: {
            description: 'Multiplication factor to determine persistent disk size for task macs2_signal_track.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of TAG-ALIGNs (BEDs) to determine required disk size of instance on GCP/AWS.'
        }
        filter_picard_java_heap: {
            description: 'Maximum Java heap (java -Xmx) in task filter.',
            group: 'resource_parameter',
            help: 'Maximum memory for Picard tools MarkDuplicates. If not defined, 90% of filter task\'s memory will be used.'
        }
        fraglen_stat_picard_java_heap: {
            description: 'Maximum Java heap (java -Xmx) in task fraglen_stat_pe (for paired end replicate only).',
            group: 'resource_parameter',
            help: 'Maximum memory for Picard tools CollectInsertSizeMetrics. If not defined, 90% of fraglen_stat tasks\'s memory will be used.'
        }
        gc_bias_picard_java_heap: {
            description: 'Maximum Java heap (java -Xmx) in task gc_bias.',
            group: 'resource_parameter',
            help: 'Maximum memory for Picard tools CollectGcBiasMetrics. If not defined, 90% of gc_bias tasks\'s memory will be used.'
        }
    }
    RuntimeEnvironment runtime_environment = {
        'docker': docker, 'singularity': singularity, 'conda': conda
    }
    RuntimeEnvironment runtime_environment_spp = {
        'docker': docker, 'singularity': singularity, 'conda': conda_spp
    }
    RuntimeEnvironment runtime_environment_macs2 = {
        'docker': docker, 'singularity': singularity, 'conda': conda_macs2
    }
    RuntimeEnvironment runtime_environment_python2 = {
        'docker': docker, 'singularity': singularity, 'conda': conda_python2
    }

    #Initalize the type of aligner, pealk caller and peak type to be used
    String aligner = 'bowtie2'
    String peak_caller = 'macs2'
    String peak_type = 'narrowPeak'
    
    #Checks if genome is defined in the genome_tsv file
    #If true, then will call the read_genome_tsv task
    if ( defined(genome_tsv) ) {
        call read_genome_tsv { input:
            genome_tsv = genome_tsv,
            runtime_environment = runtime_environment
        }
    }

    #Initalizing values from reference genome into files
    #Select first will check if the first value is defined, if it is not then the read_genome_tsv value will be used
    File ref_fa_ = select_first([ref_fa, read_genome_tsv.ref_fa])
    File bowtie2_idx_tar_ = select_first([bowtie2_idx_tar, read_genome_tsv.bowtie2_idx_tar])
    File chrsz_ = select_first([chrsz, read_genome_tsv.chrsz])
    String gensz_ = select_first([gensz, read_genome_tsv.gensz])

    #Declaring optinal file variables based on if blacklist was defined
    #If it was not defined, then the read_genome_tsv value will be used
    File? blacklist1_ = if defined(blacklist) then blacklist
        else read_genome_tsv.blacklist
    File? blacklist2_ = if defined(blacklist2) then blacklist2
        else read_genome_tsv.blacklist2        
    # Merge both balcklist files into an 
    # If the array is longer than one, call the pool_ta task
    # two blacklists can have different number of columns (3 vs 6)
    # so we limit merged blacklist's columns to 3
    Array[File] blacklists = select_all([blacklist1_, blacklist2_])
    if ( length(blacklists) > 1 ) {
        call pool_ta as pool_blacklist { input:
            tas = blacklists,
            col = 3,
            runtime_environment = runtime_environment
        }
    }
    # Assignes a blacklist file based on the length of the blacklists array
    #If the length is greater than one, then the pooled blacklist will be used
    #If the length is 1 then the first file is used
    #If the length is 0, then the second blacklist file will be used
    File? blacklist_ = if length(blacklists) > 1 then pool_blacklist.ta_pooled
        else if length(blacklists) > 0 then blacklists[0]
        else blacklist2_
    String regex_bfilt_peak_chr_name_ = select_first([regex_bfilt_peak_chr_name, read_genome_tsv.regex_bfilt_peak_chr_name])
    String genome_name_ = select_first([genome_name, read_genome_tsv.genome_name, basename(chrsz_)])

    # Read additional annotation data
    File? tss_ = if defined(tss) then tss
        else read_genome_tsv.tss
    File? dnase_ = if defined(dnase) then dnase
        else read_genome_tsv.dnase
    File? prom_ = if defined(prom) then prom
        else read_genome_tsv.prom
    File? enh_ = if defined(enh) then enh
        else read_genome_tsv.enh
    File? reg2map_ = if defined(reg2map) then reg2map
        else read_genome_tsv.reg2map
    File? reg2map_bed_ = if defined(reg2map_bed) then reg2map_bed
        else read_genome_tsv.reg2map_bed
    File? roadmap_meta_ = if defined(roadmap_meta) then roadmap_meta
        else read_genome_tsv.roadmap_meta

    ### temporary variables (do not define these)
    String aligner_ = aligner
    String peak_caller_ = peak_caller
    String peak_type_ = peak_type
    Int cap_num_peak_ = cap_num_peak
    Int mapq_thresh_ = mapq_thresh

    # # Wrap the R1 FASTQ files into a nested array (Array[Array[File]])
    Array[Array[File]] fastqs_R1 = [fastqs_rep1_R1]
    # Same thing with R2 files
    Array[Array[File]] fastqs_R2 = [fastqs_rep1_R2]
    # temporary 2-dim adapters array [rep_id][merge_id]
    Array[Array[String]] adapters_R1 = [adapters_rep1_R1]
    Array[Array[String]] adapters_R2 = [adapters_rep1_R2]
    
    # Assume no replicates: only one sample throughout the pipeline
    Int num_rep = 1

    # Verifying that program is funtionally correctly
    if ( num_rep == 0 ) {
        call raise_exception as error_input_data  { input:
            msg = 'No FASTQ/BAM/TAG-ALIGN/PEAK defined in your input JSON. Check if your FASTQs are defined as "atac.fastqs_repX_RY". DO NOT MISS suffix _R1 even for single ended FASTQ.',
            runtime_environment = runtime_environment
        }
    }

    # Check if paired_end is defined; otherwise, use the global paired_end value.
    Boolean paired_end_ = if !defined(paired_end) then select_first([paired_end])
        else paired_end
    # Check if there is an input FASTQ file for R1 and R2, and if the BAM output is missing
    Boolean has_input_of_align = length(fastqs_R1) > 0
    Boolean has_output_of_align = length(bams) > 0
    
    # If there is an input FASTQ file and no output BAM file then call the align task
    if (has_input_of_align && !has_output_of_align) {
        call align {
            input:
                fastqs_R1 = fastqs_R1[0],  # Only using the first (and only) entry
                fastqs_R2 = fastqs_R2[0],  # Same for R2
                adapter = adapter,
                adapters_R1 = adapters_R1[0],  # Use first adapter for R1
                adapters_R2 = adapters_R2[0],  # Use first adapter for R2
                paired_end = paired_end_,
                auto_detect_adapter = auto_detect_adapter,
                cutadapt_param = cutadapt_param,
            
                aligner = aligner_,
                chrsz = chrsz_,
                multimapping = multimapping,
                idx_tar = bowtie2_idx_tar_,
            
                # Resource settings
                cpu = align_cpu,
                mem_factor = align_mem_factor,
                time_hr = align_time_hr,
                disk_factor = align_disk_factor,
                runtime_environment = runtime_environment
    }
}
        #If there is an output BAM file then assign it to file
        File? bam_ = if has_output_of_align then bams[i] else align.bam

        Boolean has_input_of_filter = has_output_of_align || defined(align.bam)
        Boolean has_output_of_filter = i<length(nodup_bams)
        #If there is an input BAM file which has not been filtered the call the filter task
        #If there is an output nodup BAM file then assign it to file
        if ( has_input_of_filter && !has_output_of_filter ) {
            call filter { input :
                bam = bam_,
                paired_end = paired_end_,
                dup_marker = dup_marker,
                mapq_thresh = mapq_thresh_,
                filter_chrs = filter_chrs,
                chrsz = chrsz_,
                no_dup_removal = no_dup_removal,
                multimapping = multimapping,

                cpu = filter_cpu,
                mem_factor = filter_mem_factor,
                picard_java_heap = filter_picard_java_heap,
                time_hr = filter_time_hr,
                disk_factor = filter_disk_factor,
                runtime_environment = runtime_environment
            }
        }
        File? nodup_bam_ = if has_output_of_filter then nodup_bams[i] else filter.nodup_bam

        Boolean has_input_of_bam2ta = has_output_of_filter || defined(filter.nodup_bam)
        Boolean has_output_of_bam2ta = i<length(tas)

        # If there is a input BAM file which has not been converted to TAG-ALIGN then call the bam2ta task
        if ( has_input_of_bam2ta && !has_output_of_bam2ta ) {
            call bam2ta { input :
                bam = nodup_bam_,
                disable_tn5_shift = if pipeline_type=='atac' then false else true,
                subsample = subsample_reads,
                paired_end = paired_end_,

                cpu = bam2ta_cpu,
                mem_factor = bam2ta_mem_factor,
                time_hr = bam2ta_time_hr,
                disk_factor = bam2ta_disk_factor,
                runtime_environment = runtime_environment
            }
        }
        #If it has an output TAG-ALIGN file then assign it to file
        File? ta_ = if has_output_of_bam2ta then tas[i] else bam2ta.ta

        Boolean has_input_of_xcor = has_output_of_align || defined(align.bam)
        # If there is an output of align and cross correlation analysis is enabled then call the following tasks
        # enable_xcor is set to false by default
        if ( has_input_of_xcor && enable_xcor ) {
            call filter as filter_no_dedup { input :
                bam = bam_,
                paired_end = paired_end_,
                dup_marker = dup_marker,
                mapq_thresh = mapq_thresh_,
                filter_chrs = filter_chrs,
                chrsz = chrsz_,
                no_dup_removal = true,
                multimapping = multimapping,

                cpu = filter_cpu,
                mem_factor = filter_mem_factor,
                picard_java_heap = filter_picard_java_heap,
                time_hr = filter_time_hr,
                disk_factor = filter_disk_factor,
                runtime_environment = runtime_environment
            }
            call bam2ta as bam2ta_no_dedup { input :
                bam = filter_no_dedup.nodup_bam,  # output name is nodup but it's not deduped
                disable_tn5_shift = if pipeline_type=='atac' then false else true,
                subsample = 0,
                paired_end = paired_end_,

                cpu = bam2ta_cpu,
                mem_factor = bam2ta_mem_factor,
                time_hr = bam2ta_time_hr,
                disk_factor = bam2ta_disk_factor,
                runtime_environment = runtime_environment
            }
        Boolean has_input_of_macs2_signal_track = has_output_of_bam2ta || defined(bam2ta.ta)
        # If there is an input TAG-ALIGN file then call macs2_signal_track task
        if ( has_input_of_macs2_signal_track ) {
            # generate count signal track
            call macs2_signal_track { input :
                ta = ta_,
                gensz = gensz_,
                chrsz = chrsz_,
                pval_thresh = pval_thresh,
                smooth_win = smooth_win,

                mem_factor = macs2_signal_track_mem_factor,
                disk_factor = macs2_signal_track_disk_factor,
                time_hr = macs2_signal_track_time_hr,

                runtime_environment = runtime_environment_macs2
            }
        }

        Boolean has_input_of_call_peak = has_output_of_bam2ta || defined(bam2ta.ta)
        Boolean has_output_of_call_peak = i<length(peaks)
        # Runs if there is an input TAG-ALIGN file, no output peak file, and if the program is not set to align only
        #Align only is set to false by default
        # Assign output of call_peak to peak_ if it is defined
        if ( has_input_of_call_peak && !has_output_of_call_peak && !align_only ) {
            # call peaks on tagalign
            call call_peak { input :
                peak_caller = peak_caller_,
                peak_type = peak_type_,
                ta = ta_,
                gensz = gensz_,
                chrsz = chrsz_,
                cap_num_peak = cap_num_peak_,
                pval_thresh = pval_thresh,
                smooth_win = smooth_win,
                blacklist = blacklist_,
                regex_bfilt_peak_chr_name = regex_bfilt_peak_chr_name_,

                cpu = call_peak_cpu,
                mem_factor = call_peak_mem_factor,
                disk_factor = call_peak_disk_factor,
                time_hr = call_peak_time_hr,

                runtime_environment = if peak_caller == 'spp' then runtime_environment_spp
                    else if peak_caller == 'macs2' then runtime_environment_macs2
                    else runtime_environment
            }
        }
        File? peak_ = if has_output_of_call_peak then peaks[i] else call_peak.peak

        Boolean has_input_of_spr = has_output_of_bam2ta || defined(bam2ta.ta)
        # If there is an input TAG-ALIGN file then call the spr task
        if ( has_input_of_spr && !align_only && !true_rep_only ) {
            call spr { input :
                ta = ta_,
                paired_end = paired_end_,
                mem_factor = spr_mem_factor,
                disk_factor = spr_disk_factor,
                runtime_environment = runtime_environment
            }
        }

        Boolean has_input_of_call_peak_pr1 = defined(spr.ta_pr1)
        Boolean has_output_of_call_peak_pr1 = i<length(peaks_pr1)
        # Calls if the first replicate assigned in spr task is defined
        # Runs call peak on the 1st pseudo replicated tagalign
        # Align only and true_rep_only are set to false by default
        # Assigns output to file
        if ( has_input_of_call_peak_pr1 && !has_output_of_call_peak_pr1 &&
            !align_only && !true_rep_only ) {
            # call peaks on 1st pseudo replicated tagalign 
            call call_peak as call_peak_pr1 { input :
                peak_caller = peak_caller_,
                peak_type = peak_type_,
                ta = spr.ta_pr1,
                gensz = gensz_,
                chrsz = chrsz_,
                cap_num_peak = cap_num_peak_,
                pval_thresh = pval_thresh,
                smooth_win = smooth_win,
                blacklist = blacklist_,
                regex_bfilt_peak_chr_name = regex_bfilt_peak_chr_name_,

                cpu = call_peak_cpu,
                mem_factor = call_peak_mem_factor,
                disk_factor = call_peak_disk_factor,
                time_hr = call_peak_time_hr,

                runtime_environment = if peak_caller == 'spp' then runtime_environment_spp
                    else if peak_caller == 'macs2' then runtime_environment_macs2
                    else runtime_environment
            }
        }
        File? peak_pr1_ = if has_output_of_call_peak_pr1 then peaks_pr1[i]
            else call_peak_pr1.peak

        Boolean has_input_of_call_peak_pr2 = defined(spr.ta_pr2)
        Boolean has_output_of_call_peak_pr2 = i<length(peaks_pr2)
        # Runs if the second replicate assigned in spr task is defined
        # Rins call peak on the 2nd pseudo replicated tagalign
        # Align only and true_rep_only are set to false by default
        # Assigns output to file
        if ( has_input_of_call_peak_pr2 && !has_output_of_call_peak_pr2 &&
            !align_only && !true_rep_only ) {
            # call peaks on 2nd pseudo replicated tagalign 
            call call_peak as call_peak_pr2 { input :
                peak_caller = peak_caller_,
                peak_type = peak_type_,
                ta = spr.ta_pr2,
                gensz = gensz_,
                chrsz = chrsz_,
                cap_num_peak = cap_num_peak_,
                pval_thresh = pval_thresh,
                smooth_win = smooth_win,
                blacklist = blacklist_,
                regex_bfilt_peak_chr_name = regex_bfilt_peak_chr_name_,

                cpu = call_peak_cpu,
                mem_factor = call_peak_mem_factor,
                disk_factor = call_peak_disk_factor,
                time_hr = call_peak_time_hr,

                runtime_environment = if peak_caller == 'spp' then runtime_environment_spp
                    else if peak_caller == 'macs2' then runtime_environment_macs2
                    else runtime_environment
            }
        }
        File? peak_pr2_ = if has_output_of_call_peak_pr2 then peaks_pr2[i]
            else call_peak_pr2.peak

        Boolean has_input_of_count_signal_track = has_output_of_bam2ta || defined(bam2ta.ta)
        # enable_count_signal_track is set to false by default
        if ( has_input_of_count_signal_track && enable_count_signal_track ) {
            # generate count signal track
            call count_signal_track { input :
                ta = ta_,
                chrsz = chrsz_,
                runtime_environment = runtime_environment
            }
        }
        # tasks factored out from ATAqC
        Boolean has_input_of_tss_enrich = defined(nodup_bam_) && defined(tss_) && (
            defined(align.read_len) || i<length(read_len) )
        # Runs if tss_enrich is enabled, there is output for nodup_bam, tss is defined, and a read length is defined
        # enable_tss_enrich is set to true by default
        if ( enable_tss_enrich && has_input_of_tss_enrich ) {
            call tss_enrich { input :
                read_len = if i<length(read_len) then read_len[i]
                    else align.read_len,
                nodup_bam = nodup_bam_,
                tss = tss_,
                chrsz = chrsz_,
                runtime_environment = runtime_environment_python2
            }
        }
        # Runs if fraglen_stat is enabled, nodup_bam is defined, and paired_end is true
        # enable_fraglen_stat is set to true by default
        if ( enable_fraglen_stat && paired_end_ && defined(nodup_bam_) ) {
            call fraglen_stat_pe { input :
                nodup_bam = nodup_bam_,
                picard_java_heap = fraglen_stat_picard_java_heap,
                runtime_environment = runtime_environment
            }
        }
        #Enables preseq is set to false by default
        if ( enable_preseq && defined(bam_) ) {
            call preseq { input :
                bam = bam_,
                paired_end = paired_end_,
                mem_factor = preseq_mem_factor,
                disk_factor = preseq_disk_factor,
                picard_java_heap = preseq_picard_java_heap,
                runtime_environment = runtime_environment
            }
        }

        # Runs if gc_bias is enabled, nodup_bam is defined, and ref_fa is defined
        # enable_gc_bias is set to true by default
        if ( enable_gc_bias && defined(nodup_bam_) && defined(ref_fa_) ) {
            call gc_bias { input :
                nodup_bam = nodup_bam_,
                ref_fa = ref_fa_,
                picard_java_heap = gc_bias_picard_java_heap,
                runtime_environment = runtime_environment
            }
        }
        # Runs if annot_enrich is enabled and parameters are defined
        if ( enable_annot_enrich && defined(ta_) && defined(blacklist_) && defined(dnase_) && defined(prom_) && defined(enh_) ) {
            call annot_enrich { input :
                ta = ta_,
                blacklist = blacklist_,
                dnase = dnase_,
                prom = prom_,
                enh = enh_,
                runtime_environment = runtime_environment
            }
        }
        # Runs if enable_compare_to_roadmap is enabled and parameters are defined
        if ( enable_compare_to_roadmap && defined(macs2_signal_track.pval_bw) &&
             defined(reg2map_) && defined(roadmap_meta_) &&
             ( defined(reg2map_bed_) || defined(dnase_) ) ) {
            call compare_signal_to_roadmap { input :
                pval_bw = macs2_signal_track.pval_bw,
                dnase = dnase_,
                reg2map_bed = reg2map_bed_,
                reg2map = reg2map_,
                roadmap_meta = roadmap_meta_,
                runtime_environment = runtime_environment
            }
        }
    }

        # Since there is only one replicate, directly proceed with pooling the tagaligns (TAs)
        Boolean has_input_of_pool_ta = length(select_all(ta_)) == 1  # Only check if there is input for one replicate
        if ( has_input_of_pool_ta ) {
            # Pool tagaligns from the single replicate
            call pool_ta { 
                input :
                    tas = ta_,  # Use the tag-alignments for the single replicate
                    prefix = 'rep1',  # Prefix for the single replicate
                    runtime_environment = runtime_environment
    }
}

        # Check if there is an input for pr1 and pool it for the single replicate
        Boolean has_input_of_pool_ta_pr1 = length(spr.ta_pr1) == 1  # Check if there is only one pr1 TA
        if ( has_input_of_pool_ta_pr1 && !align_only && !true_rep_only ) {
        # Pool tagaligns from pseudo replicate 1 for the single replicate
            call pool_ta as pool_ta_pr1 { 
                input :
                    tas = spr.ta_pr1,  # Use the TA for pr1 for the single replicate
                    prefix = 'rep-pr1',  # Prefix for pseudo replicate 1
                    runtime_environment = runtime_environment
    }
}
        # Check if there is an input for pr2 and pool it for the single replicate
        Boolean has_input_of_pool_ta_pr2 = length(spr.ta_pr2) == 1  # Check if there is only one pr2 TA
        if ( has_input_of_pool_ta_pr2 && !align_only && !true_rep_only ) {
            # Pool tagaligns from pseudo replicate 2 for the single replicate
            call pool_ta as pool_ta_pr2 { 
                input :
                    tas = spr.ta_pr2,  # Use the TA for pr2 for the single replicate
                    prefix = 'rep-pr2',  # Prefix for pseudo replicate 2
                    runtime_environment = runtime_environment
    }
}

  Boolean has_input_of_jsd = defined(blacklist_) && length(select_all(nodup_bam_)) == 1
  #Will run if enable_jsd is enabled
if (has_input_of_jsd && enable_jsd) {
    # Fingerprint and JSD plot (even for one replicate)
    call jsd { input :
        nodup_bams = nodup_bam_,
        blacklist = blacklist_,
        mapq_thresh = mapq_thresh_,

        cpu = jsd_cpu,
        mem_factor = jsd_mem_factor,
        time_hr = jsd_time_hr,
        disk_factor = jsd_disk_factor,
        runtime_environment = runtime_environment
    }
}
# Reproducibility QC for overlapping peaks
#Num rep is set to 1 by default
 if ( !align_only && !true_rep_only && num_rep > 0 ) {
    # Reproducibility QC for overlapping peaks
      call reproducibility as reproducibility_overlap { input :
          prefix = 'overlap',
          peaks = select_all(overlap.bfilt_overlap_peak),
          peak_type = peak_type_,
          chrsz = chrsz_,
          runtime_environment = runtime_environment
    }
}
    # Generate final QC report and JSON
      call qc_report { input :
        pipeline_ver = pipeline_ver,
        title = title,
        description = description,
        genome = genome_name_,
        multimapping = multimapping,
        paired_ends = paired_end_,
        pipeline_type = pipeline_type,
        aligner = aligner_,
        no_dup_removal = no_dup_removal,
        peak_caller = peak_caller_,
        cap_num_peak = cap_num_peak_,
        pval_thresh = pval_thresh,
        xcor_subsample_reads = xcor_subsample_reads,

        samstat_qcs = select_all(align.samstat_qc),
        nodup_samstat_qcs = select_all(filter.samstat_qc),

        dup_qcs = select_all(filter.dup_qc),
        lib_complexity_qcs = select_all(filter.lib_complexity_qc),

        jsd_plot = jsd.plot,
        jsd_qcs = if defined(jsd.jsd_qcs) then select_first([jsd.jsd_qcs]) else [],

        frip_qcs = select_all(call_peak.frip_qc),
        frip_qcs_pr1 = select_all(call_peak_pr1.frip_qc),
        frip_qcs_pr2 = select_all(call_peak_pr2.frip_qc),

        frip_overlap_qcs = select_all(overlap.frip_qc),
        frip_overlap_qcs_pr = if defined(overlap_pr.frip_qc) then select_first([overlap_pr.frip_qc]) else [],
        frip_overlap_qc_ppr = overlap_ppr.frip_qc,
        overlap_reproducibility_qc = reproducibility_overlap.reproducibility_qc,

        annot_enrich_qcs = select_all(annot_enrich.annot_enrich_qc),
        tss_enrich_qcs = select_all(tss_enrich.tss_enrich_qc),
        tss_large_plots = select_all(tss_enrich.tss_large_plot),
        roadmap_compare_plots = select_all(compare_signal_to_roadmap.roadmap_compare_plot),
        fraglen_dist_plots = select_all(fraglen_stat_pe.fraglen_dist_plot),
        fraglen_nucleosomal_qcs = select_all(fraglen_stat_pe.nucleosomal_qc),
        gc_plots = select_all(gc_bias.gc_plot),
        preseq_plots = select_all(preseq.preseq_plot),
        picard_est_lib_size_qcs = select_all(preseq.picard_est_lib_size_qc),

        peak_region_size_qcs = select_all(call_peak.peak_region_size_qc),
        peak_region_size_plots = select_all(call_peak.peak_region_size_plot),
        num_peak_qcs = select_all(call_peak.num_peak_qc),

        overlap_opt_peak_region_size_qc = reproducibility_overlap.peak_region_size_qc,
        overlap_opt_peak_region_size_plot = reproducibility_overlap.peak_region_size_plot,
        overlap_opt_num_peak_qc = reproducibility_overlap.num_peak_qc,

        runtime_environment = runtime_environment
    }

    output {
        File report = qc_report.report
        File qc_json = qc_report.qc_json
        Boolean qc_json_ref_match = qc_report.qc_json_ref_match
    }
}

task align {
    input {        
        # for task trim_adapter
        Array[File] fastqs_R1         # [merge_id]
        Array[File] fastqs_R2

        String? adapter     # adapter for all fastqs,
                            #    this will override individual adapters in adapters_R1/R2
        Array[String] adapters_R1
        Array[String] adapters_R2
        Boolean paired_end
        Boolean auto_detect_adapter
        String cutadapt_param

        # for task align
        String aligner
        File chrsz            # 2-col chromosome sizes file
        File idx_tar        # reference index tar or tar.gz
        Int multimapping

        # resource
        Int cpu
        Float mem_factor
        Int time_hr
        Float disk_factor

        RuntimeEnvironment runtime_environment
    }
    Float input_file_size_gb = size(fastqs_R1, "G") + size(fastqs_R2, "G")
    Float mem_gb = 5.0 + size(idx_tar, "G") + mem_factor * input_file_size_gb
    Float samtools_mem_gb = 0.8 * mem_gb
    Int disk_gb = round(40.0 + disk_factor * input_file_size_gb)

    # tmp vars for task trim_adapter
    Array[Array[File]] tmp_fastqs = if paired_end then transpose([fastqs_R1, fastqs_R2])
                else transpose([fastqs_R1])
    Array[Array[String]] tmp_adapters = if paired_end then transpose([adapters_R1, adapters_R2])
                else transpose([adapters_R1])
 #shell script (inside a WDL command block) that orchestrates the steps for trimming adapters, aligning sequences, and processing alignment results.
    command {
        set -e

        # check if pipeline dependencies can be found
        if [[ -z "$(which encode_task_trim_adapter.py 2> /dev/null || true)" ]]
        then
          echo -e "\n* Error: pipeline environment (docker, singularity or conda) not found." 1>&2
          exit 3
        fi

        # trim adapter
        python3 $(which encode_task_trim_adapter.py) \
            ${write_tsv(tmp_fastqs)} \
            ${'--adapter ' + adapter} \
            --adapters ${write_tsv(tmp_adapters)} \
            ${if paired_end then '--paired-end' else ''} \
            ${if auto_detect_adapter then '--auto-detect-adapter' else ''} \
            --cutadapt-param ' ${cutadapt_param}' \
            ${'--nth ' + cpu}

        # align on trimmed/merged fastqs
        if [ '${aligner}' == 'bowtie2' ]; then
            python3 $(which encode_task_bowtie2.py) \
                ${idx_tar} \
                R1/*.fastq.gz \
                ${if paired_end then 'R2/*.fastq.gz' else ''} \
                ${if paired_end then '--paired-end' else ''} \
                ${'--multimapping ' + multimapping} \
                ${'--mem-gb ' + samtools_mem_gb} \
                ${'--nth ' + cpu}
        fi

        python3 $(which encode_task_post_align.py) \
            R1/*.fastq.gz $(ls *.bam) \
            ${'--mito-chr-name ' + mito_chr_name} \
            ${'--chrsz ' + chrsz} \
            ${'--mem-gb ' + samtools_mem_gb} \
            ${'--nth ' + cpu}
        rm -rf R1 R2
    }
    output {
        File bam = glob('*.bam')[0]
        File bai = glob('*.bai')[0]
        File samstat_qc = glob('*.samstats.qc')[0]
        File non_mito_samstat_qc = glob('non_mito/*.samstats.qc')[0]
        File read_len_log = glob('*.read_length.txt')[0]
        Int read_len = read_int(read_len_log)
    }
    runtime {
        cpu : cpu
        memory : '${mem_gb} GB'
        time : time_hr
        disks : 'local-disk ${disk_gb} SSD'
        preemptible: 0

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}
# filters a BAM file by removing:
#Low-quality alignments (based on MAPQ score)
#Specific chromosomes 
#Duplicate reads (optional, depending on no_dup_removal)
#Optionally processes paired-end or single-end data
task filter {
    input {
        File? bam
        Boolean paired_end
        Int multimapping
        String dup_marker             # picard.jar MarkDuplicates (picard) or 
                                    # sambamba markdup (sambamba)
        Int mapq_thresh                # threshold for low MAPQ reads removal
        Array[String] filter_chrs     # chrs to be removed from final (nodup/filt) BAM
        File chrsz                    # 2-col chromosome sizes file
        Boolean no_dup_removal         # no dupe reads removal when filtering BAM

        Int cpu
        Float mem_factor
        String? picard_java_heap
        Int time_hr
        Float disk_factor
        # runtime environment
        RuntimeEnvironment runtime_environment
    }
    #lets the pipeline adjust resources based on file size, 
    #avoiding crashes on large files and saving resources on small ones.
    Float input_file_size_gb = size(bam, "G")
    Float picard_java_heap_factor = 0.9
    Float mem_gb = 6.0 + mem_factor * input_file_size_gb
    Float samtools_mem_gb = 0.8 * mem_gb
    Int disk_gb = round(20.0 + disk_factor * input_file_size_gb)
    #builds a shell command to run encode_task_filter.py, 
    #configuring it based on input BAM properties and resource parameters
    command {
        set -e
        python3 $(which encode_task_filter.py) \
            ${bam} \
            ${if paired_end then '--paired-end' else ''} \
            ${'--multimapping ' + multimapping} \
            ${'--dup-marker ' + dup_marker} \
            ${'--mapq-thresh ' + mapq_thresh} \
            --filter-chrs ${sep=' ' filter_chrs} \
            ${'--chrsz ' + chrsz} \
            ${if no_dup_removal then '--no-dup-removal' else ''} \
            ${'--mem-gb ' + samtools_mem_gb} \
            ${'--nth ' + cpu} \
            ${'--picard-java-heap ' + if defined(picard_java_heap) then picard_java_heap else (round(mem_gb * picard_java_heap_factor) + 'G')}
    }
    #collects the key output files (BAM, index, and QC reports) produced by the filtering process,
    #making them available for downstream tasks in the workflow
    output {
        File nodup_bam = glob('*.bam')[0]
        File nodup_bai = glob('*.bai')[0]
        File samstat_qc = glob('*.samstats.qc')[0]
        File dup_qc = glob('*.dup.qc')[0]
        File lib_complexity_qc = glob('*.lib_complexity.qc')[0]
    }
    #specifies the computational resources and environment required to run the filter task
    runtime {
        cpu : cpu
        memory : '${mem_gb} GB'
        time : time_hr
        disks : 'local-disk ${disk_gb} SSD'

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}
#converts a BAM file into a TAGALIGN file, which is a simplified format often used in 
#ATAC-seq for peak calling
task bam2ta {
    input {
        File? bam
        Boolean paired_end
        Boolean disable_tn5_shift     # no tn5 shifting (it's for dnase-seq)
        Int subsample                 # number of reads to subsample TAGALIGN
                                    # this affects all downstream analysis
        Int cpu
        Float mem_factor
        Int time_hr
        Float disk_factor
        # runtime environment
        RuntimeEnvironment runtime_environment
    }
    Float input_file_size_gb = size(bam, "G")
    Float mem_gb = 4.0 + mem_factor * input_file_size_gb
    Float samtools_mem_gb = 0.8 * mem_gb
    Int disk_gb = round(20.0 + disk_factor * input_file_size_gb)

    command {
        set -e
        python3 $(which encode_task_bam2ta.py) \
            ${bam} \
            ${if paired_end then '--paired-end' else ''} \
            ${if disable_tn5_shift then '--disable-tn5-shift' else ''} \
            ${'--subsample ' + subsample} \
            ${'--mem-gb ' + samtools_mem_gb} \
            ${'--nth ' + cpu}
    }
    output {
        File ta = glob('*.tagAlign.gz')[0]
    }
    runtime {
        cpu : cpu
        memory : '${mem_gb} GB'
        time : time_hr
        disks : 'local-disk ${disk_gb} SSD'

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}
#erforms Self Pseudo-Replication (SPR) on a TAGALIGN (.ta) file. This is a common step in ATAC-seq quality control,
for estimating reproducibility between pseudo-replicates.
task spr {
    input {
        File? ta
        Boolean paired_end
        Int pseudoreplication_random_seed

        Float mem_factor
        Float disk_factor

        # runtime environment
        RuntimeEnvironment runtime_environment
    }
    Float input_file_size_gb = size(ta, "G")
    Float mem_gb = 4.0 + mem_factor * input_file_size_gb
    Int disk_gb = round(20.0 + disk_factor * input_file_size_gb)

    command {
        set -e
        python3 $(which encode_task_spr.py) \
            ${ta} \
            ${if paired_end then '--paired-end' else ''}
    }
    output {
        File ta_pr1 = glob('*.pr1.tagAlign.gz')[0]
        File ta_pr2 = glob('*.pr2.tagAlign.gz')[0]
    }
    runtime {
        cpu : 1
        memory : '${mem_gb} GB'
        time : 4
        disks : 'local-disk ${disk_gb} SSD'

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}
#merges multiple TAGALIGN (.ta) files into a single pooled TAGALIGN file
to increase signal depth or combine replicates.
task pool_ta {
    input {
        Array[File?] tas     # TAG-ALIGNs to be merged
        Int? col             # number of columns in pooled TA
        String? prefix         # basename prefix

        # runtime environment
        RuntimeEnvironment runtime_environment
    }
    command {
        set -e
        python3 $(which encode_task_pool_ta.py) \
            ${sep=' ' select_all(tas)} \
            ${'--prefix ' + prefix} \
            ${'--col ' + col}
    }
    output {
        File ta_pooled = glob('*.tagAlign.gz')[0]
    }
    runtime {
        cpu : 1
        memory : '8 GB'
        time : 4
        disks : 'local-disk 100 SSD'

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}
#generates a signal track from a TAGALIGN file by counting the number of mapped reads per genomic region (e.g., per base pair)
and using a chromosome sizes file to ensure the counts are scaled appropriately.
task count_signal_track {
    input {
        File? ta             # tag-align
        File chrsz            # 2-col chromosome sizes file

        # runtime environment
        RuntimeEnvironment runtime_environment
    }
    Float mem_gb = 8.0
    command {
        set -e
        python3 $(which encode_task_count_signal_track.py) \
            ${ta} \
            ${'--chrsz ' + chrsz} \
            ${'--mem-gb ' + mem_gb}
    }
    output {
        File pos_bw = glob('*.positive.bigwig')[0]
        File neg_bw = glob('*.negative.bigwig')[0]
    }
    runtime {
        cpu : 1
        memory : '${mem_gb} GB'
        time : 4
        disks : 'local-disk 50 SSD'

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}
#calling peaks from a TAGALIGN file from ATAC-seq data using a peak calling algorithm such as MACS2
task call_peak {
    input {
        String peak_caller
        String peak_type

        File? ta
        String gensz        # Genome size (sum of entries in 2nd column of 
                            # chr. sizes file, or hs for human, ms for mouse)
        File chrsz            # 2-col chromosome sizes file
        Int cap_num_peak    # cap number of raw peaks called from MACS2
        Float pval_thresh      # p.value threshold
        Int smooth_win         # size of smoothing window
        File? blacklist     # blacklist BED to filter raw peaks
        String? regex_bfilt_peak_chr_name

        Int cpu
        Float mem_factor
        Int time_hr
        Float disk_factor

        # runtime environment
        RuntimeEnvironment runtime_environment
    }
    Float input_file_size_gb = size(ta, "G")
    Float mem_gb = 4.0 + mem_factor * input_file_size_gb
    Int disk_gb = round(20.0 + disk_factor * input_file_size_gb)

    command {
        set -e

        if [ '${peak_caller}' == 'macs2' ]; then
            python3 $(which encode_task_macs2_atac.py) \
                ${ta} \
                ${'--gensz ' + gensz} \
                ${'--chrsz ' + chrsz} \
                ${'--cap-num-peak ' + cap_num_peak} \
                ${'--pval-thresh '+ pval_thresh} \
                ${'--smooth-win '+ smooth_win} \
                ${'--mem-gb ' + mem_gb}
        fi

        python3 $(which encode_task_post_call_peak_atac.py) \
            $(ls *Peak.gz) \
            ${'--ta ' + ta} \
            ${'--regex-bfilt-peak-chr-name \'' + regex_bfilt_peak_chr_name + '\''} \
            ${'--chrsz ' + chrsz} \
            ${'--peak-type ' + peak_type} \
            ${'--blacklist ' + blacklist}
    }
    output {
        File peak = glob('*[!.][!b][!f][!i][!l][!t].'+peak_type+'.gz')[0]
        # generated by post_call_peak py
        File bfilt_peak = glob('*.bfilt.'+peak_type+'.gz')[0]
        File bfilt_peak_bb = glob('*.bfilt.'+peak_type+'.bb')[0]
        File bfilt_peak_starch = glob('*.bfilt.'+peak_type+'.starch')[0]
        File bfilt_peak_hammock = glob('*.bfilt.'+peak_type+'.hammock.gz*')[0]
        File bfilt_peak_hammock_tbi = glob('*.bfilt.'+peak_type+'.hammock.gz*')[1]
        File frip_qc = glob('*.frip.qc')[0]
        File peak_region_size_qc = glob('*.peak_region_size.qc')[0]
        File peak_region_size_plot = glob('*.peak_region_size.png')[0]
        File num_peak_qc = glob('*.num_peak.qc')[0]
    }
    runtime {
        cpu : if peak_caller == 'macs2' then 1 else cpu
        memory : '${mem_gb} GB'
        time : time_hr
        disks : 'local-disk ${disk_gb} SSD'
        preemptible: 0

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}
#generates a signal track from a TAGALIGN file using the MACS2 peak calling tool, 
#and incorporates smoothing and statistical filtering. 
#This task would produce a signal track that can be visualized
task macs2_signal_track {
    input {
        File? ta
        String gensz        # Genome size (sum of entries in 2nd column of 
                            # chr. sizes file, or hs for human, ms for mouse)
        File chrsz            # 2-col chromosome sizes file
        Float pval_thresh      # p.value threshold
        Int smooth_win         # size of smoothing window

        Float mem_factor
        Int time_hr
        Float disk_factor

        # runtime environment
        RuntimeEnvironment runtime_environment
    }
    Float input_file_size_gb = size(ta, "G")
    Float mem_gb = 4.0 + mem_factor * input_file_size_gb
    Int disk_gb = round(20.0 + disk_factor * input_file_size_gb)

    command {
        set -e
        python3 $(which encode_task_macs2_signal_track_atac.py) \
            ${ta} \
            ${'--gensz '+ gensz} \
            ${'--chrsz ' + chrsz} \
            ${'--pval-thresh '+ pval_thresh} \
            ${'--smooth-win '+ smooth_win} \
            ${'--mem-gb ' + mem_gb}
    }
    output {
        File pval_bw = glob('*.pval.signal.bigwig')[0]
        File fc_bw = glob('*.fc.signal.bigwig')[0]
    }
    runtime {
        cpu : 1
        memory : '${mem_gb} GB'
        time : time_hr
        disks : 'local-disk ${disk_gb} SSD'
        preemptible: 0

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}
    input {
        String prefix         # prefix for overlap output file
        File? peak1
        File? peak2
        File? peak_pooled
        File? blacklist     # blacklist BED to filter raw peaks
        String regex_bfilt_peak_chr_name
        File? ta        # to calculate FRiP
        File chrsz            # 2-col chromosome sizes file
        String peak_type

        # runtime environment
        RuntimeEnvironment runtime_environment
}

    command {
        set -e
        touch null 
        python3 $(which encode_task_overlap.py) \
            ${peak1} ${peak2} ${peak_pooled} \
            ${'--prefix ' + prefix} \
            ${'--peak-type ' + peak_type} \
            ${'--chrsz ' + chrsz} \
            ${'--blacklist '+ blacklist} \
            --nonamecheck \
            ${'--regex-bfilt-peak-chr-name \'' + regex_bfilt_peak_chr_name + '\''} \
            ${'--ta ' + ta}
    }
    output {
        File overlap_peak = glob('*[!.][!b][!f][!i][!l][!t].'+peak_type+'.gz')[0]
        File bfilt_overlap_peak = glob('*.bfilt.'+peak_type+'.gz')[0]
        File bfilt_overlap_peak_bb = glob('*.bfilt.'+peak_type+'.bb')[0]
        File bfilt_overlap_peak_starch = glob('*.bfilt.'+peak_type+'.starch')[0]
        File bfilt_overlap_peak_hammock = glob('*.bfilt.'+peak_type+'.hammock.gz*')[0]
        File bfilt_overlap_peak_hammock_tbi = glob('*.bfilt.'+peak_type+'.hammock.gz*')[1]
        File frip_qc = if defined(ta) then glob('*.frip.qc')[0] else glob('null')[0]
    }
    runtime {
        cpu : 1
        memory : '4 GB'
        time : 4
        disks : 'local-disk 50 SSD'

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}
#calculates TSS (Transcription Start Site) enrichment, a common quality control (QC) metric for ATAC-seq
task tss_enrich {
    # based on metaseq, which is still in python2
    # python2 environment is required for this task
    input {
        Int? read_len
        File? nodup_bam
        File? tss
        File chrsz

        # runtime environment
        RuntimeEnvironment runtime_environment
    }
    command {
        set -e
        python2 $(which encode_task_tss_enrich.py) \
            ${'--read-len ' + read_len} \
            ${'--nodup-bam ' + nodup_bam} \
            ${'--chrsz ' + chrsz} \
            ${'--tss ' + tss}
    }
    output {
        File tss_plot = glob('*.tss_enrich.png')[0]
        File tss_large_plot = glob('*.large_tss_enrich.png')[0]
        File tss_enrich_qc = glob('*.tss_enrich.qc')[0]
        Float tss_enrich = read_float(tss_enrich_qc)
    }
    runtime {
        cpu : 1
        memory : '8 GB'
        time : 4
        disks : 'local-disk 150 SSD'

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}
#calculates fragment length statistics for paired-end (PE) sequencing data, using a deduplicated BAM file 
#quality control that depends on fragment size distributions.
task fraglen_stat_pe {
    # for PE only
    input {
        File? nodup_bam
        String? picard_java_heap

        # runtime environment
        RuntimeEnvironment runtime_environment
    }
    Float input_file_size_gb = size(nodup_bam, "G")
    Float mem_gb = 8.0
    Float picard_java_heap_factor = 0.9

    command {
        set -e
        python3 $(which encode_task_fraglen_stat_pe.py) \
            ${'--nodup-bam ' + nodup_bam} \
            ${'--picard-java-heap ' + if defined(picard_java_heap) then picard_java_heap else (round(mem_gb * picard_java_heap_factor) + 'G')}
    }
    output {
        File nucleosomal_qc = glob('*nucleosomal.qc')[0]
        File fraglen_dist_plot = glob('*fraglen_dist.png')[0]
    }
    runtime {
        cpu : 1
        memory : '${mem_gb} GB'
        time : 6
        disks : 'local-disk 150 SSD'

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}
#analyzes GC bias in sequencing data using a deduplicated BAM file (nodup_bam) and a reference genome FASTA file
task gc_bias {
    input {
        File? nodup_bam
        File ref_fa

        String? picard_java_heap

        # runtime environment
        RuntimeEnvironment runtime_environment
    }
    Float mem_factor = 0.3
    Float input_file_size_gb = size(nodup_bam, "G")
    Float mem_gb = 4.0 + mem_factor * input_file_size_gb
    Float picard_java_heap_factor = 0.9

    command {
        set -e
        python3 $(which encode_task_gc_bias.py) \
            ${'--nodup-bam ' + nodup_bam} \
            ${'--ref-fa ' + ref_fa} \
            ${'--picard-java-heap ' + if defined(picard_java_heap) then picard_java_heap else (round(mem_gb * picard_java_heap_factor) + 'G')}
    }
    output {
        File gc_plot = glob('*.gc_plot.png')[0]
        File gc_log = glob('*.gc.txt')[0]
    }
    runtime {
        cpu : 1
        memory : '${mem_gb} GB'
        time : 6
        disks : 'local-disk 250 SSD'

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}
#generates a comprehensive quality control (QC) report
#that summarizes the performance and quality metrics of the entire sequencing pipeline
task qc_report {
    input {
        String pipeline_ver
        String title
        String description
        String? genome
        # workflow params
        Int multimapping
        Array[Boolean] paired_ends
        String pipeline_type
        String aligner
        Boolean no_dup_removal
        String peak_caller
        Int cap_num_peak
        Float idr_thresh
        Float pval_thresh
        Int xcor_subsample_reads
        # QCs
        Array[File] frac_mito_qcs
        Array[File] samstat_qcs
        Array[File] nodup_samstat_qcs
        Array[File] dup_qcs
        Array[File] lib_complexity_qcs
        Array[File] xcor_plots
        Array[File] xcor_scores
        File? jsd_plot
        Array[File] jsd_qc
        Array[File] frip_qcs
        Array[File] frip_qcs_pr1
        Array[File] frip_qcs_pr2
        File? frip_qc_pooled
        File? frip_qc_ppr1
        File? frip_qc_ppr2
        Array[File] frip_overlap_qcs
        Array[File] frip_overlap_qcs_pr
        File? frip_overlap_qc_ppr
        File? overlap_reproducibility_qc

        Array[File] annot_enrich_qcs
        Array[File] tss_enrich_qcs
        Array[File] tss_large_plots
        Array[File] roadmap_compare_plots
        Array[File] fraglen_dist_plots
        Array[File] fraglen_nucleosomal_qcs
        Array[File] gc_plots
        Array[File] preseq_plots
        Array[File] picard_est_lib_size_qcs

        Array[File] peak_region_size_qcs
        Array[File] peak_region_size_plots
        Array[File] num_peak_qcs

        File? qc_json_ref

        # runtime environment
        RuntimeEnvironment runtime_environment
    }
    command {
        set -e
        python3 $(which encode_task_qc_report.py) \
            --pipeline-prefix atac \
            ${'--pipeline-ver ' + pipeline_ver} \
            ${"--title '" + sub(title,"'","_") + "'"} \
            ${"--desc '" + sub(description,"'","_") + "'"} \
            ${'--genome ' + genome} \
            ${'--multimapping ' + multimapping} \
            --paired-ends ${sep=' ' paired_ends} \
            --pipeline-type ${pipeline_type} \
            --aligner ${aligner} \
            ${if (no_dup_removal) then '--no-dup-removal ' else ''} \
            --peak-caller ${peak_caller} \
            ${'--cap-num-peak ' + cap_num_peak} \
            --pval-thresh ${pval_thresh} \
            --xcor-subsample-reads ${xcor_subsample_reads} \
            --samstat-qcs ${sep='_:_' samstat_qcs} \
            --nodup-samstat-qcs ${sep='_:_' nodup_samstat_qcs} \
            --dup-qcs ${sep='_:_' dup_qcs} \
            --lib-complexity-qcs ${sep='_:_' lib_complexity_qcs} \
            --xcor-plots ${sep='_:_' xcor_plots} \
            --xcor-scores ${sep='_:_' xcor_scores} \
            ${'--jsd-plot ' + jsd_plot} \
            --jsd-qcs ${sep='_:_' jsd_qcs} \
            --frip-qcs ${sep='_:_' frip_qcs} \
            --frip-qcs-pr1 ${sep='_:_' frip_qcs_pr1} \
            --frip-qcs-pr2 ${sep='_:_' frip_qcs_pr2} \
            ${'--frip-qc-ppr1 ' + frip_qc_ppr1} \
            ${'--frip-qc-ppr2 ' + frip_qc_ppr2} \
            --frip-overlap-qcs ${sep='_:_' frip_overlap_qcs} \
            --frip-overlap-qcs-pr ${sep='_:_' frip_overlap_qcs_pr} \
            --annot-enrich-qcs ${sep='_:_' annot_enrich_qcs} \
            --tss-enrich-qcs ${sep='_:_' tss_enrich_qcs} \
            --tss-large-plots ${sep='_:_' tss_large_plots} \
            --roadmap-compare-plots ${sep='_:_' roadmap_compare_plots} \
            --fraglen-dist-plots ${sep='_:_' fraglen_dist_plots} \
            --fraglen-nucleosomal-qcs ${sep='_:_' fraglen_nucleosomal_qcs} \
            --gc-plots ${sep='_:_' gc_plots} \
            --preseq-plots ${sep='_:_' preseq_plots} \
            --picard-est-lib-size-qcs ${sep='_:_' picard_est_lib_size_qcs} \
            --peak-region-size-qcs ${sep='_:_' peak_region_size_qcs} \
            --peak-region-size-plots ${sep='_:_' peak_region_size_plots} \
            --num-peak-qcs ${sep='_:_' num_peak_qcs} \
            --out-qc-html qc.html \
            --out-qc-json qc.json \
            ${'--qc-json-ref ' + qc_json_ref}
    }
    output {
        File report = glob('*qc.html')[0]
        File qc_json = glob('*qc.json')[0]
        Boolean qc_json_ref_match = read_string('qc_json_ref_match.txt')=='True'
    }
    runtime {
        cpu : 1
        memory : '4 GB'
        time : 4
        disks : 'local-disk 50 SSD'

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}
#eads and processes a genome configuration TSV file
task read_genome_tsv {
    input {
        File? genome_tsv
        String? null_s

        RuntimeEnvironment runtime_environment
    }
    command <<<
        echo "$(basename ~{genome_tsv})" > genome_name
        # create empty files for all entries
        touch ref_fa bowtie2_idx_tar chrsz gensz blacklist blacklist2
        touch bowtie2_mito_idx_tar
        touch tss tss_enrich # for backward compatibility
        touch dnase prom enh reg2map reg2map_bed roadmap_meta
        touch regex_bfilt_peak_chr_name

        python <<CODE
        import os
        with open('~{genome_tsv}','r') as fp:
            for line in fp:
                arr = line.strip('\n').split('\t')
                if arr:
                    key, val = arr
                    with open(key,'w') as fp2:
                        fp2.write(val)
        CODE
    >>>
    output {
        String? genome_name = read_string('genome_name')
        String? ref_fa = if size('ref_fa')==0 then null_s else read_string('ref_fa')
        String? ref_mito_fa = if size('ref_mito_fa')==0 then null_s else read_string('ref_mito_fa')
        String? bowtie2_idx_tar = if size('bowtie2_idx_tar')==0 then null_s else read_string('bowtie2_idx_tar')
        String? bowtie2_mito_idx_tar = if size('bowtie2_mito_idx_tar')==0 then null_s else read_string('bowtie2_mito_idx_tar')
        String? chrsz = if size('chrsz')==0 then null_s else read_string('chrsz')
        String? gensz = if size('gensz')==0 then null_s else read_string('gensz')
        String? blacklist = if size('blacklist')==0 then null_s else read_string('blacklist')
        String? blacklist2 = if size('blacklist2')==0 then null_s else read_string('blacklist2')
        String? regex_bfilt_peak_chr_name = if size('regex_bfilt_peak_chr_name')==0 then 'chr[\\dXY]+'
            else read_string('regex_bfilt_peak_chr_name')
        String? tss = if size('tss')!=0 then read_string('tss')
            else if size('tss_enrich')!=0 then read_string('tss_enrich') else null_s
        String? dnase = if size('dnase')==0 then null_s else read_string('dnase')
        String? prom = if size('prom')==0 then null_s else read_string('prom')
        String? enh = if size('enh')==0 then null_s else read_string('enh')
        String? reg2map = if size('reg2map')==0 then null_s else read_string('reg2map')
        String? reg2map_bed = if size('reg2map_bed')==0 then null_s else read_string('reg2map_bed')
        String? roadmap_meta = if size('roadmap_meta')==0 then null_s else read_string('roadmap_meta')
    }
    runtime {
        maxRetries : 0
        cpu : 1
        memory : '2 GB'
        time : 4
        disks : 'local-disk 10 SSD'

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}

task raise_exception {
    input {
        String msg

        # runtime environment
        RuntimeEnvironment runtime_environment
    }
    command {
        echo -e "\n* Error: ${msg}\n" >&2
        exit 2
    }
    output {
        String error_msg = '${msg}'
    }
    runtime {
        maxRetries : 0
        cpu : 1
        memory : '2 GB'
        time : 4
        disks : 'local-disk 10 SSD'

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}
