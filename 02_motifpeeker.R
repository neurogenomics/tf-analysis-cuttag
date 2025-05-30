#### Tweakable stuff ###########################################################
out_dir <- "data/motifpeeker-reports/CEBPG"
peak_files <- list(
    "data/peaks_chipseq/encode/CEBPG_ENCFF456USL.narrowPeak",
    "data/peaks_cuttag/5-shift_q0.01/SK322_1_R1_0.01_peaks.narrowPeak",
    "data/peaks_cuttag/5-shift_q0.01/SK322_2_R1_0.01_peaks.narrowPeak"
)
alignment_files <- list(
    "data/chipseq_reads/CEBPG_ENCFF344WOE.bam",
    "data/linear_dedup/SK322_1_R1.target.linear_dedup.sorted.bam",
    "data/linear_dedup/SK322_2_R1.target.linear_dedup.sorted.bam"
    
)
motif_files <- list(
    "data/motifs_jaspar/CEBPG_MA1636.2.jaspar"
)
exp_labels <- c(
    "ChIP_ENCFF456USL",
    "CutTag_SK322_1",
    "CutTag_SK322_2"
)
exp_type <- c(
    "chipseq",
    "cuttag",
    "cuttag"
)
motif_discovery <- FALSE # Set to TRUE if you want to run motif discovery
cpu_cores <- 8


#### Run########################################################################
# Setup
if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}

MotifPeeker::MotifPeeker(
    peak_files = peak_files,
    reference_index = 1,
    genome_build = "hg38",
    alignment_files = alignment_files,
    motif_files = motif_files,
    exp_labels = exp_labels,
    exp_type = exp_type,
    motif_discovery = motif_discovery,
    download_buttons = FALSE,
    out_dir = out_dir,
    BPPARAM = BiocParallel::MulticoreParam(workers = cpu_cores),
    verbose = TRUE,
    quiet = FALSE,
    debug = TRUE,
    meme_path = "/Users/hdash/tools/meme-5.5.8/bin"
)