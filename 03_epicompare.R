library(EpiCompare)

# Setup
data(hg38_blacklist)
peakfiles <- list(
    "CutTag_SK322_1" = "/Users/hdash/code/cuttag_Isidora_TFprofiling/data/peaks_cuttag/5-shift_q0.01/SK322_1_R1_0.01_peaks.narrowPeak",
    "CutTag_SK322_2" = "/Users/hdash/code/cuttag_Isidora_TFprofiling/data/peaks_cuttag/5-shift_q0.01/SK322_2_R1_0.01_peaks.narrowPeak"
)
reference <- list(
    "ChIP_ENCFF456USL" = "/Users/hdash/code/cuttag_Isidora_TFprofiling/data/peaks_chipseq/encode/CEBPG_ENCFF456USL.bed"
)

genome_build <- list(peakfiles="hg38", reference="hg38", blacklist="hg38")
genome_build_output <- "hg38"
out_dir <- "data/epicompare-reports/CEBPG"
cpu_cores <- 1 # Error if multi-core

# Create directory
if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}

# Read peak files, convert to GRanges object
peakfiles2 <- lapply(peakfiles, MotifPeeker::read_peak_file)
reference2 <- lapply(reference, MotifPeeker::read_peak_file)

# Run EpiCompare
EpiCompare(
    peakfiles = peakfiles2,
    reference = reference2,
    genome_build = genome_build,
    blacklist = hg38_blacklist,
    genome_build_output = genome_build_output,
    output_dir = out_dir,
    debug = TRUE,
    run_all = TRUE,
    add_download_button = FALSE,
    workers = cpu_cores
)

