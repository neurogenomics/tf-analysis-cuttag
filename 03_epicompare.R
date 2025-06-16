library(EpiCompare)

# Setup
data(hg38_blacklist)
peakfiles <- list(
    "CutTag_A303-531A_Tom" = "/Users/hdash/code/cuttag_Isidora_TFprofiling/data/peaks_cuttag/5-shift_q0.01/SK342_8_R1_0.01_peaks.narrowPeak",
    "CutTag_SAB4501855-100UG_Tom" = "/Users/hdash/code/cuttag_Isidora_TFprofiling/data/peaks_cuttag/5-shift_q0.01/SK342_9_R1_0.01_peaks.narrowPeak",
    "CutTag_A303-531A_Paulina" = "/Users/hdash/code/cuttag_Isidora_TFprofiling/data/peaks_cuttag/Paulina_q1e-05/SK342_8_R1_1e-05_peaks.narrowPeak",
    "CutTag_SAB4501855-100UG_Paulina" = "/Users/hdash/code/cuttag_Isidora_TFprofiling/data/peaks_cuttag/Paulina_q1e-05/SK342_9_R1_1e-05_peaks.narrowPeak"
)
reference <- list(
    "ChIP_ENCFF031NTF" = "/Users/hdash/code/cuttag_Isidora_TFprofiling/data/peaks_chipseq/encode/MEF2A_ENCFF031NTF.narrowPeak"
)

genome_build <- list(peakfiles="hg38", reference="hg38", blacklist="hg38")
genome_build_output <- "hg38"
out_dir <- "data/epicompare-reports/MEF2A"
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
    # run_all = TRUE,
    upset_plot = TRUE,
    stat_plot = TRUE,
    chromHMM_plot = TRUE,
    chipseeker_plot = TRUE,
    enrichment_plot = TRUE,
    tss_plot = TRUE,
    precision_recall_plot = TRUE,
    corr_plot = TRUE,
    add_download_button = FALSE,
    workers = cpu_cores,
    error = TRUE
)

