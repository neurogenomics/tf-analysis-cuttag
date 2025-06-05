# AIM: Check if peaks with lower qval are more likely to contain the motif
# Final product: A boxplot of peak qvals with motif vs without motif

#### Tweakable stuff ###########################################################
target <- "IKZF1-SK342_10"
out_dir <- paste0("data/peaks-qval-motif/", target)
peak_files <- c(
  "data/peaks_chipseq/encode/IKZF1_ENCFF637SIR.narrowPeak", # Reference
  "data/peaks_cuttag/5-shift_q0.01/SK342_10_R1_0.01_peaks.narrowPeak" # Comparison
)
motif_files <- "data/motifs_jaspar/IKZF1_MA1508.2.jaspar"
meme_path <- "/Users/hdash/tools/meme-5.5.8/bin"

#### Run #######################################################################
library(tidyverse)
library(ggplot2)
library(patchwork)
library(MotifPeeker)

dir.create(out_dir, recursive = TRUE)

# Setup
genome_build <- check_genome_build("hg38")
params <- list(
  reference_index = 1,
  BPPARAM = BiocParallel::MulticoreParam(workers = 4),
  meme_path = meme_path,
  verbose = TRUE
)

# Load peak files
peaks <- lapply(peak_files, read_peak_file)

# Load motif file
user_motifs <- lapply(motif_files, read_motif_file)

# Segregate peaks based on overlap
comparison_indices <- setdiff(seq_along(peaks), params$reference_index)
segregated_peaks <- lapply(comparison_indices, function(x) {
  segregate_seqs(peaks[[params$reference_index]], peaks[[x]])
})

# Compute motif enrichment in segregated peaks
motif_enrichment <- parallel::mclapply(
    segregated_peaks[[1]],
    function(peaks) {
        motif_enrichment(
            peaks,
            user_motifs[[1]],
            genome_build = genome_build,
            out_dir = file.path(tempdir(), sample(seq(10000,99999), 1)),
            verbose = params$verbose,
            meme_path = params$meme_path
        )
    },
    mc.preschedule = FALSE,
    mc.cores = params$BPPARAM$workers
)

# Add negative peaks
motif_enrichment2 <- lapply(names(motif_enrichment), function(name) {
    peakgroup <- motif_enrichment[[name]]
    all_peaks <- segregated_peaks[[1]][[name]]
    peakgroup$negative_peaks <- all_peaks[!all_peaks$name %in% peakgroup$positive_peaks$name]
    return(name = peakgroup)
})
names(motif_enrichment2) <- names(motif_enrichment)

# Extract q-values and combine into data.frame
qvalues_df <- lapply(names(motif_enrichment2), function(name) {
    x <- motif_enrichment2[[name]]
    data.frame(
        peak_group = name,
        peak_name = c(x$positive_peaks$name, x$negative_peaks$name),
        qvalue = c(x$positive_peaks$qValue, x$negative_peaks$qValue),
        motif_present = c(rep(TRUE, length(x$positive_peaks)), rep(FALSE, length(x$negative_peaks)))
    )
})
qvalues_df <- bind_rows(qvalues_df)

# Specify experiment technique
qvalues_df <- qvalues_df %>%
    mutate(
        experiment = factor(case_when(
            grepl("seqs1", peak_group) ~ "ChIP-Seq",
            TRUE ~ "CUT&Tag"
        ))
    )

# Transform qvalues
# MACS reports -log10(qvalue), so we need to perform inverse-operation
qvalues_df$qvalue_actual <- 10^-(qvalues_df$qvalue)

# Plot one box plot with all peak groups combined
main_boxplot <- ggplot(qvalues_df, aes(x = motif_present, y = qvalue)) +
    geom_boxplot() +
    labs(x = "Motif Present", y = "-log10(Q-value)\n Higher = Better", title = "Overall") +
    scale_y_continuous(limits = quantile(qvalues_df$qvalue, c(0.05, 0.95))) +
    theme_minimal()
main_boxplot


# Plot boxplot with peak_group facet
facet_boxplot <- ggplot(qvalues_df, aes(x = motif_present, y = qvalue)) +
    geom_boxplot() +
    facet_wrap(~ experiment, scales = "free") +
    labs(x = "Motif Present", y = "-log10(Q-value)\n Higher = Better", title = "Group-wise") +
    scale_y_continuous(limits = quantile(qvalues_df$qvalue, c(0.05, 0.95))) +
    theme_minimal()
facet_boxplot

combined_boxplots <- main_boxplot / facet_boxplot
combined_boxplots <- combined_boxplots +
    plot_annotation(
        title = "Q-values of Peaks with and without Motif",
        subtitle = paste("Target:", target, " | y-axis limits: quantile 0.05 - 0.95"),
        theme = theme(plot.title = element_text(hjust = 0.5))
    )
combined_boxplots

ggsave(file.path(out_dir, "qval_boxplot.png"), combined_boxplots, width = 10, height = 10)


# Statistical Tests
tests <- list()
tests$t_test_overall <- t.test(
    qvalues_df$qvalue[qvalues_df$motif_present == TRUE],
    qvalues_df$qvalue[qvalues_df$motif_present == FALSE]
)
tests$t_test_overall

tests$wilcox_rank_test_overall <- wilcox.test(
    qvalues_df$qvalue_actual[qvalues_df$motif_present == TRUE],
    qvalues_df$qvalue_actual[qvalues_df$motif_present == FALSE]
)
tests$wilcox_rank_test_overall

tests$t_test_chipseq <- t.test(
    qvalues_df$qvalue[qvalues_df$motif_present == TRUE & qvalues_df$experiment == "ChIP-Seq"],
    qvalues_df$qvalue[qvalues_df$motif_present == FALSE & qvalues_df$experiment == "ChIP-Seq"]
)
tests$t_test_chipseq

tests$wilcox_rank_test_chipseq <- wilcox.test(
    qvalues_df$qvalue[qvalues_df$motif_present == TRUE & qvalues_df$experiment == "ChIP-Seq"],
    qvalues_df$qvalue[qvalues_df$motif_present == FALSE & qvalues_df$experiment == "ChIP-Seq"]
)
tests$wilcox_rank_test_chipseq

tests$t_test_cuttag <- t.test(
    qvalues_df$qvalue[qvalues_df$motif_present == TRUE & qvalues_df$experiment == "CUT&Tag"],
    qvalues_df$qvalue[qvalues_df$motif_present == FALSE & qvalues_df$experiment == "CUT&Tag"]
)
tests$t_test_cuttag

tests$wilcox_rank_test_cuttag <- wilcox.test(
    qvalues_df$qvalue[qvalues_df$motif_present == TRUE & qvalues_df$experiment == "CUT&Tag"],
    qvalues_df$qvalue[qvalues_df$motif_present == FALSE & qvalues_df$experiment == "CUT&Tag"]
)
tests$wilcox_rank_test_cuttag

# Write output to file
test_results_file <- file(file.path(out_dir, "statistical_tests.txt"), open = "wt")
writeLines("Tests for comparing two independent groups (-log10(qvalue))", test_results_file)
writeLines(paste("Target:", target), test_results_file)
lapply(names(tests), function(name) {
    write(paste("\n\n", name, ":"), test_results_file, append = TRUE)
    write(capture.output(tests[[name]]), test_results_file, append = TRUE)
})
close(test_results_file)

