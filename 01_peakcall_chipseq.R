source("utils/frip.R")
source("utils/get_filename.R")
source("utils/read_peak_file.R")
source("utils/control_file.R")

#### Tweakable stuff ###########################################################
read_dir <- "data/chipseq_reads" # Directory containing the BAM files
out_dir <- "data/peaks_chipseq/default_q0.05"
qvalue <- 0.05 # Q-value threshold for peak calling

#### Run #######################################################################
# Make tempdir
temp_dir <- file.path(tempdir(), "peakcalling")
if (!dir.exists(temp_dir)) {
    dir.create(temp_dir)
}

# List all the alignment read (BAM) files in the directory
bam_files <- list.files(read_dir) |>
  grep("bam$", x = _, value = TRUE)
bam_files

# Create the output directory if it doesn't exist
if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

# Intitialise dataframe for storing metrics
metrics <- data.frame(
    read.count = integer(),
    peak.count = integer(),
    frip = numeric()
)

# Call peaks and record metrics
for (bam_file in bam_files) {
    # Skip for control files
    if (is_control(bam_file)) {
        next
    }
    
    # Get sample name from filename
    sample_name <- get_filename_noext(bam_file)
    paste("\nProcessing sample:", sample_name, "\n") |>
        cat()
    
    # Get full BAM path
    bam_path <- file.path(read_dir, bam_file)
    
    # Warn if no control file is found
    control_file <- get_control(bam_file, bam_files)
    if (is.na(control_file)) {
        warning("No control file found for sample: ", sample_name)
        control_bam_path <- NULL
    } else {
        paste("Control file found:", control_file) |>
            cat()
        control_bam_path <- file.path(read_dir, control_file)
    }
    
    # Call peaks
    peak_temp <- MACSr::callpeak(
        tfile = bam_path,
        cfile = control_bam_path,
        qvalue = qvalue,
        format = "BAMPE", # Paired-end
        name = sample_name,
        outdir = temp_dir
    )

    # Capture output name
    peak_temp <- peak_temp$outputs |>
        purrr::keep( ~ tools::file_ext(.) == "narrowPeak")
    
    # Move peak file to output directory
    peak_out_path <- file.path(out_dir, paste(sample_name, qvalue, "peaks.narrowPeak", sep="_"))
    file.rename(peak_temp, peak_out_path)
    
    # Gather metrics
    ## Alignment file
    read_count <- Rsamtools::countBam(bam_path)$records
    bam_file <- Rsamtools::BamFile(bam_path)
    
    ## Peak file
    peak_file <- read_peak_file(peak_out_path)
    peak_count <- length(peak_file)
    frip_score <- frip(bam_file, peak_file, total_reads = read_count)
    
    ## Write to file
    new_row <- data.frame(
        read.count = read_count,
        peak.count = peak_count,
        frip.score = frip_score
    )
    rownames(new_row) <- sample_name
    print(new_row)
    metrics <- rbind(metrics, new_row)
}

# Save metrics to file
metrics_out_path <- file.path(out_dir, paste0("peak_calling_qval", qvalue, "_metrics.csv"))
write.csv(metrics, metrics_out_path, row.names = TRUE)
paste("Metrics saved to:", metrics_out_path) |>
    cat()

