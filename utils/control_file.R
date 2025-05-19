source("utils/get_filename.R")

is_control <- function(filename) {
    grepl("control", get_filename_noext(filename), ignore.case = TRUE)
}

get_control <- function(filename, bam_files) {
    # Check if input is control
    if (is_control(filename)) {
        stop("Input file is a control file.")
    }
    
    # Get the control file name
    exp_code <- get_filename_noext(filename) |>
        strsplit("_") |>
        unlist() |>
        dplyr::first()
    
    control_file <- bam_files |>
        purrr::keep(is_control) |>
        purrr::keep(\(x) grepl(exp_code, x = x, ignore.case = TRUE)) |>
        dplyr::first()
}