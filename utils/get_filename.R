get_filename_noext <- function(path) {
    split_str <- basename(path) |>
        strsplit(".", fixed = TRUE) |>
        unlist()
    
    return(split_str[1])
}
