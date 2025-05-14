#' Calculate FRiP score
#' FRiP = (number of reads in peaks) / (total number of reads)
#' 
#' https://github.com/neurogenomics/motifmarker-playground/blob/master/R/frip.R
#' 
#' @param read_file A BamFile object.
#' @param peak_file A GRanges object.
#' @param single_end If TRUE, the reads classified as single-ended.
#' @param total_reads (optional) The total number of reads in the experiment.
#' @return A numeric value.
frip <- function(
        read_file,
        peak_file,
        single_end = FALSE,
        total_reads = NULL
) {
    library(GenomicAlignments)
    
    overlaps <-
        GenomicAlignments::summarizeOverlaps(peak_file,
                                             read_file,
                                             singleEnd = single_end,
                                             ignore.strand = TRUE)
    peak_reads <- sum(assay(overlaps))
    
    if (is.null(total_reads)) {
        total_reads <- countBam(read_file)$records
    }
    
    return(peak_reads / total_reads)
}