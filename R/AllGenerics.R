# AllGenerics.R -- Generic function definitions
# ------------------------------------------------------------------------------
# Following Bioconductor convention, all setGeneric() calls are collected here

#' @import methods

# ------------------------------------------------------------------------------
#  Custom generics for dbSequence functionality
# ------------------------------------------------------------------------------

#' @title Extract genomic ranges from dbSequence objects
#' @name asRanges
#'
#' @param x (required) A dbSequence object
#' @param ... (optional) Additional arguments
#' @return A view of the dbSequence with only range columns (seqnames, start, end, strand)
#' @examples
#' bed <- system.file("extdata", "example.bed", package = "dbSequence")
#' db_path <- tempfile(fileext = ".duckdb")
#' db_seq <- read_bed(bed, dest = DuckDBFile(db_path), lazy = FALSE)
#' asRanges(db_seq)
#'
#' @export
asRanges <- function(x, ...) standardGeneric("asRanges")
setGeneric("asRanges")
