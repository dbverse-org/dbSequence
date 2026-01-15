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
#' @export
setGeneric("asRanges", function(x, ...) standardGeneric("asRanges"))
