#' dbSequence: DuckDB-Backed Interface for Large-Scale Genomic Data
#'
#' dbSequence provides DuckDB-backed infrastructure for genomic interval and
#' sequence-derived data that are too large to work with comfortably in memory.
#' It imports common genomics file formats into DuckDB, keeps range operations
#' lazy, and interoperates with Bioconductor classes such as
#' [GenomicRanges::GRanges].
#'
#' @keywords internal
"_PACKAGE"
