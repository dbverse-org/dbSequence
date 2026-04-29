# compute.R -- Materialization helpers for dbSequence
# ------------------------------------------------------------------------------

#' Materialize a dbSequence in DuckDB
#'
#' @description
#' S3 method for \code{dplyr::compute()} that materializes the underlying lazy
#' query into a DuckDB table, returning a new dbSequence pointing at the
#' computed table.
#'
#' @param x A dbSequence object
#' @param name Character table name to create. Defaults to \code{tableName(x)}.
#' @param temporary Logical, create a temporary table (default: TRUE)
#' @param overwrite Logical, overwrite existing table (default: TRUE)
#' @param ... Passed to \code{dplyr::compute()}.
#'
#' @return A dbSequence object backed by the materialized table.
#'
#' @examples
#' bed <- system.file("extdata", "example.bed", package = "dbSequence")
#' db_seq <- read_bed(bed)
#' materialized <- compute.dbSequence(db_seq, name = "bed_materialized")
#' materialized
#'
#' @export
compute.dbSequence <- function(
  x,
  name = tableName(x),
  temporary = TRUE,
  overwrite = TRUE,
  ...
) {
  computed_tbl <- dplyr::compute(
    x@value,
    name = name,
    temporary = temporary,
    overwrite = overwrite,
    ...
  )

  methods::new(
    "dbSequence",
    value = computed_tbl,
    name = name,
    file_source = x@file_source
  )
}
