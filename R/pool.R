# pool.R - Aggregation (pooling) functions for dbSequence
# ==============================================================================
# Provides pool() S3 generic and methods for aggregating values by grouping
# columns. This enables database-native aggregation without loading data into
# memory.

#' Pool (aggregate) values by grouping columns
#'
#' @description
#' Aggregates values in a lazy database table by grouping columns. This is
#' the database-native equivalent of `dplyr::group_by() + summarise()` but
#' with a simpler interface optimized for common pooling operations.
#'
#' @param x A dbSequence object, tbl_duckdb_connection, or other lazy table
#' @param group_by Character vector of column names to group by
#' @param value_col Character, the column to aggregate. Default is "score"
#'   for dbSequence, "x" for tbl objects.
#' @param fun Character, aggregation function: "sum" (default), "mean",
#'   "count", "min", "max"
#' @param name Character, optional name for the resulting computed table.
#'   If NULL, returns a lazy query without materializing.
#' @param filter_zero Logical, if TRUE (default) filter out rows where
#'   aggregated value == 0
#' @param temporary Logical, if TRUE (default) create a temporary table,
#'   otherwise create a permanent table
#' @param overwrite Logical, if TRUE (default) overwrite existing table
#' @param ... Additional arguments passed to dplyr::compute()
#'
#' @return Same type as input (lazy, still in database). For dbSequence,
#'   returns dbSequence if range columns are preserved, otherwise returns tbl.
#'
#' @examples
#' bed <- system.file("extdata", "example.bed", package = "dbSequence")
#' db_seq <- read_bed(bed)
#'
#' # Pool feature scores by feature name
#' pooled <- pool(db_seq, group_by = "name", value_col = "score")
#' pooled
#'
#' @export
pool <- function(x, group_by, value_col = NULL, fun = "sum",
                 name = NULL, filter_zero = TRUE,
                 temporary = TRUE, overwrite = TRUE, ...) {
  UseMethod("pool")
}

#' @rdname pool
#' @export
pool.default <- function(x, group_by, value_col = NULL, fun = "sum",
                         name = NULL, filter_zero = TRUE,
                         temporary = TRUE, overwrite = TRUE, ...) {
  stop("pool() is not implemented for class: ", class(x)[1],
       ". Supported classes: tbl_duckdb_connection, dbSequence")
}

#' @rdname pool
#' @export
pool.tbl_duckdb_connection <- function(x, group_by, value_col = "x",
                                        fun = "sum", name = NULL,
                                        filter_zero = TRUE,
                                        temporary = TRUE, overwrite = TRUE, ...) {
  # Validate group_by columns exist
  available_cols <- colnames(x)
  if (!all(group_by %in% available_cols)) {
    missing <- setdiff(group_by, available_cols)
    stop("Column(s) not found: ", paste(missing, collapse = ", "),
         ". Available: ", paste(available_cols, collapse = ", "))
  }

  # Validate value_col exists
  if (!value_col %in% available_cols) {
    stop("Value column '", value_col, "' not found. Available: ",
         paste(available_cols, collapse = ", "))
  }

  # Build aggregation using dplyr
  result <- x |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_by)))

  # Apply aggregation function
  agg_expr <- switch(fun,
    "sum" = rlang::expr(sum(.data[[!!value_col]], na.rm = TRUE)),
    "mean" = rlang::expr(mean(.data[[!!value_col]], na.rm = TRUE)),
    "count" = rlang::expr(dplyr::n()),
    "min" = rlang::expr(min(.data[[!!value_col]], na.rm = TRUE)),
    "max" = rlang::expr(max(.data[[!!value_col]], na.rm = TRUE)),
    stop("Unsupported aggregation function: ", fun,
         ". Supported: sum, mean, count, min, max")
  )

  result <- result |>
    dplyr::summarise(!!value_col := !!agg_expr, .groups = "drop")

  # Filter zeros if requested
  if (filter_zero && fun %in% c("sum", "count")) {
    result <- result |> dplyr::filter(.data[[value_col]] > 0)
  }

  # Compute to materialized table if name provided
  if (!is.null(name)) {
    result <- dplyr::compute(
      result,
      name = name,
      temporary = temporary,
      overwrite = overwrite,
      ...
    )
  }

  result
}

#' @rdname pool
#' @export
pool.tbl_dbi <- pool.tbl_duckdb_connection

#' @rdname pool
#' @export
pool.dbSequence <- function(x, group_by, value_col = "score",
                            fun = "sum", name = NULL,
                            filter_zero = TRUE,
                            temporary = TRUE, overwrite = TRUE, ...) {
  # Use the underlying tbl
  tbl <- .dbseq_value(x)

  # Pool using tbl method
  result_tbl <- pool.tbl_duckdb_connection(
    tbl,
    group_by = group_by,
    value_col = value_col,
    fun = fun,
    name = name,
    filter_zero = filter_zero,
    temporary = temporary,
    overwrite = overwrite,
    ...
  )

  # Return new dbSequence if grouping preserves range columns
  range_cols <- c("seqnames", "start", "end", "chrom", "chr", "chromosome")
  if (any(range_cols %in% group_by)) {
    # Return dbSequence
    methods::new("dbSequence",
        value = result_tbl,
        name = name %||% tableName(x),
        file_source = fileSource(x))
  } else {
    # Just return the tbl (no longer genomic ranges)
    result_tbl
  }
}
