# plyranges-overlaps.R -- Range overlap operations for dbSequence
# ------------------------------------------------------------------------------
# Implements plyranges-style range operations that translate to efficient SQL
# queries in DuckDB.

# ------------------------------------------------------------------------------
#  filter_by_overlaps - Filter ranges overlapping another set
# ------------------------------------------------------------------------------

#' Filter ranges by overlap
#'
#' @description Filter a dbSequence object to only include ranges that overlap
#'   with a set of query ranges. This operation is performed in DuckDB using
#'   efficient SQL filtering.
#'
#' @param x A dbSequence object
#' @param y A GRanges object or dbSequence object representing query ranges
#' @param ... Additional arguments (currently unused)
#'
#' @return A dbSequence object filtered to overlapping ranges
#'
#' @details
#' Two intervals overlap if: start1 <= end2 AND start2 <= end1 AND chr1 == chr2
#'
#' The operation is pushed down to DuckDB, so only matching rows are retrieved.
#'
#' @examples
#' \dontrun{
#' # Filter fragments to a specific region
#' region <- GenomicRanges::GRanges("chr1:1000-2000")
#' filtered <- filter_by_overlaps(fragments, region)
#' }
#'
#' @export
filter_by_overlaps <- function(x, y, ...) {
  UseMethod("filter_by_overlaps")
}

#' @rdname filter_by_overlaps
#' @export
filter_by_overlaps.dbSequence <- function(x, y, ...) {
  # Validate arguments in ...
  dots <- list(...)
  if ("maxgap" %in% names(dots)) {
    if (dots$maxgap != -1L) {
      stop("Argument 'maxgap' is not supported for dbSequence objects yet.")
    }
  }
  if ("minoverlap" %in% names(dots)) {
    if (dots$minoverlap != 0L) {
      stop("Argument 'minoverlap' is not supported for dbSequence objects yet.")
    }
  }
  
  # Validate file type - BAM/CRAM not supported
  .validate_for_range_ops(x, "filter_by_overlaps")

  # Get the tbl object and connection
  tbl <- x@value
  if (is.null(tbl)) {
    cli::cli_abort("dbSequence object has no data (table may not exist)")
  }

  # Get the DuckDB connection from the tbl
  con <- dbplyr::remote_con(tbl)

  # Detect column names in x
  x_cols <- colnames(tbl)
  chr_col <- .detect_seqnames_col(x_cols)
  start_col <- .detect_start_col(x_cols)
  end_col <- .detect_end_col(x_cols)

  # Handle dbSequence y separately (table-to-table join)
  if (is(y, "dbSequence")) {
    y_tbl <- y@value
    y_cols <- colnames(y_tbl)
    y_chr_col <- .detect_seqnames_col(y_cols)
    y_start_col <- .detect_start_col(y_cols)
    y_end_col <- .detect_end_col(y_cols)

    # Rename y columns to avoid conflicts, then join
    y_renamed <- y_tbl |>
      dplyr::select(
        y_chr = !!rlang::sym(y_chr_col),
        y_start = !!rlang::sym(y_start_col),
        y_end = !!rlang::sym(y_end_col)
      )

    # Use cross join + filter for inequality join
    result_tbl <- tbl |>
      dplyr::cross_join(y_renamed) |>
      dplyr::filter(
        !!rlang::sym(chr_col) == y_chr,
        !!rlang::sym(start_col) <= y_end,
        !!rlang::sym(end_col) >= y_start
      ) |>
      dplyr::select(-y_chr, -y_start, -y_end) |>
      dplyr::distinct()

    return(new(
      "dbSequence",
      value = result_tbl,
      name = x@name,
      file_source = x@file_source
    ))
  }

  # Handle GRanges y
  if (!is(y, "GRanges")) {
    stop("y must be a GRanges or dbSequence object")
  }

  if (length(y) == 0) {
    warning("No ranges provided in y, returning empty result")
    return(x)
  }

  # Convert GRanges to data frame
  y_df <- data.frame(
    y_chr = as.character(GenomicRanges::seqnames(y)),
    y_start = as.integer(GenomicRanges::start(y)),
    y_end = as.integer(GenomicRanges::end(y))
  )

  # Write y to a temporary table in DuckDB
  temp_table_name <- paste0("temp_ranges_", gsub("[^0-9]", "", format(Sys.time(), "%H%M%S%OS3")))
  DBI::dbWriteTable(con, temp_table_name, y_df, temporary = TRUE, overwrite = TRUE)

  # Create a tbl reference to the temp table
  y_tbl <- dplyr::tbl(con, temp_table_name)

  # Perform the overlap join entirely in DuckDB using cross join + filter
  # This is more compatible with dbplyr than inequality join conditions
  result_tbl <- tbl |>
    dplyr::cross_join(y_tbl) |>
    dplyr::filter(
      !!rlang::sym(chr_col) == y_chr,
      !!rlang::sym(start_col) <= y_end,
      !!rlang::sym(end_col) >= y_start
    ) |>
    dplyr::select(-y_chr, -y_start, -y_end) |>
    dplyr::distinct()

  # Return new dbSequence with filtered tbl (still lazy!)
  new(
    "dbSequence",
    value = result_tbl,
    name = x@name,
    file_source = x@file_source
  )
}

# ------------------------------------------------------------------------------
#  compute_coverage - Compute binned coverage
# ------------------------------------------------------------------------------

#' Compute coverage over a region
#'
#' @description Compute coverage (count of overlapping fragments) in bins
#'   across a genomic region. This is useful for visualization and
#'   summarizing fragment density.
#'
#' @param x A dbSequence object
#' @param region A GRanges object or string like "chr1:1000-2000"
#' @param window Bin size in base pairs (default: 100)
#' @param ... Additional arguments
#'
#' @return A data frame with columns: bin_start, bin_end, count
#'
#' @examples
#' \dontrun{
#' region <- "chr1:1000000-2000000"
#' coverage <- compute_coverage(fragments, region, window = 1000)
#' }
#'
#' @export
compute_coverage <- function(x, region, window = 100, ...) {
  UseMethod("compute_coverage")
}

#' @rdname compute_coverage
#' @export
compute_coverage.dbSequence <- function(x, region, window = 100, ...) {
  # Parse region if it's a string
  if (is.character(region)) {
    region <- .parse_region(region)
  }

  if (!is(region, "GRanges") || length(region) != 1) {
    stop("region must be a single GRanges object or a string like 'chr1:1000-2000'")
  }

  # Extract region coordinates
  chr <- as.character(GenomicRanges::seqnames(region))
  region_start <- GenomicRanges::start(region)
  region_end <- GenomicRanges::end(region)

  # First filter to the region
  filtered <- filter_by_overlaps(x, region)
  tbl <- filtered@value

  if (is.null(tbl)) {
    # Return empty coverage
    return(data.frame(bin_start = integer(), bin_end = integer(), count = integer()))
  }

  # Detect column names
  tbl_cols <- colnames(tbl)
  start_col <- .detect_start_col(tbl_cols)

  # Compute binned counts using SQL
  # FLOOR(start / window) * window gives the bin start
  coverage_tbl <- tbl |>
    dplyr::mutate(
      bin_start = floor(!!rlang::sym(start_col) / !!window) * !!window
    ) |>
    dplyr::group_by(bin_start) |>
    dplyr::summarise(count = dplyr::n(), .groups = "drop") |>
    dplyr::arrange(bin_start) |>
    dplyr::collect()

  # Add bin_end column
  coverage_tbl$bin_end <- coverage_tbl$bin_start + window

  # Filter to only bins within the region
  coverage_tbl <- coverage_tbl |>
    dplyr::filter(bin_start >= region_start & bin_end <= region_end)

  coverage_tbl
}

# ------------------------------------------------------------------------------
#  Helper functions
# ------------------------------------------------------------------------------

#' Parse a region string to GRanges
#' @keywords internal
#' @noRd
.parse_region <- function(region) {
  if (is(region, "GRanges")) {
    return(region)
  }

  if (!is.character(region) || length(region) != 1) {
    stop("region must be a single string like 'chr1:1000-2000'")
  }

  # Remove commas from numbers
  region <- gsub(",", "", region)

  # Parse "chr1:1000-2000" format
  parts <- strsplit(region, ":")[[1]]
  if (length(parts) != 2) {
    stop("Invalid region format. Expected 'chr:start-end'")
  }

  chr <- parts[1]
  coords <- strsplit(parts[2], "-")[[1]]
  if (length(coords) != 2) {
    stop("Invalid region format. Expected 'chr:start-end'")
  }

  start <- as.integer(coords[1])
  end <- as.integer(coords[2])

  if (is.na(start) || is.na(end)) {
    stop("Could not parse coordinates from region string")
  }

  GenomicRanges::GRanges(
    seqnames = chr,
    ranges = IRanges::IRanges(start = start, end = end)
  )
}

#' Detect the seqnames/chromosome column
#' @keywords internal
#' @noRd
.detect_seqnames_col <- function(cols) {
  candidates <- c("seqnames", "chrom", "chr", "chromosome",
                  "reference_sequence_name", "ref_name", "rname")
  for (col in candidates) {
    if (col %in% cols) return(col)
  }
  stop("Could not find chromosome column. Expected one of: ",
       paste(candidates, collapse = ", "))
}

#' Detect the start position column
#' @keywords internal
#' @noRd
.detect_start_col <- function(cols) {
  candidates <- c("start", "pos", "position", "chromStart")
  for (col in candidates) {
    if (col %in% cols) return(col)
  }
  stop("Could not find start column. Expected one of: ",
       paste(candidates, collapse = ", "))
}
#' Detect the end position column
#' @keywords internal
#' @noRd
.detect_end_col <- function(cols) {
  candidates <- c("end", "chromEnd", "stop")
  for (col in candidates) {
    if (col %in% cols) return(col)
  }
  # If no end column, use start (for point data like SNPs)
  start_col <- .detect_start_col(cols)
  warning("No end column found, using start column for single-position ranges")
  return(start_col)
}
