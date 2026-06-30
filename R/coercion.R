# coercion.R -- Coercion methods between dbSequence and other `[Bioconductor]` classes
# ------------------------------------------------------------------------------
# S4 coercion methods for interoperability with `[Bioconductor]` classes

# We import methods needed for coercion
if (requireNamespace("GenomicRanges", quietly = TRUE)) {
  # Only define these if GenomicRanges is available
  #' @importFrom GenomicRanges GRanges
} else {
  # Graceful degradation when GenomicRanges not available
}

#' @importFrom methods setAs setMethod

# ------------------------------------------------------------------------------
#  dbSequence <-> GRanges coercion
# ------------------------------------------------------------------------------

# These coercion methods will be implemented when GenomicRanges is available

#' @title Coercion methods between dbSequence and GRanges
#' @name coerce-methods
#'
#' @description Methods for converting between dbSequence and GRanges objects.
#'
#' \code{as(dbSequence, "GRanges")} materializes a dbSequence object into an
#' in-memory GRanges object, triggering data collection from the DuckDB database.
#'
#' \code{as(GRanges, "dbSequence")} loads a GRanges object into a DuckDB table
#' and wraps it as dbSequence, enabling lazy evaluation of operations.
#'
#' @param from Object to convert (dbSequence or GRanges)
#' @return Converted object (GRanges or dbSequence respectively)
#'
#' @aliases coerce,dbSequence,GRanges-method coerce,GRanges,dbSequence-method
NULL

# Conversion between dbSequence and GRanges objects
# Only define these methods if GenomicRanges is available
if (requireNamespace("GenomicRanges", quietly = TRUE)) {
  setAs("dbSequence", "GRanges", function(from) {
    # Convert dbSequence (tbl) to GRanges by collecting the stored tbl.
    # This supports both table-backed objects and query-backed objects.

    # Check if GenomicRanges is available
    if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
      stop(
        "GenomicRanges package is required for this conversion. Please install it."
      )
    }

    if (is.null(.dbseq_value(from))) {
      stop("dbSequence object has no data (table/query may not exist)")
    }

    df <- tryCatch(
      dplyr::collect(.dbseq_value(from)),
      error = function(e) {
        cli::cli_abort(c(
          "Failed to collect data from dbSequence.",
          "x" = conditionMessage(e)
        ))
      }
    )

    if (nrow(df) == 0) {
      # Return empty GRanges
      return(GenomicRanges::GRanges())
    }

    # Map column names to GRanges requirements
    col_names <- names(df)

    # Find chromosome column
    seqnames_col <- NULL
    for (col in c(
      "reference_sequence_name",
      "chrom",
      "seqnames",
      "chr",
      "contig"
    )) {
      if (col %in% col_names) {
        seqnames_col <- col
        break
      }
    }

    # Find start column
    start_col <- NULL
    for (col in c("start", "pos", "position")) {
      if (col %in% col_names) {
        start_col <- col
        break
      }
    }

    # Find end column (optional for point data)
    end_col <- NULL
    for (col in c("end", "stop")) {
      if (col %in% col_names) {
        end_col <- col
        break
      }
    }

    # Find strand column (optional)
    strand_col <- NULL
    if ("strand" %in% col_names) {
      strand_col <- "strand"
    }

    # Validate required columns
    if (is.null(seqnames_col) || is.null(start_col)) {
      stop(
        "Required columns not found. Need chromosome and start position columns."
      )
    }

    # Prepare GRanges arguments
    seqnames <- df[[seqnames_col]]
    start <- df[[start_col]]

    # Handle end position
    if (!is.null(end_col)) {
      end <- df[[end_col]]
    } else {
      # For point data, end = start
      end <- start
    }

    # Handle strand
    strand <- if (!is.null(strand_col)) df[[strand_col]] else "*"

    # Create metadata (all other columns)
    metadata_cols <- setdiff(
      col_names,
      c(seqnames_col, start_col, end_col, strand_col)
    )
    mcols_df <- df[metadata_cols]

    # Create GRanges object
    gr <- GenomicRanges::GRanges(
      seqnames = seqnames,
      ranges = IRanges::IRanges(start = start, end = end),
      strand = strand
    )

    # Add metadata columns
    if (ncol(mcols_df) > 0) {
      GenomicRanges::mcols(gr) <- mcols_df
    }

    return(gr)
  })

  setAs("GRanges", "dbSequence", function(from) {
    # Load GRanges into DuckDB and wrap as dbSequence

    # Check if GenomicRanges is available
    if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
      stop(
        "GenomicRanges package is required for this conversion. Please install it."
      )
    }

    # Convert GRanges to data.frame
    df <- data.frame(
      seqnames = as.character(GenomicRanges::seqnames(from)),
      start = GenomicRanges::start(from),
      end = GenomicRanges::end(from),
      strand = as.character(GenomicRanges::strand(from))
    )

    # Add metadata columns if they exist
    if (ncol(GenomicRanges::mcols(from)) > 0) {
      mcols_df <- as.data.frame(GenomicRanges::mcols(from))
      df <- cbind(df, mcols_df)
    }

    # Create a temporary database connection
    temp_db <- tempfile(fileext = ".db")
    con <- DBI::dbConnect(duckdb::duckdb(), dbdir = temp_db)

    # Generate a unique table name
    table_name <- paste0("granges_import_", format(Sys.time(), "%Y%m%d_%H%M%S"))

    # Import the data into DuckDB
    tryCatch(
      {
        DBI::dbWriteTable(con, table_name, df)
      },
      error = function(e) {
        DBI::dbDisconnect(con)
        stop("Failed to import GRanges data into DuckDB: ", e$message)
      }
    )

    # Create and return dbSequence object
    dbSequence(table_name, temp_db)
  })
} # End of GenomicRanges conditional block

# ------------------------------------------------------------------------------
#  Convenience accessors
# ------------------------------------------------------------------------------

#' @title Extract range columns as a view
#' @name as_ranges-dbSequence-method
#'
#' @description Returns a view of the dbSequence with only the core range columns
#' (seqnames, start, end, strand). This is useful for range operations
#' that don't need metadata columns.
#'
#' @param x (required) A dbSequence object
#' @param ... (optional) Additional arguments (currently unused)
#' @return A dbSequence view with only range columns
#' @aliases asRanges,dbSequence-method
#' @examples
#' bed <- system.file("extdata", "example.bed", package = "dbSequence")
#' db_path <- tempfile(fileext = ".duckdb")
#' db_seq <- read_bed(bed, dest = DuckDBFile(db_path), lazy = FALSE)
#' asRanges(db_seq)
#'
#' @export
setMethod("asRanges", "dbSequence", function(x, ...) {
  # Return a view with chrom/start/end/strand cols only

  # Get the original table name
  original_table <- tableName(x)
  ranges_view_name <- paste0(
    original_table,
    "_ranges_",
    format(Sys.time(), "%Y%m%d_%H%M%S")
  )

  # Try to determine column names for genomic ranges
  # Different file formats use different column names
  conn_obj <- conn(x)

  # Check if we have a valid connection
  if (is.null(conn_obj)) {
    stop("Cannot extract database connection from dbSequence object")
  }

  # Get column names from the table
  column_info <- tryCatch(
    {
      DBI::dbGetQuery(conn_obj, paste0("DESCRIBE ", original_table))
    },
    error = function(e) {
      stop("Failed to get column information for table: ", original_table)
    }
  )

  col_names <- column_info$column_name

  quote_id <- function(identifier) {
    as.character(DBI::dbQuoteIdentifier(conn_obj, identifier))
  }

  # Map common genomic range column names
  range_cols <- character()

  # Chromosome/sequence name column
  if ("reference_sequence_name" %in% col_names) {
    range_cols <- c(range_cols, paste(quote_id("reference_sequence_name"), "AS", quote_id("seqnames")))
  } else if ("seqnames" %in% col_names) {
    range_cols <- c(range_cols, quote_id("seqnames"))
  } else if ("chrom" %in% col_names) {
    range_cols <- c(range_cols, paste(quote_id("chrom"), "AS", quote_id("seqnames")))
  } else if ("chr" %in% col_names) {
    range_cols <- c(range_cols, paste(quote_id("chr"), "AS", quote_id("seqnames")))
  } else if ("contig" %in% col_names) {
    range_cols <- c(range_cols, paste(quote_id("contig"), "AS", quote_id("seqnames")))
  }

  # Start position
  if ("start" %in% col_names) {
    range_cols <- c(range_cols, quote_id("start"))
  } else if ("pos" %in% col_names) {
    range_cols <- c(range_cols, paste(quote_id("pos"), "AS", quote_id("start")))
  } else if ("position" %in% col_names) {
    range_cols <- c(range_cols, paste(quote_id("position"), "AS", quote_id("start")))
  }

  # End position (for ranges like BED files)
  if ("end" %in% col_names) {
    range_cols <- c(range_cols, quote_id("end"))
  } else if ("stop" %in% col_names) {
    range_cols <- c(range_cols, paste(quote_id("stop"), "AS", quote_id("end")))
  }

  # Strand (if available)
  if ("strand" %in% col_names) {
    range_cols <- c(range_cols, quote_id("strand"))
  }

  if (length(range_cols) == 0) {
    stop(
      "No recognizable genomic range columns found in table: ",
      original_table
    )
  }

  # Create the ranges view
  select_clause <- paste(range_cols, collapse = ", ")
  create_view_query <- paste0(
    "CREATE OR REPLACE VIEW ",
    ranges_view_name,
    " AS ",
    "SELECT ",
    select_clause,
    " FROM ",
    original_table
  )

  # Execute the query
  result <- .safe_query(fileSource(x), create_view_query)

  if (is.null(result)) {
    stop("Failed to create ranges view for table: ", original_table)
  }

  # Return new dbSequence object pointing to the ranges view
  dbSequence(ranges_view_name, fileSource(x))
})
