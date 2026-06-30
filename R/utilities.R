# utilities.R -- Basic utility methods for dbSequence objects
# ------------------------------------------------------------------------------
# Methods that don't require external Bioconductor dependencies

#' @importFrom duckdb duckdb
#' @importFrom DBI dbConnect dbDisconnect dbGetQuery dbExistsTable
#' @importFrom crayon make_style

# ------------------------------------------------------------------------------
#  Internal helper functions
# ------------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.get_duckdb_connection <- function(file_path) {
  drv <- duckdb::duckdb()
  DBI::dbConnect(drv, dbdir = file_path)
}

#' @keywords internal
#' @noRd
.safe_query <- function(file_path, query) {
  tryCatch(
    {
      con <- .get_duckdb_connection(file_path)
      on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
      if (grepl("^\\s*(SELECT|PRAGMA|DESCRIBE|SHOW)\\b", query, ignore.case = TRUE)) {
        DBI::dbGetQuery(con, query)
      } else {
        DBI::dbExecute(con, query)
      }
    },
    error = function(e) {
      warning("DuckDB query failed: ", e$message)
      NULL
    }
  )
}

#' @keywords internal
#' @noRd
.table_exists <- function(file_path, table_name) {
  tryCatch(
    {
      con <- .get_duckdb_connection(file_path)
      on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
      DBI::dbExistsTable(con, table_name)
    },
    error = function(e) {
      FALSE
    }
  )
}

# ------------------------------------------------------------------------------
#  Test helper functions
# ------------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.create_test_bed <- function() {
  # Create temporary BED file
  bed_file <- tempfile(fileext = ".bed")

  # Write test BED data
  bed_data <- c(
    "chr1\t1000\t2000\tfeat1\t500\t+",
    "chr1\t3000\t4000\tfeat2\t750\t-",
    "chr2\t1500\t2500\tfeat3\t600\t+"
  )

  writeLines(bed_data, bed_file)
  return(bed_file)
}

#' @keywords internal
#' @noRd
.cleanup_test_table <- function(db_seq) {
  tryCatch(
    {
      # Get connection and table name
      con <- conn(db_seq)
      tbl_name <- tableName(db_seq)

      # Drop table if it exists
      if (DBI::dbExistsTable(con, tbl_name)) {
        DBI::dbExecute(con, paste0("DROP TABLE IF EXISTS ", tbl_name))
      }

      # Close connection if it's file-based
      if (fileSource(db_seq) != ":memory:") {
        DBI::dbDisconnect(con, shutdown = TRUE)
      }
    },
    error = function(e) {
      # Ignore cleanup errors
      invisible(NULL)
    }
  )
}

#' @title Show method for dbSequence objects
#' @name show,dbSequence-method
#'
#' @description Display a summary of the dbSequence object and preview the table data
#'
#' @param object A dbSequence object
#' @return Invisibly returns NULL.
#' @importFrom methods setMethod
#' @importFrom crayon make_style
#' @importFrom dplyr tbl
setMethod("show", "dbSequence", function(object) {
  # Style for header
  grey_color <- crayon::make_style("grey60")
  value <- .dbseq_value(object)
  file_source <- fileSource(object)

  # Print metadata header
  message(grey_color("# Class:    dbSequence"))

  # Check if object has a valid value (tbl)
  if (!is.null(value)) {
    # Object has a tbl in value slot - show it directly
    return(show(value))
  }

  # Object not initialized - try to establish connection and show table
  tryCatch(
    {
      # For file-based databases, try to create a connection
      if (file_source != ":memory:") {
        con <- .get_duckdb_connection(file_source)
        table_tbl <- dplyr::tbl(con, tableName(object))
        on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
        return(show(table_tbl))
      } else {
        # In-memory database without stored connection
        message(grey_color("# Source:   In-memory database (not initialized)"))
        message(grey_color("# Table:    ", tableName(object)))
        message(grey_color("# Note:     Connection required for data access"))
      }
    },
    error = function(e) {
      if (file_source == ":memory:") {
        message(grey_color("# Source:   In-memory database (connection lost)"))
        message(grey_color("# Table:    ", tableName(object)))
        message(grey_color(
          "# Note:     Connection may have been closed or not stored"
        ))
      } else {
        message(grey_color("# Source:   ", file_source))
        message(grey_color("# Table:    ", tableName(object)))
        message(grey_color("# Status:   Table not accessible (", e$message, ")"))
      }
      invisible(NULL)
    }
  )
})

#' @title Show method for DuckDBFile objects
#' @name show,DuckDBFile-method
#'
#' @description Display information about the DuckDBFile object.
#'
#' @param object A DuckDBFile object
#' @return Invisibly returns NULL.
setMethod("show", "DuckDBFile", function(object) {
  path <- .duckdb_file_path(object)
  message("DuckDBFile object")
  message("Path: ", path)
  if (file.exists(path)) {
    message("File exists: TRUE")
    message("File size: ", file.size(path), " bytes")
  } else {
    message("File exists: FALSE (will be created when needed)")
  }
})

#' @title Get table name from dbSequence
#' @name tableName
#'
#' @description Extract the table name from a dbSequence object.
#'
#' @param x (required) A dbSequence object
#' @return character: the table name
#' @examples
#' bed <- system.file("extdata", "example.bed", package = "dbSequence")
#' db_seq <- read_bed(bed)
#' tableName(db_seq)
#'
#' @export
tableName <- function(x) standardGeneric("tableName")
setGeneric("tableName")

#' @rdname tableName
#' @aliases tableName,dbSequence-method
#' @export
setMethod("tableName", "dbSequence", function(x) {
  .dbseq_table_name(x)
})
