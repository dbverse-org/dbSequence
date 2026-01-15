# constructors.R -- Constructor functions for dbSequence classes
# ------------------------------------------------------------------------------
# Constructor functions following Bioconductor naming conventions

#' @title Constructor for DuckDBFile objects
#' @name DuckDBFile
#'
#' @description Creates a DuckDBFile object that wraps a file path for DuckDB databases.
#'   This is a simple path wrapper following the BiocFile pattern and supports both
#'   file paths and in-memory databases.
#'
#' @param resource (required) character: path to DuckDB file. Can be a file path or
#'   ":memory:" for in-memory databases
#' @return DuckDBFile object containing the path information
#' @export
#' @importFrom methods new
#' @examples
#' \dontrun{
#' # Create DuckDBFile for existing database
#' db_file <- DuckDBFile("data.duckdb")
#'
#' # Create DuckDBFile for in-memory database
#' mem_db <- DuckDBFile(":memory:")
#' }
DuckDBFile <- function(resource) {
  # Input validation
  stopifnot(is.character(resource), length(resource) == 1, !is.na(resource))

  # DuckDBFile is now just a path wrapper following BiocFile pattern
  new("DuckDBFile", resource = resource, path = resource)
}

#' @title Constructor for dbSequence objects
#' @name dbSequence
#'
#' @description Creates a dbSequence object by wrapping a table name.
#'   This is typically called internally by import methods.
#'
#' @param table_name (required) character: name of the table in the DuckDB database
#' @param file_source (required) character or DuckDBFile: source file or database
#' @param .conn (optional) DBIConnection: existing connection for in-memory
#'   databases (default: NULL)
#' @return A dbSequence object
#' @export
#' @examples
#' \dontrun{
#' # Create a dbSequence object
#' db_file <- DuckDBFile("data.duckdb")
#' db_seq <- dbSequence("variants_table", db_file)
#' }
dbSequence <- function(
  table_name,
  file_source = NULL,
  .conn = NULL
) {
  # Allow passing a tbl object directly (e.g. from dplyr operations)
  if (inherits(table_name, "tbl_dbi")) {
    tbl_value <- table_name
    
    # Try to extract file source from connection
    if (is.null(file_source)) {
      conn <- dbplyr::remote_con(tbl_value)
      info <- tryCatch(DBI::dbGetInfo(conn), error = function(e) list())
      file_source <- info$dbname
      if (is.null(file_source)) file_source <- ":memory:"
    }
    
    # Use generic name for query results
    nm <- "custom_query"
    
    return(new(
      "dbSequence",
      value = tbl_value,
      name = nm,
      file_source = as.character(file_source)
    ))
  }

  # Standard constructor: table_name must be character
  stopifnot(
    is.character(table_name),
    length(table_name) == 1,
    !is.na(table_name)
  )
  
  if (is.null(file_source)) {
    stop("file_source is required when table_name is a character string")
  }

  # Extract file path from DuckDBFile or use as character
  if (is(file_source, "DuckDBFile")) {
    file_source_char <- file_source@path
  } else {
    file_source_char <- as.character(file_source)
  }

  # Use provided connection or create new one
  if (!is.null(.conn)) {
    # Use provided connection (typically for in-memory databases)
    conn <- .conn
    should_disconnect <- FALSE
  } else {
    # Create new connection from file path
    conn <- DBI::dbConnect(duckdb::duckdb(), dbdir = file_source_char)
    should_disconnect <- TRUE
  }

  # Create tbl object for the table if it exists
  tbl_value <- NULL
  table_exists <- DBI::dbExistsTable(conn, table_name)

  if (table_exists) {
    if (!is.null(.conn)) {
      # Use the provided connection for tbl object (don't create new one)
      tbl_value <- dplyr::tbl(conn, table_name)
    } else {
      # For file databases, create a fresh connection for the tbl object
      if (should_disconnect) {
        DBI::dbDisconnect(conn, shutdown = TRUE)
      }
      tbl_conn <- DBI::dbConnect(duckdb::duckdb(), dbdir = file_source_char)
      tbl_value <- dplyr::tbl(tbl_conn, table_name)
    }
  } else {
    # Table doesn't exist, safe to close connection if we created it
    if (should_disconnect) {
      DBI::dbDisconnect(conn, shutdown = TRUE)
    }
  }

  # Create proper S4 object with dbData architecture
  new(
    "dbSequence",
    value = tbl_value,
    name = table_name,
    file_source = file_source_char
  )
}

