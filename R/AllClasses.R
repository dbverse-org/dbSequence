#' @import methods
#' @importFrom duckdb duckdb
#' @importClassesFrom dbProject dbData
#' @importFrom dbProject conn
#' @importFrom dbProject conn<-
NULL

# ------------------------------------------------------------------------------
#  File classes for method dispatch
# ------------------------------------------------------------------------------

#' @title DuckDBFile objects
#' @name DuckDBFile-class
#'
#' @description File class for DuckDB database files, extending BiocFile
#'
#' This mirrors the *File* classes in `rtracklayer`, giving us
#' a concrete object to dispatch `import()` / `export()` methods **to**.
#' Extends BiocFile to be compatible with BiocIO import/export framework.
#'
#' @slot path Character string specifying the path to the DuckDB database file.
#'   Can be a file path or ":memory:" for in-memory databases.
#'
#' @export
setClass("DuckDBFile", contains = "BiocFile", slots = c(path = "character"))

# ------------------------------------------------------------------------------
#  Main data classes
# ------------------------------------------------------------------------------

#' @title dbSequence objects
#' @name dbSequence-class
#'
#' @description S4 class for genomic sequences stored in DuckDB
#'
#' Extends dbData to wrap genomic data stored in DuckDB tables.
#'
#' @slot value dplyr tbl representing the genomic data in the database (inherited from dbData)
#' @slot name character table name in database (inherited from dbData)
#' @slot file_source character source file path or identifier (immutable after creation)
#' @export
setClass(
  "dbSequence",
  contains = "dbData",
  slots = c(
    file_source = "character"
  )
)

# Valid file extensions for genomic data
.VALID_FILE_EXTENSIONS <- c(
  "bed",
  "vcf",
  "gff",
  "gff3",
  "gtf",
  "bam",
  "cram",
  "fasta",
  "fa",
  "fastq",
  "fq",
  "narrowpeak",
  "broadpeak"
)

#' Validity check for dbSequence objects
#' @noRd
setValidity("dbSequence", function(object) {
  errors <- character()

  # Check file_source
  file_source <- object@file_source
  if (length(file_source) != 1) {
    errors <- c(errors, "file_source must be a single character string")
  } else if (!is.na(file_source) && file_source != "") {
    # If not NA or empty, check if it's a valid file path or has valid extension
    if (!file.exists(file_source)) {
      # Not a file path - check if it looks like a valid file type indicator
      ext <- tolower(tools::file_ext(file_source))
      # Handle compressed extensions
      if (ext %in% c("gz", "bgz", "bz2")) {
        base <- tools::file_path_sans_ext(file_source)
        ext <- tolower(tools::file_ext(base))
      }
      if (!ext %in% .VALID_FILE_EXTENSIONS && ext != "") {
        errors <- c(
          errors,
          paste0(
            "file_source must be a valid file path or have a recognized extension. ",
            "Got: '",
            file_source,
            "'. ",
            "Valid extensions: ",
            paste(.VALID_FILE_EXTENSIONS, collapse = ", ")
          )
        )
      }
    }
  }

  if (length(errors) == 0) TRUE else errors
})
