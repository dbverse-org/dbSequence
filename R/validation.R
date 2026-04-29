# validation.R -- File type validation and slot protection for dbSequence
# ------------------------------------------------------------------------------

# File type definitions
.COMPRESSED_FILE_TYPES <- c("bam", "cram", "bz2", "gz", "bgz", "lz4", "zst")
.SUPPORTED_RANGE_OPS <- c("bed", "vcf", "gff", "gff3", "gtf", "narrowpeak", "broadpeak")

#' @keywords internal
#' @noRd
.is_compressed_source <- function(x) {
  if (!is(x, "dbSequence")) return(FALSE)

  file_source <- x@file_source
  if (is.null(file_source) || is.na(file_source) || file_source == "") {
    return(FALSE)
  }

  # Get file extension (handle .gz, .bgz suffixes)
  ext <- tolower(tools::file_ext(file_source))
  if (ext %in% .COMPRESSED_FILE_TYPES) {
    return(TRUE)
  }

  # Check for double extensions like .vcf.gz
  base_without_gz <- tools::file_path_sans_ext(file_source)
  inner_ext <- tolower(tools::file_ext(base_without_gz))

  # BAM and CRAM are always compressed
  if (inner_ext %in% c("bam", "cram") || ext %in% c("bam", "cram")) {
    return(TRUE)
  }

  FALSE
}

#' @keywords internal
#' @noRd
.get_file_type <- function(x) {
  if (!is(x, "dbSequence")) return(NA_character_)

  file_source <- x@file_source
  if (is.null(file_source) || is.na(file_source) || file_source == "") {
    return(NA_character_)
  }

  # Get file extension
  ext <- tolower(tools::file_ext(file_source))

  # Handle compressed files - get the inner extension
  if (ext %in% c("gz", "bgz", "bz2")) {
    base_without_gz <- tools::file_path_sans_ext(file_source)
    ext <- tolower(tools::file_ext(base_without_gz))
  }

  ext
}

#' @keywords internal
#' @noRd
.validate_for_range_ops <- function(x, operation = "range operations") {
  if (!is(x, "dbSequence")) {
    cli::cli_abort(c(
      "x" = "{.arg x} must be a {.cls dbSequence} object",
      "i" = "Got {.cls {class(x)[1]}} instead"
    ))
  }

  file_type <- .get_file_type(x)

  # BAM and CRAM are always problematic
  if (!is.na(file_type) && file_type %in% c("bam", "cram")) {
    cli::cli_abort(c(
      "x" = "{.fn {operation}} is not supported for {.field {toupper(file_type)}} files",
      "i" = "BAM/CRAM files are compressed and indexed differently",
      "i" = "For coverage on BAM files, use specialized tools like {.pkg Rsamtools}",
      "!" = "Source file: {.file {x@file_source}}"
    ))
  }

  invisible(TRUE)
}

#' @keywords internal
#' @noRd
.validate_file_type <- function(x, operation, allowed_types = .SUPPORTED_RANGE_OPS) {
  file_type <- .get_file_type(x)

  if (!is.na(file_type) && !file_type %in% allowed_types) {
    cli::cli_abort(c(
      "x" = "{.fn {operation}} is not supported for {.field {toupper(file_type)}} files",
      "i" = "Supported file types: {.field {toupper(allowed_types)}}",
      "!" = "Source file: {.file {x@file_source}}"
    ))
  }

  invisible(TRUE)
}

# ------------------------------------------------------------------------------
#  Slot protection
# ------------------------------------------------------------------------------

#' Accessor for file_source slot
#'
#' @param object A dbSequence object
#' @param value New value (will be rejected)
#' @return The file_source value
#' @examples
#' bed <- system.file("extdata", "example.bed", package = "dbSequence")
#' db_seq <- read_bed(bed)
#' fileSource(db_seq)
#' try(fileSource(db_seq) <- "other.bed")
#'
#' @export
setGeneric("fileSource", function(object) standardGeneric("fileSource"))

#' @rdname fileSource
#' @export
setMethod("fileSource", "dbSequence", function(object) {
  object@file_source
})

#' Replacement method for file_source - BLOCKED
#'
#' @description This method prevents modification of the file_source slot
#'   after object creation. The file_source is set during import and should
#'   not be changed to maintain data integrity.
#'
#' @param object A dbSequence object
#' @param value New value (will be rejected)
#' @return Always errors because file sources are immutable after creation.
#' @examples
#' bed <- system.file("extdata", "example.bed", package = "dbSequence")
#' db_seq <- read_bed(bed)
#' try(fileSource(db_seq) <- "other.bed")
#'
#' @export
setGeneric("fileSource<-", function(object, value) standardGeneric("fileSource<-"))

#' @rdname fileSource
#' @export
setReplaceMethod("fileSource", "dbSequence", function(object, value) {
  cli::cli_abort(c(
    "x" = "Cannot modify {.field file_source} after object creation",
    "i" = "The file source is set during {.fn import} or {.fn read_*} and is immutable",
    "i" = "Current value: {.file {object@file_source}}",
    ">" = "If you need a different source, create a new {.cls dbSequence} object"
  ))
})
