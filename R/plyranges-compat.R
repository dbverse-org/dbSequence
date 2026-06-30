# plyranges-compat.R -- plyranges API compatibility layer
# ------------------------------------------------------------------------------
# Provides plyranges-style function names for dbSequence users familiar with
# the plyranges package. These are thin wrappers around existing BiocIO import
# methods.
#
# Key Design Decision:
# - plyranges read_* functions return GRanges
# - dbSequence read_* functions return dbSequence (lazy, DuckDB-backed)
# - Users can chain into as_granges() to collect when needed

# ------------------------------------------------------------------------------
#  as_granges S3 generic and method
# ------------------------------------------------------------------------------

#' Convert objects to GRanges
#'
#' @description S3 generic for converting various objects to GRanges.
#'   This generic is provided for plyranges compatibility when plyranges
#'   is not loaded.
#'
#' @param .data An object to convert to GRanges
#' @param ... Additional arguments passed to methods
#' @param keep_mcols logical: whether to keep metadata columns (default: TRUE)
#'
#' @return A GRanges object
#'
#' @export
as_granges <- function(.data, ..., keep_mcols = TRUE) {
  UseMethod("as_granges")
}

#' @rdname as_granges
#' @export
as_granges.default <- function(.data, ..., keep_mcols = TRUE) {

  # If plyranges is loaded, defer to it
  if (requireNamespace("plyranges", quietly = TRUE)) {
    return(plyranges::as_granges(.data, ..., keep_mcols = keep_mcols))
  }
  stop("No as_granges method available for class: ", class(.data)[1])
}

#' Convert dbSequence to GRanges
#'
#' @description S3 method for plyranges compatibility. Converts a dbSequence
#'   object to a GRanges object by collecting data from DuckDB.
#'
#' @param .data A dbSequence object
#' @param ... Additional arguments (currently unused)
#' @param keep_mcols logical: whether to keep metadata columns (default: TRUE)
#'
#' @return A GRanges object
#'
#' @details This method enables plyranges-style workflows:
#'   \code{db_seq \%>\% as_granges()}
#'
#'   Note that this collects all data from DuckDB into memory. For large
#'   datasets, consider filtering first using dplyr verbs.
#'
#' @examples
#' # Import BED file to dbSequence
#' bed <- system.file("extdata", "example.bed", package = "dbSequence")
#' db_seq <- read_bed(bed)
#'
#' # Convert to GRanges (collects data)
#' gr <- as_granges(db_seq)
#'
#' @rdname as_granges
#' @export
as_granges.dbSequence <- function(.data, ..., keep_mcols = TRUE) {
  # Use existing S4 coercion which handles the conversion properly
  as(.data, "GRanges")
}

# ------------------------------------------------------------------------------
#  File readers (return dbSequence, not GRanges)
# ------------------------------------------------------------------------------

#' Read a BED file into DuckDB
#'
#' @description plyranges-style function to read BED files. Unlike the
#'   plyranges version which returns GRanges, this returns a dbSequence
#'   object backed by DuckDB for lazy evaluation.
#'
#' @param file Path to BED file
#' @param dest DuckDBFile or path: destination database (default: in-memory)
#' @param table_name character: name for the table (default: "bed_data")
#' @param lazy logical: if TRUE (default), do not import/copy the file into a
#'   DuckDB table. Instead, return a dbSequence backed by a DuckDB scan query
#'   (an inline SQL SELECT over the file). Call \code{dplyr::compute()} on the
#'   result to materialize when needed.
#' @param ... Additional arguments passed to import()
#'
#' @return A dbSequence object
#'
#' @examples
#' # Read BED file (lazy, stays in DuckDB)
#' bed <- system.file("extdata", "example.bed", package = "dbSequence")
#' db_seq <- read_bed(bed)
#'
#' # Collect to GRanges when needed
#' gr <- as_granges(db_seq)
#'
#' @seealso \code{\link[BiocIO]{import}}, \code{\link[rtracklayer]{BEDFile}}
#' @export
read_bed <- function(
  file,
  dest = DuckDBFile(":memory:"),
  table_name = "bed_data",
  lazy = TRUE,
  ...
) {
  if (isTRUE(lazy)) {
    return(.open_lazy_to_dbsequence(file, dest_db = dest, table_name = table_name, file_type = "bed"))
  }
  import(rtracklayer::BEDFile(file), dest = dest, table_name = table_name, ...)
}

#' Read a GFF/GTF file into DuckDB
#'
#' @description plyranges-style function to read GFF files. Returns a
#'   dbSequence object backed by DuckDB.
#'
#' @param file Path to GFF file
#' @param dest DuckDBFile or path: destination database (default: in-memory)
#' @param table_name character: name for the table (default: "gff_data")
#' @param lazy logical: if TRUE (default), do not import/copy the file into a
#'   DuckDB table. Instead, return a dbSequence backed by a DuckDB scan query
#'   (an inline SQL SELECT over the file). Call \code{dplyr::compute()} on the
#'   result to materialize when needed.
#' @param ... Additional arguments passed to import()
#'
#' @return A dbSequence object
#'
#' @examples
#' gff <- system.file("extdata", "example.gff3", package = "dbSequence")
#' db_seq <- read_gff(gff)
#' db_seq
#'
#' @seealso \code{\link[BiocIO]{import}}, \code{\link[rtracklayer]{GFFFile}}
#' @export
read_gff <- function(
  file,
  dest = DuckDBFile(":memory:"),
  table_name = "gff_data",
  lazy = TRUE,
  ...
) {
  if (isTRUE(lazy)) {
    return(.open_lazy_to_dbsequence(file, dest_db = dest, table_name = table_name, file_type = "gff"))
  }
  import(rtracklayer::GFFFile(file), dest = dest, table_name = table_name, ...)
}

#' @rdname read_gff
#' @export
read_gff3 <- function(
  file,
  dest = DuckDBFile(":memory:"),
  table_name = "gff_data",
  lazy = TRUE,
  ...
) {
  read_gff(file, dest = dest, table_name = table_name, lazy = lazy, ...)
}

#' Read a BAM file into DuckDB
#'
#' @description plyranges-style function to read BAM files. Returns a
#'   dbSequence object backed by DuckDB.
#'
#' @param file Path to BAM file
#' @param dest DuckDBFile or path: destination database (default: in-memory)
#' @param table_name character: name for the table (default: "alignments")
#' @param ... Additional arguments passed to import()
#'
#' @return A dbSequence object
#'
#' @note Unlike plyranges::read_bam which returns a DeferredGenomicRanges,
#'   this function immediately imports the BAM data into DuckDB.
#'
#' @seealso \code{\link[BiocIO]{import}}, \code{\link[Rsamtools]{BamFile}}
#' @examples
#' if (nzchar(system.file(package = "exonr"))) {
#'   bam <- system.file("extdata", "example.bam", package = "dbSequence")
#'   db_seq <- read_bam(bam, dest = DuckDBFile(tempfile(fileext = ".duckdb")))
#'   db_seq
#' }
#'
#' @export
read_bam <- function(
  file,
  dest = DuckDBFile(":memory:"),
  table_name = "alignments",
  ...
) {
  import(Rsamtools::BamFile(file), dest = dest, table_name = table_name, ...)
}

#' Read a VCF file into DuckDB
#'
#' @description Read VCF variant files into DuckDB. Returns a dbSequence
#'   object.
#'
#' @param file Path to VCF file
#' @param dest DuckDBFile or path: destination database (default: in-memory)
#' @param table_name character: name for the table (default: "variants")
#' @param lazy logical: if TRUE (default), do not import/copy the file into a
#'   DuckDB table. Instead, return a dbSequence backed by a DuckDB scan query
#'   (an inline SQL SELECT over the file). Call \code{dplyr::compute()} on the
#'   result to materialize when needed.
#' @param ... Additional arguments passed to import()
#'
#' @return A dbSequence object
#'
#' @examples
#' vcf <- system.file("extdata", "example.vcf", package = "dbSequence")
#' db_seq <- read_vcf(vcf)
#' db_seq
#'
#' @seealso \code{\link[BiocIO]{import}}, \code{\link[VariantAnnotation]{VcfFile}}
#' @export
read_vcf <- function(
  file,
  dest = DuckDBFile(":memory:"),
  table_name = "variants",
  lazy = TRUE,
  ...
) {
  if (isTRUE(lazy)) {
    return(.open_lazy_to_dbsequence(file, dest_db = dest, table_name = table_name, file_type = "vcf"))
  }
  import(VariantAnnotation::VcfFile(file), dest = dest, table_name = table_name, ...)
}
