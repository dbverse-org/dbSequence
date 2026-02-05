# import-export.R -- BiocIO-style import/export with DuckDB backend
# ------------------------------------------------------------------------------
# Implements BiocIO import() generic methods for genomic files with DuckDB storage
# using existing internal import_*_to_duckdb functions

#' @importFrom methods setMethod
#' @importFrom BiocIO import export resource path
#' @importClassesFrom BiocIO BiocFile
#' @importFrom rtracklayer BEDFile GTFFile GFFFile FastaFile
#' @importFrom Rsamtools BamFile TabixFile
#' @importFrom VariantAnnotation VcfFile
NULL

# ------------------------------------------------------------------------------
#  BiocIO import() methods - correct signature pattern
# ------------------------------------------------------------------------------

#' Import BAM files into DuckDB
setMethod(
  "import",
  signature(con = "BamFile", format = "missing", text = "ANY"),
  function(
    con,
    format,
    text,
    dest,
    table_name = "bam_data",
    .conn = NULL,
    ...
  ) {
    if (missing(dest)) {
      stop("'dest' (a DuckDBFile) is required")
    }

    .import_bam_to_duckdb(
      path(con),
      resource(dest),
      table_name = table_name,
      .conn = .conn,
      ...
    )
  }
)

#' Import VCF files into DuckDB
setMethod(
  "import",
  signature(con = "VcfFile", format = "missing", text = "ANY"),
  function(con, format, text, dest, table_name = "vcf_data", ...) {
    if (missing(dest)) {
      stop("'dest' (a DuckDBFile) is required")
    }

    .import_vcf_to_duckdb(
      path(con),
      resource(dest),
      table_name = table_name,
      ...
    )
  }
)

#' Import BED files into DuckDB
setMethod(
  "import",
  signature(con = "BEDFile", format = "missing", text = "ANY"),
  function(con, format, text, dest, table_name = "bed_data", ...) {
    if (missing(dest)) {
      stop("'dest' (a DuckDBFile) is required")
    }

    .import_bed_to_duckdb(
      resource(con),
      resource(dest),
      table_name = table_name,
      ...
    )
  }
)

#' Import GTF files into DuckDB
setMethod(
  "import",
  signature(con = "GTFFile", format = "missing", text = "ANY"),
  function(con, format, text, dest, table_name = "gtf_data", ...) {
    if (missing(dest)) {
      stop("'dest' (a DuckDBFile) is required")
    }

    .import_gtf_to_duckdb(
      resource(con),
      resource(dest),
      table_name = table_name,
      ...
    )
  }
)

#' Import GFF files into DuckDB
setMethod(
  "import",
  signature(con = "GFFFile", format = "missing", text = "ANY"),
  function(con, format, text, dest, table_name = "gff_data", ...) {
    if (missing(dest)) {
      stop("'dest' (a DuckDBFile) is required")
    }

    .import_gff_to_duckdb(
      resource(con),
      resource(dest),
      table_name = table_name,
      ...
    )
  }
)

#' Import character file paths into DuckDB (auto-detect format)
setMethod(
  "import",
  signature(con = "character", format = "missing", text = "ANY"),
  function(con, format, text, dest, table_name = NULL, ...) {
    if (missing(dest)) {
      stop("'dest' (a DuckDBFile) is required")
    }

    # Auto-detect file format from extension
    ext <- tools::file_ext(con)

    if (is.null(table_name)) {
      table_name <- paste0(ext, "_data")
    }

    switch(
      tolower(ext),
      bed = .import_bed_to_duckdb(con, resource(dest), table_name, ...),
      gtf = .import_gtf_to_duckdb(con, resource(dest), table_name, ...),
      gff = .import_gff_to_duckdb(con, resource(dest), table_name, ...),
      gff3 = .import_gff_to_duckdb(con, resource(dest), table_name, ...),
      vcf = .import_vcf_to_duckdb(con, resource(dest), table_name, ...),
      bam = .import_bam_to_duckdb(con, resource(dest), table_name, ...),
      fasta = .import_fasta_to_duckdb(con, resource(dest), table_name, ...),
      fa = .import_fasta_to_duckdb(con, resource(dest), table_name, ...),
      stop("Unsupported file format: ", ext)
    )
  }
)

#' Import FASTA files into DuckDB
setMethod(
  "import",
  signature(con = "FastaFile", format = "missing", text = "ANY"),
  function(con, format, text, dest, table_name = "fasta_data", ...) {
    if (missing(dest)) {
      stop("'dest' (a DuckDBFile) is required")
    }

    .import_fasta_to_duckdb(
      resource(con),
      resource(dest),
      table_name = table_name,
      ...
    )
  }
)

# ------------------------------------------------------------------------------
#  BiocIO import() method implementations
# ------------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.import_to_duckdb <- function(
  file_path,
  dest_db,
  table_name,
  file_type = c("bam", "vcf", "bed", "gff", "fasta", "fastq"),
  exon_options = list(),
  .conn = NULL
) {
  file_type <- match.arg(file_type)

  # Handle DuckDBFile object - extract path only
  if (is(dest_db, "DuckDBFile")) {
    db_path <- dest_db@path
  } else {
    db_path <- as.character(dest_db)
  }

  # Validate file exists
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }

  # Use provided connection or create a new one
  if (!is.null(.conn)) {
    # Use the provided connection
    con <- .conn
    should_close <- FALSE
  } else {
    # Create connection to destination database
    con <- .get_duckdb_connection(db_path)
    # For in-memory databases, we need to keep connection alive
    # For file databases, we can close after creating the table
    should_close <- (db_path != ":memory:")
  }

  if (should_close) {
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  }

  # Try to use exon if available, otherwise fall back to simple methods
  if (requireNamespace("exonr", quietly = TRUE)) {
    .import_with_exon(file_path, con, table_name, file_type, exon_options)
  } else {
    .import_with_duckdb(file_path, con, table_name, file_type, exon_options)
  }

  # Verify table was created
  if (!DBI::dbExistsTable(con, table_name)) {
    stop("Failed to create table: ", table_name)
  }

  # Return dbSequence object - for in-memory databases, pass the connection
  if (db_path == ":memory:") {
    # Pass the connection for in-memory databases to avoid creating separate databases
    if (is(dest_db, "DuckDBFile")) {
      dbSequence(table_name, dest_db, .conn = con)
    } else {
      dbSequence(table_name, db_path, .conn = con)
    }
  } else {
    # For file databases, let constructor create its own connection
    dbSequence(table_name, dest_db)
  }
}

# ------------------------------------------------------------------------------
#  Lazy file scans (no materialization)
# ------------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.open_lazy_to_dbsequence <- function(
  file_path,
  dest_db,
  table_name,
  file_type = c("bed", "gff", "vcf")
) {
  file_type <- match.arg(file_type)

  # Handle DuckDBFile object - extract path only
  if (is(dest_db, "DuckDBFile")) {
    db_path <- dest_db@path
  } else {
    db_path <- as.character(dest_db)
  }

  # Validate file exists
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }

  # Establish a connection that stays alive via the returned tbl
  con <- .get_duckdb_connection(db_path)

  file_path_norm <- normalizePath(file_path, winslash = "/", mustWork = TRUE)
  quoted_path <- DBI::dbQuoteString(con, file_path_norm)

  base_sql <- switch(
    file_type,
    "bed" = glue::glue(
      "SELECT * FROM read_csv_auto({quoted_path}, header=false, delim='\t', comment='#', ignore_errors=true, null_padding=true)"
    ),
    "gff" = glue::glue(
      "SELECT * FROM read_csv_auto({quoted_path}, header=false, delim='\t', comment='#', ignore_errors=true, null_padding=true)"
    ),
    "vcf" = glue::glue(
      "SELECT * FROM read_csv_auto({quoted_path}, header=false, delim='\t', comment='#', ignore_errors=true, null_padding=true)"
    )
  )

  tbl_value <- dplyr::tbl(con, dbplyr::sql(base_sql))

  rename_map <- switch(
    file_type,
    "bed" = c(
      seqnames = "column0",
      start = "column1",
      end = "column2",
      name = "column3",
      score = "column4",
      strand = "column5"
    ),
    "gff" = c(
      seqnames = "column0",
      source = "column1",
      type = "column2",
      start = "column3",
      end = "column4",
      score = "column5",
      strand = "column6",
      phase = "column7",
      attributes = "column8"
    ),
    "vcf" = c(
      seqnames = "column0",
      start = "column1",
      id = "column2",
      ref = "column3",
      alt = "column4",
      qual = "column5",
      filter = "column6",
      info = "column7"
    )
  )

  keep <- rename_map %in% colnames(tbl_value)
  if (any(keep)) {
    rename_syms <- rlang::set_names(
      lapply(rename_map[keep], rlang::sym),
      names(rename_map[keep])
    )
    tbl_value <- dplyr::rename(tbl_value, !!!rename_syms)
  }

  # Standardize BED coordinates: input BED is 0-based, half-open.
  # Convert to 1-based, closed intervals to match Bioconductor GRanges.
  if (identical(file_type, "bed") && all(c("start", "end") %in% colnames(tbl_value))) {
    tbl_value <- tbl_value |>
      dplyr::mutate(
        start = as.integer(start) + 1L,
        end = as.integer(end)
      )
  }

  methods::new(
    "dbSequence",
    value = tbl_value,
    name = table_name,
    file_source = db_path
  )
}

# ------------------------------------------------------------------------------
#  Convenience functions for specific file types
# ------------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.import_bam_to_duckdb <- function(
  bam_path,
  dest_db,
  table_name = "alignments",
  parse_tags = FALSE,
  mapq_min = NULL,
  chromosome = NULL,
  start = NULL,
  end = NULL,
  .conn = NULL
) {
  # Build filter clauses
  filters <- c()
  if (!is.null(mapq_min)) {
    filters <- c(filters, paste0("mapq >= ", mapq_min))
  }
  if (!is.null(chromosome)) {
    filters <- c(filters, paste0("ref_name = '", chromosome, "'"))
  }
  if (!is.null(start)) {
    filters <- c(filters, paste0("pos >= ", start))
  }
  if (!is.null(end)) {
    filters <- c(filters, paste0("pos <= ", end))
  }

  query_filter <- if (length(filters) > 0) {
    paste(filters, collapse = " AND ")
  } else {
    NULL
  }

  .import_to_duckdb(
    bam_path,
    dest_db,
    table_name,
    "bam",
    exon_options = list(parse_tags = parse_tags, query_filter = query_filter),
    .conn = .conn
  )
}

#' @keywords internal
#' @noRd
.import_vcf_to_duckdb <- function(
  vcf_path,
  dest_db,
  table_name = "variants",
  qual_min = NULL,
  chromosome = NULL
) {
  # Build filter clauses
  filters <- c()
  if (!is.null(qual_min)) {
    filters <- c(filters, paste0("qual >= ", qual_min))
  }
  if (!is.null(chromosome)) {
    # VCF files use 'chrom' in exon
    filters <- c(filters, paste0("chrom = '", chromosome, "'"))
  }

  query_filter <- if (length(filters) > 0) {
    paste(filters, collapse = " AND ")
  } else {
    NULL
  }

  .import_to_duckdb(
    vcf_path,
    dest_db,
    table_name,
    "vcf",
    exon_options = list(query_filter = query_filter)
  )
}

#' @keywords internal
#' @noRd
.import_bed_to_duckdb <- function(
  bed_path,
  dest_db,
  table_name = "intervals",
  score_min = NULL,
  chromosome = NULL
) {
  # Build filter clauses
  filters <- c()
  if (!is.null(score_min)) {
    filters <- c(filters, paste0("score >= ", score_min))
  }
  if (!is.null(chromosome)) {
    # BED files use 'reference_sequence_name' in exon
    filters <- c(
      filters,
      paste0("reference_sequence_name = '", chromosome, "'")
    )
  }

  query_filter <- if (length(filters) > 0) {
    paste(filters, collapse = " AND ")
  } else {
    NULL
  }

  res <- .import_to_duckdb(
    bed_path,
    dest_db,
    table_name,
    "bed",
    exon_options = list(query_filter = query_filter)
  )

  # Standardize BED coordinates: input BED is 0-based, half-open.
  # Convert to 1-based, closed intervals to match Bioconductor GRanges.
  if (is(res, "dbSequence")) {
    tbl_value <- res@value
    if (!is.null(tbl_value) && all(c("start", "end") %in% colnames(tbl_value))) {
      res@value <- tbl_value |>
        dplyr::mutate(
          start = as.integer(start) + 1L,
          end = as.integer(end)
        )
    }
  }

  res
}

#' @keywords internal
#' @noRd
.import_gtf_to_duckdb <- function(
  gtf_path,
  dest_db,
  table_name = "annotations"
) {
  .import_to_duckdb(gtf_path, dest_db, table_name, "gff") # GTF is handled as GFF
}

#' @keywords internal
#' @noRd
.import_gff_to_duckdb <- function(
  gff_path,
  dest_db,
  table_name = "annotations"
) {
  .import_to_duckdb(gff_path, dest_db, table_name, "gff")
}

#' @keywords internal
#' @noRd
.import_fasta_to_duckdb <- function(
  fasta_path,
  dest_db,
  table_name = "sequences"
) {
  .import_to_duckdb(fasta_path, dest_db, table_name, "fasta")
}

#' @keywords internal
#' @noRd
.import_fastq_to_duckdb <- function(
  fastq_path,
  dest_db,
  table_name = "reads"
) {
  .import_to_duckdb(fastq_path, dest_db, table_name, "fastq")
}

# ------------------------------------------------------------------------------
#  Helper functions
# ------------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.get_duckdb_connection <- function(db_path) {
  # Always create a new connection - for in-memory databases,
  # connection sharing will be handled by passing connections explicitly
  DBI::dbConnect(duckdb::duckdb(), dbdir = db_path)
}

#' @keywords internal
#' @noRd
.import_with_exon <- function(
  file_path,
  con,
  table_name,
  file_type,
  exon_options
) {
  # Create ExonRSessionContext for SQL-based processing with DataFusion
  session <- exonr::ExonRSessionContext$new()

  # Set exon options based on file type
  if (file_type == "bam" && isTRUE(exon_options$parse_tags)) {
    session$sql("SET exon.bam_parse_tags = true;")
  }

  # Build SQL query for the specific file type using exon scan functions
  scan_function <- switch(
    file_type,
    "bam" = "bam_scan",
    "vcf" = "vcf_scan",
    "bed" = "bed_scan",
    "gff" = "gff_scan",
    "fasta" = "fasta_scan",
    "fastq" = "fastq_scan",
    stop("Unsupported file type: ", file_type)
  )

  # Build the query - allow for optional filtering
  base_query <- paste0("SELECT * FROM ", scan_function, "('", file_path, "')")
  query_filter <- exon_options$query_filter
  if (!is.null(query_filter)) {
    query <- paste(base_query, "WHERE", query_filter)
  } else {
    query <- base_query
  }

  # Execute query using exon's DataFusion engine
  result <- session$sql(query)

  # Convert to RecordBatchReader for memory-efficient streaming
  batch_reader <- result$to_record_batch_reader()

  # Stream data to DuckDB using arrow if available
  # Use arrow to stream data to DuckDB - fail if arrow is not available
  tryCatch(
    {
      if (!requireNamespace("arrow", quietly = TRUE)) {
        stop("arrow package is required for streaming data to DuckDB")
      }
      arrow::to_duckdb(batch_reader, con = con, table_name = table_name)
    },
    error = function(e) {
      stop("Failed to stream genomic data to DuckDB: ", e$message)
    }
  )
}

#' @keywords internal
#' @noRd
.import_with_duckdb <- function(
  file_path,
  con,
  table_name,
  file_type,
  exon_options
) {
  # Simple fallback implementations for basic file types
  switch(
    file_type,
    "bed" = {
      # Use DuckDB's direct CSV reading for BED
      create_sql <- glue::glue(
        "
        CREATE TABLE {table_name} AS 
        SELECT * FROM read_csv_auto('{file_path}', 
          header=false, 
          delim='\t',
          comment='#',
          ignore_errors=true,
          null_padding=true
        )"
      )
      DBI::dbExecute(con, create_sql)

      # Rename standard BED columns if they have generic names
      fields <- DBI::dbListFields(con, table_name)
      bed_cols <- c("seqnames", "start", "end", "name", "score", "strand")
      
      # Rename up to the first 6 columns
      num_cols_to_rename <- min(length(fields), length(bed_cols))
      for (i in seq_len(num_cols_to_rename)) {
        current <- fields[i]
        target <- bed_cols[i]
        
        # Only rename if it looks like a generic name (column0, etc.)
        if (grepl("^column", current, ignore.case = TRUE) && current != target) {
           alter_sql <- glue::glue("ALTER TABLE {table_name} RENAME COLUMN \"{current}\" TO \"{target}\"")
           tryCatch(DBI::dbExecute(con, alter_sql), error=function(e) warning("Failed to rename BED column: ", e$message))
        }
      }
    },
    "gff" = {
      # Use DuckDB's direct CSV reading for GFF
      create_sql <- glue::glue(
        "
        CREATE TABLE {table_name} AS 
        SELECT * FROM read_csv_auto('{file_path}', 
          header=false, 
          delim='\t',
          comment='#',
          ignore_errors=true,
          null_padding=true
        )"
      )
      DBI::dbExecute(con, create_sql)

      # Rename standard GFF columns if they have generic names
      fields <- DBI::dbListFields(con, table_name)
      gff_cols <- c("seqnames", "source", "type", "start", "end",
                    "score", "strand", "phase", "attributes")

      # Rename up to the first 9 columns
      num_cols_to_rename <- min(length(fields), length(gff_cols))
      for (i in seq_len(num_cols_to_rename)) {
        current <- fields[i]
        target <- gff_cols[i]

        # Only rename if it looks like a generic name (column0, etc.)
        if (grepl("^column", current, ignore.case = TRUE) && current != target) {
          alter_sql <- glue::glue("ALTER TABLE {table_name} RENAME COLUMN \"{current}\" TO \"{target}\"")
          tryCatch(DBI::dbExecute(con, alter_sql), error = function(e) warning("Failed to rename GFF column: ", e$message))
        }
      }
    },
    "vcf" = {
      # Use DuckDB's direct CSV reading for VCF
      create_sql <- glue::glue(
        "
        CREATE TABLE {table_name} AS 
        SELECT * FROM read_csv_auto('{file_path}', 
          header=false, 
          delim='\t',
          comment='#',
          ignore_errors=true,
          null_padding=true
        )"
      )
      DBI::dbExecute(con, create_sql)

      # Rename standard VCF columns if they have generic names
      fields <- DBI::dbListFields(con, table_name)
      vcf_cols <- c("seqnames", "start", "id", "ref", "alt",
                    "qual", "filter", "info")

      # Rename up to the first 8 columns
      num_cols_to_rename <- min(length(fields), length(vcf_cols))
      for (i in seq_len(num_cols_to_rename)) {
        current <- fields[i]
        target <- vcf_cols[i]

        # Only rename if it looks like a generic name (column0, etc.)
        if (grepl("^column", current, ignore.case = TRUE) && current != target) {
          alter_sql <- glue::glue("ALTER TABLE {table_name} RENAME COLUMN \"{current}\" TO \"{target}\"")
          tryCatch(DBI::dbExecute(con, alter_sql), error = function(e) warning("Failed to rename VCF column: ", e$message))
        }
      }
    },
    {
      # For all other e.g. BAM, FASTA, FASTQ - require exonr
      stop(
        "File type '",
        file_type,
        "' requires exonr package for import. ",
        "Please install exonr for comprehensive genomic file support."
      )
    }
  )
}
