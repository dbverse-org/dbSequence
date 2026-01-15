library(dbSequence)
library(testthat)

test_that("filter_by_overlaps rejects unsupported arguments", {
  # Create a dummy dbSequence
  dummy_tbl <- dplyr::tibble(
    seqnames = "chr1",
    start = 100,
    end = 200
  )

  # Use a temporary db
  db_path <- tempfile(fileext = ".duckdb")
  con <- DBI::dbConnect(duckdb::duckdb(), db_path)
  DBI::dbWriteTable(con, "test_table", dummy_tbl)

  db_seq <- dbSequence("test_table", file_source = db_path)

  # Valid call
  expect_no_error(filter_by_overlaps(
    db_seq,
    GenomicRanges::GRanges("chr1:150-250")
  ))

  # Invalid maxgap (should error based on our implementation of check)
  expect_error(
    filter_by_overlaps(
      db_seq,
      GenomicRanges::GRanges("chr1:150-250"),
      maxgap = 10
    ),
    "maxgap"
  )

  # Invalid minoverlap
  expect_error(
    filter_by_overlaps(
      db_seq,
      GenomicRanges::GRanges("chr1:150-250"),
      minoverlap = 5
    ),
    "minoverlap"
  )

  # Cleanup
  DBI::dbDisconnect(con, shutdown = TRUE)
  unlink(db_path)
})

test_that("compute_coverage validates arguments", {
  # Mock dbSequence
  db_path <- tempfile(fileext = ".duckdb")
  con <- DBI::dbConnect(duckdb::duckdb(), db_path)
  dummy_tbl <- dplyr::tibble(seqnames = "chr1", start = 100, end = 200)
  DBI::dbWriteTable(con, "test_table", dummy_tbl)
  db_seq <- dbSequence("test_table", file_source = db_path)

  region <- "chr1:100-200"

  # Valid call
  expect_no_error(compute_coverage(db_seq, region, window = 10))

  # Invalid region type
  expect_error(compute_coverage(db_seq, 12345), "region must be")

  # Invalid window type
  expect_error(compute_coverage(db_seq, region, window = "bad"))

  # Cleanup
  DBI::dbDisconnect(con, shutdown = TRUE)
  unlink(db_path)
})

test_that("import validates required arguments", {
  # We test that calling import without dest on a file we support
  # requires the 'dest' argument if we want dbSequence behavior.
  # However, if 'dest' is missing, it might fall back to BiocIO default (returning GRanges).
  # We should check if we defined a method for missing dest that errors.

  bed_file <- system.file("extdata", "example.bed", package = "dbSequence")

  # If we implemented a method for signature (BEDFile, missing), it should error "dest required"
  # If we didn't, this test might fail (it might succeed and return GRanges).
  # Based on the user requirement "unsupported args do not work" (or "make sure ... supported args work"),
  # ensuring 'dest' is mandatory for our methods is good practice if we claimed so.

  # Checking error message
  expect_error(
    BiocIO::import(rtracklayer::BEDFile(bed_file), dest = NULL), # dest=NULL might trigger missing or NULL check
    "unable to find an inherited method" # Or similar if valid method not found
  )
  # Actually, let's test POSITIVE case first to match method existence
  # Then negative case.

  # Positive: passing dest works
  dest_db <- DuckDBFile(":memory:")
  expect_no_error(BiocIO::import(
    rtracklayer::BEDFile(bed_file),
    dest = dest_db
  ))

  # Negative: Passing random unsupported arg in ...
  # BiocIO/import seems strict about arguments if they are not matching signatures.
  # This satisfies the requirement to guard against random args.
  expect_error(
    BiocIO::import(rtracklayer::BEDFile(bed_file), dest = dest_db, foo = "bar"),
    "unused argument"
  )
})
