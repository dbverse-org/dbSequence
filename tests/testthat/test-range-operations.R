# test-range-operations.R -- Tests for range overlap and coverage functions
# ------------------------------------------------------------------------------

library(dbSequence)
library(testthat)

# Helper to create test database with fragments
.create_test_fragments <- function() {
  test_db <- tempfile(fileext = ".duckdb")
  con <- DBI::dbConnect(duckdb::duckdb(), test_db)

  # Create test data with known overlaps
  fragments <- data.frame(
    seqnames = c("chr1", "chr1", "chr1", "chr1", "chr2"),
    start = c(1000, 1500, 2000, 5000, 1000),
    end = c(1200, 1700, 2200, 5200, 1200),
    score = c(10, 20, 30, 40, 50)
  )
  DBI::dbWriteTable(con, "fragments", fragments)
  DBI::dbDisconnect(con, shutdown = TRUE)

  dbSequence("fragments", test_db)
}

# ==============================================================================
#  filter_by_overlaps tests
# ==============================================================================

test_that("filter_by_overlaps returns correct overlapping ranges", {
  db_seq <- .create_test_fragments()
  on.exit(file.remove(fileSource(db_seq)))

  # Region overlaps ranges at 1500-1700 only
  region <- GenomicRanges::GRanges("chr1:1400-1600")
  filtered <- filter_by_overlaps(db_seq, region)

  result <- dplyr::collect(dbSequence:::.dbseq_value(filtered))
  expect_equal(nrow(result), 1)
  expect_equal(result$start, 1500)
  expect_equal(result$end, 1700)
})

test_that("filter_by_overlaps handles empty results", {
  db_seq <- .create_test_fragments()
  on.exit(file.remove(fileSource(db_seq)))

  # Region does not overlap any ranges
  region <- GenomicRanges::GRanges("chr3:1-1000")
  filtered <- filter_by_overlaps(db_seq, region)

  result <- dplyr::collect(dbSequence:::.dbseq_value(filtered))
  expect_equal(nrow(result), 0)
})

test_that("filter_by_overlaps handles multiple overlaps", {
  db_seq <- .create_test_fragments()
  on.exit(file.remove(fileSource(db_seq)))

  # Region overlaps first three ranges on chr1
  region <- GenomicRanges::GRanges("chr1:1000-2200")
  filtered <- filter_by_overlaps(db_seq, region)

  result <- dplyr::collect(dbSequence:::.dbseq_value(filtered))
  expect_equal(nrow(result), 3)
})

test_that("filter_by_overlaps result is lazy (not collected)", {
  db_seq <- .create_test_fragments()
  on.exit(file.remove(fileSource(db_seq)))

  region <- GenomicRanges::GRanges("chr1:1400-1600")
  filtered <- filter_by_overlaps(db_seq, region)

  expect_true(inherits(dbSequence:::.dbseq_value(filtered), "tbl_lazy"))
  expect_s4_class(filtered, "dbSequence")
})

test_that("filter_by_overlaps rejects BAM files", {
  db_seq <- new("dbSequence", file_source = "test.bam")

  expect_error(
    filter_by_overlaps(db_seq, GenomicRanges::GRanges("chr1:1-100")),
    "not supported for BAM"
  )
})

# ==============================================================================
#  compute_coverage tests
# ==============================================================================

test_that("compute_coverage returns binned counts", {
  db_seq <- .create_test_fragments()
  on.exit(file.remove(fileSource(db_seq)))

  coverage <- compute_coverage(db_seq, "chr1:1000-2500", window = 500)

  expect_s3_class(coverage, "data.frame")
  expect_true("bin_start" %in% names(coverage))
  expect_true("count" %in% names(coverage))
  expect_true("bin_end" %in% names(coverage))
})

test_that("compute_coverage parses region strings", {
  db_seq <- .create_test_fragments()
  on.exit(file.remove(fileSource(db_seq)))

  # Both formats should work
  cov1 <- compute_coverage(db_seq, "chr1:1000-2500", window = 500)
  cov2 <- compute_coverage(db_seq, GenomicRanges::GRanges("chr1:1000-2500"), window = 500)

  expect_equal(nrow(cov1), nrow(cov2))
})

test_that("compute_coverage respects window size", {
  db_seq <- .create_test_fragments()
  on.exit(file.remove(fileSource(db_seq)))

  cov_500 <- compute_coverage(db_seq, "chr1:1000-3000", window = 500)
  cov_1000 <- compute_coverage(db_seq, "chr1:1000-3000", window = 1000)

  # Smaller window = more bins

  expect_gte(nrow(cov_500), nrow(cov_1000))
})

