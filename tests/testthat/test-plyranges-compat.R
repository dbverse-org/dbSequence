# test-plyranges-compat.R -- Tests for plyranges compatibility layer
# ------------------------------------------------------------------------------

library(dbSequence)
library(testthat)

# ==============================================================================
#  as_granges tests
# ==============================================================================

test_that("as_granges is exported", {
  expect_true(exists("as_granges", where = "package:dbSequence"))
})

test_that("as_granges.dbSequence returns GRanges", {
  skip_if_not_installed("GenomicRanges")

  # Create test database
  test_db <- tempfile(fileext = ".duckdb")
  con <- DBI::dbConnect(duckdb::duckdb(), test_db)
  test_data <- data.frame(
    seqnames = c("chr1", "chr1", "chr2"),
    start = c(100, 200, 100),
    end = c(150, 250, 150),
    strand = c("+", "-", "+")
  )
  DBI::dbWriteTable(con, "ranges", test_data)
  DBI::dbDisconnect(con, shutdown = TRUE)

  db_seq <- dbSequence("ranges", test_db)
  on.exit(file.remove(test_db))

  gr <- as_granges(db_seq)

  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), 3)
})

# ==============================================================================
#  read_* function tests
# ==============================================================================

test_that("read_bed is exported", {
  expect_true(exists("read_bed", where = "package:dbSequence"))
})

test_that("read_gff is exported", {
  expect_true(exists("read_gff", where = "package:dbSequence"))
})

test_that("read_bam is exported", {
  expect_true(exists("read_bam", where = "package:dbSequence"))
})

test_that("read_vcf is exported", {
  expect_true(exists("read_vcf", where = "package:dbSequence"))
})

test_that("read_bed returns dbSequence object", {
  bed_file <- system.file("extdata", "example.bed", package = "dbSequence")
  skip_if(!file.exists(bed_file), "No example BED file in package")

  db_seq <- read_bed(bed_file)

  expect_s4_class(db_seq, "dbSequence")
  expect_true(inherits(db_seq@value, "tbl_lazy"))
})

test_that("read_bed is lazy by default and can be materialized", {
  bed_file <- system.file("extdata", "example.bed", package = "dbSequence")
  skip_if(!file.exists(bed_file), "No example BED file in package")

  temp_db <- tempfile(fileext = ".duckdb")
  on.exit({
    if (file.exists(temp_db)) file.remove(temp_db)
  })

  db_seq <- read_bed(
    bed_file,
    dest = DuckDBFile(temp_db),
    table_name = "bed_data"
  )

  expect_s4_class(db_seq, "dbSequence")
  expect_true(inherits(db_seq@value, "tbl_lazy"))

  # Lazy scan should not create a physical table by default
  con <- dbplyr::remote_con(db_seq@value)
  expect_false(DBI::dbExistsTable(con, "bed_data"))

  collected <- dplyr::collect(db_seq@value)
  expect_true(all(c("seqnames", "start", "end") %in% colnames(collected)))
  expect_gt(nrow(collected), 0)

  # Materialize explicitly
  db_seq2 <- dplyr::compute(db_seq, name = "bed_data", temporary = FALSE)
  con2 <- dbplyr::remote_con(db_seq2@value)
  expect_true(DBI::dbExistsTable(con2, "bed_data"))
})

test_that("read_vcf is lazy by default and can be materialized", {
  vcf_file <- tempfile(fileext = ".vcf")
  on.exit({
    if (file.exists(vcf_file)) file.remove(vcf_file)
  }, add = TRUE)

  vcf_lines <- c(
    "##fileformat=VCFv4.2",
    "##source=dbSequence-test",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\trs1\tA\tC\t60\tPASS\tDP=10"
  )
  writeLines(vcf_lines, vcf_file)

  temp_db <- tempfile(fileext = ".duckdb")
  on.exit({
    if (file.exists(temp_db)) file.remove(temp_db)
  }, add = TRUE)

  db_seq <- read_vcf(
    vcf_file,
    dest = DuckDBFile(temp_db),
    table_name = "vcf_data"
  )

  expect_s4_class(db_seq, "dbSequence")
  expect_true(inherits(db_seq@value, "tbl_lazy"))

  con <- dbplyr::remote_con(db_seq@value)
  expect_false(DBI::dbExistsTable(con, "vcf_data"))

  df <- dplyr::collect(db_seq@value)
  expect_true(all(c("seqnames", "start") %in% colnames(df)))
  expect_gt(nrow(df), 0)

  db_seq2 <- dplyr::compute(db_seq, name = "vcf_data", temporary = FALSE)
  con2 <- dbplyr::remote_con(db_seq2@value)
  expect_true(DBI::dbExistsTable(con2, "vcf_data"))
})
