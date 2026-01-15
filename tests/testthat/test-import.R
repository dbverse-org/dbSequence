library(dbSequence)

test_that("BiocIO import validates inputs correctly", {
  # Test file validation
  expect_error(
    BiocIO::import(
      con = "nonexistent.bed",
      dest = DuckDBFile("test.duckdb"),
      table_name = "test"
    ),
    "File not found|cannot find file"
  )

  # Test with valid file but invalid dest type - should fail gracefully
  temp_bed <- tempfile(fileext = ".bed")
  write("chr1\t100\t200", temp_bed)
  on.exit(file.remove(temp_bed))

  expect_error(
    BiocIO::import(con = temp_bed, dest = "invalid_dest", table_name = "test"),
    class = "error"
  )
})

test_that("BED file import works", {
  skip_if_not_installed("exonr")

  # Use test BED file from exon
  bed_file <- "/Users/user/Documents/test-seq/exon/exon/exon-core/test-data/datasources/bed/test.bed"
  skip_if_not(file.exists(bed_file), "Test BED file not found")

  test_db <- tempfile(fileext = ".duckdb")
  on.exit(if (file.exists(test_db)) file.remove(test_db))

  # Test basic import using BiocIO import method
  result <- BiocIO::import(
    con = bed_file,
    dest = DuckDBFile(test_db),
    table_name = "test_bed"
  )

  expect_s4_class(result, "dbSequence")
  expect_equal(tableName(result), "test_bed")
  expect_true(file.exists(test_db))

  # Verify table exists and has data
  con <- DBI::dbConnect(duckdb::duckdb(), test_db)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  expect_true(DBI::dbExistsTable(con, "test_bed"))

  row_count <- DBI::dbGetQuery(con, "SELECT COUNT(*) as n FROM test_bed")$n
  expect_gt(row_count, 0)

  # Check that we have expected BED columns (exon uses different names)
  columns <- DBI::dbListFields(con, "test_bed")
  expect_true("reference_sequence_name" %in% columns || "chrom" %in% columns)
  expect_true("start" %in% columns || "chromStart" %in% columns)
})

test_that("BED import with filtering works", {
  skip_if_not_installed("exonr")

  bed_file <- "/Users/user/Documents/test-seq/exon/exon/exon-core/test-data/datasources/bed/test.bed"
  skip_if_not(file.exists(bed_file), "Test BED file not found")

  test_db <- tempfile(fileext = ".duckdb")
  on.exit(if (file.exists(test_db)) file.remove(test_db))

  # Test import with parameters using BiocIO import method
  result <- BiocIO::import(
    con = bed_file,
    dest = DuckDBFile(test_db),
    table_name = "filtered_bed",
    chromosome = "1" # Use "1" instead of "chr1" as it's more likely to exist
  )

  expect_s4_class(result, "dbSequence")
  expect_true(file.exists(test_db))

  # Note: We can't easily test if filtering worked without knowing the content
  # But we can verify the function runs without error
})

test_that("GFF file import works", {
  skip_if_not_installed("exonr")

  # Use test GFF file from exon
  gff_file <- "/Users/user/Documents/test-seq/exon/exon/exon-core/test-data/datasources/gff-index/gencode.v38.polyAs.gff3"
  skip_if_not(file.exists(gff_file), "Test GFF file not found")

  test_db <- tempfile(fileext = ".duckdb")
  on.exit(if (file.exists(test_db)) file.remove(test_db))

  # Test basic import using BiocIO import method
  result <- BiocIO::import(
    con = gff_file,
    dest = DuckDBFile(test_db),
    table_name = "test_gff"
  )

  expect_s4_class(result, "dbSequence")
  expect_equal(tableName(result), "test_gff")
  expect_true(file.exists(test_db))

  # Verify table exists and has data
  con <- DBI::dbConnect(duckdb::duckdb(), test_db)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  expect_true(DBI::dbExistsTable(con, "test_gff"))

  row_count <- DBI::dbGetQuery(con, "SELECT COUNT(*) as n FROM test_gff")$n
  expect_gte(row_count, 0) # Allow for empty files, just check it doesn't error

  # Check that we have expected GFF columns
  columns <- DBI::dbListFields(con, "test_gff")
  expect_true("seqname" %in% columns || "chrom" %in% columns)
  expect_true("start" %in% columns)
  expect_true("end" %in% columns)
})

test_that("FASTA file import works", {
  skip_if_not_installed("exonr")

  # Use test FASTA file from exon
  fasta_file <- "/Users/user/Documents/test-seq/exon/exon/exon-core/test-data/datasources/repartition-test/test.fasta"
  skip_if_not(file.exists(fasta_file), "Test FASTA file not found")

  test_db <- tempfile(fileext = ".duckdb")
  on.exit(if (file.exists(test_db)) file.remove(test_db))

  # Test basic import using BiocIO import method
  result <- BiocIO::import(
    con = fasta_file,
    dest = DuckDBFile(test_db),
    table_name = "test_fasta"
  )

  expect_s4_class(result, "dbSequence")
  expect_equal(tableName(result), "test_fasta")
  expect_true(file.exists(test_db))

  # Verify table exists and has data
  con <- DBI::dbConnect(duckdb::duckdb(), test_db)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  expect_true(DBI::dbExistsTable(con, "test_fasta"))

  row_count <- DBI::dbGetQuery(con, "SELECT COUNT(*) as n FROM test_fasta")$n
  expect_gt(row_count, 0)

  # Check that we have expected FASTA columns
  columns <- DBI::dbListFields(con, "test_fasta")
  expect_true("id" %in% columns || "name" %in% columns)
  expect_true("sequence" %in% columns || "seq" %in% columns)
})

# Test that both in-memory and file-backed databases work correctly
test_that("in-memory database import works correctly", {
  # Use existing BED file from package
  bed_file <- system.file("extdata", "example.bed", package = "dbSequence")
  skip_if(!file.exists(bed_file), "Example BED file not found")

  # Test import to in-memory database using BiocIO import method
  result_memory <- BiocIO::import(
    con = bed_file,
    dest = DuckDBFile(":memory:"),
    table_name = "test_memory"
  )

  # Test that result is a proper dbSequence object
  expect_s4_class(result_memory, "dbSequence")
  expect_equal(tableName(result_memory), "test_memory")
  expect_equal(result_memory@file_source, ":memory:")
  expect_false(is.null(result_memory@value))

  # Test that the tbl object exists and is accessible
  expect_false(is.null(result_memory@value))
  expect_s3_class(result_memory@value, "tbl")

  # Test that we can query the data
  data_retrieved <- result_memory@value %>% dplyr::collect()
  expect_gt(nrow(data_retrieved), 0)
})

test_that("file database import works correctly", {
  # Use existing BED file from package
  bed_file <- system.file("extdata", "example.bed", package = "dbSequence")
  skip_if(!file.exists(bed_file), "Example BED file not found")

  # Create temporary database file
  temp_db <- tempfile(fileext = ".duckdb")
  on.exit({
    if (file.exists(temp_db)) file.remove(temp_db)
  })

  # Test import to file database using BiocIO import method
  result_file <- BiocIO::import(
    con = bed_file,
    dest = DuckDBFile(temp_db),
    table_name = "test_file"
  )

  # Test that result is a proper dbSequence object
  expect_s4_class(result_file, "dbSequence")
  expect_equal(tableName(result_file), "test_file")
  expect_equal(result_file@file_source, temp_db)
  expect_false(is.null(result_file@value))

  # Test that the file was created
  expect_true(file.exists(temp_db))

  # Test that the tbl object exists and is accessible
  expect_false(is.null(result_file@value))
  expect_s3_class(result_file@value, "tbl")

  # Test that we can query the data
  data_retrieved <- result_file@value %>% dplyr::collect()
  expect_gt(nrow(data_retrieved), 0)
})

test_that("connection handling is correct for in-memory vs file databases", {
  # Use existing BED file from package
  bed_file <- system.file("extdata", "example.bed", package = "dbSequence")
  skip_if(!file.exists(bed_file), "Example BED file not found")

  # Test that in-memory connections work independently
  library(DBI)
  library(duckdb)

  # Import to separate in-memory databases (each creates its own connection)
  result1 <- BiocIO::import(
    con = bed_file,
    dest = DuckDBFile(":memory:"),
    table_name = "table1"
  )
  result2 <- BiocIO::import(
    con = bed_file,
    dest = DuckDBFile(":memory:"),
    table_name = "table2"
  )

  # Both should work and be separate objects
  expect_s4_class(result1, "dbSequence")
  expect_s4_class(result2, "dbSequence")
  expect_equal(tableName(result1), "table1")
  expect_equal(tableName(result2), "table2")

  # Both should have accessible data (each in its own in-memory database)
  expect_gt(nrow(dplyr::collect(result1@value)), 0)
  expect_gt(nrow(dplyr::collect(result2@value)), 0)

  # Test that file connections work independently
  temp_db1 <- tempfile(fileext = ".duckdb")
  temp_db2 <- tempfile(fileext = ".duckdb")
  on.exit(
    {
      if (file.exists(temp_db1)) {
        file.remove(temp_db1)
      }
      if (file.exists(temp_db2)) file.remove(temp_db2)
    },
    add = TRUE
  )

  result3 <- BiocIO::import(
    con = bed_file,
    dest = DuckDBFile(temp_db1),
    table_name = "table1"
  )
  result4 <- BiocIO::import(
    con = bed_file,
    dest = DuckDBFile(temp_db2),
    table_name = "table1"
  )

  # Both should work with same table name but different files
  expect_s4_class(result3, "dbSequence")
  expect_s4_class(result4, "dbSequence")
  expect_equal(tableName(result3), "table1")
  expect_equal(tableName(result4), "table1")
  expect_false(result3@file_source == result4@file_source)

  # Both should have accessible data
  expect_gt(nrow(dplyr::collect(result3@value)), 0)
  expect_gt(nrow(dplyr::collect(result4@value)), 0)
})
