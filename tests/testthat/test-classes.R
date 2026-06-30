library(dbSequence)

test_that("S4 classes are defined correctly", {
  # Test DuckDBFile class
  expect_s4_class(new("DuckDBFile", path = "test.duckdb"), "DuckDBFile")
  expect_s4_class(new("DuckDBFile", path = "test.duckdb"), "BiocFile")

  # Test dbSequence class
  db_seq <- new("dbSequence", file_source = "test.duckdb")
  expect_s4_class(db_seq, "dbSequence")
  expect_s4_class(db_seq, "dbData")
  expect_equal(fileSource(db_seq), "test.duckdb")
})

test_that("Constructor functions work", {
  # Test DuckDBFile constructor
  db_file <- DuckDBFile("test.duckdb")
  expect_s4_class(db_file, "DuckDBFile")
  expect_equal(dbSequence:::.duckdb_file_path(db_file), "test.duckdb")

  # Test dbSequence constructor
  db_seq <- dbSequence("test_table", "test.duckdb")
  expect_s4_class(db_seq, "dbSequence")
  expect_equal(tableName(db_seq), "test_table")
  expect_equal(fileSource(db_seq), "test.duckdb")
})

test_that("Accessor methods work", {
  db_seq <- dbSequence("test_table", "test.duckdb")

  expect_equal(tableName(db_seq), "test_table")
  expect_equal(fileSource(db_seq), "test.duckdb")
})

test_that("Show methods work without errors", {
  db_file <- DuckDBFile("test.duckdb")
  db_seq <- dbSequence("test_table", "test.duckdb")

  # These should not throw errors
  expect_message(show(db_file), "DuckDBFile")
  expect_message(show(db_seq), "dbSequence")
})
