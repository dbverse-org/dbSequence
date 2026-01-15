library(dbSequence)

test_that("Database connection helper works", {
  test_db <- tempfile(fileext = ".duckdb")
  on.exit(if (file.exists(test_db)) file.remove(test_db))

  # Test connection creation
  con <- dbSequence:::.get_duckdb_connection(test_db)
  expect_s4_class(con, "duckdb_connection")

  # Test basic operation
  DBI::dbExecute(con, "CREATE TABLE test (id INTEGER, name TEXT)")
  expect_true(DBI::dbExistsTable(con, "test"))

  DBI::dbDisconnect(con, shutdown = TRUE)
})

test_that("Safe query helper handles errors gracefully", {
  test_db <- tempfile(fileext = ".duckdb")

  # Test with invalid SQL query to force an error
  expect_warning(
    result <- dbSequence:::.safe_query(
      test_db,
      "SELECT * FROM definitely_nonexistent_table_12345"
    ),
    "DuckDB query failed"
  )
  expect_null(result)

  # Clean up
  if (file.exists(test_db)) file.remove(test_db)
})

test_that("Basic utility functions work", {
  # Test basic package loading and access to exported functions
  expect_true(exists("tableName"))
  expect_true(exists("fileSource"))

  # Test with a simple dbSequence object
  db_seq <- dbSequence(
    file_source = "test.duckdb",
    table_name = "test_table"
  )
  expect_equal(tableName(db_seq), "test_table")
  expect_equal(fileSource(db_seq), "test.duckdb")
})
