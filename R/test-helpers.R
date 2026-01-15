# test-helpers.R -- Test utilities and example data creation
# ------------------------------------------------------------------------------
# Functions to create test datasets for demonstrating dbSequence functionality

#' @keywords internal
#' @noRd
#' @importFrom duckdb duckdb
#' @importFrom DBI dbConnect dbDisconnect dbWriteTable
#' @importFrom stats rpois runif
.create_test_dataset <- function(
  db_path,
  table_name = "test_cells",
  n_cells = 1000,
  n_features = 2000
) {
  # Create sample cell metadata
  set.seed(42) # For reproducible results

  cell_data <- data.frame(
    cell_id = paste0("Cell_", seq_len(n_cells)),
    n_features = rpois(n_cells, lambda = 800) + 200, # 200-1500 features per cell
    total_counts = rpois(n_cells, lambda = 5000) + 1000, # 1000-10000 counts per cell
    mito_pct = runif(n_cells, min = 0, max = 25), # 0-25% mitochondrial
    cell_type = sample(
      c("T_cell", "B_cell", "Monocyte", "NK_cell"),
      n_cells,
      replace = TRUE
    ),
    treatment = sample(c("Control", "Treated"), n_cells, replace = TRUE)
  )

  # Add some quality filtering scenarios
  # Make some cells fail QC
  poor_quality_idx <- sample(n_cells, size = n_cells * 0.1)
  cell_data$n_features[poor_quality_idx] <- rpois(
    length(poor_quality_idx),
    lambda = 100
  )
  cell_data$mito_pct[poor_quality_idx] <- runif(
    length(poor_quality_idx),
    min = 30,
    max = 50
  )

  # Connect to DuckDB and create table
  con <- .get_duckdb_connection(db_path)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))

  # Write the data
  DBI::dbWriteTable(con, table_name, cell_data, overwrite = TRUE)

  cat("Created test dataset:\n")
  cat("  Database:", db_path, "\n")
  cat("  Table:", table_name, "\n")
  cat("  Cells:", n_cells, "\n")
  cat("  Features simulated:", n_features, "\n")

  # Return dbSequence object
  dbSequence(table_name, db_path)
}

#' @keywords internal
#' @noRd
.import_from_dataframe <- function(
  data,
  db_path,
  table_name,
  overwrite = TRUE
) {
  # Connect to DuckDB
  con <- .get_duckdb_connection(db_path)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))

  # Write the table
  DBI::dbWriteTable(con, table_name, data, overwrite = overwrite)

  # Return dbSequence object
  dbSequence(table_name, db_path)
}

#' @keywords internal
#' @noRd
.cleanup_test_files <- function(...) {
  files <- list(...)
  for (file in files) {
    if (file.exists(file)) {
      unlink(file)
    }
  }
}
