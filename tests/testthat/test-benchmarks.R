library(testthat)
library(dbSequence)
library(duckdb)
library(DBI)
library(dplyr)

skip_if_not_installed("microbenchmark")
library(microbenchmark)

test_that("BAM file benchmarks - dbSequence vs samtools", {
  # Use the test BAM file from inst/extdata
  bam_file <- system.file(
    "extdata",
    "CytAssist_FFPE_Human_Colon_Post_Xenium_Rep1_possorted_genome_bam.bam",
    package = "dbSequence"
  )

  skip_if_not(file.exists(bam_file), "No BAM test file found")
  skip_if(file.size(bam_file) == 0, "BAM file is empty")
  skip_if(Sys.which("samtools") == "", "samtools not available")

  # Check if BAM index exists, if not create it
  bai_file <- paste0(bam_file, ".bai")
  if (!file.exists(bai_file)) {
    message("Creating BAM index...")
    system2("samtools", c("index", bam_file))
    skip_if_not(file.exists(bai_file), "Failed to create BAM index")
  }

  # Test approach 1: Use dbSequence package's BAM import with region filtering
  # This should leverage exonr's capabilities under the hood
  test_db_path <- tempfile(fileext = ".duckdb")
  on.exit(if (file.exists(test_db_path)) file.remove(test_db_path), add = TRUE)

  # Import BAM using dbSequence - this should use exonr if available
  message("Testing dbSequence BAM import...")
  bam_result <- import(
    BamFile(bam_file),
    dest = DuckDBFile(test_db_path),
    table_name = "bam_data"
  )

  # Get connection to query the imported data
  con <- bam_result@value$src$con

  # Verify the import worked
  total_reads <- DBI::dbGetQuery(con, "SELECT COUNT(*) as n FROM bam_data")$n
  message(sprintf("Successfully imported %d reads", total_reads))
  expect_gt(total_reads, 0)

  # Test 1: Single region query - comparing dbSequence/DuckDB vs samtools
  message("\n=== Test 1: Single Region Query ===")

  # First, get available chromosomes to use for testing
  chrom_query <- "SELECT DISTINCT reference as chrom, COUNT(*) as count FROM bam_data WHERE reference IS NOT NULL GROUP BY reference ORDER BY count DESC LIMIT 3"
  available_chroms <- DBI::dbGetQuery(con, chrom_query)
  skip_if(nrow(available_chroms) == 0, "No chromosomes found in BAM data")

  test_chrom <- available_chroms$chrom[1]
  message(sprintf("Using chromosome %s for testing", test_chrom))

  # DuckDB/dbSequence query for region
  region_start <- 1000000
  region_end <- 2000000

  duckdb_time <- system.time({
    duckdb_result <- DBI::dbGetQuery(
      con,
      sprintf(
        "SELECT COUNT(*) as count FROM bam_data WHERE reference = '%s' AND start >= %d AND start <= %d",
        test_chrom,
        region_start,
        region_end
      )
    )
  })

  # samtools equivalent query
  samtools_region <- sprintf("%s:%d-%d", test_chrom, region_start, region_end)
  samtools_time <- system.time({
    samtools_result <- as.numeric(system2(
      "samtools",
      c("view", "-c", bam_file, samtools_region),
      stdout = TRUE
    ))
  })

  message(sprintf(
    "DuckDB result: %d reads in %.3f seconds",
    duckdb_result$count,
    duckdb_time[["elapsed"]]
  ))
  message(sprintf(
    "samtools result: %d reads in %.3f seconds",
    samtools_result,
    samtools_time[["elapsed"]]
  ))

  # Test 2: Multiple regions - showing DuckDB's strength
  message("\n=== Test 2: Multiple Region Query ===")

  # DuckDB: Single query for multiple regions
  multi_region_query <- sprintf(
    "
    SELECT COUNT(*) as total_count FROM bam_data r
    WHERE (r.reference = '%s' AND r.start >= %d AND r.start <= %d)
       OR (r.reference = '%s' AND r.start >= %d AND r.start <= %d)
       OR (r.reference = '%s' AND r.start >= %d AND r.start <= %d)",
    test_chrom,
    1000000,
    2000000,
    test_chrom,
    3000000,
    4000000,
    test_chrom,
    5000000,
    6000000
  )

  duckdb_multi_time <- system.time({
    duckdb_multi_result <- DBI::dbGetQuery(con, multi_region_query)
  })

  # samtools: requires multiple separate calls
  samtools_multi_time <- system.time({
    region1_count <- as.numeric(system2(
      "samtools",
      c(
        "view",
        "-c",
        bam_file,
        sprintf("%s:%d-%d", test_chrom, 1000000, 2000000)
      ),
      stdout = TRUE
    ))
    region2_count <- as.numeric(system2(
      "samtools",
      c(
        "view",
        "-c",
        bam_file,
        sprintf("%s:%d-%d", test_chrom, 3000000, 4000000)
      ),
      stdout = TRUE
    ))
    region3_count <- as.numeric(system2(
      "samtools",
      c(
        "view",
        "-c",
        bam_file,
        sprintf("%s:%d-%d", test_chrom, 5000000, 6000000)
      ),
      stdout = TRUE
    ))
    samtools_multi_result <- region1_count + region2_count + region3_count
  })

  message(sprintf(
    "DuckDB multi-region: %d reads in %.3f seconds",
    duckdb_multi_result$total_count,
    duckdb_multi_time[["elapsed"]]
  ))
  message(sprintf(
    "samtools multi-region: %d reads in %.3f seconds (3 separate calls)",
    samtools_multi_result,
    samtools_multi_time[["elapsed"]]
  ))

  # Test 3: Complex analytical query - DuckDB's true strength
  message("\n=== Test 3: Complex Analytical Query ===")

  analytical_query <- "
    SELECT 
      reference as chromosome,
      FLOOR(start / 1000000) * 1000000 as megabase_bin,
      COUNT(*) as read_count,
      AVG(mapping_quality::INTEGER) as avg_mapq
    FROM bam_data 
    WHERE reference IS NOT NULL 
      AND mapping_quality IS NOT NULL
    GROUP BY reference, FLOOR(start / 1000000)
    ORDER BY read_count DESC
    LIMIT 10"

  analytical_time <- system.time({
    analytical_result <- DBI::dbGetQuery(con, analytical_query)
  })

  message(sprintf(
    "DuckDB analytical query: %d result rows in %.3f seconds",
    nrow(analytical_result),
    analytical_time[["elapsed"]]
  ))
  message(
    "(This type of analytical query would be very difficult with samtools alone)"
  )

  # Summary
  message("\n=== BENCHMARK SUMMARY ===")
  message(sprintf(
    "Single region - DuckDB: %.3fs, samtools: %.3fs",
    duckdb_time[["elapsed"]],
    samtools_time[["elapsed"]]
  ))
  message(sprintf(
    "Multi region - DuckDB: %.3fs, samtools: %.3fs",
    duckdb_multi_time[["elapsed"]],
    samtools_multi_time[["elapsed"]]
  ))
  message(sprintf(
    "Analytical query - DuckDB: %.3fs (samtools: not feasible)",
    analytical_time[["elapsed"]]
  ))

  expect_true(TRUE)
})
