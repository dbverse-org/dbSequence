# dbSequence

`dbSequence` provides a DuckDB-backed interface for genomic interval and
sequence-derived data that are too large to work with comfortably in memory. It
imports common genomics file formats into DuckDB, keeps range operations lazy,
and interoperates with Bioconductor classes such as `GRanges` through
`BiocIO`, `GenomicRanges`, `IRanges`, `Rsamtools`, `rtracklayer`, and
`VariantAnnotation`.

## Installation

Once `dbSequence` is available in Bioconductor, install it with:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("dbSequence")
```

Development versions can be installed from GitHub with:

```r
# install.packages("pak")
pak::pak("dbverse-org/dbsequence-r")
```

## Example

```r
library(dbSequence)
library(GenomicRanges)

bed_file <- system.file("extdata", "example.bed", package = "dbSequence")

# Import BED data into DuckDB; the result is lazy and database-backed
fragments <- read_bed(bed_file)

# Restrict to a region of interest without collecting into memory
promoter <- GRanges("chr1:100-500")
promoter_fragments <- filter_by_overlaps(fragments, promoter)

# Materialize only when Bioconductor-native ranges are needed
as_granges(promoter_fragments)
```

## Features

- Import genomic ranges from BED, GFF/GTF, VCF, and BAM-backed sources
- Keep interval filtering and overlap operations lazy in DuckDB
- Convert filtered results back to `GRanges` only when collection is needed
- Use a `plyranges`-compatible API for common interval workflows

## Documentation

The package includes executable vignettes covering:

- data ingestion with `BiocIO` and `DuckDBFile`
- overlap filtering and lazy `plyranges`-style workflows
