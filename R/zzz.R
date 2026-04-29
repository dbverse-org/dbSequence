# zzz.R -- Package initialization and SQL translator setup
# ------------------------------------------------------------------------------
# Package hooks for loading and SQL translation registration

# Suppress R CMD check NOTEs for NSE variables used in dplyr/ggplot2
utils::globalVariables(c(
  # dplyr NSE variables
  "bin_start", "bin_end", "count", "group",
  "start", "end", "y", "y_chr", "y_start", "y_end",
  # ggplot2 aes() variables
  "xmin", "xmax", ".data", ":=",
  # base R
  "setNames"
))

#' @keywords internal
#' @noRd
#' @importFrom dbplyr sql_variant sql_infix
.onLoad <- function(libname, pkgname) {
  # Register DuckDB-specific SQL translations for genomic operations
  # This allows dplyr verbs to generate optimal SQL for range overlaps
  # register_sql_translator(
  #   variant = sql_variant(
  #     scalar = sql_translate_env(
  #       # Custom function for genomic range overlaps
  #       # overlaps(start1, end1, start2, end2) -> boolean
  #       overlaps = sql_infix("OVERLAPS"),
  #
  #       # Distance between ranges
  #       # range_distance(start1, end1, start2, end2) -> integer
  #       range_distance = function(start1, end1, start2, end2) {
  #         sql(paste0("GREATEST(0, GREATEST(", start1, ", ", start2,
  #                    ") - LEAST(", end1, ", ", end2, "))"))
  #       },
  #
  #       # Range width calculation
  #       # range_width(start, end) -> integer
  #       range_width = function(start, end) {
  #         sql(paste0("GREATEST(0, ", end, " - ", start, " + 1)"))
  #       }
  #     )
  #   ),
  #   con_class = "duckdb_connection"
  # )

  invisible()
}
