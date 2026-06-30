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
