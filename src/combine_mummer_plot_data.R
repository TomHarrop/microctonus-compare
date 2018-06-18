#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)

###########
# GLOBALS #
###########

plot_data_files <- snakemake@input[["plot_data"]]
plot_data_file <- snakemake@output[["plot_data"]]

########
# MAIN #
########

plot_data <- lapply(plot_data_files, readRDS)
all_plot_data <- rbindlist(plot_data)

# write to R binary file
saveRDS(all_plot_data, plot_data_file)

# write log
sessionInfo()

