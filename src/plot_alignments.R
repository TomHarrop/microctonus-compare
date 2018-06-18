#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)
library(ggplot2)

###########
# GLOBALS #
###########

plot_data_file <- snakemake@input[["plot_data"]]
plot_file <- snakemake@output[["plot_file"]]

########
# MAIN #
########

# plot the alignments
YlOrRd <- RColorBrewer::brewer.pal(6, "YlOrRd")
gp <- ggplot(readRDS(plot_data_file),
             aes(x = ref_coord / 1e6,
                 y = query_coord / 1e6,
                 colour = `%IDY`)) +
    theme(strip.background = element_blank(),
          strip.placement = "outside") +
    facet_grid(query_label ~ ref_label,
               as.table = FALSE,
               switch = "both",
               labeller = label_parsed) +
    xlab(NULL) + ylab(NULL) +
    scale_colour_gradientn(colours = YlOrRd,
                           guide = guide_colourbar(title = "Identity (%)")) +
    geom_point(shape = 18, alpha = 0.5)

# write to R binary file
ggsave(plot_file, gp, width = 10, height = 7.5, units = "in")

# write log
sessionInfo()

