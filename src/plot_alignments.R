#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)
library(ggplot2)
library(grid)

###########
# GLOBALS #
###########

plot_data_file <- snakemake@input[["plot_data"]]
plot_file <- snakemake@output[["plot_file"]]

########
# MAIN #
########

pd <- readRDS(plot_data_file)
gp <- ggplot(pd,
             aes(x = ref_coord / 1e6,
                 y = query_coord / 1e6,
                 colour = `%IDY`)) +
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          legend.position = c(3.5/4, 1.5/4),
          legend.key.size = unit(8, "pt")) +
    facet_grid(query_label ~ ref_label,
               as.table = TRUE,
               switch = "both",
               labeller = label_parsed) +
    xlab("Position in reference (MB)") + 
    ylab("Position in query (MB)") +
    scale_colour_gradientn(colours = YlOrRd,
                           guide = guide_colourbar(title = "Identity (%)")) +
    geom_point(shape = 18, alpha = 0.5)

# hide the empty panels
hidden_panels <- c("panel-1-2",
                   "panel-1-3",
                   "panel-1-4",
                   "panel-2-3",
                   "panel-2-4",
                   "panel-3-4")
g <- ggplotGrob(gp)
for(panel in hidden_panels){
    g$grobs[[which(g$layout$name == panel)]] <- nullGrob()
}

# write to pdf
cairo_pdf(plot_file, width = 10, height = 7.5, bg = "transparent")
grid.newpage()
grid.draw(g)
dev.off()

# write log
sessionInfo()

