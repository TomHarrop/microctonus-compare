#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)
library(ggplot2)
library(extrafont)
library(grid)

###########
# GLOBALS #
###########

plot_data_file <- snakemake@input[["plot_data"]]
jpeg_file <- snakemake@output[["jpeg_file"]]
pdf_file <- snakemake@output[["pdf_file"]]

# dev
# plot_data_file <- "output/plot_data/mummer_test_data.Rds"
# jpeg_file <- "alignment_subset.jpg"
# pdf_file <- "alignment_subset.pdf"

########
# MAIN #
########

# setup ggplot theme
theme_poster <- ggplot2::theme_grey(base_size = 16, 
                                       base_family = "Lato") +
    ggplot2::theme(plot.background =
                       ggplot2::element_rect(fill = "transparent",
                                             colour = NA),
                   legend.background =
                       ggplot2::element_rect(fill = "transparent",
                                             colour = NA))

# setup colours
set1 <- RColorBrewer::brewer.pal(9, "Set1")
heatscale <- RColorBrewer::brewer.pal(6, "YlOrRd")

# read data
pd <- readRDS(plot_data_file)

# make plot
gp <- ggplot(pd,
             aes(x = ref_coord / 1e6,
                 y = query_coord / 1e6,
                 colour = `%IDY`)) +
    theme_poster +
    theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = c(3.5/4, 1.5/4)) +
    facet_grid(query_label ~ ref_label,
               scales = "free",
               as.table = TRUE,
               switch = "both",
               labeller = label_parsed) +
    xlab("Position in reference (MB)") + 
    ylab("Position in query (MB)") +
    scale_colour_gradientn(colours = heatscale,
                           guide = guide_colourbar(title = "Identity (%)")) +
    geom_point(shape = 18, alpha = 0.5)

# hide the empty panels
hidden_panels <- c(
    "panel-1-2",
    "panel-1-3",
    "panel-1-4",
    "panel-2-3",
    "panel-2-4",
    "panel-3-4"
)
g <- ggplotGrob(gp)
g$widths[[4]] <- unit(0.75, "cm")
for(panel in hidden_panels){
    g$grobs[[which(g$layout$name == panel)]] <- nullGrob()
}

# write to jpeg since PDF is too big
jpeg(jpeg_file,
     width = 208,
     height = 205,
     units = "mm",
     res = 600,
     type = "cairo")
grid.newpage()
grid.draw(g)
dev.off()

# try pdf anyway
cairo_pdf(pdf_file,
          width = grid::convertUnit(grid::unit(208, "mm"),
                                    "in",
                                    valueOnly = TRUE),
          height = grid::convertUnit(grid::unit(205, "mm"),
                                     "in",
                                     valueOnly = TRUE),
          family = "Lato",
          bg = "transparent")
grid.newpage()
grid.draw(g)
dev.off()

# write log
sessionInfo()

