library(data.table)
library(ggplot2)
library(grid)

# plot the alignments
YlOrRd <- RColorBrewer::brewer.pal(6, "YlOrRd")

pd <- readRDS("output/plot_data/mummer_test_data.Rds")

gp <- ggplot(pd,
       aes(x = ref_coord / 1e6,
           y = query_coord / 1e6,
           colour = `%IDY`)) +
    coord_fixed() +
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          legend.position = c(3.5/4, 1.5/4),
          legend.key.size = unit(8, "pt")) +
    facet_grid(query_label ~ ref_label,
               scales = "free",
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
plot_file <- "test.pdf"
pdf(plot_file, width = 10, height = 7.5, bg = "transparent")
grid.newpage()
grid.draw(g)
dev.off()

# write to jpeg
jpeg_file <- "test.jpeg"
jpeg(jpeg_file, width = 10, height = 7.5, units = "in", res = 300, bg = "transparent")
grid.newpage()
grid.draw(g)
dev.off()

grid.newpage()
grid.draw(g)

# you can't see fopius on this version, because the matches are so short, but that's ok -
# fopius is only for genetic distance

p <- ggplot(filter(coord_data, LEN2 > 2000),
            aes(x = S1_coord / 1e6,
                xend = E1_coord / 1e6,
                y = S2_coord / 1e6,
                yend = E2_coord / 1e6,
                colour = `%IDY`)) +
    theme(strip.background = element_blank(),
          strip.placement = "outside") +
    facet_grid(query_assembly ~ ref_assembly,
               as.table = FALSE,
               switch = "both") +
    xlab("Position in REF (MB)") + ylab("Position in QUERY (MB)") +
    scale_colour_gradientn(colours = YlOrRd,
                           guide = guide_colourbar(title = "Identity (%)")) +
    #    coord_fixed() +
    geom_segment(size = 1)

