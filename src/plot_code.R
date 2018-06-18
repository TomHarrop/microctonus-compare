library(data.table)
library(ggplot2)

# plot the alignments
YlOrRd <- RColorBrewer::brewer.pal(6, "YlOrRd")

pd <- readRDS("pd.Rds")[LEN2 > 250]

ggplot(readRDS("pd.Rds")[LEN2 > 250],
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

g <- ggplotGrob(p)
grid.newpage()
#g$grobs[[which(g$layout$name == "panel-2-2")]] <- nullGrob()
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

