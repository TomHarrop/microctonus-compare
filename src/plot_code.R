library(data.table)
library(ggplot2)
library(ggtree)
library(grid)
library(extrafont)

#############
# FUNCTIONS #
#############


MakeLabels <- function(x) {
    if(x == "fopius_arisanus"){
        return("italic('Fopius arisanus')")
    } else {
        my_name_data <- unlist(strsplit(x, "_"))
        names(my_name_data) <- c("species", "strain", "processing", "kmer", "diplo")
        my_species <- ifelse(my_name_data["species"] == "ma",
                             "M. aethiopoides",
                             "M. hyperodae")
        my_spec_strain <- ifelse(my_name_data["strain"] == "UNK",
                                 paste0("italic('",
                                        my_species,
                                        "')"),
                                 paste0("italic('",
                                        my_species,
                                        "')~'",
                                        my_name_data["strain"],
                                        "'"))
        return(my_spec_strain)
    }
}

#########
# SETUP #
#########

# ggplot theme
theme_poster <- ggplot2::theme_grey(base_size = 18,
                                    base_family = "Lato") +
    ggplot2::theme(plot.background =
                       ggplot2::element_rect(fill = "transparent",
                                             colour = NA),
                   legend.background =
                       ggplot2::element_rect(fill = "transparent",
                                             colour = NA))

# colours
set1 <- RColorBrewer::brewer.pal(9, "Set1")
heatscale <- RColorBrewer::brewer.pal(6, "YlOrRd")

##############
# alignments #
##############

# (have to do on BCC because of file size)
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

# write to jpeg
jpeg_file <- "test.jpeg"
jpeg(jpeg_file, width = 10, height = 7.5, units = "in", res = 300, bg = "transparent")
grid.newpage()
grid.draw(g)
dev.off()


#################
# distance tree #
#################

distance_matrix <- readRDS("output/plot_data/distance_matrix.Rds")

nj <- ape::root(ape::nj(distance_matrix), "mh_UNK_trim-decon_k41_diplo1")
negbranch <- which(nj$edge.length < 0)
branchdiff <- abs(nj$edge.length[negbranch])
nj$edge.length[negbranch] <- 0
nj$edge.length[negbranch + 1] <- nj$edge.length[negbranch + 1] + branchdiff

nj2 <- ape::drop.tip(nj, "fopius_arisanus")

# tree annotation
dd <- data.table(tip.label = nj2$tip.label)
dd[, name := MakeLabels(tip.label), by = tip.label]
dd[, c("species", "strain", "processing", "k", "mer_mode") :=
       tstrsplit(tip.label, "_")]
dd[species == "ma" & strain %in% c("FR", "MA"), sexual := "Sexual"]
dd[species == "ma" & strain == "IE", sexual := "Asexual"]
dd[species == "mh", sexual := "Asexual"]
dd[species == "ma" & sexual == "Sexual", social := "Solitary"]
dd[species == "ma" & sexual == "Asexual", social := "Gregarious"]
dd[species == "mh", social := "Solitary"]
dd[strain == "MA", img_url := "circ_img/lucerne.jpg"]
dd[species == "mh", 
   img_url := "circ_img/ASW.jpg"]
dd[is.na(img_url), img_url := "circ_img/CRW.jpg"]


# tree. scale is average # of SNPs per base
gt <- ggtree(nj2, ladderize = FALSE, size = 1) +
    xlim(0, 0.055) +
    theme_poster +
    theme(legend.position = "right",
          panel.background = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(hjust = 1),
          axis.line.x = element_line(size = 0.5)) +
    scale_color_brewer(palette = "Set1",
                       guide = guide_legend(title = NULL)) +
    scale_shape(guide = guide_legend(title = NULL)) +
    xlab(expression("SNPs"~"per"~"aligned"~"base"))
gt2 <- gt  %<+%
    dd +
    geom_tiplab(aes(label = name,
                    colour = sexual),
                parse = TRUE,
                align = FALSE,
                offset = 0.0075,
                hjust = 0,
                family = "Lato",
                size = 14)+
    geom_tiplab(aes(label = name,
                    colour = sexual,
                    image = img_url),
                offset = 0.0005,
                hjust = 0,
                align = FALSE,
                geom = "image",
                size = 0.1) +
    geom_tippoint(aes(shape = social), size = 3)

# set aspect ratio
gr <- ggplotGrob(gt2)
for (k in 1:length(gr$grobs[[6]]$children[[6]]$children)) {
    gr$grobs[[6]]$children[[6]]$children[[k]]$height <- unit(0.1, "native")
    gr$grobs[[6]]$children[[6]]$children[[k]]$width <- unit(0.1, "native")
}

 pdf(file = "test.pdf",
     width = convertUnit(unit(216, "mm"), "in", valueOnly = TRUE),
     height = convertUnit(unit(216, "mm"), "in", valueOnly = TRUE))
grid.newpage()
grid.draw(gr)
dev.off()

