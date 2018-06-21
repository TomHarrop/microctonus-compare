#!/usr/bin/env Rscript

library(extrafont)
library(tidyverse)
library(scales)

###########
# GLOBALS #
###########

# dev
hist_files <- list.files("data/microctonus_assemblies/norm",
                         pattern = "hist.txt",
                         full.names = TRUE)

hist_out_files <- list.files("data/microctonus_assemblies/norm",
                             pattern = "hist_out.txt",
                             full.names = TRUE)

peaks_files <- list.files("data/microctonus_assemblies/norm",
                          pattern = "peaks.txt",
                          full.names = TRUE)

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

# sizes
width.out <- 210.25
height.out <- 194.564

# colours
set1 <- RColorBrewer::brewer.pal(9, "Set1")

########
# MAIN #
########

# set up names
names(hist_files) <- sub("_hist.txt", "", basename(hist_files))
names(hist_out_files) <- sub("_hist_out.txt", "", basename(hist_out_files))
names(peaks_files) <- sub("_peaks.txt", "", basename(peaks_files))


# read histogram data
hist_data <- bind_rows(lapply(hist_files, read_tsv),
                       .id = "assembly_name")
hist_out_data <- bind_rows(lapply(hist_out_files, read_tsv),
                           .id = "assembly_name")
all_hist_data <- bind_rows(list("Raw" = hist_data,
                                "Normalised" = hist_out_data),
                           .id = "type") %>%
    separate(assembly_name, c("species", "strain"), "-") %>%
    mutate(species = if_else(species == "Ma",
                             "M. aethiopoides",
                             "M. hyperodae"),
           spec_strain = if_else(strain == "UNK",
                                 paste0(
                                     "italic('",
                                     species,
                                     "')"),
                                 paste0(
                                     "italic('",
                                     species,
                                     "')~'",
                                     gsub("[[:digit:]]+", "", strain),
                                     "'")))


# read metadata
bbnorm_metadata <- bind_rows(
    lapply(peaks_files, read_tsv, col_names = c("key", "value"), n_max = 13),
    .id = "assembly_name") %>% 
    mutate(key = sub("#", "", key)) %>%
    separate(assembly_name, c("species", "strain"), "-") %>%
    mutate(species = if_else(species == "Ma",
                         "M. aethiopoides",
                         "M. hyperodae"),
       spec_strain = if_else(strain == "UNK",
                             paste0(
                                 "italic('",
                                 species,
                                 "')"),
                             paste0(
                                 "italic('",
                                 species,
                                 "')~'",
                                 gsub("[[:digit:]]+", "", strain),
                                 "'")))
       

# read peaks
main_peaks <- lapply(peaks_files,
       read_tsv,
       comment = "#",
       col_names = c("start", "centre", "stop", "max", "volume")) %>% 
    bind_rows(.id = "assembly_name") %>% 
    group_by(assembly_name) %>% 
    summarise(main_peak = centre[which.max(volume)],
              peak_start = start[which.max(volume)],
              peak_stop = stop[which.max(volume)]) %>%
    separate(assembly_name, c("species", "strain"), "-") %>%
    mutate(species = if_else(species == "Ma",
                             "M. aethiopoides",
                             "M. hyperodae"),
           spec_strain = if_else(strain == "UNK",
                                 paste0(
                                     "italic('",
                                     species,
                                     "')"),
                                 paste0(
                                     "italic('",
                                     species,
                                     "')~'",
                                     gsub("[[:digit:]]+", "", strain),
                                     "'")))


# test plots (move to poster)
Set1 <- RColorBrewer::brewer.pal(9, "Set1")
kmer_plot <- ggplot(filter(all_hist_data, type == "Raw"),
       aes(x = `#Depth`, y = Unique_Kmers, colour = type)) +
    theme_poster + 
    xlab("31-mer depth") + ylab("Number of unique 31-mers") +
    scale_y_continuous(
        trans = "log10",
        labels = trans_format("log10", math_format(10^.x)),
        breaks = trans_breaks("log10", function(x) 10^x)) +
    scale_x_continuous(trans = log_trans(base = 4),
                       breaks = trans_breaks(function(x) log(x, 4),
                                             function(x) 4^x)) +
    scale_colour_brewer(palette = "Set1",
                        guide = FALSE) +
    facet_wrap(~spec_strain, labeller = label_parsed) +
    geom_vline(data = main_peaks,
               mapping = aes(xintercept = main_peak),
               linetype = 2,
               colour = alpha(Set1[2], 0.5)) +
    geom_text(data = main_peaks,
               mapping = aes(x = main_peak,
                             y = 10^8,
                             label = paste0(main_peak, "×")),
               colour = Set1[2], hjust = -0.1) +
    geom_path()



ggsave("kmer_plot.pdf",
       kmer_plot,
       width = width.out,
       height = height.out,
       units = "mm",
       device = cairo_pdf)

