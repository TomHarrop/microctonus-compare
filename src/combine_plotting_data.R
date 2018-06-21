#!/usr/bin/env Rscript

library(tidyverse)
library(extrafont)

###########
# GLOBALS #
###########

# dev
stats_files <- list.files("output/stats", pattern = ".tsv", full.names = TRUE)
busco_files <- unlist(sapply(list.dirs("output/busco", recursive = FALSE),
                             list.files,
                             recursive = FALSE,
                             pattern = "full_table*",
                             full.names = TRUE))

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

# set names
names(stats_files) <- sub(".tsv", "", basename(stats_files))
names(busco_files) <- gsub("^full_table_(.+)\\.tsv$",
                           "\\1",
                           basename(busco_files))

# read files
stats_data <- bind_rows(lapply(stats_files, read_tsv),
                        .id = "assembly_name")
busco_data <- bind_rows(lapply(busco_files, read_tsv, skip = 4),
                        .id = "assembly_name")

# count buscos
n_total_buscos <- length(unique(busco_data$`# Busco id`))
busco_results <- busco_data %>%
    group_by(assembly_name, Status) %>% 
    summarise(busco_percent = length(unique(`# Busco id`)) * 100 /
                  n_total_buscos)
busco_cols <- unique(busco_results$Status)

# merge stats
plot_data <- spread(busco_results, Status, busco_percent) %>% 
    full_join(stats_data, by = "assembly_name") %>% 
    gather(key, value, -assembly_name) %>% 
    filter(assembly_name != "fopius_arisanus") %>% 
    arrange(assembly_name, key) %>% 
    separate(col = assembly_name,
             into = c("species", "strain", "processing", "kmer", "diplo"),
             sep = "_") %>% 
    mutate(species = if_else(species == "ma",
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
                                     strain,
                                     "'")),
           kmer = as.numeric(gsub("[[:alpha:]]+", "", kmer))) %>% 
    select(-strain)

label_vector <- plot_data %>% 
    select(spec_strain) %>% 
    distinct() %>% 
    unlist(c())
names(label_vector) <- label_vector
label_vector <- sapply(label_vector, function(x) parse(text = bquote(.(x))))

# busco plot, move to poster
busco_plot <- ggplot(filter(plot_data, key %in% busco_cols),
       aes(x = spec_strain, y = value, fill = key)) +
    theme_poster +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
    xlab(NULL) + ylab("Percent") +
    scale_x_discrete(labels = label_vector) +
    scale_fill_brewer(palette = "Set1",
                      guide = guide_legend(title = NULL)) +
    geom_col(position = position_dodge())

ggsave("busco_plot.pdf",
       busco_plot,
       width = width.out,
       height = height.out,
       units = "mm",
       device = cairo_pdf)

# stats plot, move to poster
plot_keys <- c("contig_bp", "n_scaffolds",
               "scaf_L50", "scaf_max")

pd2 <- filter(plot_data, key %in% plot_keys) %>% 
    spread(key, value) %>% 
    transmute(
        spec_strain,
        "'Assembly length (MB)'" = contig_bp / 1e6,
        "Scaffolds" = n_scaffolds,
        "italic(L)[50]~'(KB)'" = scaf_L50 / 1e3,
        "'Longest scaffold (KB)'" = scaf_max / 1e3) %>% 
    gather(key, value, -spec_strain)

ggplot(pd2,
       aes(x = spec_strain, y = value, fill = spec_strain)) +
    theme(strip.background = element_blank(),
          strip.placement = "outside") + 
    xlab(NULL) + ylab(NULL) +
    scale_fill_brewer(palette = "Set1",
                      guide = FALSE) +
    facet_wrap(~ key,
               scales = "free_y",
               labeller = label_parsed,
               strip.position = "left") +
    geom_col()


