#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(tidyverse)

#############
# FUNCTIONS #
#############

ReadMcoordAndAddCoordinates <- function(coord_file){
    # define columns
    mcoords_cols <- c(
        "S1", "E1",
        "S2", "E2",
        "LEN1", "LEN2",
        "%IDY",
        "LENR", "LENQ",
        "COVR", "COVQ",
        "R", "Q")
    # get assembly names
    result_name <- sub("output/mummer/", "", dirname(coord_file))
    # read data
    coords <- read_tsv(coord_file,
                       col_names = mcoords_cols)
    # sort by reference postion
    sorted_coords <- coords %>% 
        arrange(-LENR, S1) %>% 
        mutate(R = factor(R, levels = unique(R))) %>% 
        arrange(R) %>% 
        mutate(hit_id = 1:n())
    # calculate ref coordinates
    hits_with_ref_coords <- sorted_coords %>% 
        distinct(R, LENR, .keep_all = TRUE) %>% 
        transmute(R, ref_scaf_start_coord = cumsum(LENR) + 1 - LENR) %>% 
        full_join(sorted_coords, by = "R") %>% 
        mutate(S1_coord = ref_scaf_start_coord + S1 - 1,
               E1_coord = ref_scaf_start_coord + E1 - 1)
    # calculate query coordinates          
    hits_with_query_coords <- hits_with_ref_coords %>% 
        arrange(R, S1_coord) %>% 
        mutate(Q = factor(Q, levels = unique(Q))) %>% 
        distinct(Q, LENQ, .keep_all = TRUE) %>% 
        arrange(Q) %>% 
        transmute(Q = as.character(Q),
                  query_scaf_start_coord = cumsum(LENQ) + 1 - LENQ) %>% 
        full_join(hits_with_ref_coords, by = "Q") %>% 
        mutate(S2_coord = query_scaf_start_coord + S2 - 1,
               E2_coord = query_scaf_start_coord + E2 - 1)
    # split the result name
    return(    
        hits_with_query_coords %>% 
            mutate(R = as.character(R),
                   result_name = result_name) %>% 
            separate(result_name,
                     c("ref_assembly", "query_assembly"),
                     "-vs-"))
}

MakeLabels <- function(x) {
    if(x == "fopius_arisanus"){
        return("'Position in'~italic('Fopius arisanus')~'(MB)'")
    }
    my_name_data <- tibble(assembly_name = x) %>% 
        separate(col = assembly_name,
                 into = c("species", "strain", "processing", "kmer", "diplo"),
                 sep = "_") %>% 
        mutate(
            species = if_else(
                species == "ma",
                "M. aethiopoides",
                "M. hyperodae"),
            spec_strain = if_else(
                strain == "UNK",
                paste0("'Position in'~",
                       "italic(",
                       species,
                       ")~'(MB)'"),
                paste0("'Position in'~",
                       "italic('",
                       species,
                       "')~'",
                       strain,
                       "'~'(MB)'")))
    return(my_name_data$spec_strain)
}

MakeMatchTibble <- function(coord_tbl){
    ref_coord <- seq(coord_tbl$S1_coord, coord_tbl$E1_coord)
    query_coord <- seq(coord_tbl$S2_coord, coord_tbl$E2_coord)
    if (length(ref_coord) > length(query_coord)){
        length_diff <- length(ref_coord) - length(query_coord)
        query_coord <- c(query_coord, rep(NA, length_diff))
    } else if (length(ref_coord) < length(query_coord)){
        length_diff <- length(query_coord) - length(ref_coord)
        ref_coord <- c(ref_coord, rep(NA, length_diff))
    }
    tibble(ref_assembly = coord_tbl$ref_assembly[[1]],
           query_assembly = coord_tbl$query_assembly[[1]],
           ref_coord, 
           query_coord,
           `%IDY` = coord_tbl$`%IDY`[[1]])
}


###########
# GLOBALS #
###########

coord_file <- snakemake@input[["coords_file"]]
plot_data_file <- snakemake@output[["plot_data"]]

# dev
# coord_files <- list.files("output/mummer",
#                           recursive = TRUE,
#                           pattern = "out.1coords",
#                           full.names = TRUE)
# coord_file <- coord_files[[1]]

coord_file <- "output/mummer/ma_FR_norm_k71_diplo1-vs-ma_IE_trim-decon_k71_diplo0//out.1coords"

########
# MAIN #
########

# read the mummer file
coord_data <- ReadMcoordAndAddCoordinates(coord_file)

# make into a matrix of matches
plot_data <- filter(coord_data, LEN2 > 250) %>%
    group_by(ref_assembly, query_assembly, hit_id) %>% 
    do(MakeMatchTibble(.)) %>% 
    ungroup(x) %>% 
    group_by(ref_assembly, query_assembly, ref_coord, query_coord) %>% 
    slice(which.max(`%IDY`)) %>% 
    ungroup() %>% 
    mutate(ref_label = MakeLabels(ref_assembly[[1]]),
           query_label = MakeLabels(query_assembly[[1]]))

# write to R binary file
saveRDS(plot_data, plot_data_file)
saveRDS(plot_data, "pd.Rds")

# write log
sessionInfo()

