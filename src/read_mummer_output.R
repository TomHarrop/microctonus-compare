#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)

#############
# FUNCTIONS #
#############

FillWithNAs <- function(x1, x2, y1, y2) {
    xs <- seq(x1, x2)
    ys <- seq(y1, y2)
    if(length(xs) == length(ys)) {
        return(list(ref_coord = xs,
                    query_coord = ys))
    }else if(length(xs) > length(ys)){
        ld <- length(xs) - length(ys)
        return(list(ref_coord = xs,
                    query_coord = c(ys, rep(NA, ld))))
    } else if(length(xs) > length(ys)) {
        ld <- length(ys) - length(xs)
        return(list(ref_coord = c(xs, rep(NA, ld)),
                    query_coord = ys))
    }
}

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
                                        my_name_data["strain"]))
        return(my_spec_strain)
    }
}

###########
# GLOBALS #
###########

coord_file <- snakemake@input[["coords_file"]]
plot_data_file <- snakemake@output[["plot_data"]]

mcoords_cols <- c(
    "S1", "E1",
    "S2", "E2",
    "LEN1", "LEN2",
    "%IDY",
    "LENR", "LENQ",
    "COVR", "COVQ",
    "R", "Q")

# dev
# coord_files <- list.files("output/mummer",
#                           recursive = TRUE,
#                           pattern = "out.1coords",
#                           full.names = TRUE)
# coord_file <- coord_files[[1]]
# coord_file <- "output/mummer/ma_FR_norm_k71_diplo1-vs-ma_IE_trim-decon_k71_diplo0//out.1coords"

########
# MAIN #
########

# read the mummer file
coord_data <- fread(coord_file,
                    col.names = mcoords_cols)

# get assembly names
result_name <- sub("output/mummer/", "", dirname(coord_file))

# sort by reference postion
setorder(coord_data, -LENR, S1)
coord_data[, R := factor(R, levels = unique(R))]
setorder(coord_data, R)
coord_data[, hit_id := seq(1, .N)]

# calculate ref coordinates
ref_start_coords <- unique(coord_data, by = "R")[
    , .(R, ref_scaf_start_coord = cumsum(LENR) + 1 - LENR)]
coord_data <- merge(coord_data,
                    ref_start_coords,
                    all.x = TRUE,
                    by = "R")
coord_data[, c("S1_coord", "E1_coord") :=
               .(ref_scaf_start_coord + S1 - 1,
                 ref_scaf_start_coord + E1 - 1)]

# calculate query coordinates
setorder(coord_data, R, S1_coord)
coord_data[, Q := factor(Q, levels = unique(Q))]
setorder(coord_data, Q)
query_start_coords <- unique(coord_data, by = "Q")[
    , .(Q, query_scaf_start_coord = cumsum(LENQ) + 1 - LENQ)]
coord_data <- merge(coord_data,
                    query_start_coords,
                    all.x = TRUE,
                    by = "Q")
coord_data[, c("S2_coord", "E2_coord") :=
               .(query_scaf_start_coord + S2 - 1,
                 query_scaf_start_coord + E2 - 1)]
coord_data[, c("ref_assembly", "query_assembly") := 
               tstrsplit(result_name, split = "-vs-")]

# make full coordinates
plot_data <- merge(
    coord_data,
    coord_data[, FillWithNAs(S1_coord, E1_coord, S2_coord, E2_coord),
               by = hit_id],
    all = TRUE,
    by = "hit_id")

# add labels
plot_data[, ref_label := MakeLabels(ref_assembly), by = ref_assembly]
plot_data[, query_label := MakeLabels(query_assembly), by = query_assembly]

# write to R binary file
saveRDS(plot_data[LEN2 > 250], plot_data_file)
# saveRDS(plot_data, "pd.Rds")

# write log
sessionInfo()

