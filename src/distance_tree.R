#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)

#############
# FUNCTIONS #
#############

ParseTotalLength <- function(x) {
    fread(paste('grep "TotalLength"', x),
          header = FALSE,
          nrows = 1,
          col.names = "TotalLength",
          drop = c("V1", "V2"))
}

ParseTotalSnps <- function(x) {
    fread(paste('grep "TotalSNPs"', x),
          header = FALSE,
          nrows = 1,
          col.names = "TotalSNPs",
          drop = c("V1", "V2"))
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
                                        my_name_data["strain"],
                                        "'"))
        return(my_spec_strain)
    }
}

###########
# GLOBALS #
###########

report_files <- snakemake@input[["report_files"]]
distance_matrix_file <- snakemake@output[["distance_matrix"]]

# dev
# report_files <- list.files("output/mummer",
#                            recursive = TRUE,
#                            pattern = "out.report",
#                            full.names = TRUE)
# distance_matrix_file <- "output/plot_data/distance_matrix.Rds"

########
# MAIN #
########

# get assembly names
names(report_files) <- sub("output/mummer/", "", dirname(report_files))

# stats
total_snps <- rbindlist(lapply(report_files, ParseTotalSnps),
                        idcol = "result_name")
total_length <- rbindlist(lapply(report_files, ParseTotalLength),
                          idcol = "result_name")
mummer_stats <- merge(total_length, total_snps, by = "result_name")
mummer_stats[, distance := TotalSNPs / TotalLength]
identity_data <- copy(mummer_stats)

# read the mummer file
# identity_data <- rbindlist(lapply(report_files, ReadMummerReport),
#                            idcol = "result_name")
# identity_data[, aligned_pct := gsub("^[[:digit:]]+\\(([^%]+).*", "\\1", AlignedBases)]
# identity_data[, aligned_no := gsub("\\(.*$", "", AlignedBases)]
# identity_data[, scaled_bases := scales::rescale(as.numeric(aligned_no)]

# make a distance matrix
identity_data[, c("ref_assembly", "query_assembly") := 
                  tstrsplit(result_name, split = "-vs-")]
reverse_comparisons <- rbind(
    identity_data[, .(query_assembly, ref_assembly, distance)],
    identity_data[, .(query_assembly = ref_assembly,
                      ref_assembly = query_assembly,
                      distance)])
filtered_comparisons <- reverse_comparisons[!(query_assembly == "fopius_arisanus" |
                                                  ref_assembly == "fopius_arisanus")]
distance_dt <- dcast(reverse_comparisons,
                     query_assembly ~ ref_assembly,
                     value.var = "distance",
                     fill = 0)
distance_matrix <- as.matrix(data.frame(
    distance_dt, row.names = "query_assembly"))

# save output
saveRDS(distance_matrix, distance_matrix_file)

# write log
sessionInfo()
