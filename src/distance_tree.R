library(data.table)
library(ggtree)

#############
# FUNCTIONS #
#############

ReadMummerReport <- function(x){
    fread(paste('grep "AlignedBases"', x),
          header = FALSE,
          nrows = 1,
          col.names = "AlignedBases",
          drop = c("V1", "V2"))
}

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

###########
# GLOBALS #
###########


# dev
report_files <- list.files("output/mummer",
                           recursive = TRUE,
                           pattern = "out.report",
                           full.names = TRUE)
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
# identity_data[, distance := (100 - as.numeric(aligned_pct)) / 100]

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

nj <- ape::root(ape::nj(distance_matrix), "mh_UNK_trim-decon_k41_diplo1")
negbranch <- which(nj$edge.length < 0)
branchdiff <- abs(nj$edge.length[negbranch])
nj$edge.length[negbranch] <- 0
nj$edge.length[negbranch + 1] <- nj$edge.length[negbranch + 1] + branchdiff

nj2 <- ape::drop.tip(nj, "fopius_arisanus")

# tree. scale is average # of SNPs per base
gt <- ggplot(nj2, aes(x, y)) +
    geom_tree() +
    theme_tree() +
    geom_treescale(offset = -0.2) + # offset moves the scale bar label up and down
    geom_tiplab()
gt

ggsave("test.pdf", gt, width = 10, height = 7.5, units = "in")
