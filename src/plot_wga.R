#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log, type = "message")
sink(log, type = "output", append = TRUE)

library(circlize)
library(data.table)

fai_file <- snakemake@input[["ref_fai"]]
query_fai_file <- snakemake@input[["query_fai"]]
paf_file <- snakemake@input[["paf"]]
plot_file <- snakemake@output[["plot"]]

# dev
# fai_file <- "wga/tmp/m_aethiopoides.fa.fai"
# query_fai_file <- "wga/020_oriented/m_hyperodae.fa.fai"
# paf_file <- "wga/030_wga/wga.paf"
# plot_file <- "wga.pdf"
# indiv_name <- "BB31"

# need the names for the FAI
fai_names <- c("name", "length", "offset", "line_bases", "line_bytes")

# only get the main chromosomes from the ref assembly
fai <- fread(fai_file, col.names = fai_names)[
    startsWith(name, "PGA_scaffold")][
        order(length, decreasing = TRUE)]
ref_order <- fai[, unique(name)]

# read the query genome fai
query_fai <- fread(query_fai_file,
                   col.names = fai_names)[
                       startsWith(name, "PGA_scaffold")][
                           order(length, decreasing = TRUE)]


# read the paf
pafnames <- c(
    "query",
    "query_length",
    "query_start",
    "query_end",
    "strand",
    "target",
    "target_length",
    "target_start",
    "target_end",
    "residue_matches",
    "alignment_block_length",
    "mapping_quality")
mypaf <- fread(paf_file, select = 1:12, col.names = pafnames, fill = TRUE)

# we only want hits > x kb against the ref chroms
block_size_min <- 1e3
score_min <- 30
filtered_paf <- mypaf[alignment_block_length > block_size_min &
                          mapping_quality > score_min &
                          startsWith(query, "PGA_scaffold")]

# only plot one target per query - this is to tidy up the weird cross-links
n_per_query <- filtered_paf[, .N, by = .(query, target)]
keep_rows <- n_per_query[, .I[which.max(N)], by = query][, V1]
n_per_query[keep_rows]

filtered_paf <- filtered_paf[n_per_query[keep_rows], on = .(query, target)]

# order by ref contig which they have their biggest hit on?
keep <- filtered_paf[, .I[which.max(alignment_block_length)], by = query][, V1]
biggest_hits <- filtered_paf[keep]
biggest_hits[, ref_fact := factor(target, levels = ref_order)]
setorder(biggest_hits, ref_fact)
query_order <- biggest_hits[, unique(query)]

# links
x1 <- filtered_paf[query %in% query_order &
                       target %in% ref_order,
                   .(target, target_start, target_end)]

# subtract the length for query positions to reverse the coordinates
x2 <- filtered_paf[query %in% query_order &
                       target %in% ref_order,
                   .(query,
                     query_length - query_start,
                     query_length - query_end)]


# join the 2 fais
query_fai_subset <- query_fai[name %in% query_order]
joined_fai <- rbind(fai,
                    query_fai_subset,
                    fill = TRUE)

myfai <- joined_fai[, .(name, 0, length)]


# generate chr colours
ncols <- fai[, .N + 2]
col_pal <- viridis::viridis(ncols)
refcol <- col_pal[1]
qcol <- col_pal[length(col_pal)]
linkcols <- col_pal[-length(col_pal)][-1]
fai[, chr_col := linkcols[.I], by = name]
linkcol <- fai[x1, on=c(name="target"), chr_col]

# draw the plot
wo <- grid::convertUnit(grid::unit(398, "pt"), "in", valueOnly = TRUE)
ho <- grid::convertUnit(grid::unit(227, "pt"), "in", valueOnly = TRUE)

cairo_pdf(plot_file,
          family = "Lato",
          width = wo,
          height = ho)
par(xpd = NA)
circos.initializeWithIdeogram(cytoband = myfai,
                              plotType = c("axis"),
                              chromosome.index = c(ref_order, rev(query_order)),
                              axis.labels.cex = 0.3)

text(0, 1.1, "M. hyperodae", cex = 0.75, col = qcol, font = 3)
text(0, -1.1, "M. aethiopoides", cex = 0.75, col = refcol, font = 3)

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    # named chromosomes - don't fit
    # chr = CELL_META$sector.index
    chr = ""
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    circos.rect(xlim[1], 0,
                xlim[2], 1,
                col = ifelse(CELL_META$sector.index %in% query_order,
                             qcol,
                             refcol),
                lwd = 0.3)
    circos.text(mean(xlim), mean(ylim), chr, cex = 0.7, col = "white",
                facing = "inside", niceFacing = TRUE)},
    track.height = 0.1, bg.border = NA)

circos.genomicLink(x1, x2,
                   col = NA,
                   border = linkcol,
                   lwd = 0.3)
dev.off()

sessionInfo()
