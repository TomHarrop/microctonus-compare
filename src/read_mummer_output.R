library(tidyverse)


mcoords_cols <- c(
    "S1", "E1",
    "S2", "E2",
    "LEN1", "LEN2",
    "%IDY",
    "LENR", "LENQ",
    "COVR", "COVQ",
    "R", "Q"
)


mcoords <- read_tsv("output/mummer/ma_FR_norm_k71_diplo1-vs-ma_IE_trim-decon_k71_diplo0/out.1coords",
                    col_names = mcoords_cols)

# sort by reference postion
sorted_mcoords <- mcoords %>% 
    arrange(-LENR, S1) %>% 
    mutate(R = factor(R, levels = unique(R))) %>% 
    arrange(R) %>% 
    mutate(
        hit_id = 1:n())

# calculate ref coordinates
hits_with_ref_coords <- sorted_mcoords %>% 
    distinct(R, LENR, .keep_all = TRUE) %>% 
    transmute(R, ref_scaf_start_coord = cumsum(LENR) + 1 - LENR) %>% 
    full_join(sorted_mcoords, by = "R") %>% 
    mutate(S1_coord = ref_scaf_start_coord + S1 - 1,
           E1_coord = ref_scaf_start_coord + E1 - 1)

# calculate query coordinates          
hits_with_query_coords <- hits_with_ref_coords %>% 
    arrange(R, S1_coord) %>% 
    mutate(Q = factor(Q)) %>% 
    distinct(Q, LENQ, .keep_all = TRUE) %>% 
    transmute(Q = as.character(Q),
              query_scaf_start_coord = cumsum(LENQ) + 1 - LENQ) %>% 
    full_join(hits_with_ref_coords, by = "Q") %>% 
    mutate(S2_coord = query_scaf_start_coord + S2 - 1,
           E2_coord = query_scaf_start_coord + E2 - 1)
    
YlOrRd <- RColorBrewer::brewer.pal(6, "YlOrRd")
ggplot(hits_with_query_coords,
       aes(x = S1_coord / 1e6,
           xend = E1_coord / 1e6,
           y = S2_coord / 1e6,
           yend = E2_coord / 1e6,
           colour = `%IDY`)) +
    xlab("Position in REF (MB)") + ylab("Position in QUERY (MB)") +
    scale_colour_gradientn(colours = YlOrRd,
                           guide = guide_colourbar(title = "Identity (%)")) +
    coord_fixed() +
    geom_segment(size = 1)

