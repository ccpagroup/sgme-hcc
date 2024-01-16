### "========================================================================="
### Plot the localization of the COMET's PATH genes using STx
### Note: Step 4 must be run first
###
### Copyright (c) 2021-2024. Bioinformatics Institute, A*STAR, Singapore
### License: MIT License. See LICENSE for the full license.
### Authors: Lit-Hsin Loo (loolh@bii.a-star.edu.sg)
### "========================================================================="
library(Seurat)
library(tidyverse)
library(here)
library(org.Hs.eg.db)
library(fs)

### Load the configuration file
source(here("conf", "study_conf.R"))

#### Load the gene sets
CEG_list_vip <- readRDS(here("results", "RNA", "CEG_list_vip.rds"))

gene_sets <- list()
for (name in names(CEG_list_vip)) {
    ens_ids <- CEG_list_vip[[name]]$CEG_ensembls
    annots <- dplyr::select(
        org.Hs.eg.db, keys=ens_ids, 
        columns="SYMBOL", keytype="ENSEMBL"
    )
    gene_sets[name] <- unique(list(annots$SYMBOL))
}

names(gene_sets) <- c(
    "ME-normal", "ME-low-grade", "ME-high-grade",
    "ME-necrotic", "ME-fibrotic", "ME-steatotic"
)
sapply(gene_sets, length)

### Function to save visium output
ExportVisium <- function(suerat_obj, img_filename) {

    ## suerat_obj <- F008
    ## suerat_obj <- F011

    ### Get the count data genes x spots
    cnts <- GetAssayData(object = suerat_obj, layer = "count")

    ### Convert to matrix for faster calculation
    cnts <- as.matrix(cnts)

    ### Remove genes that have very low counts
    gene_cnts <- rowSums(cnts)
    min_cnt <- 100
    cnts <- cnts[gene_cnts >= min_cnt, ]

    # mean_cnts <- apply(cnts + 1, 2, median)
    cnts <- cnts + 1.0
    spot_mean_cnts <- colMeans(cnts, na.rm = TRUE)
    # spot_mean_cnts <- colSums(cnts, na.rm=TRUE)

    coeff <- matrix(rep(spot_mean_cnts, each = nrow(cnts)), nrow = nrow(cnts))

    ### Normal all genes with respect to the mean counts per spot
    cnts_norm <- cnts / coeff

    gs_exp_mat <- NULL

    for (gs_name in names(gene_sets)) {
        message(paste0("Processing ", gs_name, "..."))
        genes_found <- rownames(cnts) %in% gene_sets[[gs_name]]

        message(
            paste0(
                " * Found ",
                sum(genes_found), "/",
                length(gene_sets[[gs_name]]), " genes"
            )
        )

        ### Get the counts for the gene sets
        gs_cnts_norm <- cnts_norm[genes_found, ]

        ### Normal all genes across all the counts
        for (row_idx in 1:nrow(gs_cnts_norm)) {
            gs_cnts_norm[row_idx, ] <-
                gs_cnts_norm[row_idx, ] * 100 /
                    sum(gs_cnts_norm[row_idx, ], na.rm = TRUE)
        }

        gs_exp_mat <- rbind(
            gs_exp_mat,
            #colMeans(gs_cnts_norm, na.rm = TRUE)
            #apply(gs_cnts_norm, 2, quantile, 0.5)
            apply(gs_cnts_norm, 2, quantile, 0.75)
        )
    }

    # gs_exp_mat[!is.finite(gs_exp_mat)] <- 0

    rownames(gs_exp_mat) <- names(gene_sets)

    suerat_obj <- AddMetaData(suerat_obj, t(gs_exp_mat))

    p <- SpatialFeaturePlot(
        object = suerat_obj,
        pt.size.factor = 2,
        image.alpha = 0, #0.8,
        # keep.scale = "all",
        stroke = NA,
        features = names(gene_sets)
    ) & scale_fill_gradientn(
        colours = SpatialColors(n = 100),
        guide = guide_colorbar(
            frame.colour = "black", ticks.colour = "black"
        )
    ) & theme(
        panel.border = element_rect(
            colour = "black", fill = NA, linewidth = 0.75
        ),
        legend.position = "right"
    )

    fs::dir_create(here("figures", "RNA"))
    
    ggsave(
        width = 11, height = 6, units = "in",
        filename = here("figures", "RNA", paste0(img_filename, ".pdf")),
        plot = p, dpi = 600
    )
    ggsave(
        width = 11, height = 6, units = "in",
        filename = here("figures", "RNA", paste0(img_filename, ".png")),
        plot = p, dpi = 600
    )

    p <- SpatialFeaturePlot(
        object = suerat_obj,
        pt.size.factor = 2,
        image.alpha = 0,
        # keep.scale = "all",
        stroke = NA,
        features = names(gene_sets)
    ) & scale_fill_gradientn(
        colours = SpatialColors(n = 100),
        guide = guide_colorbar(
            frame.colour = "black", ticks.colour = "black"
        )
    ) & theme(
        panel.border = element_rect(
            colour = "black", fill = NA, linewidth = 0.75
        ),
        legend.position = "right"
    )

    ggsave(
        width = 11, height = 6, units = "in",
        filename = here("figures", "RNA", paste0(img_filename, "_clear.pdf")),
        plot = p, dpi = 600
    )
    ggsave(
        width = 11, height = 6, units = "in",
        filename = here("figures", "RNA", paste0(img_filename, "_clear.png")),
        plot = p, dpi = 600
    )

}

### Load the Space ranger output
F008 <- Load10X_Spatial(paste0(visium_raw_dir, "/F008_CA/outs"))
F011 <- Load10X_Spatial(paste0(visium_raw_dir, "/F011_CA/outs"))

### Normalize the data
### (After this, counts are corrected)
F008 <- SCTransform(F008, assay = "Spatial", verbose = FALSE)
F011 <- SCTransform(F011, assay = "Spatial", verbose = FALSE)

### Define the colormap
SpatialColors <- colorRampPalette(
    colors = rev(x = RColorBrewer::brewer.pal(n = 11, name = "Spectral"))
)

### Generate the figures (Fig. 4l)
ExportVisium(suerat_obj = F008, img_filename = "F008_visium")
ExportVisium(suerat_obj = F011, img_filename = "F011_visium")