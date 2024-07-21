### "========================================================================="
### Analyze the spatial patterns of PC metabolism enzymes in HCC usinbg STx
### Note: Step 4 must be run first
###
### Copyright (c) 2024. Bioinformatics Institute, A*STAR, Singapore
### License: MIT License. See LICENSE for the full license.
### Authors: Lit-Hsin Loo (loolh@bii.a-star.edu.sg)
### "========================================================================="
library(Seurat)
library(tidyverse)
library(here)
library(org.Hs.eg.db)
library(fs)
library(DESeq2)
library(clusterProfiler)
library(enrichplot)

source(here("lib", "MSI_lib.R"))
source(here("lib", "Study_lib.R"))

### Load the configuration file
source(here("conf", "study_conf.R"))

### "================================================================="
### Load RNA data
### "================================================================="
### * RNA_vst_mat (normalized data for visualization or ML)
### * dds
### * gene_names
tmp <- LoadRNAData()
dds <- tmp$dds
RNA_vst_mat <- tmp$RNA_vst_mat
all_gene_info <- tmp$all_gene_info
rm(tmp)

GetRNARes <- function(dds, stage_name) {
    dds_res <- results(
        dds,
        contrast = c("condition", stage_name, "Normal"),
        lfcThreshold = RNA_log2FC_thres,
        alpha = RNA_padj_thres
    )

    as_tibble(
        dds_res[, c("log2FoldChange", "pvalue", "padj")],
        rownames = "ENSEMBL"
    )
}

F008 <- Load10X_Spatial(paste0(visium_raw_dir, "/F008_CA/outs"))
F011 <- Load10X_Spatial(paste0(visium_raw_dir, "/F011_CA/outs"))


### Normalize the data
F008 <- SCTransform(F008, assay = "Spatial", verbose = FALSE)
F011 <- SCTransform(F011, assay = "Spatial", verbose = FALSE)

### '--------------------------------------------------------------------------'
### PC metabolism pathways
### * Kennedy pathway, PEMT pathway, Lands cycle
### '--------------------------------------------------------------------------'
### Find all the glycerophospholipid metabolism eynzmes
pc_idx <- grepl("^PLA2G|^LPCAT|^CEPT1$|^CHPT1$|^PEMT$", all_gene_info$SYMBOL)
pc_gene_info <- all_gene_info[pc_idx, ]
pc_dds <- dds[rownames(dds) %in% pc_gene_info$ENSEMBL,]

pc_res <- GetRNARes(pc_dds, "StageII") |>
    left_join(pc_gene_info)


pc_plot <- pc_res |>
    mutate(
        is_sig = if_else(padj < 0.05 & abs(log2FoldChange) > 1, "Y", "N"),
        SYMBOL = fct_relevel(
            SYMBOL,
            "CEPT1", "CHPT1",
            "PEMT",
            "LPCAT1", "LPCAT2", "LPCAT3", "LPCAT4",
            "PLA2G1B", "PLA2G2A", "PLA2G2D", "PLA2G4A", "PLA2G4B", "PLA2G4C",
            "PLA2G5", "PLA2G6", "PLA2G7", "PLA2G12A", "PLA2G12B", "PLA2G15"
        )
    ) |>
    ggplot(aes(x = log2FoldChange, y = SYMBOL, fill = is_sig)) +
    geom_bar(stat = "identity", color = "black", width = 0.8) +
    xlab("log2-FoldChange") +
    ylab("Genes coding key PC enzymes") +
    scale_fill_manual(
        name = "",
        values = c("#377eb8", "#e41a1c"),
        labels = c("Padj>0.05 or |log2FC|<1", "Padj<0.05 and |log2FC|>1")
    ) +
    coord_flip() +
    theme_classic() +
    theme(
        legend.position = "top",
        legend.key.size = unit(0.4, 'cm'),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

SavePlot(
    pc_plot, width = 4, height = 2.4,
    dir_name = "RNA", file_name = "RNA_PC_genes"
)


pc_degs <- pc_res |> arrange(padj) |> filter(padj < 0.05) |> pull(SYMBOL)

## We focus on the following genes
pc_examples <- c(
    "CHPT1",
    "PEMT",
    "LPCAT1", ### Land Cycle
    "PLA2G5"  ### Land Cycle (specific to PC)
)

### Define the colormap
SpatialColors <- colorRampPalette(
    colors = rev(x = RColorBrewer::brewer.pal(n = 11, name = "Spectral"))
)

### View selected enzymes
### F011 = Grade 1
F011_plot <- SpatialFeaturePlot(
        F011,
        features = pc_examples, 
        pt.size.factor = 2,
        image.alpha = 0, #0.8,
        # keep.scale = "all",
        stroke = NA,
        #features = names(gene_sets)
    ) & scale_fill_gradientn(
        limits = c(0, 2),
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


### F008 = Grade 4
F008_plot <- SpatialFeaturePlot(
        F008,
        features = pc_examples, 
        pt.size.factor = 2,
        image.alpha = 0, #0.8,
        # keep.scale = "all",
        stroke = NA,
        #features = names(gene_sets)
    ) & scale_fill_gradientn(
        limits = c(0, 2),
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

SavePlot(
    F011_plot, width = 8, height = 5,
    dir_name = "RNA", file_name = "RNA_PC_STx_F011"
)

SavePlot(
    F008_plot, width = 8, height = 5,
    dir_name = "RNA", file_name = "RNA_PC_STx_F008"
)
