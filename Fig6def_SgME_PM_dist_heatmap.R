### "========================================================================="
### Plot SgME heatmaps and other statistics
###
### Copyright (c) 2021-2024. Bioinformatics Institute, A*STAR, Singapore
### License: MIT License. See LICENSE for the full license.
### Authors: Lit-Hsin Loo (loolh@bii.a-star.edu.sg)
### "========================================================================="
library(here)
source(here("lib", "MSI_lib.R"))
source(here("lib", "Study_lib.R"))
library(ropls)

library(ComplexHeatmap)
library(dendsort)
library(ggupset)

### Load the classifiers =====================================================
plsda_multi <- LoadSgMEClassifier("plsda", "plsda_multi")
plsda_steatotic <- LoadSgMEClassifier("plsda", "plsda_steatotic")
plsda_fibrotic <- LoadSgMEClassifier("plsda", "plsda_fibrotic")

### Load the SgME maps 
multi_pred_labels <- plsda_multi$pred_labels
steatotic_pred_labels <- plsda_steatotic$pred_labels
fibrotic_pred_labels <- plsda_fibrotic$pred_labels

### Fig. 6e =====================================================
### Load the raw data as matrix
msi_tissue_mat <- LoadAllRawMSIMat(
    slide_ids = test_slide_ids,
    is_normalized = is_normalized_study,
    max_na_fraction = 0.99,
    ### Since our purpose is to study the abundance levels and distributions,
    ### we will fill the undetected peaks
    is_undetected_min = TRUE
)

### Rename the rows
peak_ids <- sub("^(.+):.+", "\\1", rownames(msi_tissue_mat))
rownames(msi_tissue_mat) <- peak_ids

### Find the average level at the normal regions
norm_mean <- rowMeans(
    msi_tissue_mat[, multi_pred_labels == "Normal"],
    na.rm = TRUE
)
g2_mean <- rowMeans(
    msi_tissue_mat[, multi_pred_labels == "G2"],
    na.rm = TRUE
)
g3_mean <- rowMeans(
    msi_tissue_mat[, multi_pred_labels == "G3"],
    na.rm = TRUE
)
necrotic_mean <- rowMeans(
    msi_tissue_mat[, multi_pred_labels == "Necrotic"],
    na.rm = TRUE
)
steatotic_mean <- rowMeans(
    msi_tissue_mat[, steatotic_pred_labels == "Steatotic"],
    na.rm = TRUE
)
fibrotic_mean <- rowMeans(
    msi_tissue_mat[, fibrotic_pred_labels == "Fibrotic"],
    na.rm = TRUE
)

msi_mat <- cbind(
    "ME-steatotic" = steatotic_mean - norm_mean,
    "ME-low-grade" = g2_mean - norm_mean,
    "ME-fibrotic" = fibrotic_mean - norm_mean,
    "ME-high-grade" = g3_mean - norm_mean,
    "ME-necrotic" = necrotic_mean - norm_mean
)

### Remove those with norm_mean is zero
msi_mat <- msi_mat[norm_mean != 0, ]

### Normalize
msi_mat <- t(apply(msi_mat, 1, function(x) x / sum(abs(x))))

### 
load(here("data", "lcms", "lcms_peaks.Rdata"))
lcms_hit_peaks %>%
    filter(peak_id == "P35") %>%
    pull(MSMS_annotation)

### Sample marking
show_names <- c(
    "P135"="P135 (ME-transform PLS-DA)", #"SAPC", ### Transformed 
    "P120"="P120 (ME-low-grade PLS-DA)", # "Digitoxin" ### Transformed/G2  - PC(36:4)
    "N71"="N71 (ME-high-grade PLS-DA)", # "PI(19:0/18:3)", ### G3 - N71
    "N14"="N14 (ME-fibrotic PLS-DA)", # "Cholesterol sulfate",  ### Fibrotic - N14
    "P73"="P73 (ME-steatotic PLS-DA)", # "DG(16:0/0:0/18:1n9)",  ### Steatotic - P73
    "P58"="P58 (ME-normal PLS-DA)", # PC(18:1/0:0)[U], "Normal"
    "P29"="P29 (ME-normal PLS-DA)", # "CHEBI:34173", ### P29 =Normal - Laudanosine
    "P37"="P37 (ME-fibrotic/necrotic PLS-DA)", # "Oleoylcarnitine",  ### Fibrotic - P37
    #"P60"="P60 (ME-low-grade PLS-DA)", # "LysoPC(18:0/0:0)", ### G2 - P60    
    #"P38"="P38 (ME-necrotic PLS-DA)",  # "13-OH-alpha-Tocopherol" ### P38 - Necrotic
    
    P70 = "P70", # "DG(16:0/0:0/18:3)", ### ME-low-grade+/ME-necrotic- PMs
    N54 = "N54", # d-Tocotrienol", ### ME_necrotic+ PMs
    P108 = "P108" # "PC(15:0/18:0)" ## ME-steatotic+ PMs
)

sample_marking <- PrepareSampleMarking(msi_mat, show_names, is_peak_id = TRUE)

sample_annotation_list <- PrepareMetaboliteAnnotation(
    msi_mat, metabolite_classes_colors
)

row_dend <- hclust(
    dist(msi_mat),
    method = "ward.D2"
)

### Enable this to determine optimum cluster cut, then use the custom ordering
row_splits <- factor(
    cutree(row_dend, 6)
)

row_splits <- factor(
    cutree(row_dend, 6),
    levels = c(3, 5, 1, 6, 2, 4),
    labels = c(
        "ME-\nsteatotic+",           # 3
        "ME-low-grade+/\nME-necrotic-",    # 5
        "ME-trans+/\nME-fibrotic+", # 1
        "ME-high-grade+/\nME-necrotic+",  # 6
        "ME-\nnecrotic+",          # 2
        "ME-\nnormal+"            # 4
    )
)

### Save the results
PM_groups <- row_splits
levels(PM_groups) <- gsub("\n", "", levels(PM_groups))
saveRDS(PM_groups, here("results", "PM_groups.rds"))

### Get the labels
clust_sizes <- table(row_splits)
size_labels <- paste0("\n(", clust_sizes, ", ", signif(clust_sizes / sum(clust_sizes) * 100, 2), "%)")
levels(row_splits) <- paste0(levels(row_splits), size_labels)

msi_heatmap <- ComplexHeatmap::Heatmap(
    msi_mat,
    cluster_row_slices = FALSE,
    left_annotation = sample_annotation_list[["metabolite_annotation"]],
    right_annotation = sample_marking,
    cluster_columns = FALSE,
    #cluster_rows = row_dend,
    row_split = row_splits,
    clustering_method_rows = "ward.D2",
    border = TRUE,
    rect_gp = gpar(lty = 0),
    #row_title = paste0(
    #    "Highly-abundant putative metabolites (", nrow(msi_mat), ")"
    #),
    heatmap_legend_param = list(
        title = "Normalized fold-change\nto ME-normal region",
        #title_position = "leftcenter-rot",
        #legend_height = unit(4, "cm"),
        title_position = "topcenter",
        direction = "horizontal",
        legend_width = unit(4, "cm"),
        border = "black"
    ),
    show_row_names = FALSE,
    height = unit(8, "inch"),
    width = unit(1, "inch"),
    use_raster = TRUE,
    raster_quality = 10
)

fs::dir_create(here("figures", "msi"))

width <- 8
height <- 12
for (img_type in c(".png", ".pdf")) {

    file_name <- here(
        "figures", "msi",
        paste0("msi_peak_dist_heatmap", img_type)
    )

    if (img_type == ".png") {
        png(
            file_name,
            width = width, height = height, units = "in", res = 300
        )
    } else {
        pdf(
            file_name,
            width = width, height = height
        )
    }

    ht <- draw(
        msi_heatmap,
        annotation_legend_list = list(
            sample_annotation_list[["metabolite_classes_legend"]]
            # roi_annot_legend,
            # case_annot_legend
        ),
        # merge_legend = TRUE,
        legend_gap = unit(0.3, "inch"),
        heatmap_legend_side = "right",
        annotation_legend_side = "bottom"
    )

    dev.off()
}

### "=========================================================================="
### Fig. 6d. Plot individual PMs
### "=========================================================================="
pm_list <- c("P73", "P120", "P135", "N14", "N71", "N54", "P29")
pm_labels <- lcms_hit_peaks %>%
    filter(peak_id %in% pm_list) %>%
    pull(MSMS_annotation, name = peak_id)
pm_labels <- pm_labels[pm_list]
pm_labels <- paste0(names(pm_labels), "\n", pm_labels)

plot_data <- as_tibble(
        msi_mat[pm_list, ],
        rownames = "peak_id"
    ) |> pivot_longer(
        cols = !peak_id, 
        names_to = "MER", 
        values_to = "abd_norm"
    ) |> mutate(
        peak_id = factor(
            peak_id,
            levels = pm_list,
            labels = pm_labels
        ),
        MER = factor(
            MER,
            levels = c(
                "ME-steatotic",
                "ME-low-grade",
                "ME-fibrotic",
                "ME-high-grade",
                "ME-necrotic"
            )
        )
    )

MER_annot_colors <- c(
    "ME-steatotic" = roi_annot_colors[["Steatotic"]],
    "ME-low-grade" = roi_annot_colors[["G2"]],
    "ME-fibrotic"  = roi_annot_colors[["Fibrotic"]],
    "ME-high-grade"= roi_annot_colors[["G3"]],
    "ME-necrotic"  = roi_annot_colors[["Necrotic"]]
)


cur_p <- plot_data |> ggplot(aes(x = MER, y = abd_norm, fill = MER)) +
    geom_bar(stat = "identity", color = "black", width = 0.7) +
    geom_hline(yintercept = 0) +
    xlab("") +
    ylab(expression(paste("Normalized fold-change (", tilde(Delta[MER[i]]), ")"))) +
    scale_fill_manual(
        name = "Predicted MERs",
        values = MER_annot_colors
    ) +
    scale_y_continuous(
        limits = c(-0.6, 0.6),
        breaks = c(-0.6, -0.3, 0, 0.3, 0.6)
    ) +
    facet_wrap(vars(peak_id), nrow = 1) +
    theme_classic() +
    theme(
        legend.position = "bottom",
        legend.key.size = unit(0.4, 'cm'),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 9),
        panel.grid.major.y = element_line(color = "gray", linewidth = 0.5, linetype = 1),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x=element_blank(),
        strip.text = element_text(face = "bold", margin = margin(0,0,2,0), size=7),
        strip.background = element_blank()
    )

SavePlot(
    cur_p, width = 7.6, height = 2.4,
    dir_name = "msi", file_name = "msi_norm_fc"
)
    
### "=========================================================================="
### Plot Fig. 6f
### "=========================================================================="
PlotDistUpset <- function(
    plot_mat, y_title, file_name, width, y_lim = c(0, 20)
){
    ### Update the colnames with the size
    colnames(plot_mat) <- paste0(colnames(plot_mat), " (", colSums(plot_mat), ")")

    plot_data <- as_tibble(plot_mat, rownames = "peak_id") %>%
        pivot_longer(
            cols = !peak_id,
            names_to = "ROI",
            values_to = "is_member"
        ) %>%
        filter(is_member) %>%
        dplyr::select(-is_member) %>%
        group_by(peak_id) %>%
        summarize(ROIs = list(ROI))

    inc_p <- plot_data %>%
        ggplot(aes(x = ROIs)) +
        geom_bar(
            fill = "black",
            width = 0.7
        ) +
        xlab("MER combinations") +
        scale_y_continuous(
            name = y_title,
            limits = y_lim,
            expand = c(0, 0)
        ) +
        #    breaks = c(0, 10, 20, 30),
        #    name = paste0(
        #        "Number of significantly changed\n",
        #        "and highly-abundant putative\n",
        #        "metabolites (Padj < ", p_thres,")"
        #    )
        # ) +
        theme_classic() +
        theme(
            axis.text = element_text(colour = "black")
        ) +
        scale_x_upset() +
        theme_combmatrix(
            combmatrix.label.text = element_text(colour = "black")
        )

    SavePlot(
        inc_p,
        width = width, height = 2.2,
        dir_name = "msi",
        file_name = file_name
    )
}

####
PlotDistUpset(
    msi_mat > log2(1.2),
    y_title = "Number of increased PM\n(FC to ME-normal > 120%)",
    file_name = "msi_metabolities_upset",
    width = 3.7,
    y_lim = c(0, 60)
)

PlotDistUpset(
    msi_mat < -log2(1.2),
    y_title = "Number of decreased PM\n(FC to ME-normal < 80%)",
    file_name = "msi_metabolities_downset",
    width = 3,
    y_lim = c(0, 60)
)
