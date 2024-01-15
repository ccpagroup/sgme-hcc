### "========================================================================="
### Plot DESI-MSI ROI heatmap
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

### "-------------------------------------------------------------------------"
### Load common data
### "-------------------------------------------------------------------------"
### Get clinical info
clinical_info <- LoadClinicalInfo()

### Load data
msi_mat <- LoadAllMSIMat(
    slide_ids = ref_slide_ids,
    is_normalized = is_normalized_study,
    max_na_fraction = 0.50
)

### Load the reference peaks
load(here("data", "lcms", "lcms_peaks.Rdata"))

### Add class information
lcms_hit_peaks <- AddClassAnnotation(lcms_hit_peaks)

### "-------------------------------------------------------------------------"
### Get all the annotations
### "-------------------------------------------------------------------------"
### Prepare annotations
roi_types <- FindROITypes(msi_mat)

### ROI annotations
roi_annots <- rep(NA, length(roi_types$is_transformed))
roi_annots[roi_types$is_necrotic] <- "Necrotic"
roi_annots[roi_types$is_diff_grade3] <- "G3"
roi_annots[roi_types$is_diff_grade2] <- "G2"
roi_annots[roi_types$is_nontransformed] <- "Normal"
roi_annots[roi_types$is_fibrotic] <- "Fibrotic"
roi_annots <- factor(
    roi_annots,
    levels = c("Necrotic", "G3", "G2", "Fibrotic", "Normal")
)

### Steatotic annotations
is_steatotic <- rep(NA, length(roi_types$is_transformed))
is_steatotic[roi_types$is_steatotic] <- "Steatotic"
is_steatotic <- factor(is_steatotic)

### Case ids
slide_ids <- roi_types$slide_ids
case_ids <- clinical_info$case_id[match(slide_ids, clinical_info$slide_id)]
case_ids <- factor(case_ids)

### "-------------------------------------------------------------------------"
### Prepare sample (row) annotations
### "-------------------------------------------------------------------------"
### Sample marking
show_names <- c(    
    ### Some special case
    P135 = "P135", #"SAPC", ### Transformed 
    P120 = "P120", #"PC(20:4/16:0)", ### Low grade
    N71 = "N71", #"PI(20:2/17:1)", ### High grade    
    N14 = "N14", # "Cholesterol sulfate",    
    P73 = "P73", # "DG(16:0/0:0/18:1)",    
    P58 = "P58", #"LysoPC(18:1/0:0)", ### Normal
    P29 = "P29", #"Laudanosine", ### Normal
    P37 = "P37", #"Oleoylcarnitine", ### Necrotic/fibrotic

    P70 = "P70", # "DG(16:0/0:0/18:3)", ### ME-low-grade+/ME-necrotic- PMs
    N54 = "N54", #"d-Tocotrienol", ### ME_necrotic+ PMs
    P108 = "P108" # "PC(15:0/18:0)" ## ME-steatotic+ PMs
)

sample_marking <- PrepareSampleMarking(msi_mat, show_names, is_peak_id = TRUE)

sample_annotation_list <- PrepareMetaboliteAnnotation(
    msi_mat, metabolite_classes_colors
)

### "-------------------------------------------------------------------------"
### Prepare ROI (column) annotations
### "-------------------------------------------------------------------------"
### Create the legend
roi_annot_legend <- Legend(
    title = "ROI types",
    title_position = "topcenter",
    border = "black",
    legend_gp = gpar(fill = roi_annot_colors),
    labels = c(
        paste0("Normal (", sum(roi_types$is_nontransformed), ")"),
        paste0("Fibrotic (", sum(roi_types$is_fibrotic), ")"),
        paste0("G3 (", sum(roi_types$is_diff_grade3), ")"),
        paste0("G2 (", sum(roi_types$is_diff_grade2), ")"),
        paste0("Necrotic (", sum(roi_types$is_necrotic), ")"),
        paste0("Steatotic (", sum(roi_types$is_steatotic), ")")
    ),
    column_gap = unit(0.4, "inch"),
    row_gap = unit(0.1, "inch"),
    ncol = 1
)

case_dist <- table(case_ids)
case_labels <- paste0(names(case_dist), " (", case_dist, ")")
names(case_labels) <- names(case_dist)
case_labels <- case_labels[names(case_annot_colors)]

case_annot_legend <- Legend(
    title = "Patients",
    title_position = "topcenter",
    border = "black",
    gap = unit(0.1, "inch"),
    legend_gp = gpar(fill = case_annot_colors),
    labels = case_labels,
    column_gap = unit(0.1, "inch"),
    row_gap = unit(0.1, "inch"),
    ncol = 1
)

roi_annotation <- HeatmapAnnotation(
    roi_annots = anno_simple(
        roi_annots,
        col = roi_annot_colors,
        border = TRUE,
        gp = gpar(col = "black"),
        simple_anno_size = unit(0.08, "inch"),
        na_col = "white"
    ),
    steatotic_annots = anno_simple(
        is_steatotic,
        col = steatotic_annot_colors,
        border = TRUE,
        gp = gpar(col = "black"),
        simple_anno_size = unit(0.08, "inch"),
        na_col = "white"
    ),
    case_annots = anno_simple(
        case_ids,
        col = case_annot_colors,
        border = TRUE,
        gp = gpar(col = "black"),
        simple_anno_size = unit(0.08, "inch"),
        na_col = "white"
    ),
    show_annotation_name = FALSE
)

### "-------------------------------------------------------------------------"
### Plot the heatmap
### "-------------------------------------------------------------------------"
roi_order_weight <-
    ### Transformed has the highest
    if_else(roi_types$is_necrotic, 20, 0) +
    if_else(roi_types$is_transformed, 10, 0) +
    if_else(roi_types$is_diff_grade3, 10, 0) +
    if_else(roi_types$is_fibrotic, 9, 0) +
    if_else(roi_types$is_steatotic, 0, 0)

# Pairwise correlation between samples (columns)
cols.cor <- cor(msi_mat, use = "pairwise.complete.obs", method = "pearson")
# Pairwise correlation between rows (genes)
rows.cor <- cor(t(msi_mat), use = "pairwise.complete.obs", method = "pearson")

#row_dend <- dendsort(hclust(as.dist(1 - rows.cor), method = "ward.D2"))
#col_dend <- dendsort(hclust(as.dist(1 - cols.cor), method = "ward.D2"))

row_dend <- dendsort(hclust(dist(msi_mat), method="ward.D2"))
col_dend <- dendsort(hclust(dist(t(msi_mat)), method="ward.D2"))

col_fun <- circlize::colorRamp2(c(-8, 0, 5), c("blue", "white", "red"))

msi_heatmap <- Heatmap(
    msi_mat,
    col = col_fun,
    column_split = 4,
    top_annotation = roi_annotation,
    left_annotation = sample_annotation_list[["metabolite_annotation"]],
    right_annotation = sample_marking,
    #cluster_rows = row_dend,
    cluster_columns = col_dend,
    column_dend_reorder = roi_order_weight,
    #clustering_method_rows = "average", # "ward.D2",
    #clustering_method_columns = "average", # "ward.D2",
    border = TRUE,
    rect_gp = gpar(lty = 0),
    row_title = paste0(
        "Highly-abundant PMs (", nrow(msi_mat), ")"
    ),
    heatmap_legend_param = list(
        title = "Mean normalized\nabundance level (log2)",
        title_position = "leftcenter-rot",
        legend_height = unit(4, "cm"),
        #title_position = "topcenter",
        #direction = "horizontal",
        #legend_width = unit(3, "cm"),
        border = "black"
    ),
    column_title = paste0("ROIs (", ncol(msi_mat), ")"),
    #column_title_side = "bottom",
    show_row_names = FALSE,
    show_column_names = FALSE,
    height = unit(4, "inch"),
    width = unit(4, "inch"),
    use_raster = TRUE,
    raster_quality = 10
)

### Draw png
fs::dir_create(here("figures", "msi"))

width <- 10
height <- 9
png(
    here("figures", "msi", "msi_heatmap.png"),
    width = width, height = height, units = "in", res = 300
)
draw(
    msi_heatmap,
    annotation_legend_list = list(
        sample_annotation_list[["metabolite_classes_legend"]],
        roi_annot_legend,
        case_annot_legend
    ),
    #merge_legend = TRUE,
    legend_gap = unit(0.3, "inch"),
    heatmap_legend_side = "right",
    annotation_legend_side = "bottom"
)
dev.off()

