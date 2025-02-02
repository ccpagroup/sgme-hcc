### "========================================================================="
### Plot DESI-MSI ROI PLS-DA classifier discriminative PMs
### Note: Step 3 must be run first
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

### Load all MSI data
### Note: NA values are due to missing neg samples for A-03-01
msi_mat <- LoadAllMSIMat(
    slide_ids = ref_slide_ids,
    is_normalized = is_normalized_study,
    max_na_fraction = 0.50
)
roi_types <- FindROITypes(msi_mat)

### Load the raw data
load(here("data", "lcms", "lcms_peaks.Rdata"))

### Add class information
lcms_hit_peaks <- AddClassAnnotation(lcms_hit_peaks)

### "--------------------------------------------------------------------------"
### Load all the models
### "--------------------------------------------------------------------------"
msi_multi <- readRDS(here("results", "plsda", "plsda_multi.rds"))
msi_fibrotic <- readRDS(here("results", "plsda", "plsda_fibrotic.rds"))
msi_steatotic <- readRDS(here("results", "plsda", "plsda_steatotic.rds"))

### "-------------------------------------------------------------------------"
### Get the list of all markers
### "-------------------------------------------------------------------------"
all_markers <- LoadROITopMarkers()

### "-------------------------------------------------------------------------"
### Plot the transformation features
### "-------------------------------------------------------------------------"
PlotVIPvsCoeff(
    model_name = "multi",
    histopath_type = "Normal",
    up_col = roi_annot_colors[["Normal"]],
    down_col = roi_annot_colors[["G3"]],
    x_label = "ME-normal-region coefficient",
    x_lim = c(-0.3, 0.3),
    x_breaks = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3),
    y_lim = c(0.5, 2),
    y_breaks = c(0.5, 1, 1.5, 2),
    is_show_annotation = TRUE,
    top_num_show = 5,
    file_name = "msi_normal_markers"
)

PlotVIPvsCoeff(
    model_name = "multi",
    histopath_type = "G2",
    up_col = roi_annot_colors[["G2"]],
    down_col = roi_annot_colors[["Normal"]],
    x_label = "ME-low-grade-region coefficient",
    x_lim = c(-0.3, 0.3),
    x_breaks = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3),
    y_lim = c(0.5, 2),
    y_breaks = c(0.5, 1, 1.5, 2),
    is_show_annotation = TRUE,
    top_num_show = 5,
    file_name = "msi_G2_markers"
)

PlotVIPvsCoeff(
    model_name = "multi",
    histopath_type = "G3",
    up_col = roi_annot_colors[["G3"]],
    down_col = roi_annot_colors[["Normal"]],
    x_label = "ME-high-grade coefficient",
    x_lim = c(-0.3, 0.3),
    x_breaks = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3),
    y_lim = c(0.5, 2),
    y_breaks = c(0.5, 1, 1.5, 2),
    is_show_annotation = TRUE,
    top_num_show = 5,
    file_name = "msi_G3_markers"
)

PlotVIPvsCoeff(
    model_name = "multi",
    histopath_type = "Necrotic",
    up_col = roi_annot_colors[["Necrotic"]],
    down_col = roi_annot_colors[["Normal"]],
    x_label = "ME-necrotic-region coefficient",
    x_lim = c(-0.3, 0.3),
    x_breaks = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3),
    y_lim = c(0.5, 2),
    y_breaks = c(0.5, 1, 1.5, 2),
    is_show_annotation = TRUE,
    top_num_show = 5,
    file_name = "msi_necrotic_markers"
)

PlotVIPvsCoeff(
    model_name = "fibrotic",
    histopath_type = "Fibrotic",
    up_col = roi_annot_colors[["Fibrotic"]],
    down_col = roi_annot_colors[["Normal"]],
    x_label = "ME-fibrotic-region coefficient",
    x_lim = c(-0.6, 0.6),
    x_breaks = c(-0.6, -0.3, 0, 0.3, 0.6),
    y_lim = c(0.5, 2),
    y_breaks = c(0.5, 1, 1.5, 2),
    is_show_annotation = TRUE,
    top_num_show = 5,
    file_name = "msi_fibrotic_markers"
)

PlotVIPvsCoeff(
    model_name = "steatotic",
    histopath_type = "Steatotic",
    up_col = roi_annot_colors[["Steatotic"]],
    down_col = roi_annot_colors[["Normal"]],
    x_label = "ME-steatotic-region coefficient",
    x_lim = c(-0.3, 0.3),
    x_breaks = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3),
    y_lim = c(0.5, 2),
    y_breaks = c(0.5, 1, 1.5, 2),
    is_show_annotation = TRUE,
    top_num_show = 5,
    file_name = "msi_steatotic_markers"
)
