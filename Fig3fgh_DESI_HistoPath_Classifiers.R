### "========================================================================="
### Plot DESI-MSI ROI PLS-DA classifiers
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
msi_multi <- readRDS(here("results", "plsda", "msi_multi.rds"))
msi_fibrotic <- readRDS(here("results", "plsda", "msi_fibrotic.rds"))
msi_steatotic <- readRDS(here("results", "plsda", "msi_steatotic.rds"))

### "-------------------------------------------------------------------------"
### Get the list of all markers
### "-------------------------------------------------------------------------"
all_markers <- LoadROITopMarkers()

### "-------------------------------------------------------------------------"
### Plot the transformation region classifier
### "-------------------------------------------------------------------------"
SavePLSDAPlot(
    width = 3.7, height = 3.7,
    msi_multi,
    fct_relevel(
        msi_multi$grps,
        "G2", "G3", "Necrotic", "Normal"
    ),
    plot_components = c(1, 2),
    grp_labels = c(
        "Normal" = "0",
        "G2" = "2",
        "G3" = "3",
        "Necrotic" = "T"
    ),
    grp_colors = c(
        "Normal" = roi_annot_colors[["Normal"]],
        "G3" = roi_annot_colors[["G3"]],
        "G2" = roi_annot_colors[["G2"]],
        "Necrotic" = roi_annot_colors[["Necrotic"]]
    ),
    "DESI-MSI profiles",
    "msi_multi_12.pdf"
)

SavePLSDAPlot(
    width = 3.7, height = 3.7,
    msi_multi,
    fct_relevel(
        msi_multi$grps,
        "G2", "G3", "Necrotic", "Normal"
    ),
    plot_components = c(3, 2),
    grp_labels = c(
        "Normal" = "0",
        "G2" = "2",
        "G3" = "3",
        "Necrotic" = "T"
    ),
    grp_colors = c(
        "Normal" = roi_annot_colors[["Normal"]],
        "G3" = roi_annot_colors[["G3"]],
        "G2" = roi_annot_colors[["G2"]],
        "Necrotic" = roi_annot_colors[["Necrotic"]]
    ),
    "DESI-MSI profiles",
    "msi_multi_23.pdf"
)

SavePLSDAPlot(
    width = 3.7, height = 3.7,
    msi_fibrotic,
    fct_relevel(msi_fibrotic$grps, "Fibrotic", "Nonfibrotic"),
    plot_components = c(1, 2),
    grp_labels = c(
        "Nonfibrotic" = "0",
        "Fibrotic" = "F"
    ),
    grp_colors = c(
        "Nonfibrotic" = roi_annot_colors[["Normal"]],
        "Fibrotic" = roi_annot_colors[["Fibrotic"]]
    ),
    "DESI-MSI profiles",
    "msi_fibrotic.pdf"
)

SavePLSDAPlot(
    width = 3.7, height = 3.7,
    msi_steatotic,
    fct_relevel(
        msi_steatotic$grps,
        "Nonsteatotic", "Steatotic"
    ),
    plot_components = c(1, 2),
    grp_labels = c(
        "Nonsteatotic" = "0",
        "Steatotic" = "S"
    ),
    grp_colors = c(
        "Nonsteatotic" = roi_annot_colors[["Normal"]],
        "Steatotic" = roi_annot_colors[["Steatotic"]]
    ),
    "DESI-MSI profiles",
    "msi_steatotic.pdf"
)