### "========================================================================="
### Plot MER distributions based on SgME maps
###
### Copyright (c) 2021-2024. Bioinformatics Institute, A*STAR, Singapore
### License: MIT License. See LICENSE for the full license.
### Authors: Lit-Hsin Loo (loolh@bii.a-star.edu.sg)
### "========================================================================="
library(here)
source(here("lib", "MSI_lib.R"))
source(here("lib", "Study_lib.R"))
library(ropls)

### Load the results ==========================================================
plsda_multi <- LoadSgMEClassifier("plsda", "plsda_multi")
plsda_steatotic <- LoadSgMEClassifier("plsda", "plsda_steatotic")
plsda_fibrotic <- LoadSgMEClassifier("plsda", "plsda_fibrotic")

### "=========================================================================="
### Fig. 5c. Prepare multi composition plot_data
### "=========================================================================="
### Get clinical info
clinical_info <- LoadClinicalInfo()

### Load the multi predictor
multi_dist_data <- GetGrpDist(plsda_multi$pred_labels, is_by_case = FALSE)

### Put back group IIIA and IIIB
multi_dist_data <- multi_dist_data |>
    mutate(
        tumor_stage_AJCC_V8 = if_else(
            tumor_stage_AJCC_V8 == "TNM Stage III",
            "TNM Stage IIIB", tumor_stage_AJCC_V8
        ),
        tumor_stage_AJCC_V8 = if_else(
            case_id == "HEP0152",
            "TNM Stage IIIA", tumor_stage_AJCC_V8
        )
    )

### We first order the case_id by total G3
ordered_case_ids <- multi_dist_data %>%
    filter(pred_label == "G3") %>%
    group_by(case_id, tumor_stage_AJCC_V8) %>%
    summarize(
        G3_total = mean(percent),
        .groups = "drop"
    ) %>%
    arrange(tumor_stage_AJCC_V8, G3_total) %>%
    pull(case_id)

### Manually overwrite the order
ordered_case_ids <- c(
    "HEP0321",
    "HEP0214",
    "HEP0194",
    "HEP0277",
    "HEP0152",
    "HEP0319",
    "B003",
    "B008",
    "HEP0268"
)

multi_dist_data$case_id <- factor(
    multi_dist_data$case_id,
    levels = ordered_case_ids
)

ordered_sample_ids <- multi_dist_data %>%
    filter(pred_label == "G3") %>%
    select(sample_id, case_id, section_id, tumor_stage_AJCC_V8, percent) %>%
    arrange(case_id, percent) %>%
    pull(sample_id)

multi_dist_data$sample_id <- factor(
    multi_dist_data$sample_id,
    levels = ordered_sample_ids
)

multi_dist_data$pred_label <- fct_relevel(
    multi_dist_data$pred_label, "Necrotic", "G3", "G2", "Normal"
)

### Prepare steatotic data
pred_labels <- plsda_steatotic$pred_labels
#pred_labels[plsda_multi$pred_labels == "Necrotic"] <- NA

steatotic_dist_data <- GetGrpDist(pred_labels, is_by_case = FALSE)
    
steatotic_dist_data$case_id <- factor(
    steatotic_dist_data$case_id,
    levels = ordered_case_ids
)
steatotic_dist_data$sample_id <- factor(
    steatotic_dist_data$sample_id,
    levels = ordered_sample_ids
)
steatotic_dist_data$pred_label <- fct_relevel(
    steatotic_dist_data$pred_label, "Steatotic", "Nonsteatotic"
)

### Prepare fibrotic data
pred_labels <- plsda_fibrotic$pred_labels
#pred_labels[plsda_multi$pred_labels == "Necrotic"] <- NA
fibrotic_dist_data <- GetGrpDist(pred_labels, is_by_case = FALSE)

fibrotic_dist_data$case_id <- factor(
    fibrotic_dist_data$case_id,
    levels = ordered_case_ids
)
fibrotic_dist_data$sample_id <- factor(
    fibrotic_dist_data$sample_id,
    levels = ordered_sample_ids
)
fibrotic_dist_data$pred_label <- fct_relevel(
    fibrotic_dist_data$pred_label, "Fibrotic", "Nonfibrotic"
)

### "=========================================================================="
### Plot the multi model composition
### "=========================================================================="
PlotPredDist2(
    pred_dist = multi_dist_data,
    y_title = "% of tissue area",
    annot_colors = roi_annot_colors,
    width = 11,
    height = 2,
    file_name = "msi_multi_pred_dist"
)

### "=========================================================================="
### Plot the steatosis model composition
### "=========================================================================="
PlotPredDist2(
    pred_dist = steatotic_dist_data,
    y_title = "% of tissue area",
    annot_colors = steatotic_annot_colors,
    width = 11,
    height = 2,
    file_name = "msi_steatotic_pred_dist"
)

### "=========================================================================="
### Plot the steatosis model composition
### "=========================================================================="
PlotPredDist2(
    pred_dist = fibrotic_dist_data,
    y_title = "% of tissue area",
    annot_colors = fibrotic_annot_colors,
    width = 11,
    height = 2,
    file_name = "msi_fibrotic_pred_dist"
)
