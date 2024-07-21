### "========================================================================="
### Perform SgMERdeconv and plot all Figure 6
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
library(glmnet)

### Load the classifiers =====================================================
plsda_multi <- LoadSgMEClassifier("plsda", "plsda_multi")
plsda_steatotic <- LoadSgMEClassifier("plsda", "plsda_steatotic")
plsda_fibrotic <- LoadSgMEClassifier("plsda", "plsda_fibrotic")

### Load the SgME maps 
multi_pred_labels <- plsda_multi$pred_labels

steatotic_pred_labels <- plsda_steatotic$pred_labels
#steatotic_pred_labels[multi_pred_labels == "Necrotic"] <- NA

fibrotic_pred_labels <- plsda_fibrotic$pred_labels
#fibrotic_pred_labels[multi_pred_labels == "Necrotic"] <- NA

### Find the percentages of all tissues
tissue_labels <- sub("^(.+)\\.xy.+", "\\1", names(multi_pred_labels))
slide_id <- sub("^([^\\_]+)_.+", "\\1", tissue_labels)
section_id <- sub("^.+\\.([A-Za-z][0-9]*).*-.+$", "\\1", tissue_labels)

clinical_info <- LoadClinicalInfo()

### Generate a tibble with all the ME_region labels
ME_region_info <- tibble(
    tissue_label = tissue_labels,
    slide_id = slide_id,
    section_id = section_id,
    multi_pred_labels = multi_pred_labels,
    steatotic_pred_labels = steatotic_pred_labels,
    fibrotic_pred_labels = fibrotic_pred_labels
)

ME_region_stats <- ME_region_info %>%
    group_by(slide_id, section_id) %>%
    summarise(
        total_area = n(),
        normal_area = sum(multi_pred_labels == "Normal"),
        G2_area = sum(multi_pred_labels == "G2"),
        G3_area = sum(multi_pred_labels == "G3"),
        necrotic_area = sum(multi_pred_labels == "Necrotic"),
        fibro_area = sum(fibrotic_pred_labels == "Fibrotic", na.rm = TRUE),
        steatotic_area = sum(steatotic_pred_labels == "Steatotic", na.rm = TRUE),
        .groups = "drop"
    ) %>%
    left_join(
        clinical_info %>% dplyr::select(
            slide_id, case_id, tumor_stage_AJCC_V8, edmonson_grade,
            metavir_fibrosis_score, steatosis_score, necrosis_score
        ) %>% filter(!is.na(slide_id))
    ) %>%
    mutate(
        sample_id = paste0(case_id, "_", section_id)
    )

#ME_region_stats %>% arrange(tumor_stage_AJCC_V8) %>% print(n=100)


### Load all the LCMS data
lcms_mat <- LoadLCMSMat(is_remove_CCA = TRUE)

### Define the models
model_names <- c("glmnet", "plsda")

#### Test G3
for (model_name in model_names) {
    for (is_hit_peaks_only in c(TRUE, FALSE)) {
        TestSgMEClinicalModel(
            lcms_mat, ME_region_stats,
            model_name = model_name,
            clinical_phenotype = "edmonson_grade",
            ME_region_name = "G3",
            test_xlabel = "Edmonson-Steiner's grade",
            train_xlabel = "Tissue coverage based on\nSgME map (log10[%])",
            is_log10 = TRUE,
            alternative = "less",
            is_hit_peaks_only = is_hit_peaks_only,
            train_dim = c(3, 3),
            test_dim = c(2.5, 4),
            train_ylim = c(-1, 2),
            test_ylim = c(-1, 3)
        )
    }
}

#### Test G2
for (model_name in model_names) {
    for (is_hit_peaks_only in c(TRUE, FALSE)) {
        TestSgMEClinicalModel(
            lcms_mat, ME_region_stats,
            model_name = model_name,
            clinical_phenotype = "edmonson_grade",
            ME_region_name = "G2",
            test_xlabel = "Edmonson-Steiner's grade",
            train_xlabel = "Tissue coverage based on\nSgME map (log10[%])",
            is_log10 = TRUE,
            alternative = "less",
            is_hit_peaks_only = is_hit_peaks_only,
            train_dim = c(3, 3),
            test_dim = c(2.5, 4),
            train_ylim = c(-1, 2),
            train_xlim = c(-1, 2),
            test_ylim = c(-1, 3)
        )
    }
}

#### Test Normal
for (model_name in model_names) {
    for (is_hit_peaks_only in c(TRUE, FALSE)) {
        TestSgMEClinicalModel(
            lcms_mat, ME_region_stats,
            model_name = model_name,
            clinical_phenotype = "edmonson_grade",
            ME_region_name = "normal",
            test_xlabel = "Edmonson-Steiner's grade",
            train_xlabel = "Tissue coverage based on\nSgME map (log10[%])",
            is_log10 = TRUE,
            alternative = "greater",
            is_hit_peaks_only = is_hit_peaks_only,
            train_dim = c(3, 3),
            test_dim = c(2.5, 4),
            train_ylim = c(-1, 2),
            test_ylim = c(-1, 3),
            cat_step_size = 0.08
        )
    }
}

### Test steaotosis
for (model_name in model_names) {
    for (is_hit_peaks_only in c(TRUE, FALSE)) {
        TestSgMEClinicalModel(
            lcms_mat, ME_region_stats,
            model_name = model_name,
            clinical_phenotype = "steatosis_score",
            ME_region_name = "steatotic",
            test_xlabel = "Steatosis score (%)",
            train_xlabel = "Tissue coverage based on\nSgME map (log10[%])",
            is_log10 = TRUE,
            alternative = "greater",
            is_hit_peaks_only = is_hit_peaks_only,
            train_dim = c(3, 3),
            test_dim = c(2.5, 4),
            train_xlim = c(-1, 2),
            train_ylim = c(-1, 2),
            test_ylim = c(-1, 2),
            rnd_seed = 4
        )
    }
}

### Test fibrosis (only look at F0-F3)
for (model_name in model_names) {
    for (is_hit_peaks_only in c(TRUE, FALSE)) {
        TestSgMEClinicalModel(
            lcms_mat,
            ME_region_stats |> filter(metavir_fibrosis_score < 4),
            model_name = model_name,
            clinical_phenotype = "metavir_fibrosis_score",
            ME_region_name = "fibro",
            test_xlabel = "Metavir fibrosis score",
            train_xlabel = "Tissue coverage based on\nSgME Map (log10[%])",
            is_log10 = FALSE,
            alternative = "less",
            is_hit_peaks_only = is_hit_peaks_only,
            train_dim = c(3, 3),
            test_dim = c(3, 4),
            train_xlim = c(0, 100),
            train_ylim = c(0, 100),
            cat_step_size = 0.08
        )
    }
}

### Test necrosis
for (model_name in model_names) {
    for (is_hit_peaks_only in c(TRUE, FALSE)) {
        TestSgMEClinicalModel(
            lcms_mat, ME_region_stats,
            model_name = model_name,
            clinical_phenotype = "necrosis_score",
            ME_region_name = "necrotic",
            test_xlabel = "Necrosis score (%)",
            train_xlabel = "Tissue coverage based on\nSgME Map (log10[%])",
            is_log10 = TRUE,
            alternative = "greater",
            is_hit_peaks_only = is_hit_peaks_only,
            train_dim = c(3, 3),
            test_dim = c(2.5, 4),
            train_xlim = c(-1.5, 2),
            train_ylim = c(-1.5, 2),
            test_xlim = c(0, 100),
            test_ylim = c(-1, 1.5)
        )
    }
}
