### "========================================================================="
### Build SgME map classifiers
###
### Copyright (c) 2021-2025. Bioinformatics Institute, A*STAR, Singapore
### License: MIT License. See LICENSE for the full license.
### Authors: Lit-Hsin Loo (loolh@bii.a-star.edu.sg)
### "========================================================================="
library(here)
source(here("lib", "MSI_lib.R"))
source(here("lib", "Study_lib.R"))
library(ropls)
library(doParallel)
library(caret)

## Register parallel clusters for caret
cl <- makePSOCKcluster(20, autoStop = TRUE)
registerDoParallel(cl)

## When you are done:
#stopCluster(cl)

## Do we want to post proces the SgME Map (such as median filter)?
is_post_processing_study <- TRUE

### "--------------------------------------------------------------------------"
### Load common data
### "--------------------------------------------------------------------------"
### Only need to run once to save all the ROI data
SaveStudyMSIData(
    is_normalized = is_normalized_study,
    peak_summary_mode = peak_summary_mode_study,
    slide_ids = ref_slide_ids
)

### Load DESI-MSI ROI data for the whole study
study_data <- LoadStudyMSIData()

### Load the summarized MSI data
msi_roi_mat <- LoadAllMSIMat(
    slide_ids = ref_slide_ids,
    is_normalized = is_normalized_study,
    max_na_fraction = 0.50
)
roi_types <- FindROITypes(msi_roi_mat)

### Load the raw data as matrix to be applied to construct SgME map
msi_tissue_mat <- LoadAllRawMSIMat(
    slide_ids = test_slide_ids,
    is_normalized = is_normalized_study,
    #max_na_fraction = 0.99,
    is_undetected_min = TRUE
)

### Rename the rows
peak_ids <- sub("^(.+):.+", "\\1", rownames(msi_tissue_mat))
rownames(msi_tissue_mat) <- peak_ids

### Only consider peaks that were found in ROIs
msi_tissue_mat <- msi_tissue_mat[rownames(msi_roi_mat), , drop = FALSE]

### Common configs
impute_nn <- 5
smote_nn <- 3
smote_over_ratio <- 0.5

### "-------------------------------------------------------------------------"
### Build the SVM model for multiple regions
### "-------------------------------------------------------------------------"
### Do the tuning first
svm_multi_tune <- RunSgMEClassifier(
    msi_roi_mat = msi_roi_mat,
    msi_tissue_mat = NULL,
    pos_types = c(
        G2 = "is_diff_grade2",
        G3 = "is_diff_grade3",
        Necrotic = "is_necrotic"
    ),
    neg_types = c(
        Normal = "is_nontransformed"
    ),
    model_type = "svm",
    model_opts = list(
        kernel = "radial"
    ),
    is_post_processing = is_post_processing_study,
    is_impute = TRUE,
    impute_nn = impute_nn,
    is_balancing = FALSE,
    rnd_seed = 0,
    file_name = "svmrbf_multi_tune"
)

### sigma = 0.01, and C = 32 is optimum
svm_multi <- RunSgMEClassifier(
    msi_roi_mat = msi_roi_mat,
    msi_tissue_mat = msi_tissue_mat,
    pos_types = c(
        G2 = "is_diff_grade2",
        G3 = "is_diff_grade3",
        Necrotic = "is_necrotic"
    ),
    neg_types = c(
        Normal = "is_nontransformed" ## Could be F or S
    ),
    model_type = "svm",
    model_opts = list(
        kernel = "radial",
        sigma = 0.01, C = 32
    ),
    is_post_processing = is_post_processing_study,
    is_impute = TRUE,
    impute_nn = impute_nn,
    is_balancing = FALSE,
    smote_nn = smote_nn, smote_over_ratio = smote_over_ratio,
    rnd_seed = 0,
    file_name = "svmrbf_multi"
)

### "-------------------------------------------------------------------------"
### Build the SVM model for multiple regions
### "-------------------------------------------------------------------------"
### Do the tuning first
svmlinear_multi_tune <- RunSgMEClassifier(
    msi_roi_mat = msi_roi_mat,
    msi_tissue_mat = NULL,
    pos_types = c(
        G2 = "is_diff_grade2",
        G3 = "is_diff_grade3",
        Necrotic = "is_necrotic"
    ),
    neg_types = c(
        Normal = "is_nontransformed"
    ),
    model_type = "svm",
    model_opts = list(
        kernel = "linear"
    ),
    is_post_processing = is_post_processing_study,
    is_impute = TRUE,
    impute_nn = impute_nn,
    is_balancing = FALSE,
    rnd_seed = 0,
    file_name = "svmlinear_multi_tune"
)

### C = 0.5 is optimum
svmlnear_multi <- RunSgMEClassifier(
    msi_roi_mat = msi_roi_mat,
    msi_tissue_mat = msi_tissue_mat,
    pos_types = c(
        G2 = "is_diff_grade2",
        G3 = "is_diff_grade3",
        Necrotic = "is_necrotic"
    ),
    neg_types = c(
        Normal = "is_nontransformed" ## Could be F or S
    ),
    model_type = "svm",
    model_opts = list(
        kernel = "linear",
        C = 0.5
    ),
    is_post_processing = is_post_processing_study,
    is_impute = TRUE,
    impute_nn = impute_nn,
    is_balancing = FALSE,
    smote_nn = smote_nn, smote_over_ratio = smote_over_ratio,
    rnd_seed = 0,
    file_name = "svmlinear_multi"
)

### "-------------------------------------------------------------------------"
### Build the PLS-DA model for multiple regions
### "-------------------------------------------------------------------------"
### Do the tuning first
plsda_multi_tune <- RunSgMEClassifier(
    msi_roi_mat = msi_roi_mat,
    msi_tissue_mat = NULL,
    pos_types = c(
        G2 = "is_diff_grade2",
        G3 = "is_diff_grade3",
        Necrotic = "is_necrotic"
    ),
    neg_types = c(
        Normal = "is_nontransformed"
    ),
    model_type = "plsda",
    model_opts = list(orthoI = 0, predI = 10, permI = 0),
    is_post_processing = is_post_processing_study,
    is_impute = FALSE,
    is_balancing = FALSE,
    file_name = "plsda_multi_tune"
)

### 5 is optimum
plsda_multi <- RunSgMEClassifier(
    msi_roi_mat = msi_roi_mat,
    msi_tissue_mat = msi_tissue_mat,
    pos_types = c(
        G2 = "is_diff_grade2",
        G3 = "is_diff_grade3",
        Necrotic = "is_necrotic"
    ),
    neg_types = c(
        Normal = "is_nontransformed" ## Could be F or S
    ),
    model_type = "plsda",
    model_opts = list(orthoI = 0, predI = 5, permI = 0),
    is_post_processing = is_post_processing_study,
    is_impute = FALSE,
    is_balancing = FALSE,
    rnd_seed = 0,
    file_name = "plsda_multi"
)

### "-------------------------------------------------------------------------"
### Build the steatotic peaks (based on ROI p-values)
### "-------------------------------------------------------------------------"
steatotic_peak_ids <- study_data %>%
    filter(
        steatotic_padj < 0.05,
        steatotic_log2FC > 1,
        transformed_padj > 0.1
    ) %>%
    dplyr::pull(
        peak_id
    )

### "-------------------------------------------------------------------------"
### Build the SVM model for steatotic region
### "-------------------------------------------------------------------------"
svmlinear_steatotic_tune <- RunSgMEClassifier(
    msi_roi_mat = msi_roi_mat,
    msi_tissue_mat = NULL,
    selected_peak_ids = steatotic_peak_ids,
    pos_types = c(Steatotic = "is_steatotic"),
    neg_types = list(
        Nonsteatotic = c(
            "is_diff_grade2",
            "is_diff_grade3",
            "is_fibrotic",
            "is_necrotic",
            "is_normal"
        )
    ),
    model_type = "svm",
    model_opts = list(kernel = "linear"),
    is_post_processing = is_post_processing_study,
    is_impute = TRUE,
    impute_nn = impute_nn,
    is_balancing = FALSE,
    rnd_seed = 0,
    file_name = "svmlinear_steatotic_tune"
)

### C = 0.5
svmlinear_steatotic <- RunSgMEClassifier(
    msi_roi_mat = msi_roi_mat,
    msi_tissue_mat = msi_tissue_mat,
    selected_peak_ids = steatotic_peak_ids,
    pos_types = c(Steatotic = "is_steatotic"),
    neg_types = list(
        Nonsteatotic = c(
            "is_diff_grade2",
            "is_diff_grade3",
            "is_fibrotic",
            "is_necrotic",
            "is_normal"
        )
    ),
    model_type = "svm",
    model_opts = list(
        kernel = "linear", 
        C = 0.5
    ),
    is_post_processing = is_post_processing_study,
    is_impute = TRUE,
    impute_nn = impute_nn,
    is_balancing = FALSE,
    rnd_seed = 0,
    file_name = "svmlinear_steatotic"
)

### "-------------------------------------------------------------------------"
### Build the SVM model for steatotic region
### "-------------------------------------------------------------------------"
svmrbf_steatotic_tune <- RunSgMEClassifier(
    msi_roi_mat = msi_roi_mat,
    msi_tissue_mat = NULL,
    selected_peak_ids = steatotic_peak_ids,
    pos_types = c(Steatotic = "is_steatotic"),
    neg_types = list(
        Nonsteatotic = c(
            "is_diff_grade2",
            "is_diff_grade3",
            "is_fibrotic",
            "is_necrotic",
            "is_normal"
        )
    ),
    model_type = "svm",
    model_opts = list(kernel = "radial"),
    is_post_processing = is_post_processing_study,
    is_impute = TRUE,
    impute_nn = impute_nn,
    is_balancing = FALSE,
    rnd_seed = 0,
    file_name = "svmrbf_steatotic_tune"
)

### sigma = 0.1, C = 8
svmrbf_steatotic <- RunSgMEClassifier(
    msi_roi_mat = msi_roi_mat,
    msi_tissue_mat = msi_tissue_mat,
    selected_peak_ids = steatotic_peak_ids,
    pos_types = c(Steatotic = "is_steatotic"),
    neg_types = list(
        Nonsteatotic = c(
            "is_diff_grade2",
            "is_diff_grade3",
            "is_fibrotic",
            "is_necrotic",
            "is_normal"
        )
    ),
    model_type = "svm",
    model_opts = list(
        kernel = "radial",
        sigma = 0.1, C = 8
    ),
    is_post_processing = is_post_processing_study,
    is_impute = TRUE,
    impute_nn = impute_nn,
    is_balancing = FALSE,
    rnd_seed = 0,
    file_name = "svmrbf_steatotic"
)

### "-------------------------------------------------------------------------"
### Build the PLS-DA model for steatotic region
### "-------------------------------------------------------------------------"
plsda_steatotic_tune <- RunSgMEClassifier(
    msi_roi_mat = msi_roi_mat,
    msi_tissue_mat = NULL,
    selected_peak_ids = steatotic_peak_ids,
    pos_types = c(Steatotic = "is_steatotic"),
    neg_types = list(
        Nonsteatotic = c(
            "is_diff_grade2",
            "is_diff_grade3",
            "is_fibrotic",
            "is_necrotic",
            "is_normal"
        )
    ),
    model_type = "plsda",
    model_opts = list(orthoI = 0, predI = 10, permI = 0),
    is_post_processing = is_post_processing_study,
    is_impute = FALSE,
    is_balancing = FALSE,
    file_name = "plsda_steatotic_tune"
)

### 3 is optimum
plsda_steatotic <- RunSgMEClassifier(
    msi_roi_mat = msi_roi_mat,
    msi_tissue_mat = msi_tissue_mat,
    selected_peak_ids = steatotic_peak_ids,
    pos_types = c(Steatotic = "is_steatotic"),
    neg_types = list(
        Nonsteatotic = c(
            "is_diff_grade2",
            "is_diff_grade3",
            "is_fibrotic",
            "is_necrotic",
            "is_normal"
        )
    ),
    model_type = "plsda",
    model_opts = list(orthoI = 0, predI = 3, permI = 0),    
    is_post_processing = is_post_processing_study,
    is_impute = FALSE,
    is_balancing = FALSE,
    file_name = "plsda_steatotic"
)

### "-------------------------------------------------------------------------"
### Prepare the fibrotic region peaks
### "-------------------------------------------------------------------------"
fibrotic_peak_ids <- study_data %>%
    filter(
        fibrotic_padj < 0.05
    ) %>%
    dplyr::pull(
            peak_id
    )

### "-------------------------------------------------------------------------"
### Build the SVM model for fibrotic region
### "-------------------------------------------------------------------------"
svmlinear_fibrotic_tune <- RunSgMEClassifier(
    msi_roi_mat = msi_roi_mat,
    msi_tissue_mat = NULL,
    selected_peak_ids = fibrotic_peak_ids,
    pos_types = c(Fibrotic = "is_fibrotic"),
    neg_types = list(
        Nonfibrotic = c(
            "is_diff_grade2",
            "is_diff_grade3",
            "is_steatotic",
            "is_normal"
        )
    ),
    model_type = "svm",
    model_opts = list(kernel = "linear"),
    is_post_processing = is_post_processing_study,
    is_impute = TRUE,
    is_balancing = TRUE,
    smote_nn = 3, smote_over_ratio = 0.25,
    rnd_seed = 0,
    file_name = "svmlinear_fibrotic_tune"
)

## C = 0.25
svmlinear_fibrotic <- RunSgMEClassifier(
    msi_roi_mat = msi_roi_mat,
    msi_tissue_mat = msi_tissue_mat,
    selected_peak_ids = fibrotic_peak_ids,
    pos_types = c(Fibrotic = "is_fibrotic"),
    neg_types = list(
        Nonfibrotic = c(
            "is_diff_grade2",
            "is_diff_grade3",
            "is_steatotic",
            "is_normal"
        )
    ),
    model_type = "svm",
    model_opts = list(
        kernel = "linear",
        C = 0.25
    ),
    is_post_processing = is_post_processing_study,
    is_impute = TRUE,
    is_balancing = TRUE,
    smote_nn = 3, smote_over_ratio = 0.25,
    rnd_seed = 0,
    file_name = "svmlinear_fibrotic"
)

### "-------------------------------------------------------------------------"
### Build the SVM model for fibrotic region
### "-------------------------------------------------------------------------"
svmrbf_fibrotic_tune <- RunSgMEClassifier(
    msi_roi_mat = msi_roi_mat,
    msi_tissue_mat = NULL,
    selected_peak_ids = fibrotic_peak_ids,
    pos_types = c(Fibrotic = "is_fibrotic"),
    neg_types = list(
        Nonfibrotic = c(
            "is_diff_grade2",
            "is_diff_grade3",
            "is_steatotic",
            "is_normal"
        )
    ),
    model_type = "svm",
    model_opts = list(kernel = "radial"),
    is_post_processing = is_post_processing_study,
    is_impute = TRUE,
    is_balancing = TRUE,
    smote_nn = 3, smote_over_ratio = 0.25,
    rnd_seed = 0,
    file_name = "svmrbf_fibrotic_tune"
)

## sigma = 0.1, C = 2
svmrbf_fibrotic <- RunSgMEClassifier(
    msi_roi_mat = msi_roi_mat,
    msi_tissue_mat = msi_tissue_mat,
    selected_peak_ids = fibrotic_peak_ids,
    pos_types = c(Fibrotic = "is_fibrotic"),
    neg_types = list(
        Nonfibrotic = c(
            "is_diff_grade2",
            "is_diff_grade3",
            "is_steatotic",
            #"is_necrotic",
            "is_normal"
        )
    ),
    model_type = "svm",
    model_opts = list(
        kernel = "radial",
        sigma = 0.1, C = 2
    ),
    is_post_processing = is_post_processing_study,
    is_impute = TRUE,
    is_balancing = TRUE,
    smote_nn = 3, smote_over_ratio = 0.25,
    rnd_seed = 0,
    file_name = "svmrbf_fibrotic"
)

### "-------------------------------------------------------------------------"
### Build the PLS-DA model for fibrotic region
### "-------------------------------------------------------------------------"
plsda_fibrotic_tune <- RunSgMEClassifier(
    msi_roi_mat = msi_roi_mat,
    msi_tissue_mat = NULL,
    selected_peak_ids = fibrotic_peak_ids,
    pos_types = c(Fibrotic = "is_fibrotic"),
    neg_types = list(
        Nonfibrotic = c(
            "is_diff_grade2",
            "is_diff_grade3",
            "is_steatotic",
            "is_normal"
        )
    ),
    model_type = "plsda",
    model_opts = list(orthoI = 0, predI = 8, permI = 0),
    is_post_processing = is_post_processing_study,
    is_impute = FALSE,
    is_balancing = TRUE,
    smote_nn = 3, smote_over_ratio = 0.25,
    file_name = "plsda_fibrotic_tune"
)

### 3 is optimum
plsda_fibrotic <- RunSgMEClassifier(
    msi_roi_mat = msi_roi_mat,
    msi_tissue_mat = msi_tissue_mat,
    selected_peak_ids = fibrotic_peak_ids,
    pos_types = c(Fibrotic = "is_fibrotic"),
    neg_types = list(
        Nonfibrotic = c(
            "is_diff_grade2",
            "is_diff_grade3",
            "is_steatotic",
            #"is_necrotic",
            "is_normal"
        )
    ),
    model_type = "plsda",
    model_opts = list(orthoI = 0, predI = 3, permI = 0),
    is_post_processing = is_post_processing_study,
    is_impute = FALSE,
    is_balancing = TRUE,
    smote_nn = 3, smote_over_ratio = 0.25,
    file_name = "plsda_fibrotic"
)

### "========================================================================="
### Fig. 6ab and S9. Save all the prediction as images
### "========================================================================="
plsda_multi <- LoadSgMEClassifier("plsda", "plsda_multi")
saveRDS(plsda_multi, here("results", "plsda", "msi_multi.rds"))

SaveSgMEMapAsImage(
    pred_model = plsda_multi,
    model_name = "plsda_multi"
)

plsda_steatotic <- LoadSgMEClassifier("plsda", "plsda_steatotic")
saveRDS(plsda_steatotic, here("results", "plsda", "msi_steatotic.rds"))

SaveSgMEMapAsImage(
    pred_model = plsda_steatotic,
    model_name = "plsda_steatotic"
)

plsda_fibrotic <- LoadSgMEClassifier("plsda", "plsda_fibrotic")
saveRDS(plsda_fibrotic, here("results", "plsda", "msi_fibrotic.rds"))

SaveSgMEMapAsImage(
    pred_model = plsda_fibrotic,
    model_name = "plsda_fibrotic"
)



