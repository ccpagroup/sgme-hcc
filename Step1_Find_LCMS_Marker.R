### "========================================================================="
### Train LC-MS PLS-DA models
###
### Copyright (c) 2021-2024. Bioinformatics Institute, A*STAR, Singapore
### License: MIT License. See LICENSE for the full license.
### Authors: Lit-Hsin Loo (loolh@bii.a-star.edu.sg)
### "========================================================================="
library("here")
source(here("lib", "MSI_lib.R"))
source(here("lib", "Study_lib.R"))
library(ropls)

### Define the dataset
msi_conf <- read_json(here("conf", "msi_conf.json"))

### "--------------------------------------------------------------------------"
### Compare HCC and CCA
### "--------------------------------------------------------------------------"
### Load the LCMS mat
lcms_mat <- LoadLCMSMat(is_remove_CCA = FALSE)
lcms_grps <- GetSampleGroups(lcms_mat)

### Find the optimum component to separate tumors
lcms_tumor_tune <- opls(
    lcms_mat, lcms_grps$tumor_grps,
    crossvalI = 5, predI = 10
)
 
fs::dir_create(here("results", "plsda"))
saveRDS(
    lcms_tumor_tune,
    file = here("results", "plsda", "lcms_tumor_tune.rds"),
    compress = "xz"
)

### 9 is optimum. So, retrain the model
lcms_tumor <- opls(
    lcms_mat, lcms_grps$tumor_grps,
    crossvalI = 5, predI = 9
)
save(
    lcms_tumor, lcms_grps,
    file = here("results", "plsda", "lcms_tumor.Rdata"),
    compress = "xz"
)

### "--------------------------------------------------------------------------"
### Determine the optimum component for LCMS (12)
### "--------------------------------------------------------------------------"
### Reload the LCMS mat without CCA
lcms_mat <- LoadLCMSMat(is_remove_CCA = TRUE)
lcms_grps <- GetSampleGroups(lcms_mat)

### We determine the optimum components by trying opls up to 12 components.
lcms_hcc_stage_tune <- opls(
    lcms_mat, lcms_grps$hcc_stage_grps,
    crossvalI = 5, predI = 15
)
saveRDS(
    lcms_hcc_stage_tune,
    file = here("results", "plsda", "lcms_hcc_stage_tune.rds"),
    compress = "xz"
)

### Optimum is 12
### We rerun the final model.
lcms_hcc_stage <- opls(
    lcms_mat, lcms_grps$hcc_stage_grps,
    crossvalI = 5, predI = 12
)
save(
    lcms_hcc_stage,
    lcms_grps,
    file = here("results", "plsda", "lcms_hcc_stage.Rdata"),
    compress = "xz"
)

### "--------------------------------------------------------------------------"
### Determine the optimum component for RNAseq (10)
### "--------------------------------------------------------------------------"
### Load RNA data (RNA_vst_mat)
tmp <- LoadRNAData()
RNA_mat <- tmp$RNA_vst_mat
all_gene_info <- tmp$all_gene_info
rm(tmp)

### Get the group info
RNA_grps <- GetSampleGroups(RNA_mat)

### Try up to 15
RNA_hcc_stage_tune <- opls(
    RNA_mat, RNA_grps$hcc_stage_grps, crossvalI = 5, predI = 15
)
saveRDS(
    RNA_hcc_stage_tune,
    file = here("results", "plsda", "RNA_hcc_stage_tune.rds"),
    compress = "xz"
)

### Optimum is 10. ### Rerun the final model
RNA_hcc_stage <- opls(
    RNA_mat, RNA_grps$hcc_stage_grps,
    crossvalI = 5, predI = 10
)
save(
    RNA_hcc_stage, RNA_grps,
    file = here("results", "plsda", "RNA_hcc_stage.Rdata"),
    compress = "xz"
)

### "--------------------------------------------------------------------------"
### Run and load the stats for bulk LCMS and RNA data
### "--------------------------------------------------------------------------"
### Load the bulk sample stats
lcms_stats <- readRDS(here("data", "lcms", "lcms_stats.rds"))
RNA_stats <- readRDS(here("data", "RNA", "RNA_stats.rds"))

### "-------------------------------------------------------------------------"
### Find the final raw lcms data
### "-------------------------------------------------------------------------"
### Load the MS/MS data
msms_peaks <- LoadMSMSData()

### We only consider msms_peaks that have high abundances
pct90_thres <- 19

### Create a table for all the peaks
lcms_raw_peaks <- tibble(
        mz_id = lcms_stats$mz_id,
        lcms_mode = str_sub(mz_id, start = -4),
        ion_type = tolower(str_sub(mz_id, start = -3)),
        iqr = apply(lcms_mat[, mz_id], 2, IQR),
        mean = apply(lcms_mat[, mz_id], 2, mean),
        med = apply(lcms_mat[, mz_id], 2, median),
        pct90 = apply(lcms_mat[, mz_id], 2, quantile, 0.9),
        is_tested = mz_id %in% msms_peaks$mz_id,
        is_msms = mz_id %in% msms_peaks$mz_id[msms_peaks$is_confirmed]
    ) %>%
    mutate(
        is_tested = if_else(pct90 < pct90_thres, FALSE, is_tested),
        is_msms   = if_else(pct90 < pct90_thres, FALSE, is_msms)
    )

### Add the annotations
lcms_raw_peaks <- lcms_raw_peaks %>%
    left_join(
        msms_peaks %>%
            dplyr::select(mz_id, MSMS_annotation, Comment, mz, exact_mz),
        by = "mz_id"
    )

### Create a table for all the genes
RNA_raw_data <- tibble(
    gene_id = RNA_stats$gene_id
)

### Add the P-values for IETH and IRTH
lcms_raw_peaks <- AddStats(lcms_raw_peaks, lcms_stats, "mz_id")
RNA_raw_data  <- AddStats(RNA_raw_data, RNA_stats, "gene_id")

### We select highly abundant peaks to study their spatial
### distributions and microenvironment heterogeneity. Note: not all of them
### may be the markers of the early stage
lcms_hit_peaks <- lcms_raw_peaks %>%
    filter(
        is_msms == TRUE,
        ### Only select DESI peaks
        mz > 40, mz < 1210
    )

### Add ID based on unique mz
lcms_hit_peaks <- AddPeakID(lcms_hit_peaks)

fs::dir_create(here("data", "lcms"))
save(
    lcms_hit_peaks, lcms_raw_peaks,
    file = here("data", "lcms", "lcms_peaks.Rdata"),
    compress = "xz"
)

fs::dir_create(here("data", "RNA"))
save(
    RNA_raw_data,
    file = here("data", "RNA", "RNA_data.Rdata"),
    compress = "xz"
)