### "========================================================================="
### Exporting the Supp Data
###
### Copyright (c) 2021-2024. Bioinformatics Institute, A*STAR, Singapore
### License: MIT License. See LICENSE for the full license.
### Authors: Lit-Hsin Loo (loolh@bii.a-star.edu.sg)
### "========================================================================="

library("here")
source(here("lib", "MSI_lib.R"))
source(here("lib", "Study_lib.R"))

### "==========================================================================="
### Export Data S1 LCMS raw data (107 cases x 8742 MFs)
### "==========================================================================="
lcms_mat <- LoadLCMSMat(is_remove_CCA = TRUE)
export_data <- as_tibble(lcms_mat, rownames = "sample_id") %>%
    separate_wider_delim(cols="sample_id", delim = "_", names = c("case_id", "section_id"))

fs::dir_create(here("outputs"))
readr::write_excel_csv(
    export_data,
    here("outputs", "Supp_Data_S1_Raw_LCMS_data.csv")
)

### "==========================================================================="
### Export Data S2 DESI MSI training data (117 PMs x 116 ROIs)
### "==========================================================================="
### Load the summarized MSI ROI data
msi_roi_mat <- LoadAllMSIMat(
    slide_ids = ref_slide_ids,
    is_normalized = is_normalized_study,
    max_na_fraction = 0.50
)

export_data <- as_tibble(t(msi_roi_mat), rownames = "sample_id") %>%
    separate_wider_delim(
        cols = "sample_id", delim = ".", names = c("sample_id", "section_id")
    ) %>%
    separate_wider_delim(
        cols = "section_id", delim = "_", names = c("section_id", "roi_grade")
    ) %>%
    separate_wider_delim(
        cols = "section_id", delim = "-", names = c("section_id", "roi_id")
    ) %>%
    separate_wider_regex(
        cols = "section_id", patterns = c(roi_type = ".", section_id = ".*")
    ) %>%
    mutate(
        section_prefix = case_when(
            roi_type == "N" ~ "N", TRUE ~ "T"
        )
    ) %>%
    unite(
        "section_id", c("section_prefix", "section_id"),
        sep = ""
    ) %>%
    relocate(
        section_id, .after="sample_id"
    )

readr::write_excel_csv(
    export_data,
    here("outputs", "Supp_Data_S2_DESI_ROI_data.csv")
)

### "==========================================================================="
### Export Data S3 Correlated gene sets
### "==========================================================================="
CEG_list_vip <- readRDS(here("results", "RNA", "CEG_list_vip.rds"))

### Save as csv
for (reg_name in names(CEG_list_vip)) {
    readr::write_excel_csv(
        CEG_list_vip[[reg_name]][["export_data"]],
        here("outputs", paste0("RNA_corr_", reg_name, ".csv"))
    )
}

### "==========================================================================="
### Export Table S2 Highly abundant PMs
### "==========================================================================="
### Load the peaks
load(here("data", "lcms", "lcms_peaks.Rdata"))
export_data <- lcms_hit_peaks %>%
    dplyr::select(
        peak_id,
        mz_id, lcms_mode, ion_type,
        iqr, mean, med, pct90,
        MSMS_annotation,
        mz,
        exact_mz
    ) %>% arrange(peak_id)

readr::write_excel_csv(
    export_data,
    here("outputs", "Supp_Table_S2_Highly_abundant_PMs.csv")
)