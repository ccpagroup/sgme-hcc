### "========================================================================="
### Peform COMET's analysis
### Note: Step 3 must be run first
###
### Copyright (c) 2021-2024. Bioinformatics Institute, A*STAR, Singapore
### License: MIT License. See LICENSE for the full license.
### Authors: Lit-Hsin Loo (loolh@bii.a-star.edu.sg)
### "========================================================================="
library(here)
source(here("lib", "MSI_lib.R"))
source(here("lib", "Study_lib.R"))
library(DESeq2)
library(clusterProfiler)
library(enrichplot)

### Load the LCMS mat
lcms_mat <- LoadLCMSMat(is_remove_CCA = TRUE)

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

### "================================================================="
### Make the lcms and RNA data are the same
### "================================================================="
### We only perform correlation analyis on highly variable transcripts
RNA_iqr <- apply(RNA_vst_mat, 2, IQR)
RNA_vst_mat <- RNA_vst_mat[, RNA_iqr >= log2(2)]

### Get the overlapping samples
common_sample_ids <- intersect(rownames(lcms_mat), rownames(RNA_vst_mat))
lcms_mat <- lcms_mat[common_sample_ids, ]
RNA_vst_mat <- RNA_vst_mat[common_sample_ids, ]

### "========================================================================="
### Find Correlation Genes for top VIP genes
### "========================================================================="
all_markers <- LoadROITopMarkers()

### Find number of dicrminative PMs
all_markers %>%
    filter(vip >= vip_thres_study, abs(coeff) >= coeff_thres_study) %>%
    pull(peak_id) %>%
    unique() %>%
    length()

### Load the raw data
load(here("data", "lcms", "lcms_peaks.Rdata"))

GetPMs <- function(x, dir) {

    selected_markers <- if (dir == "pos") {
        all_markers %>%
            filter(
                histopath_type == x,
                vip >= vip_thres_study,
                coeff >= coeff_thres_study
            )
    } else if (dir == "neg") {
        all_markers %>%
            filter(
                histopath_type == x,
                vip >= vip_thres_study,
                coeff <= -coeff_thres_study
            )

    } else if (dir == "both") {
        all_markers %>%
            filter(
                histopath_type == x,
                vip >= vip_thres_study,
                abs(coeff) >= coeff_thres_study
            )

    } else {
        stop("Unknown direction")
    }

    selected_markers %>%
        left_join(
            lcms_hit_peaks %>% dplyr::select(peak_id, mz_id),
            by = "peak_id"
        ) %>%
        dplyr::pull(mz_id) %>%
        unique()
}

top_vip_list <- list(
    Normal = GetPMs("Normal", "pos"),
    G2 = GetPMs("G2", "pos"),
    G3 = GetPMs("G3", "pos"),
    Necrotic = GetPMs("Necrotic", "pos"),
    Fibrotic = GetPMs("Fibrotic", "pos"),
    Steatotic = GetPMs("Steatotic", "pos")
)

sapply(top_vip_list, length)

top_lcms_list <- lapply(
    top_vip_list, function(x) {
        lcms_mat[, x, drop = FALSE]
    }
)

min_metabolite_coverages <- c(
    "Transformed" = 0.5,
    "Normal" = 0.5,
    "G2" = 0.5,
    "G3" = 0.5,
    "Necrotic" = 0.5,
    "Fibrotic" = 0.5,
    "Steatotic" = 0.5
)

CEG_list_vip <- lapply(
    names(top_lcms_list),
    GetCEGs,
    top_lcms_list,
    RNA_vst_mat,
    corr_padj_thres = 0.10,
    min_metabolite_coverages = min_metabolite_coverages
)
names(CEG_list_vip) <- names(top_lcms_list)

### Uni-directions
# Normal, total metabolites = 16, passed RNAs = 67 (0.5)
# G2, total metabolites = 15, passed RNAs = 225 (0.5)
# G3, total metabolites = 12, passed RNAs = 59 (0.5)
# Necrotic, total metabolites = 13, passed RNAs = 470 (0.5)
# Fibrotic, total metabolites = 3, passed RNAs = 461 (0.5)
# Steatotic, total metabolites = 7, passed RNAs = 185 (0.5)

output_pathname <- here("results", "RNA", "CEG_list_vip.rds")

fs::dir_create(path_dir(output_pathname))
saveRDS(CEG_list_vip, file = output_pathname, compress = "xz")
