### "========================================================================="
### Plot InterTH and IntraTH
###
### Copyright (c) 2021-2024. Bioinformatics Institute, A*STAR, Singapore
### License: MIT License. See LICENSE for the full license.
### Authors: Lit-Hsin Loo (loolh@bii.a-star.edu.sg)
### "========================================================================="
library("here")
source(here("lib", "MSI_lib.R"))
source(here("lib", "Study_lib.R"))
library(ropls)
library(ComplexHeatmap)

### Load the raw data
load(here("data", "lcms", "lcms_peaks.Rdata"))
load(here("data", "RNA", "RNA_data.Rdata"))

### Summarize the heterogeneity
lcms_InterTH_mat <- lcms_raw_peaks %>%
    filter(
        (S1_vs_N_padj < LCMS_padj_thres & abs(S1_log2FC) > LCMS_log2FC_thres)|
        (S2_vs_N_padj < LCMS_padj_thres & abs(S2_log2FC) > LCMS_log2FC_thres)|
        (S3_vs_N_padj < LCMS_padj_thres & abs(S3_log2FC) > LCMS_log2FC_thres)
    ) %>%
    dplyr::select(mz_id, S1_IETH_type, S2_IETH_type, S3_IETH_type) %>%
    arrange(-S1_IETH_type, -S2_IETH_type, -S3_IETH_type) %>%
    column_to_rownames("mz_id") %>%
    as.matrix()

lcms_IntraTH_mat <- lcms_raw_peaks %>%
    filter(
        (S1_vs_N_padj < LCMS_padj_thres & abs(S1_log2FC) > LCMS_log2FC_thres)|
        (S2_vs_N_padj < LCMS_padj_thres & abs(S2_log2FC) > LCMS_log2FC_thres)|
        (S3_vs_N_padj < LCMS_padj_thres & abs(S3_log2FC) > LCMS_log2FC_thres)
    ) %>%
    dplyr::select(mz_id, S1_IRTH_type, S2_IRTH_type, S3_IRTH_type) %>%
    arrange(-S1_IRTH_type, -S2_IRTH_type, -S3_IRTH_type) %>%
    column_to_rownames("mz_id") %>%
    as.matrix()

lcms_hit_InterTH_mat <- lcms_hit_peaks %>%
    filter(
        (S1_vs_N_padj < LCMS_padj_thres & abs(S1_log2FC) > LCMS_log2FC_thres)|
        (S2_vs_N_padj < LCMS_padj_thres & abs(S2_log2FC) > LCMS_log2FC_thres)|
        (S3_vs_N_padj < LCMS_padj_thres & abs(S3_log2FC) > LCMS_log2FC_thres)
    ) %>%
    dplyr::select(mz_id, S1_IETH_type, S2_IETH_type, S3_IETH_type) %>%
    arrange(-S1_IETH_type, -S2_IETH_type, -S3_IETH_type) %>%
    column_to_rownames("mz_id") %>%
    as.matrix()

lcms_hit_IntraTH_mat <- lcms_hit_peaks %>%
    filter(
        (S1_vs_N_padj < LCMS_padj_thres & abs(S1_log2FC) > LCMS_log2FC_thres)|
        (S2_vs_N_padj < LCMS_padj_thres & abs(S2_log2FC) > LCMS_log2FC_thres)|
        (S3_vs_N_padj < LCMS_padj_thres & abs(S3_log2FC) > LCMS_log2FC_thres)
    ) %>%
    dplyr::select(mz_id, S1_IRTH_type, S2_IRTH_type, S3_IRTH_type) %>%
    arrange(-S1_IRTH_type, -S2_IRTH_type, -S3_IRTH_type) %>%
    column_to_rownames("mz_id") %>%
    as.matrix()

RNA_InterTH_mat <- RNA_raw_data %>%
    filter(
        (S1_vs_N_padj < RNA_padj_thres & abs(S1_log2FC) > RNA_log2FC_thres)|
        (S2_vs_N_padj < RNA_padj_thres & abs(S2_log2FC) > RNA_log2FC_thres)|
        (S3_vs_N_padj < RNA_padj_thres & abs(S3_log2FC) > RNA_log2FC_thres)
    ) %>%
    dplyr::select(gene_id, S1_IETH_type, S2_IETH_type, S3_IETH_type) %>%
    arrange(-S1_IETH_type, -S2_IETH_type, -S3_IETH_type) %>%
    column_to_rownames("gene_id") %>%
    as.matrix()

RNA_IntraTH_mat <- RNA_raw_data %>%
    filter(
        (S1_vs_N_padj < RNA_padj_thres & abs(S1_log2FC) > RNA_log2FC_thres)|
        (S2_vs_N_padj < RNA_padj_thres & abs(S2_log2FC) > RNA_log2FC_thres)|
        (S3_vs_N_padj < RNA_padj_thres & abs(S3_log2FC) > RNA_log2FC_thres)
    ) %>%
    dplyr::select(gene_id, S1_IRTH_type, S2_IRTH_type, S3_IRTH_type) %>%
    arrange(-S1_IRTH_type, -S2_IRTH_type, -S3_IRTH_type) %>%
    column_to_rownames("gene_id") %>%
    as.matrix()

RNA_row_title <- "Changed RNAs"
RNA_column_labels <- c("Stage IB", "Stage II", "Stage III")

lcms_row_title <- "Changed MFs"
lcms_column_labels <- c("Stage IB", "Stage II", "Stage III")

InterTH_legend_title <- "Compared to InterTH-normal sections"
IntraTH_legend_title <- "Compared to InterTH-normal sections"

PlotTHHeatmap(
    TH_mat = RNA_InterTH_mat,
    row_title = RNA_row_title,
    column_labels = RNA_column_labels,
    legend_title = InterTH_legend_title,
    filename = "RNA_InterTH"
)

PlotTHHeatmap(
    TH_mat = RNA_IntraTH_mat,
    row_title = RNA_row_title,
    column_labels = RNA_column_labels,
    legend_title = IntraTH_legend_title,
    filename = "RNA_IntraTH"
)

PlotTHHeatmap(
    TH_mat = lcms_InterTH_mat,
    row_title = lcms_row_title,
    column_labels = lcms_column_labels,
    legend_title = InterTH_legend_title,
    filename = "lcms_InterTH"
)

PlotTHHeatmap(
    TH_mat = lcms_IntraTH_mat,
    row_title = lcms_row_title,
    column_labels = lcms_column_labels,
    legend_title = IntraTH_legend_title,
    filename = "lcms_IntraTH"
)

PlotTHHeatmap(
    TH_mat = lcms_hit_InterTH_mat,
    row_title = lcms_row_title,
    column_labels = lcms_column_labels,
    legend_title = InterTH_legend_title,
    filename = "lcms_hit_InterTH"
)

PlotTHHeatmap(
    TH_mat = lcms_hit_IntraTH_mat,
    row_title = lcms_row_title,
    column_labels = lcms_column_labels,
    legend_title = IntraTH_legend_title,
    filename = "lcms_hit_IntraTH"
)

