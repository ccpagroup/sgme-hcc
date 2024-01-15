### "========================================================================="
### Study-specific common functions
###
### Copyright (c) 2021-2024. Bioinformatics Institute, A*STAR, Singapore
### License: MIT License. See LICENSE for the full license.
### Authors: Lit-Hsin Loo (loolh@bii.a-star.edu.sg)
### "========================================================================="

### Data management ==========================================================
### "========================================================================="
### Add stats
### "========================================================================="
AddStats <- function(
    raw_data, stats, id_name
) {
    left_join(
        raw_data, stats,
        by = id_name
    ) %>%
    mutate(
        ### This is now computed by Step 1
        #TvsN_padj = p.adjust(TvsN_pval, method = "BH"),
        #S1_vs_N_padj = p.adjust(S1_vs_N_pval, method = "BH"),
        #S2_vs_N_padj = p.adjust(S2_vs_N_pval, method = "BH"),
        #S3_vs_N_padj = p.adjust(S3_vs_N_pval, method = "BH"),

        S1_IETH_padj = p.adjust(S1_IETH_pval, method = "BH"),
        S2_IETH_padj = p.adjust(S2_IETH_pval, method = "BH"),
        S3_IETH_padj = p.adjust(S3_IETH_pval, method = "BH"),
        S1_IRTH_padj = p.adjust(S1_IRTH_pval, method = "BH"),
        S2_IRTH_padj = p.adjust(S2_IRTH_pval, method = "BH"),
        S3_IRTH_padj = p.adjust(S3_IRTH_pval, method = "BH"),

        S1_IETH_log2FC = S1_IETH_median - normal_IETH_median,
        S2_IETH_log2FC = S2_IETH_median - normal_IETH_median,
        S3_IETH_log2FC = S3_IETH_median - normal_IETH_median,
        S1_IRTH_log2FC = S1_IRTH_median - normal_IETH_median,
        S2_IRTH_log2FC = S2_IRTH_median - normal_IETH_median,
        S3_IRTH_log2FC = S3_IRTH_median - normal_IETH_median,

        S1_IETH_type = case_when(
            S1_IETH_padj < 0.05 & S1_IETH_log2FC > 0 ~ 1,
            S1_IETH_padj < 0.05 & S1_IETH_log2FC < 0 ~ -1,
            TRUE ~ 0
        ),
        S2_IETH_type = case_when(
            S2_IETH_padj < 0.05 & S2_IETH_log2FC > 0 ~ 1,
            S2_IETH_padj < 0.05 & S2_IETH_log2FC < 0 ~ -1,
            TRUE ~ 0
        ),
        S3_IETH_type = case_when(
            S3_IETH_padj < 0.05 & S3_IETH_log2FC > 0 ~ 1,
            S3_IETH_padj < 0.05 & S3_IETH_log2FC < 0 ~ -1,
            TRUE ~ 0
        ),
        S1_IRTH_type = case_when(
            S1_IRTH_padj < 0.05 & S1_IRTH_log2FC > 0 ~ 1,
            S1_IRTH_padj < 0.05 & S1_IRTH_log2FC < 0 ~ -1,
            TRUE ~ 0
        ),
        S2_IRTH_type = case_when(
            S2_IRTH_padj < 0.05 & S2_IRTH_log2FC > 0 ~ 1,
            S2_IRTH_padj < 0.05 & S2_IRTH_log2FC < 0 ~ -1,
            TRUE ~ 0
        ),
        S3_IRTH_type = case_when(
            S3_IRTH_padj < 0.05 & S3_IRTH_log2FC > 0 ~ 1,
            S3_IRTH_padj < 0.05 & S3_IRTH_log2FC < 0 ~ -1,
            TRUE ~ 0
        )
    )
}

### "========================================================================="
### Load RNA Mat
### (including preprocessing)
### "========================================================================="
LoadRNAData <- function() {
    load(here("data", "RNA", "RNA_normalized_data.Rdata"))

    ### Return
    list(
        dds = dds,
        RNA_vst_mat = RNA_vst_mat,
        all_gene_info = all_gene_info
    )
}

### "========================================================================="
### Load LCMS Mat
### (including preprocessing)
### "========================================================================="
LoadLCMSMat <- function(
    is_remove_CCA = TRUE,
    is_averaging_tumor = FALSE
) {
    ### Load all the raw data
    lcms_raw <- readRDS(here("data", "lcms", "lcms_raw.rds"))

    ### Remove the two outliers
    lcms_raw <- lcms_raw %>%
        filter(
            !(case_id == "HEP0321" & section_id == "N"), ## LNEG data is missing!
            !(case_id == "B014" & section_id == "T3")
        )

    ### Do we need to average all the tumors?
    if (is_averaging_tumor) {
        lcms_raw <- lcms_raw %>%
            mutate(
                section_type = if_else(
                    startsWith(section_id, "T"), "T", "N", NA
                )
            ) %>%
            group_by(
                case_id, mz_id, section_type
            ) %>%
            summarize(
                lcms_abd = mean(lcms_abd, na.rm = TRUE),
                .groups = "drop"
            ) %>%
            unite(
                "sample_id",
                c("case_id", "section_type")
            )

    } else {
        lcms_raw <- lcms_raw %>%
            unite("sample_id", case_id:section_id)
    }

    ### Convert the data to a matrix
    lcms_mat <- lcms_raw %>%
        select(sample_id, mz_id, lcms_abd) %>%
        pivot_wider(
            names_from = mz_id, values_from = lcms_abd, values_fill = 0
        ) %>%
        column_to_rownames("sample_id") %>%
        as.matrix()

    ### Check the distribution of abundance levels to find out where is the
    ### detection limits
    # plot(density(lcms_mat))

    ### Remove very low-level features
    # is_low <- (lcms_mat > 0) & (lcms_mat < 5)
    # low_num <- sum(is_low)
    # lcms_mat[is_low] <- 0
    # if (low_num > 0) {
    #    message("* Setting ", low_num, " low-level features to zero.")
    # }

    ### Remove constant features
    mz_all_sds <- apply(lcms_mat, 2, sd)

    num_before <- ncol(lcms_mat)
    lcms_mat <- lcms_mat[, mz_all_sds > 0.1]
    num_after <- ncol(lcms_mat)

    if (num_after != num_before) {
        message(
            "* Removed ", num_before - num_after,
            " near-constant features from ", num_before, " features."
        )
    }

    ### Check for outliers
    #lcms_pca <- opls(lcms_mat)
    #plot(lcms_pca, typeVc = "x-score")

    ### Do some filtering
    iqrs <- apply(lcms_mat, 2, IQR)
    sums <- apply(lcms_mat, 2, sum)

    ### View the plot to decide on the thresholds
    # plot(x = iqrs, y = sums)
    lcms_mat <- lcms_mat[, (sums > 500) | (iqrs > 5)]

    ### Do we need to remove CCA?
    if (is_remove_CCA) {
        ### Remove CCA
        lcms_mat <- lcms_mat[
            !startsWith(rownames(lcms_mat), "HEP0275_T") &
            !startsWith(rownames(lcms_mat), "HEP0290_T"),
        ]
    }

    ### Return
    lcms_mat
}

### "========================================================================="
### Load the raw MSMS data
### "========================================================================="
LoadRawMSMSData <- function() {
    #### Load Set A
    SetA <- read_excel(
        here("data", "lcms", "20230314_HCC sig lcms features+ID_updated6.xlsx"),
        sheet = "Set A",
        range = cell_cols("A:H")
    )

    #### Load Set B
    SetB <- read_excel(
        here("data", "lcms", "20230314_HCC sig lcms features+ID_updated6.xlsx"),
        sheet = "Set B",
        range = cell_cols("A:H")
    )

    ### return
    rbind(SetA, SetB)
}

### "========================================================================="
### Load MSMS dataset
### "========================================================================="
LoadMSMSData <- function() {

    ### Combine the two datasets
    msms_data <- LoadRawMSMSData()

    ### Determine the mz based on MSMS peaks
    msms_data$mz <- sapply(
            msms_data$Result,
            function(tmp) {
                mean(
                    as.double(
                        sapply(
                            str_split(
                                str_split(tmp, ", ")[[1]],
                                "_"
                            ),
                            function(x) {
                                if (length(x) > 1) {
                                    x[[2]]
                                } else {
                                    NA
                                }
                            }
                        )
                    )
                )
            }
        )


    msms_data <- msms_data %>%
        rename(
            #`Searched ID` = "MSMS_annotation",
            `short_name` = "MSMS_annotation",
            `LCMS mode` = "lcms_mode"
        ) %>%
        mutate(
            mz_id = paste0(Compound, "_", lcms_mode),
            ion_type = tolower(substr(lcms_mode, 2, 4))
        ) %>%
        replace_na(
            list(MSMS_annotation = "Unannotated")
        ) %>%
        mutate(
            is_confirmed = case_when(
                `Result` %in% c("x", "X") ~ FALSE,
                is.na(mz) ~ FALSE,
                TRUE ~ TRUE
            )
        ) %>%
        select(
            mz_id, ion_type, mz, is_confirmed,
            MSMS_annotation, exact_mz, Result, Comment
        ) %>%
        distinct() %>%
        arrange(mz)

    msms_data
}

### "========================================================================"
### Determine the fold change and IETH and IATH
### data_mat = sample_ids x feature_ids
### "========================================================================"
### Delta function to get heterogeneity measurements
QuantifyHeterogeneity <- function(
    tumor_samples, prepend_name, normal_inter_pdist
) {
    ### Get the tumor inter-tumor heterogeneity
    tumor_means <- sapply(tumor_samples, mean, na.rm = TRUE)
    tumor_inter_pdist <- as.numeric(dist(tumor_means))

    ### Get the tumor intra-tumor heterogeneity
    tumor_intra_pdists <- unlist(lapply(
        tumor_samples, function(y) as.numeric(dist(y))
    ))

    ### Find the difference to normal ITETH
    tumor_inter_pdist_pval <- wilcox.test(
        tumor_inter_pdist, normal_inter_pdist,
        na.rm = TRUE
    )$p.val

    tumor_intra_pdist_pval <- wilcox.test(
        tumor_intra_pdists, normal_inter_pdist,
        na.rm = TRUE
    )$p.val

    ### Give the output
    res <- c(
        IETH_median = median(tumor_inter_pdist),
        IETH_pval   = tumor_inter_pdist_pval,
        IRTH_median = median(tumor_intra_pdists),
        IRTH_pval   = tumor_intra_pdist_pval
    )
    names(res) <- paste0(prepend_name, "_", names(res))

    res
}

### "========================================================================="
### Save pls-da plot
### "========================================================================="
SavePLSDAPlot <- function(
    width, height, final_model, grp_ids,
    plot_components = c(1, 2),
    grp_labels, grp_colors,
    title_txt, filename
) {
    ### Get the training label
    training_true_labels <- final_model$model@suppLs$y[
        final_model$data_types == "balanced"
    ]

    fs::dir_create(here("figures", "plsda"))

    pdf(
        file = here("figures", "plsda", filename),
        width = width, height = height
     )
    par(pty="s")
    ### Get all labels
    plot(
        final_model$model,
        typeVc = "x-score",
        parCompVi = plot_components,
        parPaletteVc = grp_colors[levels(grp_ids)],
        parTitleL = FALSE
    )

    for (grp_cur in levels(training_true_labels)) {
        sample_idx <- training_true_labels == grp_cur
        sample_col <- grp_colors[grp_cur]

        points(
            x = final_model$model@scoreMN[sample_idx, plot_components[1]],
            y = final_model$model@scoreMN[sample_idx, plot_components[2]],
            pch = 21,
            cex = 1,
            bg = sample_col
        )
    }
    #plot(
    #    final_model$model,
    #    parCexN = 0.8,
    #    parCompVi = plot_components,
    #    typeVc = "x-score",
    #    parLabVc = all_labels,
    #    parPaletteVc = grp_colors[levels(grp_ids)],
    #    parTitleL = FALSE
    #)
    title(title_txt)
    dev.off()

    ### Print the number of samples
    print(table(training_true_labels))
}

### "========================================================================="
### Find the sample info
### "========================================================================="
GetSampleGroups <- function(data_mat) {
    ### Load the clinical information
    clinical_info <- LoadClinicalInfo()

    ### Get the case ids
    case_ids <- gsub("^([^_]+)_.*$", "\\1", rownames(data_mat))
    sample_idx <- match(case_ids, clinical_info$case_id)

    ### Check if any missing
    if (any(is.na(sample_idx))) {
        stop("Data matrix has missing samples")
    }

    ### "----------------------------------------------------------------------"
    ### Get the key classes
    ### "----------------------------------------------------------------------"
    is_tumor <- (gsub("^[^_]+_(.).*$", "\\1", rownames(data_mat)) == "T")
    is_hcc <- (clinical_info$histological_diagnosis[sample_idx] == "HCC")
    is_cca <- (clinical_info$histological_diagnosis[sample_idx] == "CCA")

    ### "----------------------------------------------------------------------"
    ### Tumor groups
    ### "----------------------------------------------------------------------"
    tumor_grps <- rep(NA, length(case_ids))
    tumor_grps[!is_tumor] <- paste0(
        "Normal (", length(case_ids[!is_tumor]), ")"
    )
    tumor_grps[is_tumor & is_cca] <- paste0(
        "CCA Tumor (", length(case_ids[is_tumor & is_cca]), ")"
    )
    tumor_grps[is_tumor & is_hcc] <- paste0(
        "HCC Tumor (", length(case_ids[is_tumor & is_hcc]), ")"
    )
    tumor_grps <- factor(tumor_grps)

    ### "----------------------------------------------------------------------"
    ### Divide all HCC tumors according grades
    ### "----------------------------------------------------------------------"
    hcc_grade_grps <- rep("Normal", length(case_ids))
    hcc_grade_grps[is_tumor & is_cca] <- "CCA"

    is_hcc_grade_2 <- is_tumor & is_hcc &
        (clinical_info$edmonson_grade[sample_idx] == "2")

    is_hcc_grade_3_4 <- is_tumor & is_hcc &
        (clinical_info$edmonson_grade[sample_idx] %in% c("3", "4"))

    hcc_grade_grps[is_hcc_grade_2] <- "Grade 2"
    hcc_grade_grps[is_hcc_grade_3_4] <- "Grade 3/4"
    hcc_grade_grps <- factor(hcc_grade_grps)

    ### "----------------------------------------------------------------------"
    ### Divide all HCC tumors according stages
    ### "----------------------------------------------------------------------"
    hcc_stage_grps <- rep("Normal", length(case_ids))
    hcc_stage_grps[is_tumor & is_cca] <- "CCA"

    is_hcc_stage_I <- is_tumor & is_hcc &
        (clinical_info$tumor_stage_AJCC_V8[sample_idx] == "TNM Stage IB")
    is_hcc_stage_II <- is_tumor & is_hcc &
        (clinical_info$tumor_stage_AJCC_V8[sample_idx] == "TNM Stage II")
    is_hcc_stage_III <- is_tumor & is_hcc &
        ((clinical_info$tumor_stage_AJCC_V8[sample_idx] == "TNM Stage IIIA") |
        (clinical_info$tumor_stage_AJCC_V8[sample_idx] == "TNM Stage IIIB"))

    hcc_stage_grps[is_hcc_stage_I] <- "HCC Stage IB"
    hcc_stage_grps[is_hcc_stage_II] <- "HCC Stage II"
    hcc_stage_grps[is_hcc_stage_III] <- "HCC Stage III"

    hcc_stage_grps <- factor(
        hcc_stage_grps,
        levels = c("Normal", "HCC Stage IB", "HCC Stage II", "HCC Stage III")
    )

    ### "----------------------------------------------------------------------"
    ### Stages vs normal samples
    ### "----------------------------------------------------------------------"
    S1_vs_N_grps <- rep(NA, length(case_ids))
    S1_vs_N_grps[!is_tumor] <- "Normal"
    S1_vs_N_grps[is_hcc_stage_I] <- "Stage IB"
    S1_vs_N_grps <- factor(S1_vs_N_grps, levels = c("Normal", "Stage IB"))

    S2_vs_N_grps <- rep(NA, length(case_ids))
    S2_vs_N_grps[!is_tumor] <- "Normal"
    S2_vs_N_grps[is_hcc_stage_II] <- "Stage II"
    S2_vs_N_grps <- factor(S2_vs_N_grps, levels = c("Normal", "Stage II"))

    S2_vs_S1_grps <- rep(NA, length(case_ids))
    S2_vs_S1_grps[is_hcc_stage_I] <- "Stage IB"
    S2_vs_S1_grps[is_hcc_stage_II] <- "Stage II"
    S2_vs_S1_grps <- factor(S2_vs_S1_grps, levels = c("Stage IB", "Stage II"))

    S3_vs_N_grps <- rep(NA, length(case_ids))
    S3_vs_N_grps[!is_tumor] <- "Normal"
    S3_vs_N_grps[is_hcc_stage_III] <- "Stage III"
    S3_vs_N_grps <- factor(S3_vs_N_grps, levels = c("Normal", "Stage III"))

    S3_vs_S2_grps <- rep(NA, length(case_ids))
    S3_vs_S2_grps[is_hcc_stage_II] <- "Stage II"
    S3_vs_S2_grps[is_hcc_stage_III] <- "Stage III"
    S3_vs_S2_grps <- factor(S3_vs_S2_grps, levels = c("Stage II", "Stage III"))

    ### "----------------------------------------------------------------------"
    ### Tumor vs normal samples
    ### "----------------------------------------------------------------------"
    T_vs_N_groups <- rep("Normal", length(case_ids))
    T_vs_N_groups[is_tumor] <- "Tumor"
    T_vs_N_groups <- factor(T_vs_N_groups)

    ### Return
    list(
        case_ids = case_ids,
        hcc_stage_grps = hcc_stage_grps,
        hcc_grade_grps = hcc_grade_grps,
        T_vs_N_groups = T_vs_N_groups,
        S1_vs_N_grps = S1_vs_N_grps,
        S2_vs_N_grps = S2_vs_N_grps,
        S2_vs_S1_grps = S2_vs_S1_grps,
        S3_vs_N_grps = S3_vs_N_grps,
        S3_vs_S2_grps = S3_vs_S2_grps,
        tumor_grps = tumor_grps,
        steatosis_score = clinical_info$steatosis_score[sample_idx],
        is_steatosis = factor(
            clinical_info$steatosis_score[sample_idx] >= 5,
            levels = c(FALSE, TRUE),
            labels = c("<5%", ">=5%")
        ),
        fibrosis_score = clinical_info$metavir_fibrosis_score[sample_idx],
        is_fibrotic = factor(
            clinical_info$metavir_fibrosis_score[sample_idx] > 2,
            levels = c(FALSE, TRUE),
            labels = c("F0-2", "F3-4")
        ),
        necrosis_score = clinical_info$necrosis_score[sample_idx]*100
    )
}

### "========================================================================="
### Load clinical info
### "========================================================================="
LoadClinicalInfo <- function() {
    read_excel(
        here("data", "lcms", "clinical_info.xlsx"),
        sheet = "Sample Key"
    ) #%>%
    #filter(
    #    `Histological diagnosis` == "HCC"
    #)
}

### "========================================================================="
### Load metabolite classes
### "========================================================================="
AddClassAnnotation <- function(lcms_hit_peaks) {
    ### Load the raw MSMS data
    raw_msms_data <- LoadRawMSMSData() %>%
        dplyr::rename(
            MSMS_annotation = `Searched ID`
        ) %>%
        replace_na(
            list(
                MSMS_annotation = "Unannotated",
                short_name = "Unannotated"
            )
        ) %>%
        select(MSMS_annotation, short_name) %>%
        distinct() %>%
        arrange(MSMS_annotation)

    ### Export the MSMS names
    #readr::write_excel_csv(
    #    raw_msms_data, here("data", "lcms", "lcms_msms_names.csv")
    #)
    
    ### Load the data
    metabolite_class_data <- readr::read_csv(
            file = here(
                "data", "lcms", "lcms_hit_peaks_annotated_updated_v2.csv"
            ),
            locale = readr::locale(encoding = "latin1"),
            show_col_types = FALSE
        ) %>%
        select(short_name, Sub_Class, Class) %>%
        distinct()

    ### Check to see if there is any missing labels
    missing_raw_data <- sort(
        metabolite_class_data |>
            filter(!(short_name %in% raw_msms_data$short_name)) |>
            pull(short_name)
    )
    if (length(missing_raw_data) > 1) {
        message("Missing the following raw_data:")
        print(missing_raw_data)
        stop()
    }

    missing_annotation <- sort(
        raw_msms_data |>
            filter(!(short_name %in% metabolite_class_data$short_name)) |>
            pull(short_name)
    )
    if (length(missing_annotation) > 1) {
        message("Missing annotations for the following raw_data:")
        print(missing_annotation)
        stop()
    }

    ### Create the lookup tibble
    class_lookup_tibble <- left_join(
        raw_msms_data,
        metabolite_class_data,
        by = "short_name"
    ) %>%
    mutate(
        MSMS_annotation = NULL
    ) %>%
    rename(
        short_name = "MSMS_annotation"
    ) %>%
    distinct()

    ### Merge the data
    left_join(
        lcms_hit_peaks,
        class_lookup_tibble,
        by = "MSMS_annotation"
    )
}

### Graph plotting ============================================================
### "========================================================================="
### Plot tumor heterogeneity heatmap
### "========================================================================="
PlotTHHeatmap <- function(
    TH_mat, row_title, column_labels, legend_title, filename
) {
    col_fun <- circlize::colorRamp2(
        c(-1, 0, 1),
        c("#4dac26", "#dddddd", "#d01c8b")
    )

    TH_dist <- apply(
        TH_mat, 2, function(x) table(factor(x, levels = c(-1, 0, 1)))
    ) / nrow(TH_mat) * 100

    rownames(TH_dist) <- c("Lower", "N.S.", "Higher")
    TH_anno <- anno_barplot(
        t(TH_dist),
        bar_width = 1,
        extend = 0,
        #height = unit(0.25, "inch"),
        annotation_height = unit(0.3, "inch"),
        gp = gpar(fill = col_fun(c(-1, 0, 1)))
    )

    p <- Heatmap(
        TH_mat,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = FALSE,
        border = TRUE,
        rect_gp = gpar(lty = 0),
        col = col_fun,
        height = unit(2.5, "inch"),
        width = unit(0.7, "inch"),
        row_title = paste0(row_title, " (", nrow(TH_mat), ")"),
        column_labels = column_labels,
        column_split = c("I", "II", "III"),
        column_title = " ",
        show_heatmap_legend = TRUE,
        heatmap_legend_param = list(
            title = legend_title,
            color_bar = "discrete",
            border = "black",
            at = c(-1, 0, 1),
            labels = c(
                "Lower (Padj < 0.05)",
                "N.S. (Padj > 0.05)",
                "Higher (Padj < 0.05)"
            )
        ),
        top_annotation = HeatmapAnnotation(
            "%" = TH_anno
        ),
        use_raster = TRUE,
        raster_quality = 10
    )

    fs::dir_create(here("figures", "heatmaps"))

    ### Draw png
    png(
        here("figures", "heatmaps", paste0(filename, ".png")),
        width = 4, height = 4, units = "in", res = 300
    )
    draw(p)
    dev.off()

    pdf(
        here("figures", "heatmaps", paste0(filename, ".pdf")),
        width = 4, height = 4
    )
    draw(p)
    dev.off()
}

### "========================================================================="
### Plot the predicted ROI size vs pathological grades
### "========================================================================="
PlotPredvsPath <- function(
    plot_data, test_data, y_title, y_lim, y_breaks, comparisons,
    step_increase, file_name
) {
    stat_test <- test_data %>%
        wilcox_test(
            formula = pred_score ~ plot_grp,
            comparisons = comparisons,
            p.adjust.method = "fdr",
            alternative = "less"
        ) %>%
        add_y_position(
            fun = "max",
            step.increase = step_increase
        )

    print(stat_test)

    cur_p <- plot_data %>%
        ggplot(aes(x = plot_grp, y = pred_score)) +
        stat_summary(
            fun = median,
            geom = "bar",
            fill = "gray",
            width = 0.75
            # shape = 18,
            # size = 4
        ) + 
        # geom_point() +
        ggbeeswarm::geom_beeswarm(
            size = 1,
            cex = 6 ## point spacing
        ) +
        stat_pvalue_manual(
            stat_test, 
            label = "Padj={p.adj}, n={n1 + n2}",
            tip.length = 0.02,
            vjust = -0.15,
            size = 2.5
        ) +
        scale_y_continuous(
            name = y_title,
            limits = y_lim,
            breaks = y_breaks,
            expand = expansion(mult = c(0, 0.01))
        ) +
        coord_cartesian(clip = "off") +
        theme_classic(
            base_family = "Arial"
        ) +
        theme(
            # plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
            #legend.key.size = unit(0.5, 'cm'),
            #legend.title = element_text(size = 8, face = "bold"),
            #legend.text = element_text(size = 8),
            #legend.position = "top",
            #legend.position = c(0.8, 0.8),
            #legend.direction="horizontal",
            axis.title.x = element_blank(),
            axis.text = element_text(colour = "black"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
        )

    ### Save the plot
    for (img_type in c(".png", ".pdf")) {
        ggsave(
            filename = here(
                "figures", "msi",
                paste0("pred_vs_path_", file_name, img_type)
            ),
            plot = cur_p,
            width = 1.6,
            height = 2.7,
            units = "in",
            dpi = 300,
            device = ifelse(img_type == ".png", png, cairo_pdf)
        )
    }
}

### "========================================================================="
### Plot the correlation between predicted ROI size and pathological scores
### "========================================================================="
PlotPredvsPathCorr <- function(
    plot_data, x_title, y_title, stat_type, file_name
) {
    ### Find the sample size
    section_dist <- plot_data %>%
        group_by(section_type) %>%
        count() %>%
        dplyr::pull(n, name = section_type)
    section_labels <- paste0(names(section_dist), " (n = ", section_dist, ")")
    names(section_labels) <- names(section_dist)

    cur_p <- plot_data %>%
        ggscatter(
            x = "path_score", y = "pred_score",
            add = "reg.line",
            #conf.int = TRUE,
            add.params = list(color = "blue", fill = "gray"),
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(
                method = stat_type,
                cor.coef.name = ifelse(stat_type == "kendall", "tau", "R"),
                alternative = "greater",
                label.x.npc = 0.5,
                hjust = 0.5,
                size = 3
            )
        ) +
        #ggplot(aes(x = steatosis_score, y = tumor_percent)) +
        #geom_smooth(method = "lm", formula = y ~ x) +
        #geom_point() +
        facet_wrap(
            section_type ~ .,
            ncol = 1,
            scales = "free",
            labeller = as_labeller(section_labels),
            strip.position = "top"
        ) +
        scale_x_continuous(
            name = x_title
        ) +
        scale_y_continuous(
            name = y_title,
            #expand = expansion(mult = c(0.1, 0.1)),
            limits = c(0, 100),
            breaks = c(0, 20, 40, 60, 80, 100)
            # limits = y_lim
        ) +
        theme_classic(
            base_family = "Arial"
        ) +
        theme(
            # plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
            #legend.key.size = unit(0.5, 'cm'),
            #legend.title = element_text(size = 8, face = "bold"),
            #legend.text = element_text(size = 8),
            #legend.position = "top",
            #legend.position = c(0.8, 0.8),
            #legend.direction="horizontal",
            axis.text = element_text(colour = "black"),
            #axis.line = element_line(),
            strip.text = element_text(face = "bold"),
            strip.background = element_blank()
            #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        )

    ### Save the plot
    for (img_type in c(".png", ".pdf")) {
        ggsave(
            filename = here(
                "figures", "msi",
                paste0("pred_vs_path_corr_", file_name, img_type)
            ),
            plot = cur_p,
            width = 2,
            height = 3,
            units = "in",
            dpi = 300,
            device = ifelse(img_type == ".png", png, cairo_pdf)
        )
    }
}

### "========================================================================="
### Plot the composition of prediction
### "========================================================================="
PlotPredDist2 <- function(
    pred_dist, y_title, annot_colors, width, height, file_name
) {
    cur_p <- pred_dist %>%
        ggplot(
            aes(x = section_id, y = percent, fill = pred_label)
        ) +
        # geom_bar(stat="identity", position="fill",  colour="black") +
        geom_col(
            width = 0.8,
            color = "black"
        ) +
        xlab("") +
        scale_y_continuous(
            name = y_title, #"Section area (%)",
            # limits = y_lim,
            expand = expansion(mult = c(0, .01))
        ) +
        scale_fill_manual(
            name = "Predicted MERs",
            values = annot_colors,
        ) +
        facet_grid(. ~ case_id, scales = "free", space = "free") +
        # facet_wrap(~ case_id, scales = 'free', space = "free", nrow = 2) +
        theme_classic() +
        theme(
            legend.position = "top",
            legend.key.size = unit(0.5, 'cm'),
            legend.title = element_text(size = 8, face = "bold"),
            legend.text = element_text(size = 8),
            axis.text = element_text(colour = "black"),
            strip.text = element_text(face = "bold", margin = margin(0,0,2,0)),
            strip.background = element_blank()
        )
    
    ### save the plots
    SavePlot(
        cur_p, width = width, height = height,
        dir_name = "msi", file_name = file_name
    )    
}

PlotPredDist <- function(
    pred_dist, y_title, y_lim, annot_colors,
    width, height, file_name
) {
    section_dist <- pred_dist %>%
        group_by(section_type) %>%
        select(section_type, slide_id) %>%
        distinct() %>%
        count() %>%
        dplyr::pull(n, name = section_type)

    section_labels <- paste0(names(section_dist), "(n=", section_dist, ")")
    names(section_labels) <- names(section_dist)

    cur_p <- pred_dist %>%
        ggplot(
            aes(x = case_id, y = tumor_percent, fill = pred_label)
        ) +
        geom_col(
            width = 0.8,
            color = "black"
        ) +
        scale_y_continuous(
            name = y_title,
            limits = y_lim,
            expand = expansion(mult = c(0, .01))
        ) +
        facet_wrap(
            section_type ~ .,
            ncol = 1,
            #scales = "fixed"
            labeller = as_labeller(section_labels),
            strip.position = "top"
        ) +
        scale_fill_manual(
            name = "Histopathology",
            values = annot_colors,
        ) +
        geom_hline(aes(yintercept=-Inf)) +
        coord_cartesian(clip = "off") +
        #xlab("Case ID") +
        theme_classic() +
        theme(
            legend.key.size = unit(0.5, 'cm'),
            legend.title = element_text(size = 8, face = "bold"),
            legend.text = element_text(size = 8),
            legend.position = "none", #"top",
            #legend.position = c(0.8, 0.8),
            #legend.direction="horizontal",
            axis.text = element_text(colour = "black"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            axis.title.x = element_blank(),
            strip.text = element_text(face = "bold", margin = margin(0,0,2,0)),
            strip.background = element_blank()
            #panel.border = element_rect(colour = "black", fill = NA)
        )

    ### save the plots
    SavePlot(
        cur_p, width = width, height = height,
        dir_name = "msi", file_name = file_name
    )
}

PrepareSampleMarking <- function(
    msi_mat, show_names, is_peak_id = FALSE
) {
    ### Load the reference peaks
    load(here("data", "lcms", "lcms_peaks.Rdata"))

    ### Get the peak_ids
    show_peak_ids <- names(show_names)
    
    ### Sample marking
    if (!is_peak_id) {
        ### Show the annotation
        show_peak_labels <- lcms_hit_peaks$peak_id[
            match(show_names, lcms_hit_peaks$MSMS_annotation)
        ]
    } else {
        ### Show peak_id
        show_peak_labels <- show_names
    }

    show_row_idx <- match(show_peak_ids, rownames(msi_mat))

    sample_marking <- rowAnnotation(
        selected_samples = anno_mark(
            at = show_row_idx,
            labels = show_peak_labels
        )
    )

    sample_marking
}

PrepareMetaboliteAnnotation <- function(
    msi_mat, metabolite_classes_colors, ncol = 2
) {
    ### Load the reference peaks
    load(here("data", "lcms", "lcms_peaks.Rdata"))

    ### Add class information
    lcms_hit_peaks <- AddClassAnnotation(lcms_hit_peaks)

    ### Metabolite annotations
    idx <- match(rownames(msi_mat), lcms_hit_peaks$peak_id)
    metabolite_classes <- lcms_hit_peaks$Class[idx]
    metabolite_classes <- factor(metabolite_classes)

    metabolite_dist <- table(metabolite_classes)
    metabolite_classes_labels <- paste0(
        names(metabolite_dist), " (", metabolite_dist, ")"
    )
    names(metabolite_classes_labels) <- names(metabolite_dist)

    ### No need to reorder
    # metabolite_classes_labels <- metabolite_classes_labels[
    #    names(metabolite_classes_colors)
    # ]
    metabolite_classes_colors_cur <- metabolite_classes_colors[
        names(metabolite_classes_labels)
    ]

    metabolite_classes_legend <- Legend(
        title = "Metabolite types",
        title_position = "topcenter",
        border = "black",
        gap = unit(0.1, "inch"),
        legend_gp = gpar(fill = metabolite_classes_colors_cur),
        labels = metabolite_classes_labels,
        column_gap = unit(0.15, "inch"),
        row_gap = unit(0.1, "inch"),
        ncol = ncol
    )

    ### Metabolite class
    metabolite_annotation <- rowAnnotation(
        metabolite_classes = anno_simple(
            metabolite_classes,
            col = metabolite_classes_colors_cur,
            border = TRUE,
            gp = gpar(col = "black"),
            na_col = "white"
        ),
        show_annotation_name = FALSE
    )

    list(
        metabolite_annotation = metabolite_annotation,
        metabolite_classes_legend = metabolite_classes_legend
    )
}

### "-------------------------------------------------------------------------"
### Load the top markers for each classifiers
### "-------------------------------------------------------------------------"
LoadROITopMarkers <- function() {

    ### Load study data
    study_data <- LoadStudyMSIData()

    ### Load all the models
    msi_multi <- readRDS(here("results", "plsda", "plsda_multi.rds"))
    msi_fibrotic <- readRDS(here("results", "plsda", "plsda_fibrotic.rds"))
    msi_steatotic <- readRDS(here("results", "plsda", "plsda_steatotic.rds"))

    ### Get the list of all markers
    model_names <- c("multi", "fibrotic", "steatotic")

    all_markers <- NULL
    for (model_name in model_names) {
        ### Get the vip of the model
        cur_model <- get(paste0("msi_", model_name))
        vips <- ropls::getVipVn(cur_model$model)

        if (nlevels(cur_model$model@suppLs$y) > 2) {
            histopath_types <- colnames(cur_model$model@coefficientMN)
        } else {
            histopath_types <- levels(cur_model$model@suppLs$y)[2]
        }

        for (histopath_type in histopath_types) {
            ### metabolie_num x class_num
            if (nlevels(cur_model$model@suppLs$y) > 2) {
                coeffs <- cur_model$model@coefficientMN[, histopath_type]
            } else {
                coeffs <- cur_model$model@coefficientMN[, 1]
            }

            all_markers <- bind_rows(
                all_markers,
                study_data %>%
                    mutate(
                        vip = vips[peak_id],
                        coeff = coeffs[peak_id],
                        model_name = !!model_name,
                        histopath_type = !!histopath_type
                    ) %>%
                    # filter(vip >= 1.0) %>%
                    arrange(-vip) %>%
                    select(
                        model_name, histopath_type, peak_id, MSMS_annotations,
                        vip, coeff
                    )
            )
        }
    }

    all_markers
}

### "========================================================================="
### Train a regression model for different region and plot the regression plot
### "========================================================================="
TestSgMEClinicalModel <- function(
    lcms_mat, ME_region_stats, 
    model_name, clinical_phenotype, ME_region_name,
    test_xlabel, train_xlabel, is_log10,
    train_dim, test_dim,
    alternative,
    test_type = "wilcox_test",
    is_hit_peaks_only = FALSE,
    train_xlim = NULL, train_ylim = NULL,
    test_xlim = NULL, test_ylim = NULL,
    predI = 2, cat_step_size = 0.04, rnd_seed = 1
) {
    set.seed(rnd_seed)

    ### Only hit peaks?
    if (is_hit_peaks_only) {
        ### We want to focus on the most abundant peaks
        load(here("data", "lcms", "lcms_peaks.Rdata"))
        lcms_mat <- lcms_mat[, lcms_hit_peaks$mz_id]
    }
    
    ### Only subset ME-regions to those with LCMS data
    ME_region_stats <- ME_region_stats %>%
        filter(sample_id %in% rownames(lcms_mat))

    ### Convert the ME_region_name
    ME_region_name_disp <- if (ME_region_name == "G3") {
        "high-grade"
    } else if (ME_region_name == "G2") {
        "low-grade"
    } else if (ME_region_name == "fibro") {
        "fibrotic"
    } else {
        ME_region_name
    }

    ylabel <- paste0("Predicted tissue coverage\nbased on LC-MS ")
    if (is_log10) {
        ylabel <- paste0(ylabel, "(log10[%])")
    } else {
        ylabel <- paste0(ylabel, "(%)")
    }

    #is_log10 <- TRUE

    ### Get the ME_region_areas for the training samples
    ME_region_areas <- ME_region_stats[[paste0(ME_region_name, "_area")]]
    ME_region_areas <- if_else(ME_region_areas == 0, 1, ME_region_areas)

    ### Get the training data
    lcms_train_mat <- lcms_mat[ME_region_stats$sample_id, ]

    train_values <- if (is_log10) {
        log10(ME_region_areas / ME_region_stats$total_area * 100)
    } else {
        ME_region_areas / ME_region_stats$total_area * 100
    }

    ### Train the model
    if (model_name == "glmnet") {
        trained_model <- cv.glmnet(
            lcms_train_mat,
            ### We must do log10 because the predictions are about low level classifications
            train_values,
            # family = "binomial",
            # type.measure = "mse",
            # alpha = 0.7,
            keep = TRUE
        )

        train_preds <- predict(
            trained_model,
            lcms_train_mat #, s = "lambda.min"
        ) # lambda.1se or lambda.min
    
        ### Convert from a matrix to numeric
        train_preds <- train_preds[, 1]
        #assess.glmnet(cfit$fit.preval, newy = y, family = "binomial")
        #a <- ME_hit_region_stats %>% select(total_area, normal_area) %>% as.matrix()
        #assess.glmnet(cv_fit, newx = lcms_mat, newy = a)


    } else if (model_name == "plsda") {

        trained_model <- ropls::opls(
            lcms_train_mat,    
            train_values,
            crossvalI = 5, predI = predI
        )
        train_preds <- predict(
            trained_model, lcms_train_mat
        )

    } else if (model_name == "svm") {
        trained_model <- e1071::svm(
            x = lcms_train_mat,
            y = train_values,
            cost = 1,
            gamma = 1,
            scale = TRUE
        )
        train_preds <- predict(
            trained_model, lcms_train_mat
        )

    } else {
        stop("Unknown model provided")
    }

    ### Compute the RMSE
    train_tibble <- tibble(
        train_values = train_values,
        train_preds  = train_preds
    )
    
    train_lm <- lm(train_preds ~ train_values, train_tibble)

    lm_n <- sum(!is.na(train_values))
    lm_RMSE <- sqrt(mean((train_preds - train_values)^2))
    lm_R2 <- summary(train_lm)$r.squared
    lm_pval <- summary(train_lm)$coefficients["train_values", "Pr(>|t|)"]
    # stat_text <- bquote(atop(textstyle(R^2 == .(signif(lm_R2, 3)) * ", " ~ P == .(signif(lm_pval, 3))), textstyle(RMSE == .(signif(lm_RMSE, 3)))))
    stat_text <- bquote(atop(NA, atop(textstyle(R^2 == .(signif(lm_R2, 3)) * ", " ~ P == .(signif(lm_pval, 3))), textstyle(RMSE == .(signif(lm_RMSE, 3)) * ", " ~ s == .(lm_n)))))
    # stat_text <- c(
    #    as.expression(bquote(R^2 == .(signif(lm_R2, 3)) * ", " ~ P == .(signif(lm_pval, 3)))),
    #    as.expression(bquote(RMSE == .(signif(lm_RMSE, 3))))
    # )

    ### Prepare the plot
    train_p <- train_tibble %>%
        ggplot(aes(x = train_values, y = train_preds)) +
        labs(
            title = paste0("ME-", ME_region_name_disp, " regions"),
            # subtitle = stat_text,
            x = train_xlabel,
            y = ylabel
        ) +
        geom_point() +
        coord_fixed() +
        stat_smooth(
            method = "lm",
            formula = y ~ x,
            geom = "smooth"
        ) +
        annotate(
            "text",
            x = -Inf, y = Inf,
            hjust = -0.05, vjust = 0.8,
            label = stat_text, size = 3.3
        ) +
        theme_classic() +
        theme(
            axis.text = element_text(colour = "black"),
            plot.title = element_text(
                hjust = 0.5,
                margin=margin(0, 0, 0, 0)
            ),
            plot.subtitle = element_text(                
                hjust = 0.5,
                lineheight = 0
            )
        )
    
    if (!is.null(train_xlim)) {
        train_p <- train_p + xlim(train_xlim[1], train_xlim[2])
    }
    if (!is.null(train_ylim)) {
        train_p <- train_p + ylim(train_ylim[1], train_ylim[2])
    }

    ### Apply the model to all the samples
    all_preds <- predict(trained_model, lcms_mat) # , s = "lambda.min")
    if (is.matrix(all_preds)) {
        all_preds <- all_preds[, 1]
    }

    ### Plot it compare to clinical phenotypes
    plot_data <- as_tibble(all_preds, rownames = "sample_id") %>%
        separate(sample_id, c("case_id", "section_id")) %>%
        left_join(
            clinical_info %>% select(
                case_id, tumor_stage_AJCC_V8, edmonson_grade,
                metavir_fibrosis_score, steatosis_score, necrosis_score
            )
        ) %>%
        mutate(
            is_tumor = startsWith(section_id, "T"),

            ### Re-arrange some of the columns
            edmonson_grade = case_when(
                edmonson_grade == 4 ~ 3,
                !is_tumor ~ 0,
                TRUE ~ edmonson_grade
            ),
            edmonson_grade = factor(
                edmonson_grade
            ),
            metavir_fibrosis_score = factor(
                metavir_fibrosis_score,
                levels = 0:3, labels = paste0("F", 0:3)
            ),
            necrosis_score = necrosis_score * 100
            #necrosis_score = case_when(
            #    necrosis_score == 0 ~ 1,
            #    TRUE ~ necrosis_score
            #),
            #steatosis_score = case_when(
            #    steatosis_score == 0 ~ 1,
            #    TRUE ~ steatosis_score
            #)
        ) 

    plot_data$y <- plot_data[["value"]]
    plot_data$x <- plot_data[[clinical_phenotype]]

    if (is.factor(plot_data[[clinical_phenotype]])) {
        ### "-----------------------------------------------------------------"
        ### For categorical data
        ### "-----------------------------------------------------------------"
        plot_data <- plot_data %>%
            filter(!is.na(x))
        
        stat_test <- if (test_type == "t_test") {
            plot_data %>%
                t_test(
                    y ~ x,
                    alternative = alternative, p.adjust.method = "BH"
                ) %>%
                add_y_position()

        } else if (test_type == "wilcox_test") {
            plot_data %>%
                wilcox_test(
                    y ~ x,
                    alternative = alternative, p.adjust.method = "BH"
                ) %>%
                add_y_position()

        } else {
            stop("Unknown test")
        }

        print(stat_test)

        median_range <- function(x) {            
            qu <- quantile(x, probs = c(0.25, 0.5, .75), na.rm=T)
            data.frame(y=qu[2], ymin=qu[1], ymax=qu[3])
        }

        test_p <- plot_data %>%
            ggplot(aes(x = .data[[clinical_phenotype]], y = value, group = 1)) +
            ggbeeswarm::geom_quasirandom(
                col = "darkgray"
            )

        median_cl_boot <- function(x, conf = 0.95) {
            lconf <- (1 - conf)/2
            uconf <- 1 - lconf
            require(boot)
            bmedian <- function(x, ind) median(x[ind])
            bt <- boot(x, bmedian, 1000)
            bb <- boot.ci(bt, type = "perc")
            data.frame(y = median(x), ymin = quantile(bt$t, lconf), ymax = quantile(bt$t, 
                uconf))
        }

        if (test_type == "t_test") {
            test_p <- test_p +
                stat_summary(
                    #fun.data = median_q1q3, # mean_cl_boot,
                    fun.data = mean_cl_boot, #mean_se,
                    geom = "errorbar", width = 0.25
                ) +
                stat_summary(fun = mean, geom = "line") +
                stat_summary(fun = mean, geom = "point", size = 2)
                #stat_summary(fun = median, geom = "line") +
                #stat_summary(fun = median, geom = "point", size = 2) +

        } else {
            test_p <- test_p +
                stat_summary(
                    fun.data = median_cl_boot,
                    geom = "errorbar", width = 0.25
                ) +                
                stat_summary(fun = median, geom = "line") +
                stat_summary(fun = median, geom = "point", size = 2)
        }
        
        test_p <- test_p +
            stat_pvalue_manual(
                stat_test,
                step.increase = cat_step_size
                # hide.ns = TRUE
            ) +
            labs(
                title = paste0("ME-", ME_region_name_disp, " regions"),
                subtitle = paste0("s = ", sum(!is.na(plot_data$x))),
                x = test_xlabel,
                y = ylabel
            ) +
            theme_classic() +
            theme(
                axis.text = element_text(colour = "black"),
                plot.title = element_text(hjust = 0.5),
                plot.subtitle = element_text(hjust = 0.5)
            )
    
    } else {
        ### "-----------------------------------------------------------------"
        ### For numerical data
        ### "-----------------------------------------------------------------"
        test_p <- plot_data %>%
            ggplot(aes(x = x, y = y)) +
            geom_point(color = "darkgray") +
            # coord_fixed() +
            stat_cor(
                method = "pearson",
                alternative = alternative
            ) +
            stat_smooth(
                method = "lm",
                formula = y ~ x,
                geom = "smooth",
                colour="black"
            ) +
            labs(
                title = paste0("ME-", ME_region_name_disp, " regions"),
                subtitle = paste0("s = ", sum(!is.na(plot_data$x))),
                x = test_xlabel,
                y = ylabel
            ) +
            theme_classic() +
            theme(
                axis.text = element_text(colour = "black"),
                plot.title = element_text(hjust = 0.5),
                plot.subtitle = element_text(hjust = 0.5)
            ) 
    }

    if (!is.null(test_xlim)) {
        test_p <- test_p + xlim(test_xlim[1], test_xlim[2])
    }
    if (!is.null(test_ylim)) {
        test_p <- test_p + ylim(test_ylim[1], test_ylim[2])
    }

    ### For the filename
    filename <- paste0(
        ME_region_name, "_", clinical_phenotype, "_", model_name, "_",
        ifelse(is_hit_peaks_only, "hit", "raw")
    )


    ### Save the results
    fs::dir_create(here("results", "SgMERdeconv"))
    saveRDS(
        all_preds,
        file = here("results", "SgMERdeconv", paste0(filename, ".RDS"))
    )

    ### save the plots
    SavePlot(
        train_p,
        width = train_dim[1],
        height = train_dim[2],
        dir_name = "clinical",
        file_name = paste0(filename, "_train")
    )

    SavePlot(
        test_p,
        width = test_dim[1],
        height = test_dim[2],
        dir_name = "clinical",
        file_name = paste0(filename, "_test")
    )
}

### "-------------------------------------------------------------------------"
### Plot the transformation features
### http://dx.doi.org/10.1016/j.jhydrol.2012.05.035
### "-------------------------------------------------------------------------"
PlotVIPvsCoeff <- function(
    model_name, histopath_type, up_col, down_col,
    x_label, x_lim, x_breaks,
    y_lim, y_breaks,
    is_show_annotation = TRUE,
    top_num_show, file_name
) {
    control_col <- "white"

    plot_data <- all_markers %>%
        filter(
            model_name == !!model_name,
            histopath_type == !!histopath_type
        ) %>%
        ### remove all peaks without VIP because they were not used to
        ### build the model
        filter(!is.na(vip)) %>%
        arrange(-vip) %>%
        mutate(vip_order_id = row_number()) %>%
        arrange(-coeff) %>%
        mutate(coeff_dec_id = row_number()) %>%
        arrange(coeff) %>%
        mutate(
            coeff_inc_id = row_number(),
            importance = factor(
                case_when(
                    #vip < vip_thres ~ "Unimportant",
                    (coeff >= coeff_thres_study) & (vip >= vip_thres_study) ~ "Increasing",
                    (coeff <= -coeff_thres_study) & (vip >= vip_thres_study) ~ "Decreasing",
                    TRUE ~ "NS"
                ),
                levels = c("Increasing", "Decreasing", "NS")
            ),
            is_label =
                (
                    vip_order_id <= top_num_show |
                    coeff_dec_id <= 2 |
                    coeff_inc_id <= 2
                ) &
                (vip >= vip_thres_study) & 
                (abs(coeff) >= coeff_thres_study),
            MSMS_annotations = case_when(
                is_label & is_show_annotation ~
                        paste0(peak_id, " - ", MSMS_annotations),
                is_label & !is_show_annotation ~ peak_id,
                TRUE ~ NA
            ),
            MSMS_annotations = str_wrap(
                MSMS_annotations, 15,
                whitespace_only = FALSE
            )
        )

    #if (is_Pval) {
    #    vip_vs_fc_p <- plot_data %>%
    #        ggplot(aes(x = uni_log2FC, y = -log10(uni_padj)))

    #    y_label <- "Padj (-log10)"

    #} else {
        vip_vs_fc_p <- plot_data %>%
            ggplot(aes(x = coeff, y = vip))

        y_label <- "VIP"
    # }

    vip_vs_fc_p <- vip_vs_fc_p +
        geom_vline(xintercept = 0, linetype = "solid", color = "gray") +
        geom_hline(yintercept = vip_thres_study, linetype = "dotted", color = "gray") +
        geom_vline(xintercept = coeff_thres_study, linetype = "dotted", color = "gray") +
        geom_vline(xintercept = -coeff_thres_study, linetype = "dotted", color = "gray") +
        # geom_hline(yintercept = 1, linetype = "dashed", color = "darkgray") +
        geom_point(
            data = ~ subset(., importance == "NS"),
            aes(fill = importance),
            shape = 21
        ) +
        geom_point(
            data = ~ subset(., importance != "NS"),
            aes(fill = importance),
            shape = 21
        ) +
        geom_text_repel(
            aes(label = MSMS_annotations),
            lineheight = 0.8,
            col = "black",
            box.padding = 0.5,
            force = 3,
            min.segment.length = 0,
            max.overlaps = Inf,
            nudge_x = 0.02,
            nudge_y = 0.01,
            size = 3,
            segment.size = 0.2,
            show.legend = FALSE
        ) +
        scale_x_continuous(
            name = x_label,
            breaks = x_breaks,
            limits = x_lim
        ) +
        scale_y_continuous(
            name = y_label,
            breaks = y_breaks,
            limits = y_lim
        ) +
        scale_fill_manual(
            name = "Projection",
            labels = c(
                Increasing = "Positive",
                Decreasing = "Negative",
                NS = "Unimportant"
            ),
            values = c(
                Increasing = up_col, 
                Decreasing = down_col,
                NS = control_col
            )
        ) +
        theme_classic() +
        theme(
            # plot.title = element_text(hjust = 0.5, size = 11),
            legend.title = element_text(hjust = 0.5, size = 9),
            legend.position = c(0.2, 0.15),
            legend.text = element_text(size = 8),
            legend.spacing.y = unit(0.01, 'inch'),
            legend.key.size = unit(0.15, "inch"),
            legend.background = element_blank(),
            axis.text = element_text(colour = "black")
        )

    ### save the plots
    SavePlot(
        vip_vs_fc_p,
        width = 2.8, height = 2.6,
        dir_name = "msi",
        file_name = file_name
    )
}
