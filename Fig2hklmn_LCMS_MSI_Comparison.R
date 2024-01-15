### "========================================================================="
### Plot LC-MS and MSI comparisons
###
### Copyright (c) 2021-2024. Bioinformatics Institute, A*STAR, Singapore
### License: MIT License. See LICENSE for the full license.
### Authors: Lit-Hsin Loo (loolh@bii.a-star.edu.sg)
### "========================================================================="

library(here)
source(here("lib", "MSI_lib.R"))
source(here("lib", "Study_lib.R"))
library(ggpmisc)

### "========================================================================"
### Common things that need to be run
### "========================================================================"
### Load the LCMS hits
load(here("data", "lcms", "lcms_peaks.Rdata"))

sample_list <- read_excel(
        here("data", "lcms", "Supp_Tables_S1 to S4_v3.xlsx"),
        sheet = "Table_S1b Sample_Usages"
    ) %>%
    slice_head(n=117) %>%
    filter(!is.na(`DESI-MSI`)) %>%
    select(`Case ID`, `Section ID`) %>%
    rename(`Case ID` = "case_id", `Section ID` = "section_id")

### "========================================================================"
### Fig. 2h. Distribution of the DESI-MSI samples
### "========================================================================"
### Get clinical info
clinical_info <- LoadClinicalInfo()

sample_list <- sample_list |>
    mutate(
        section_type = if_else(
            startsWith(section_id, "T"), "T", "N"
        ),
        .after = "section_id"
    ) |>
    group_by(case_id, section_type) |>
    summarize(
        section_num = n(),
        .groups = "drop"
    )
    
sample_dist <- left_join(
        sample_list,
        clinical_info %>%
            select(
                case_id,
                histological_diagnosis, tumor_stage_AJCC_V8,
                edmonson_grade
            ),
        by = c("case_id")
    ) %>%
    mutate(
        sample_type = factor(
            case_when(
                section_type == "T" & histological_diagnosis == "HCC" ~ "HCC tumor",
                section_type == "T" & histological_diagnosis == "CCA" ~ "iCCA tumor",
                TRUE ~ "Normal"
            ),
            levels = c("HCC tumor", "iCCA tumor", "Normal")
        ),
        tumor_stage = factor(
            case_when(
                tumor_stage_AJCC_V8 == "TNM Stage IB" ~ "Stage IB",
                tumor_stage_AJCC_V8 == "TNM Stage II" ~ "Stage II",
                tumor_stage_AJCC_V8 == "TNM Stage IIIA" ~ "Stage III",
                tumor_stage_AJCC_V8 == "TNM Stage IIIB" ~ "Stage III",
                TRUE ~ "Unknown"
            ),
            levels = c("Stage IB", "Stage II", "Stage III")
        ),
        tumor_grade = factor(
            case_when(
                edmonson_grade == 2 ~ "G2",
                edmonson_grade == 3 ~ "G3",
                edmonson_grade == 4 ~ "G4",
                TRUE ~ "Unknown"
            ),
            levels = c("G2", "G3", "G4")
        )
        #case_id = ifelse(hbv_status == "B", paste0("*", case_id), case_id)
    ) %>%
    arrange(
        tumor_stage, case_id
    ) %>%
    rename(section_num = "sample_num")

#ordered_case_ids <- unique(sample_dist$case_id)
ordered_case_ids <- sample_dist %>%
    group_by(tumor_stage, case_id) %>%
    summarize(n = sum(sample_num)) %>%
    arrange(tumor_stage, n) %>%
    pull(case_id)

sample_dist$case_id <- factor(sample_dist$case_id, levels = ordered_case_ids)

sample_dist_plot <- sample_dist %>%
    ggplot(aes(x = case_id, y = sample_num, fill = sample_type)) +
    geom_col(colour = "black", width = 0.75) +
    scale_fill_manual(
        name = "Tissue type",
        values = c(
            roi_annot_colors[["G3"]],
            "#4daf4a",
            roi_annot_colors[["Normal"]]
        ),
        guide = guide_legend(
            direction = "horizontal",
            title.position = "top"
        )
    ) +
    ylab("Section number") +
    xlab("Case ID") +
    theme_classic() +
    scale_y_continuous(expand = c(0, 0), breaks=c(0, 2, 4, 6)) +
    theme(
        # plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_text(hjust = 0.5, size = 8, face = "bold"),
        legend.text = element_text(size = 8),
        legend.position = "top",
        legend.box = "horizontal",
        #legend.position = c(0.8, 0.8),
        #legend.direction="horizontal",
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
    )

### save the plots
fs::dir_create(here("figures", "msi"))

for (img_type in c(".png", ".pdf")) {
    ggsave(
        filename = here(
            "figures", "msi",
            paste0("msi_sample_dist", img_type)
        ),
        plot = sample_dist_plot,
        width = 2,
        height = 2,
        units = "in",
        dpi = 300
    )
}

### "========================================================================"
### Fig. 2k Plot the overlap with lcms and MSI
### "========================================================================"
PlotMSOverlap(peak_summary_mode_study)

### "========================================================================"
### Fig. 2l Plot the correlation of all the peaks
### "========================================================================"
for (ion_type in c("pos", "neg")) {
    x_lim <- if(ion_type == "pos") {c(0, 30)} else {c(15, 30)}
    is_legend <- ifelse(ion_type == "pos", TRUE, FALSE)

    PlotMSIvsLCMS(
        ion_type = ion_type,
        peak_summary_mode = peak_summary_mode_study,
        x_lim = x_lim,
        is_legend = is_legend,
        is_normalized = TRUE,
        slide_ids = test_slide_ids
    )    
}


### "========================================================================"
### LCMS and MSI correlation
### "========================================================================"
annotations <- lcms_hit_peaks$MSMS_annotation
names(annotations) <- lcms_hit_peaks$peak_id

min_sample_num <- 4

mean_stats <- bind_rows(
    LoadMSILCMSTMeanTibble(
        ion_type = "pos",
        peak_summary_mode = peak_summary_mode_study,
        is_normalized = is_normalized_study,
        slide_ids = test_slide_ids
    ),
    LoadMSILCMSTMeanTibble(
        ion_type = "neg",
        peak_summary_mode = peak_summary_mode_study,
        is_normalized = is_normalized_study,
        slide_ids = test_slide_ids
    )
)

test_res <- mean_stats %>%
    unite("sample_id", c("case_id", "section_type")) %>%
    select(peak_id, ave_msi_mean, ave_lcms_mean) %>%
    group_by(peak_id) %>%
    mutate(sample_num = n()) %>%
    filter(sample_num > min_sample_num) %>%
    cor_test(
        ave_msi_mean, ave_lcms_mean,
        alternative = "greater", method = "spearman"
    ) %>%
    adjust_pvalue(method = "fdr") %>%
    mutate(
        p.adj.signif = case_when(
            p.adj > 0.10 ~ "ns",
            p.adj > 0.05 ~ "*",
            TRUE ~ "**"
        ),
        # MSMS_annotation = annotations[peak_id],
        MSMS_annotation = if_else(
            p.adj < 0.1, paste0(peak_id, " - ", annotations[peak_id]), NA
        )
    ) %>%
    arrange(cor)

### "========================================================================"
### Fig. 2m and n T_vs_N fold change for selected metabolites
### "========================================================================"
PlotMetaboliteLCMSvsMSI <- function(
    peak_name, mean_stats, test_res, xlim, ylim
) {
    ### Get the peak_id and metabolite name
    plot_data <- mean_stats %>%
        filter(MSMS_annotation == !!peak_name)

    peak_id <- plot_data$peak_id[1]
    MSMS_annotation <- plot_data$MSMS_annotation[1]
    mz_val <- plot_data$mz[1]
    rho <- test_res %>% filter(peak_id == !!peak_id) %>% dplyr::pull(cor)
    p_val <- test_res %>% filter(peak_id == !!peak_id) %>% dplyr::pull(p)

    cur_plot <- plot_data %>%
        ggplot(aes(x = ave_lcms_mean, y = ave_msi_mean)) +
        stat_poly_line(col = "black", lty = "dashed", lwd = 0.75) +
        stat_correlation(
            mapping = use_label(c("R", "P", "s")),
            method = "spearman", alternative = "greater",
            label.y = 1
        ) +
        stat_poly_eq(
        ) +
        geom_point(aes(col = section_type)) +
        geom_text_repel(
            aes(label = case_id, col = section_type),
            size = 3,
            show.legend = FALSE
        ) +
        scale_color_manual(
            name = "Tissue type",
            values = c(
                T = roi_annot_colors[["G3"]], N = "black"
            ),
            labels = c(T = "HCC tumor", N = "Normal")
        ) +
        scale_y_continuous(
            name = "Mean normalized DESI-MSI\nabundance level (log2)",
            limits = ylim
        ) +
        scale_x_continuous(
            name = "Mean normalized LC-MS\nabundance level (log2)",
            limits = xlim
        ) +
        ggtitle(
            paste0(peak_id, " - ", MSMS_annotation) # "\n(m/z = ", mz_val, ")")
        ) +
        theme_classic() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 11),
            legend.position = "bottom", #c(0.8, 0.15),
            #legend.key.size = unit(0.25, "inch"),
            axis.text = element_text(colour = "black")
        )

    ### save the plots
    SavePlot(
        cur_plot,
        width = 3,
        height = 4,
        dir_name = "lcms",
        file_name = paste0("corr_lcms_msi_", peak_id)
    )
}

### Plot some examples
### Fig. 2m
PlotMetaboliteLCMSvsMSI(
    "SAPC", mean_stats, test_res,
    #xlim = NULL, ylim = NULL,
    xlim = c(21.8, 26), ylim = c(-2, 2)
)

### Fig. 2n
### N14 Cholesterol Sulfate
PlotMetaboliteLCMSvsMSI(
    "Cholesterol sulfate", mean_stats, test_res,
    xlim = NULL, ylim = NULL
)

