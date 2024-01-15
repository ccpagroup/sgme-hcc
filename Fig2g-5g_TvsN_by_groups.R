### "========================================================================="
### Plot tissue-averaged tumor vs normal comparisons
### Note: this can only be run after the PM group has been saved in
###       Fig5def_SgME_PM_dist_heatmap.R
### Copyright (c) 2021-2024. Bioinformatics Institute, A*STAR, Singapore
### License: MIT License. See LICENSE for the full license.
### Authors: Lit-Hsin Loo (loolh@bii.a-star.edu.sg)
### "========================================================================="
library(here)
source(here("lib", "MSI_lib.R"))
source(here("lib", "Study_lib.R"))

### Load the data
load(here("data", "lcms", "lcms_peaks.Rdata"))

mean_MSI_data <- bind_rows(
    LoadMSIData(
        ion_type = "pos",
        peak_type = "ref_peaks",
        peak_summary_mode = peak_summary_mode_study,
        stat_type = "tissue",
        msi_conf = msi_conf
    ),
    LoadMSIData(
        ion_type = "neg",
        peak_type = "ref_peaks",
        peak_summary_mode = peak_summary_mode_study,
        stat_type = "tissue",
        msi_conf = msi_conf
    )
) |> mutate(
    section_type = factor(substring(section_id, 1, 1)),
    .after = "section_id"
) |> filter(
    slide_id %in% test_slide_ids,
    !is.nan(msi_norm_mean)
)

### Must have at least two samples per group
good_ids <- mean_MSI_data |>
    group_by(peak_id) |>
    count(section_type, .drop = FALSE) |>    
    summarize(n_min = min(n), .groups = "drop") |> 
    filter(n_min > 1) |>
    pull(peak_id)

test_res <- mean_MSI_data |>
    filter(peak_id %in% good_ids) |>
    group_by(peak_id) |>
    t_test(formula = msi_norm_mean ~ section_type, ref.group = "N", detailed = TRUE, p.adjust.method="none") %>%
    adjust_pvalue(method = "BH") %>%
    mutate(log2FC = estimate2 - estimate1) %>%
    arrange(p)

PM_groups <- readRDS(here("results", "PM_groups.rds"))

label_peak_annotations <- c(
    ### Some special case
    P135 = "SAPC", ### Transformed 
    P120 = "PC(20:4/16:0)", ### Low grade
    N71 = "PI(20:2/17:1)", ### High grade    
    N14 = "Cholesterol sulfate",    
    P73 = "DG(16:0/0:0/18:1)",    
    P58 = "LysoPC(18:1/0:0)", ### Normal
    P29 = "Laudanosine", ### Normal
    P37 = "Oleoylcarnitine", ### Necrotic/fibrotic

    P70 = "DG(16:0/0:0/18:3)", ### ME-low-grade+/ME-necrotic- PMs
    N54 = "d-Tocotrienol", ### ME_necrotic+ PMs
    P108 = "PC(15:0/18:0)" ## ME-steatotic+ PMs
    #P35 = "L-Palmitoylcarnitine",    
    
)

x_label <- "Difference in the mean normalized\nabundance levels between\nT and N sections (log2FC)"
x_breaks <- c(-2, -1, 0, 1, 2)
x_lim <- c(-2.5, 2.5)
y_breaks <- c(0, 2, 4, 6, 8)
y_lim <- c(0, 8)
text_col <- "#e41a1c"
text2_col <- "#377eb8"
control_col <- "darkgray"

PlotVolcano <- function(
    study_data, x_label, PM_group,
    x_breaks, x_lim, y_breaks, y_lim,
    text_col, 
    control_col = "darkgray"
) {
    peaks_ids_cur <- if (PM_group == "all") {
        names(PM_groups)
    } else {
        names(PM_groups)[(PM_groups == PM_group)]
    }

    plot_data <- study_data %>%
        filter(peak_id %in% peaks_ids_cur) %>%
        dplyr::rename(
            uni_pval = "p",
            uni_padj = "p.adj",
            uni_log2FC = "log2FC"
        ) %>%
        mutate(
            is_signif = factor(
                case_when(
                    uni_padj < 0.1 ~ "<0.1",
                    uni_padj < 0.2 ~ "<0.2",
                    TRUE ~ ">0.2"
                ),
                levels = c("<0.1", "<0.2", ">0.2")
            ),
            MSMS_annotations = if_else(
                peak_id %in% names(label_peak_annotations), peak_id, NA
                # mz_ids %in% label_mz_ids, MSMS_annotations, NA
            )
        ) %>%
        select(uni_log2FC, uni_pval, uni_padj, is_signif, MSMS_annotations) 
    
    up_num <- plot_data %>%
        filter(uni_padj < 0.1, uni_log2FC > 0) %>%
        nrow()
    down_num <- plot_data %>%
        filter(uni_padj < 0.1, uni_log2FC < 0) %>%
        nrow()
    total_num <- plot_data %>% nrow()

    ### Annotations
    num_annots <- data.frame(
        xpos = c(-1.5, 1.5),
        ypos =  c(Inf, Inf),
        annotatedText = c(
            paste0(down_num, "/", total_num, "\n(", signif(down_num/total_num*100, 3),"%)"),
            paste0(up_num, "/", total_num, "\n(", signif(up_num/total_num*100, 3),"%)")
        ),
        hjustvar = c(0.5, 0.5), #c(0, 1),
        vjustvar = c("inward", "inward") #,c(2, 2)
    )

    msi_TvN_p <- plot_data %>%
        ggplot(aes(x = uni_log2FC, y = -log10(uni_pval))) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
        # geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
        geom_point(
            data = ~ subset(., is_signif == ">0.2"), # , shape = 21
            aes(col = is_signif)
        ) +
        geom_point(
            data = ~ subset(., is_signif != ">0.2"), # , shape = 21
            aes(col = is_signif)
        ) +
        geom_text_repel(
            aes(label = MSMS_annotations),
            col = "black",
            box.padding = 0.4,
            force = 3,
            min.segment.length = 0,
            # nudge_x = -0.1,
            nudge_y = 0.5,
            size = 3,
            segment.size = 0.4,
            show.legend = FALSE
        ) +
        scale_x_continuous(
            name = x_label,
            breaks = x_breaks,
            limits = x_lim
        ) +
        scale_y_continuous(
            name = "-log10(P value)",
            breaks = y_breaks,
            limits = y_lim
        ) +
        scale_color_manual(
            name = "Padj",
            values = c(
                "<0.1" = text_col,
                "<0.2" = text2_col,
                ">0.2" = control_col
            )
        ) +
        geom_text(
            data = num_annots,
            aes(
                x = xpos, y = ypos,
                hjust = hjustvar, vjust = vjustvar,
                label = annotatedText
            ),
            size = 3,
            lineheight = 0.75
        ) +
        ggtitle(paste0(PM_group, " PMs")) +
        theme_classic() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
            legend.title.align = 0.5,
            legend.position = c(0.87, 0.15),
            legend.title = element_text(size = 10, face = "bold"),
            legend.text = element_text(size = 10),
            legend.direction = "vertical",
            legend.key.size = unit(0.15, "inch"),
            legend.background = element_blank(),
            axis.text = element_text(colour = "black")
        )

    ### save the plots
    SavePlot(
        msi_TvN_p,
        width = 2.5,
        height = 3.5,
        dir_name = "msi",
        file_name = paste0("msi_TvN_", str_replace_all(PM_group, "/", ""))
    )
}

for (PM_group in levels(PM_groups)) {
    PlotVolcano(
        test_res, x_label,
        PM_group = PM_group,
        x_breaks, x_lim, y_breaks, y_lim,
        text_col
    )
}

PlotVolcano(
    test_res, x_label,
    PM_group = "all",
    x_breaks, x_lim, y_breaks, y_lim,
    text_col
)