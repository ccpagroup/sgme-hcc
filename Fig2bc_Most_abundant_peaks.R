### "========================================================================="
### Plot highly-abundant peaks
###
### Copyright (c) 2021-2024. Bioinformatics Institute, A*STAR, Singapore
### License: MIT License. See LICENSE for the full license.
### Authors: Lit-Hsin Loo (loolh@bii.a-star.edu.sg)
### "========================================================================="
library("here")
source(here("lib", "MSI_lib.R"))
source(here("lib", "Study_lib.R"))

### Load the peaks
load(here("data", "lcms", "lcms_peaks.Rdata"))

### "--------------------------------------------------------------------------"
### Plot Fig. 2b
### "--------------------------------------------------------------------------"
label_peak_annotations <- c(
    ### Some special case
    "Vitamin A",
    "Cholesterol sulfate",
    "SAPC",
    "LysoPC(18:1/0:0)",
    "PC(20:4/16:0)",
    "PI(20:2/17:1)",
    "DG(16:0/0:0/18:1)"
)
label_mz_ids <- lcms_hit_peaks$mz_id[
    lcms_hit_peaks$MSMS_annotation %in% label_peak_annotations
]

### Get some numbers
total_peaks <- nrow(lcms_raw_peaks)
total_highly_abundant_selected <- sum(lcms_raw_peaks$is_tested)

total_MSMS_confirmed <- sum(lcms_hit_peaks$is_msms)
total_database_matched <- sum(
    lcms_hit_peaks$is_msms & lcms_hit_peaks$MSMS_annotation != "Unannotated"
)

### A mini graph without the labels
msms_p_mini <- lcms_raw_peaks %>%
    arrange(-mean) %>%
    mutate(
        MSMS_annotation = str_wrap(MSMS_annotation, width = 20),
        peak_type = case_when(
            is_msms & MSMS_annotation != "Unannotated" ~ "matched",
            is_msms ~ "confirmed",
            TRUE ~ "all"
        )
    ) %>%
    dplyr::select(
        mz_id, mean, pct90, peak_type, is_tested, is_msms, MSMS_annotation
    ) %>%
    ggplot(
        aes(x = mean, y = pct90, col = peak_type)
    ) +
    geom_point(data = ~ subset(., peak_type == "all"), size = 0.5) +
    geom_point(data = ~ subset(., peak_type == "confirmed"), size = 0.5) +
    geom_point(data = ~ subset(., peak_type == "matched"), size = 0.5) +
    geom_text_repel(
        aes(label = MSMS_annotation),
        col = "black",
        data = ~ subset(., mz_id %in% label_mz_ids),
        force = 1,
        size = 2.5,
        segment.size = 0.3,
        box.padding = 0.5,
        max.overlaps = Inf,
        max.time = 5, max.iter = 1e6, 
        min.segment.length = 0, # draw all line segments
        nudge_x = -0.2,
        #nudge_y = 0.2,
        lineheight = 1,
        #xlim  = c(30, NA),
        #ylim = c(NA, 30),
        #direction = "y",
        hjust = "center"
    ) +
    scale_x_continuous(
        name = "Mean abundance level (log2)",
        breaks = c(0, 10, 20, 30),
        limits = c(0, 30)
    ) +
    scale_y_continuous(
        name = "90th-%tile abundance level (log2)",
        breaks = c(10, 20, 30),
        limits = c(7, 30)
    ) +
    scale_color_manual(
        name = "LC/ToF-MS features",
        breaks = c("all", "confirmed", "matched"),
        labels = c(
            all = paste0(
                "All (", total_peaks, " mass features)"
            ),
            confirmed = paste0(
                "MS/MS confirmed\n(", total_MSMS_confirmed,
                " putative metabolites)"
            ),
            matched = paste0(
                "Database matched\n(", total_database_matched,
                " annotated metabolites)"
            )
        ),
        values = c(
            all = "gray", confirmed = "#377eb8", matched = "#e41a1c"
        )
    ) +
    theme_classic() +
    theme(
        axis.text = element_text(colour="black"),
        legend.key.size = unit(0.5, 'cm'),
        legend.title.align = 0.5,
        legend.position = c(0.73, 0.18),
        legend.title = element_text(size = 7), # , face = "bold"),
        legend.text = element_text(size = 7),
        legend.background = element_blank()
    )

### Rasterize the plot
# library(ggrastr)
# msms_p_mini <- rasterize(msms_p_mini, layers = 'Point', dpi = 600)

### save the plots
fs::dir_create(here("figures", "lcms"))

for (img_type in c(".png", ".pdf")) {
    ggsave(
        filename = here(
            "figures", "lcms",
            paste0("MSMS_selection_mini", img_type)
        ),
        plot = msms_p_mini,
        width = 3,
        height = 3,
        units = "in",
        dpi = 600
    )
}

### "=========================================================================="
### Fig. 2c) Plot Treemap
### "=========================================================================="
lcms_hit_peaks <- AddClassAnnotation(lcms_hit_peaks)

compounds_summary <- lcms_hit_peaks %>%
    group_by(Class) %>%
    summarise(num_compounds = n(), .groups = "drop") %>%
    mutate(
        class_color = metabolite_classes_colors[Class]
    )

fs::dir_create(here("figures", "lcms"))
pdf(
    here("figures", "lcms", "lcms_metabolite_treemap.pdf"),
    width = 5, height = 4
)
treemap::treemap(
    compounds_summary,
    index = c("Class"),
    vSize = "num_compounds",
    vColor = "class_color",
    type = "color",
    # algorithm = "pivotSize",
    sortID = "-num_compounds",
    title = "",
    fontsize.labels = c(12, 8),
    border.lwds = c(2, 0.75),
    aspRatio = 1.2,
    force.print.labels = TRUE
)
dev.off()