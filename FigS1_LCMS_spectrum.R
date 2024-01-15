### "========================================================================="
### Plot Fig S1
###
### Copyright (c) 2021-2024. Bioinformatics Institute, A*STAR, Singapore
### License: MIT License. See LICENSE for the full license.
### Authors: Lit-Hsin Loo (loolh@bii.a-star.edu.sg)
### "========================================================================="
library("here")
source(here("lib", "MSI_lib.R"))
source(here("lib", "Study_lib.R"))

### Load all the raw data
lcms_raw <- readRDS(here("data", "lcms", "lcms_raw.rds"))

### Get clinical info
clinical_info <- LoadClinicalInfo()

### Load the peaks
load(here("data", "lcms", "lcms_peaks.Rdata"))

### "=========================================================================="
### Fig. S1 The LC/MS raw data plot
### "=========================================================================="
### Plot a heatmap for all the tumors
GetLCMSPlotData <- function(lcms_mode, is_diff) {
    plot_data <- left_join(
        lcms_raw,
        clinical_info %>% select(case_id, histological_diagnosis),
        by = "case_id"
    ) %>%
        ### Filter for a partical case and section
        filter(
            #case_id == !!case_id,
            #histological_diagnosis == "HCC",
            endsWith(mz_id, lcms_mode)
        )

    ### Get case number
    sample_ids_all <- unique(paste0(lcms_raw$case_id, "-", lcms_raw$section_id))
    sample_ids_cur <- unique(paste0(plot_data$case_id, "-", plot_data$section_id))

    case_num <- length(sample_ids_cur)
    message("Case number = ", case_num)
    message(
        "Missing = ", 
        paste(sample_ids_all[!(sample_ids_all %in% sample_ids_cur)], collapse = ", ")
    )

    if (is_diff) {
        plot_data <- plot_data %>%
            ### Take the average of all tissues for a case
            group_by(case_id, mz_id) %>%
            summarize(
                lcms_abd = mean(
                    lcms_abd[startsWith(section_id, "T")],
                    na.rm = TRUE
                ) - mean(
                    lcms_abd[section_id == "N"],
                    na.rm = TRUE
                ),
                .groups = "drop"
            ) %>%
            ### Take the average of all cases
            group_by(mz_id) %>%
            summarize(
                lcms_abd = mean(lcms_abd, na.rm = TRUE),
                .groups = "drop"
            )

    } else {
        plot_data <- plot_data %>%
            ### Take the average of all tissues for a case
            group_by(case_id, mz_id) %>%
            summarize(
                lcms_abd = mean(
                    lcms_abd,  #[startsWith(section_id, section_type)],
                    na.rm = TRUE
                ),
                .groups = "drop"
            ) %>%
            ### Take the average of all cases
            group_by(mz_id) %>%
            summarize(
                lcms_abd = mean(lcms_abd, na.rm = TRUE),
                .groups = "drop"
            )
    }

    ### Get the rt and mz
    plot_data <- plot_data %>%
        extract(
            mz_id,
            c("rt", "mz"),
            "([[:digit:]\\.]+)_([[:digit:]\\.]+)[[:alpha:]\\/]+_[[:alpha:]]+"
        ) %>%
        mutate(
            rt = as.double(rt),
            mz = as.double(mz)
        )

    ### Find min_max
    message("m/z 1% and 99%tile = ")
    print(round(quantile(plot_data$mz, c(0.01, 0.99)), 2))

    plot_data %>%
        ggplot(
            aes(x = mz, y = rt, col = lcms_abd)
        ) +
        geom_point(shape = 20, size = 0.75) +
        scale_color_gradient(
            limits = c(0, 30),
            name = "Mean ion\nabundance (log2)",
            low = "white", high = "black",
            guide = guide_colorbar(
                frame.colour = "black", ticks.colour = "black",
                barheight = 0.5
            )
        ) +
        ggtitle(paste0(lcms_mode, " ions (s = ", case_num, ")")) +
        scale_x_continuous(limits = c(0, 2000), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, 12), expand = c(0, 0)) +
        ylab("Retention time [min]") +
        xlab("m/z") +
        theme_classic() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 8),
            #legend.direction="horizontal",
            legend.position = "top",
            axis.text = element_text(colour="black")
        )
}

### Plot the LC/MS heatmap for all tissues
# case_id <- "HEP0152"
# lcms_mode <- "LPOS"
# section_type <- "T"

lcms_heatmaps <- list(
    LNEG = GetLCMSPlotData(
        lcms_mode = "LNEG", is_diff = FALSE
    ),
    LPOS = GetLCMSPlotData(
        lcms_mode = "LPOS", is_diff = FALSE
    ),
    HPOS = GetLCMSPlotData(
        lcms_mode = "HPOS", is_diff = FALSE
    )
)

### save the plots
fs::dir_create(here("figures", "lcms"))

for (heatmap_name in names(lcms_heatmaps)) {
    for (img_type in c(".png", ".pdf")) {
        ggsave(
            filename = here(
                "figures", "lcms",
                paste0("lcms_heatmap_", heatmap_name, img_type)
            ),
            plot = lcms_heatmaps[[heatmap_name]],
            width = 3,
            height = 4,
            units = "in",
            dpi = 300
        )
    }
}