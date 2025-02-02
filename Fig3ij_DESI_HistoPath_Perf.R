### "========================================================================="
### Plot DESI-MSI ROI PLS-DA classifier performance
### Note: Step 3 must be run first
###
### Copyright (c) 2021-2024. Bioinformatics Institute, A*STAR, Singapore
### License: MIT License. See LICENSE for the full license.
### Authors: Lit-Hsin Loo (loolh@bii.a-star.edu.sg)
### "========================================================================="
library(here)
source(here("lib", "MSI_lib.R"))
source(here("lib", "Study_lib.R"))
library(pROC)

### "=========================================================================="
### Fig. 3i) Plot the ROC
### "=========================================================================="
data_type <- "roi"

roc_list <- list()
roc_list <- LoadROC(data_type, model_name = "plsda_multi")

roc_list_tmp <- LoadROC(data_type, model_name = "plsda_steatotic")
roc_list$training <- c(roc_list$training, roc_list_tmp$training)
roc_list$test <- c(roc_list$test, roc_list_tmp$test)

roc_list_tmp <- LoadROC(data_type, model_name = "plsda_fibrotic")
roc_list$training <- c(roc_list$training, roc_list_tmp$training)
roc_list$test <- c(roc_list$test, roc_list_tmp$test)

for (set_name in c("test", "training")) {

    set_name_upper <- paste0(
        toupper(substr(set_name, 1, 1)), substr(set_name, 2, nchar(set_name))
    )

    test_roc_p <- ggroc(roc_list[[set_name]]) +
        ggtitle("10 x 10-fold CV") +
        xlab(paste0(set_name_upper, " specificity (%)")) +
        ylab(paste0(set_name_upper, " sensitivity (%)")) +
        scale_color_manual(
            name = "ME regions",
            values = roi_annot_colors,
            labels = c(
                Normal="Normal", G2="Low-grade", G3="High-grade",
                Necrotic="Necrotic", Steatotic="Steatotic", Fibrotic="Fibrotic"
            )
        ) +
        geom_segment(
            aes(x = 100, xend = 0, y = 0, yend = 100),
            color = "darkgrey", linetype = "dashed"
        ) +
        coord_fixed() +
        theme_classic() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 9, face = "bold", margin=margin(0,0,0,0)),
            legend.title = element_text(hjust = 0.5, size = 7),
            # legend.title.align = 0.5,
            legend.position = "right", #c(0.8, 0.3),
            legend.text = element_text(size = 7),
            legend.spacing.y = unit(0.01, "inch"),
            legend.key.size = unit(0.12, "inch"),
            legend.background = element_blank(),
            axis.text = element_text(colour = "black")
        )

    #test_roc_p <- rasterize(test_roc_p, layers = 'Point', dpi = 600)

    SavePlot(
        test_roc_p,
        width = 3,
        height = 1.8,
        dir_name = "plsda",
        file_name = paste0("msi_roi_", set_name, "_roc")
    )
}

### "=========================================================================="
### Fig. 3j) Table of performance values
### "=========================================================================="
ExportSgMEPerf <- function(name) {
    data_type <- "roi"

    roc_list <- list()
    roc_list <- LoadROC(data_type, model_name = paste0(name, "_multi"))

    roc_list_tmp <- LoadROC(data_type, model_name = paste0(name, "_steatotic"))
    roc_list$training <- c(roc_list$training, roc_list_tmp$training)
    roc_list$test <- c(roc_list$test, roc_list_tmp$test)

    roc_list_tmp <- LoadROC(data_type, model_name = paste0(name, "_fibrotic"))
    roc_list$training <- c(roc_list$training, roc_list_tmp$training)
    roc_list$test <- c(roc_list$test, roc_list_tmp$test)

    perf_list <- list()
    perf_list <- LoadCVPerf(data_type, model_name = paste0(name, "_multi"))

    perf_list_tmp <- LoadCVPerf(data_type, model_name = paste0(name, "_steatotic"))
    perf_list$training <- rbind(perf_list$training, perf_list_tmp$training)
    perf_list$test <- rbind(perf_list$test, perf_list_tmp$test)

    perf_list_tmp <- LoadCVPerf(data_type, model_name = paste0(name, "_fibrotic"))
    perf_list$training <- rbind(perf_list$training, perf_list_tmp$training)
    perf_list$test <- rbind(perf_list$test, perf_list_tmp$test)

    auc_list <- lapply(roc_list, sapply, auc)
    perf_list$training <- cbind(AUC = auc_list$training, perf_list$training * 100)
    perf_list$test <- cbind(AUC = auc_list$test, perf_list$test * 100)

    mean(perf_list$test[, "Balanced Accuracy"])

    fs::dir_create(here("outputs"))

    readr::write_excel_csv(
        as_tibble(perf_list$training, rownames = "Region"),
        here("outputs", paste0("roi_", name, "_training_perf.csv"))
    )

    readr::write_excel_csv(
        as_tibble(perf_list$test, rownames = "Region"),
        here("outputs", paste0("roi_", name, "_test_perf.csv"))
    )
}

ExportSgMEPerf("svmlinear")
ExportSgMEPerf("svmrbf")
ExportSgMEPerf("plsda")
