### "========================================================================="
### General SgME Profiling library
###
### Copyright (c) 2021-2024. Bioinformatics Institute, A*STAR, Singapore
### License: MIT License. See LICENSE for the full license.
### Authors: Lit-Hsin Loo (loolh@bii.a-star.edu.sg)
### "========================================================================="

library("tidyverse")
library("ggrepel")
library(ggbeeswarm)
library(ggpubr)
library(rstatix)
library("Cardinal")
library("MsCoreUtils")
library("ijtiff")
library("RImageJROI")
library("jsonlite")
library("here")
library("parallel")
library("readxl")
library("RColorBrewer")
library("png")
library("fs")
library("MetaboCoreUtils")
RNGkind("L'Ecuyer-CMRG")
library("e1071")
library(glmnet)

set.seed(1)
setCardinalBPPARAM(MulticoreParam(workers = 20))
mc.cores <- 20

#setCardinalBPPARAM(SerialParam())

### Load the study configurations
source(here("conf", "study_conf.R"))

### Constants ================================================================
valid_ROI_types <- c("T", "F", "C", "S", "N")
valid_section_types <- c("T", "N")
valid_reg_types <- unique(c(valid_ROI_types, valid_section_types))

### General Utilities ========================================================
### "========================================================================="
### Filter the samples
### "========================================================================="
GetSamples <- function(msi_conf, ref_slide_ids) {

    ### Do we have reference ids?
    if (missing(ref_slide_ids)) {
        ### Identify all the samples
        list(
            pos = Filter(
                function(x) x$ion_type == "pos",
                msi_conf$samples
            ),
            neg = Filter(
                function(x) x$ion_type == "neg",
                msi_conf$samples
            )
        )
    } else {
        ### Identify all the samples
        list(
            pos = Filter(
                function(x) x$ion_type == "pos",
                msi_conf$samples
            ),
            pos_ref = Filter(
                function(x) x$ion_type == "pos" && x$slide_id %in% ref_slide_ids,
                msi_conf$samples
            ),
            neg = Filter(
                function(x) x$ion_type == "neg",
                msi_conf$samples
            ),
            neg_ref = Filter(
                function(x) x$ion_type == "neg" && x$slide_id %in% ref_slide_ids,
                msi_conf$samples
            )
        )
    }
}

### "========================================================================"
### Load all MSI data for the whole study (only for ROI!!!)
### "========================================================================"
LoadStudyMSIData <- function() {
    readRDS(file = here("results", "study_data.rds"))
}

### "========================================================================"
### Save all MSI ROI data for the whole study and calculate the p.adj
### "========================================================================"
SaveStudyMSIData <- function(
    is_normalized, peak_summary_mode, slide_ids
) {

    ### Define the dataset
    msi_conf <- read_json(here("conf", "msi_conf.json"))

    ### Load the LCMS hits
    load(here("data", "lcms", "lcms_peaks.Rdata"))

    ### Some of LC/MS peaks may have the same exact mass, and thus the same
    ### peak ID. So, we need to consolidate them here
    consol_lcms_peaks <- lcms_hit_peaks %>%
        group_by(peak_id, mz) %>%
        summarize(
            MSMS_annotations = paste0(MSMS_annotation, collapse = ";"),
            mz_ids = paste0(mz_id, collapse = ";"),
            .groups = "drop"
        )

    pos_metabolites_ROI <- AnalyzeMetabolites(
        ion_type = "pos",
        slide_ids = slide_ids,
        msi_conf = msi_conf,
        peak_summary_mode = peak_summary_mode,
        is_normalized = is_normalized
    )

    export_pos_hits <-
        left_join(
            pos_metabolites_ROI,
            consol_lcms_peaks,
            by = "peak_id"
        ) %>% relocate(
            `mz`, `mz_ids`, `MSMS_annotations`,
            .after = `peak_id`
        )

    neg_metabolites_ROI <- AnalyzeMetabolites(
        ion_type = "neg", 
        slide_ids = slide_ids,
        msi_conf = msi_conf,
        peak_summary_mode = peak_summary_mode,
        is_normalized = is_normalized
    )

    export_neg_hits <-
        left_join(
            neg_metabolites_ROI,
            consol_lcms_peaks,
            by = "peak_id"
        ) %>% relocate(
            `mz`, `mz_ids`, `MSMS_annotations`,
            .after = `peak_id`
        )

    ### Bind and adjust p-values
    study_data <- bind_rows(
        export_pos_hits,
        export_neg_hits
    ) %>%
    mutate(
        transformed_padj = p.adjust(transformed_pval, method = "fdr"),
        fibrotic_padj = p.adjust(fibrotic_pval, method = "fdr"),
        necrotic_padj = p.adjust(necrotic_pval, method = "fdr"),
        steatotic_padj = p.adjust(steatotic_pval, method = "fdr"),
        diff_grade_padj = p.adjust(diff_grade_pval, method = "fdr"),
        diff_grade2_padj = p.adjust(diff_grade2_pval, method = "fdr"),
        diff_grade3_padj = p.adjust(diff_grade3_pval, method = "fdr")
    ) %>% relocate(
        `transformed_padj`,
        `fibrotic_padj`,
        `necrotic_padj`,
        `steatotic_padj`,
        `diff_grade_padj`, `diff_grade2_padj`, `diff_grade3_padj`,
        .after = `diff_grade3_pval`
    ) %>% arrange(
        transformed_padj
    )

    ### Add VIP value
    #study_data <- left_join(
    #    study_data,
    #    lcms_hit_peaks %>% dplyr::select(peak_id, vip, TvsN_padj, TvsN_log2FC),
    #    by = "peak_id"
    #) %>% relocate(
    #    vip, .after = "mz_ids"
    #)

    ### Save the data
    saveRDS(
        study_data,
        file = here("results", "study_data.rds"),
        compress = "xz"
    )
}

### "========================================================================"
### Find ROI type
### * The format of the region label is:
###   - "Region type""Section index"-"Region index"_"differentiation grade"
###   - For example, S4a-4_3
###     Region type = ROI because Region index != 0
###     ROI type = S
###     Section index = T4a because section label index is not 0
###     Region index = ROI #4
###     Diff grade = Edmonson Grade 3
### "========================================================================"
FindROITypes <- function(
    data_all
) {
    ### Load the clinical information
    clinical_info <- LoadClinicalInfo()

    ### Get all slide ids
    slide_ids <- gsub("^([^\\_]*)\\_.*$", "\\1", colnames(data_all))
    #slide_ids <- unique(slide_ids)

    region_ids <- gsub("^[^\\.]*\\.(.*)$", "\\1", colnames(data_all))

    ### Determine the tumor grade level from T, F, S, C sections
    ### * Valid entry should be U, 0, 1, 2, 3, 4
    ### * 0 or missing = This is not a tumor section or ROI.
    diff_grades <- gsub("^[TFCS].*-._(.)$", "\\1", region_ids)
    is_no_grade <- (diff_grades == region_ids)

    ### Regions from tumor sections must have defined diff grades
    if (any((!startsWith(region_ids, "N")) & is_no_grade)) {
        stop("Differentiation grade undefined")
    }

    ### Those undefined are from non-tumor sections
    diff_grades[is_no_grade] <- "0"

    ### Check again to make sure all grades are valid
    if (any(!(diff_grades %in% c("U", "0", "1", "2", "3", "4")))) {
        stop("Differentiation grade invalid")
    }

    ### Find out the categories for each regions.
    ### Some of them may be completely empty!

    ### Usually on the N regions, but could be F or S regions too!
    is_nontransformed <- (diff_grades == "0") #& startsWith(region_ids, "N")

    ### Non-transformed and non-steatotic regions
    is_normal <- (diff_grades == "0") & startsWith(region_ids, "N")

    ### Must be from the T ROIs
    is_transformed <-
        (diff_grades %in% c("1", "2", "3", "4")) & startsWith(region_ids, "T")

    ### Make sure there is no overlap
    if (any(is_transformed & is_nontransformed)) {
        stop("is_tramsformed overlaps with is_nontransformed")
    }

    is_diff_grade2 <- (diff_grades == "2") & startsWith(region_ids, "T")
    is_diff_grade3 <- (diff_grades == "3") & startsWith(region_ids, "T")

    is_fibrotic  <- startsWith(region_ids, "F")
    is_necrotic  <- startsWith(region_ids, "C")

    is_steatotic <-
        startsWith(region_ids, "S") &
        (diff_grades %in% c("0", "1", "2"))

    ### Return
    list(
        slide_ids = slide_ids,
        diff_grades = diff_grades,
        is_nontransformed = is_nontransformed,
        is_normal = is_normal,
        is_transformed = is_transformed,
        is_fibrotic = is_fibrotic,
        is_necrotic = is_necrotic,
        is_steatotic = is_steatotic,
        is_diff_grade2 = is_diff_grade2,
        is_diff_grade3 = is_diff_grade3
    )
}

### "========================================================================"
### Load all MSI data as a matrix
### * Please use this function to load all MSI data, instead of GetPeakStats()
###   because it will add the additional metadata and normalization
### "========================================================================"
LoadAllMSIMat <- function(
    slide_ids, ## ref_slide_ids
    is_normalized, ## is_normalized_study
    max_na_fraction # 0.50
) {

    ### Define the dataset
    msi_conf <- read_json(here("conf", "msi_conf.json"))

    ### Load and combine two MSI datasets
    msi_mat_pos <- LoadMSIMat(
        ion_type = "pos",
        slide_ids = slide_ids,
        peak_type = "ref_peaks",
        peak_summary_mode = peak_summary_mode_study,
        stat_type = "roi",
        msi_conf = msi_conf,
        is_normalized = is_normalized,
        max_na_fraction = max_na_fraction
    )

    colnames(msi_mat_pos) <- sub("_pos.", ".", colnames(msi_mat_pos))

    msi_mat_neg <- LoadMSIMat(
        ion_type = "neg",
        slide_ids = slide_ids,
        peak_type = "ref_peaks",
        peak_summary_mode = peak_summary_mode_study,
        stat_type = "roi",
        msi_conf = msi_conf,
        is_normalized = is_normalized,
        max_na_fraction = max_na_fraction
    )

    colnames(msi_mat_neg) <- sub("_neg.", ".", colnames(msi_mat_neg))

    ROI_names <- sort(unique(c(colnames(msi_mat_pos), colnames(msi_mat_neg))))
    peak_ids <- c(rownames(msi_mat_pos), rownames(msi_mat_neg))

    msi_mat <- matrix(
        NA,
        nrow = length(peak_ids),
        ncol = length(ROI_names),
        dimnames = list(peak_ids, ROI_names)
    )

    ### Fill up the pos
    pos_roi_names <- colnames(msi_mat_pos)
    pos_peak_ids <- rownames(msi_mat_pos)
    for (peak_id in pos_peak_ids) {
        msi_mat[peak_id, pos_roi_names] <- msi_mat_pos[peak_id, pos_roi_names]
    }

    ### Fill up the neg
    neg_roi_names <- colnames(msi_mat_neg)
    neg_peak_ids <- rownames(msi_mat_neg)
    for (peak_id in neg_peak_ids) {
        msi_mat[peak_id, neg_roi_names] <- msi_mat_neg[peak_id, neg_roi_names]
    }

    msi_mat
}

### "========================================================================"
### Load raw MSI data for an ion type as a matrix
### "========================================================================"
LoadAllRawMSIMat <- function(
    slide_ids, is_normalized, 
    max_na_fraction = 0.99,
    ### Do we want to assign zero to the min values?
    ### We will do it for classifiers, because some features may not be detected
    ### But for fold-change analysis, we will remove all the missing samples
    is_undetected_min = TRUE
) {

    ### Define the dataset
    msi_conf <- read_json(here("conf", "msi_conf.json"))

    msi_raw <- c(
        GetRawPeaks(
            ion_type = "pos",
            msi_conf,
            peak_type = "ref_peaks",
            peak_summary_mode = "lcms_only",
            stat_type = "tissue"
        ),
        GetRawPeaks(
            ion_type = "neg",
            msi_conf,
            peak_type = "ref_peaks",
            peak_summary_mode = "lcms_only",
            stat_type = "tissue"
        )
    )

    ### Find all the common peak_names
    peak_names <- NULL
    for (sample_id in names(msi_raw)) {
        for (reg_name in names(msi_raw[[sample_id]])) {
            peak_names <- c(
                peak_names,
                rownames(msi_raw[[sample_id]][[reg_name]])
            )
        }
    }
    peak_names <- unique(peak_names)

    ### Find out all sample ids
    sample_ids <- unique(sub("^(.+)_[^\\_]+", "\\1", names(msi_raw)))

    msi_raw_combined <- list()

    for (sample_id in sample_ids) {
        ### Get the pos and neg sample Id
        pos_sample_id <- paste0(sample_id, "_pos")
        neg_sample_id <- paste0(sample_id, "_neg")

        ### Check to make sure the sections matches
        pos_section_names <- names(msi_raw[[pos_sample_id]])
        neg_section_names <- names(msi_raw[[neg_sample_id]])

        common_section_names <- unique(pos_section_names, neg_section_names)

        for (section_name in common_section_names) {
            ### The pixels may not be the same, because some pixels may be
            ### skipped or missing in some of the samples
            pos_xy <- colnames(msi_raw[[pos_sample_id]][[section_name]])
            neg_xy <- colnames(msi_raw[[neg_sample_id]][[section_name]])
            common_xy <- unique(c(pos_xy, neg_xy))

            ### Get the pos data
            if (length(pos_xy) > 0) {
                pos_data <- matrix(
                    NA,
                    nrow = nrow(msi_raw[[pos_sample_id]][[section_name]]),
                    ncol = length(common_xy),
                    dimnames = list(
                        rownames(msi_raw[[pos_sample_id]][[section_name]]),
                        common_xy
                    )
                )
                pos_data[, pos_xy] <- msi_raw[[pos_sample_id]][[section_name]]

            } else {
                pos_data <- NULL
            }

            ### Get the neg data
            if (length(neg_xy) > 0) {
                neg_data <- matrix(
                    NA,
                    nrow = nrow(msi_raw[[neg_sample_id]][[section_name]]),
                    ncol = length(common_xy),
                    dimnames = list(
                        rownames(msi_raw[[neg_sample_id]][[section_name]]),
                        common_xy
                    )
                )
                neg_data[, neg_xy] <- msi_raw[[neg_sample_id]][[section_name]]

            } else {
                neg_data <- NULL
            }

            ### Note: One or both of the matrices may be empty!!
            cur_data <- rbind(pos_data, neg_data)

            ### Init an empty matrix
            new_data <- matrix(
                NA,
                nrow = length(peak_names), ncol = ncol(cur_data),
                dimnames = list(
                    peak_names,
                    paste0(
                        sample_id, ".", section_name, ".", colnames(cur_data)
                    )
                )
            )

            ### Copy the data
            new_data[rownames(cur_data), ] <- cur_data

            ### Update the original
            msi_raw_combined[[sample_id]][[section_name]] <- new_data
        }
    }

    ### Merge and return
    msi_raw_mat <- do.call(
        cbind, lapply(msi_raw_combined, function(x) do.call(cbind, x))
    )

    ### Only select those in the ref_slide_ids
    good_idx <- sub(
        "^([^\\_]+)\\_.+$", "\\1", colnames(msi_raw_mat)
    ) %in% slide_ids
    msi_raw_mat <- msi_raw_mat[,good_idx, drop = FALSE]

    ### Find out where are the zeros before normalization
    is_zero_mat <- !is.na(msi_raw_mat) & (msi_raw_mat == 0)

    ### Perform normalization
    if (is_normalized) {
        ### We must replace zero with NA before normalization to avoid the 
        ### normalization scaling the zeros. The users must decide what to do 
        ### with the NA values. The recommended way is to set NAs with an 
        ### arbitrary low value
        msi_raw_mat[is_zero_mat] <- NA

        ### Perform normalization
        norm_coeffs <- c(
            LoadNormCoeff(ion_type = "pos"),
            LoadNormCoeff(ion_type = "neg")
        )

        sample_ids <- sub("^([^\\.]+)\\..+$", "\\1", colnames(msi_raw_mat))

        pos_row_norm_coeffs <- norm_coeffs[paste0(sample_ids, "_pos")]
        neg_row_norm_coeffs <- norm_coeffs[paste0(sample_ids, "_neg")]

        for (row_idx in seq_len(nrow(msi_raw_mat))) {
            ### Is the row positive or negative?
            if (startsWith(peak_names[row_idx], "P")) {
                msi_raw_mat[row_idx, ] <-
                    msi_raw_mat[row_idx, ] - pos_row_norm_coeffs

            } else {
                msi_raw_mat[row_idx, ] <-
                    msi_raw_mat[row_idx, ] - neg_row_norm_coeffs
            }
        }

        ### Put the zero back
        if (is_undetected_min) {
            msi_raw_mat[is_zero_mat] <- 0

            ### Replace zeros with min value
            msi_raw_mat <- ReplaceZeros(msi_raw_mat, is_zero_mat)
        }
    }

    ### We will only retain MSI peaks detected in > 50% of the samples
    msi_raw_mat <- RemoveRarePeaks(
        msi_raw_mat, is_zero_mat, max_na_fraction
    )

    ### return
    msi_raw_mat
}

### "========================================================================"
### Load the summarized MSI data for an ion type as a matrix
### * Please use this function to load all MSI data, instead of GetPeakStats()
###   because it will add the additional metadata and normalization
### "========================================================================"
LoadMSIMat <- function(
    ion_type, slide_ids, peak_type, peak_summary_mode, stat_type, msi_conf,
    is_normalized,
    max_na_fraction = 0.5
) {

    ### Load the data
    roi_stats <- LoadMSIData(
        ion_type = ion_type,
        peak_type = peak_type,
        peak_summary_mode = peak_summary_mode,
        stat_type = stat_type,
        msi_conf = msi_conf
    ) %>%
    filter(
        slide_id %in% slide_ids
    ) %>%
    mutate(
        sample_label = paste0(sample_id, ".", reg_name),
        ## We replace NA or NaN back with zeros
        msi_mean = replace_na(msi_mean, 0),
        ### msi_norm_mean is msi_mean divided by the norm_coeff.
        ### So, if msi_norm_mean is NA, that means msi_mean is NA
        msi_norm_mean = replace_na(msi_norm_mean, 0)
    )

    ### Convert the data to a matrix
    res <- if (is_normalized) {
        roi_stats %>%
            dplyr::select(
                sample_label, peak_id, msi_norm_mean
            ) %>%
            pivot_wider(
                names_from = peak_id,
                values_from = msi_norm_mean,
                values_fill = NA
            )
    } else {
        roi_stats %>%
            dplyr::select(
                sample_label, peak_id, msi_mean
            ) %>%
            pivot_wider(
                names_from = peak_id,
                values_from = msi_mean,
                values_fill = NA
            )
    }

    ### Convert to a matrix
    sample_labels <- res$sample_label
    data_all <- t(as.matrix(res[, -1]))
    colnames(data_all) <- sample_labels

    ### Replace NA with the min value
    is_zero_mat <- data_all == 0
    data_all <- ReplaceZeros(data_all, is_zero_mat)

    ### We will only retain MSI peaks detected in > 50% of the samples
    data_all <- RemoveRarePeaks(data_all, is_zero_mat, max_na_fraction)

    ### Return data_all
    data_all
}

### "========================================================================"
### Remove columns with mostly NAs or zeros
### "========================================================================"
RemoveRarePeaks <- function(data_mat, is_zero_mat, max_na_fraction) {

    zero_nums <- apply(is_zero_mat, 1, sum, na.rm = TRUE)
    non_na_nums <- apply(is_zero_mat, 1, function(x) sum(!is.na(x)))

    is_retained <- zero_nums <= (max_na_fraction * non_na_nums)

    ### Is there anything to remove?
    if (sum(is_retained) < nrow(is_zero_mat)) {
        message("Old peak num = ", nrow(data_mat))
        data_mat <- data_mat[is_retained, , drop = FALSE]
        message("New peak num = ", nrow(data_mat))
    }

    ### return
    data_mat
}

### "========================================================================"
### Replace zeros of a matrix with min value
### * data_mat may have NA values
### "========================================================================"
ReplaceZeros <- function(data_mat, is_zero_mat) {

    if (any(dim(data_mat) != dim(is_zero_mat))) {
        stop("data_mat and is_zero_mat have different dimensions")
    }

    for (i in 1:nrow(is_zero_mat)) {
        ### Find out where the NAs and zeros are
        is_zero <- is_zero_mat[i, ]
        is_na   <- is.na(data_mat[i, ])

        ### It cannot be all zeros
        if (sum(!is_na & is_zero) != sum(!is_na)) {
            ### Also, there must be zeros and nonzeros
            if (
                (sum(!is_na & is_zero) > 0) &&
                (sum(!is_na & !is_zero) > 0)
            ) {
                data_mat[i, !is_na & is_zero] <- min(
                    data_mat[i, !is_na & !is_zero],
                    na.rm = TRUE
                ) - 1.0
            }
        }
    }

    data_mat
}

### "========================================================================"
### Load all MSI data for an ion type
### * Note: Almost all of the time, you should be using LoadMSIMat()
###   which include all the preprocessing of the data. This function should
###   only be used to load the raw un-processed data
### "========================================================================"
LoadMSIData <- function(
    ion_type, peak_type, peak_summary_mode, stat_type, msi_conf
) {
    readRDS(
        file = here(
            "data", "msi", peak_type, ion_type,
            paste0("peak_", stat_type, "_", peak_summary_mode, ".rds")
        )
    )
}

### "========================================================================"
### Get reference roi names
### "========================================================================"
GetRefROINames <- function(case_conf) {
    ### Make sure type is valid
    #if (!(roi_type %in% c("T", "F", "C"))) {
    #    stop("Unknown roi_type")
    #}

    ### Get all the ROIs
    regs <- read.ijzip(here(case_conf$roi_path), names = TRUE)

    ### Get all the reference rois
    ref_roi_ids <- NULL

    for (reg_name in names(regs)){
        ### We assume the references are always rect
        if (
            regs[[reg_name]]$strType == "rect" &&
            any(startsWith(reg_name, valid_ROI_types))
        ) {
            ref_roi_ids <- c(ref_roi_ids, reg_name)
        }
    }

    ref_roi_ids
}

### "========================================================================="
### Get the names of the given region type
### * if reg_name_in is provided, only that reg_name is returned
### "========================================================================="
GetRegNames <- function(
    regs, reg_name_in = NULL, reg_type = "freehand"
) {
    ### Init the reg_names
    reg_names <- NULL

    for (reg_name in names(regs)) {
        ### Do we want to process this ROI?
        if (!is.null(reg_name_in) && (reg_name != reg_name_in)) {
            next
        }

        ### Define the tissue region name
        reg <- regs[[reg_name]]

        ### We only process freehand
        if (reg$strType != reg_type) {
            next
        }

        reg_names <- c(reg_names, reg_name)
    }

    reg_names
}

### "========================================================================="
### Query the size of a MSI image
### "========================================================================="
GetMSIDim <- function(
    sample_id, ion_mode
) {
    ### Define the dataset
    msi_conf <- read_json(here("conf", "msi_conf.json"))
    full_sample_id <- paste0(sample_id, "_", ion_mode)
    # sample_id <- "A-01-01_1_pos"

    ### Get the case conf
    case_conf <- Filter(
        function(x) x$sample_id == full_sample_id,
        msi_conf$samples
    )[[1]]

    ### Return
    c(case_conf[["width"]], case_conf[["height"]])
}

### Get the prediction tibble
GetPredTibble <- function(
    tissue_pred_labels
) {
    ### Form a tibble
    tibble(
            id = names(tissue_pred_labels),
            pred_label = tissue_pred_labels
        ) %>%
        separate(
            col = id, into = c("sample_id", "section_id", "pos"), sep  = "\\."
        ) %>%
        separate(
            col = pos, into = c("tmp", "x", "y"), sep = "_"
        ) %>%
        mutate(
            tmp = NULL
        )
}

### Save the PLDA predcitons as images
SaveSgMEMapAsImage <- function(
    pred_model, model_name
) {
    ### Get the tibble
    tissue_pred_tibble <- GetPredTibble(pred_model$pred_labels)
    pred_levels <- levels(tissue_pred_tibble$pred_label)

    for (sample_id in unique(tissue_pred_tibble$sample_id)) {

        ### Prepare the full pred img, so that we can do a median filter
        sample_pred_img <- pred_model$pred_img_list[[sample_id]]
        level_values <- pred_model$level_values

        ### Loop through all the pred_levels
        for (pred_level in pred_levels) {
            pred_img <- matrix(
                0,
                nrow = nrow(sample_pred_img),
                ncol = ncol(sample_pred_img)
                #nrow = nrow(sample_pred$pred_img),
                #ncol = ncol(sample_pred$pred_img)
            )

            pred_img[
                sample_pred_img == level_values[pred_level]
                #sample_pred$pred_img == sample_pred$level_values[pred_level]
            ] <- 255

            ## We need to flip the data
            pred_img <- flip(flop(pred_img))

            fs::dir_create(here("images", "SgME_Maps"))
            write_tif(
                pred_img,
                here(
                    "images", "SgME_Maps",
                    paste0(
                        model_name, "_", sample_id, "_", pred_level, ".tiff"
                    )
                ),
                bits_per_sample = 8L, # ifelse(is_16bit, 16L, 8L),
                compression = "LZW", overwrite = TRUE
            )
        }
    }
}

### "========================================================================="
### Save a gg plot
### "========================================================================="
SavePlot <- function(
    cur_p, width, height, dir_name, file_name
) {
    fs::dir_create(here("figures", dir_name))

    for (img_type in c(".png", ".pdf")) {
        ggsave(
            filename = here(
                "figures", dir_name,
                paste0(file_name, img_type)
            ),
            plot = cur_p,
            width = width,
            height = height,
            units = "in",
            dpi = 300,
            device = ifelse(img_type == ".png", png, pdf)
        )
    }
}

### "========================================================================="
### Get the sample conf
### "========================================================================="
GetSampleConf <- function(sample_id) {
    ### Define the dataset
    msi_conf <- read_json("conf/msi_conf.json")

    ### Identify all the samples
    samples <- GetSamples(msi_conf, ref_slide_ids)

    ### Find out cur_sample
    if (endsWith(sample_id, "pos")) {
        ion_type <- "pos"

    } else if (endsWith(sample_id, "neg")) {
        ion_type <- "neg"

    } else {
        stop("Invalid sample_id provided")

    }

    cur_samples <- samples[[ion_type]]

    ### Find cur_sample
    idx <- which(sapply(cur_samples, function(x) x$sample_id) == sample_id)
    if (length(idx) == 0) {
        stop("sample_id cannot be found")
    }

    ### Return
    cur_samples[[idx]]
}

### General peak processing ==================================================
### "========================================================================="
### Build ID Lookup table based on unique mz
### "========================================================================="
AddPeakID <- function(lcms_peaks) {
    ### Add ID based on unique mz
    id_lookup_table <- NULL
    for (ion_type in c("pos", "neg")) {
        peaks <- lcms_peaks %>%
            filter(`ion_type` == !!ion_type) %>%
            arrange(mz) %>%
            dplyr::pull(mz) %>%
            unique()

        ### Add the labels, which is sorted by mz
        front_id <- ifelse(ion_type == "pos", "P", "N")

        id_lookup_table <- bind_rows(
            id_lookup_table,
            tibble(
                ion_type = ion_type,
                mz = peaks,
                peak_id = paste0(front_id, seq_along(peaks)),
            )
        )
    }

    left_join(
        lcms_peaks,
        id_lookup_table,
        by = c("ion_type", "mz")
    )
}

### "========================================================================="
### Get LC/MS ions
### "========================================================================="
GetLCMSIons <- function(ion_type) {
    ### According to https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9529534/
    ### and https://pubs.acs.org/doi/10.1021/acs.analchem.0c04720#
    ### [M+H]+ 74.0%
    ### [M−H]- 80.7
    if (ion_type == "pos") {
        c(
            "[M]+", "[M+H]+"
            # "[M+Na]+", "[M+K]+"
            #"[M+NH4]+", "[M+H-H2O]+"
            # "[2M+H]+", "[M+H+Na]2+"
        )
    } else {
        c(
            "[M]-", "[M-H]-"
            #"[2M-H]-", "[M+Na-2H]-", "[M+K-2H]-"
            # "[M+Cl]-", "[M+CHO2]-", "[M+C2H3O2]-",
        )
    }
}

### "========================================================================="
### Get DESI ions
### "========================================================================="
GetDESIIons <- function(ion_type) {
    if (ion_type == "pos") {
        c(
           "[M]+", "[M+H]+"
           # "[M+Na]+", "[M+K]+"
        )
    } else {
        c("[M]-", "[M-H]-")
    }
}

### "========================================================================="
### Generate possible ions based on the mass types
### * peak_list = a list with lcms_mzs and (optional) exact_mzs
### "========================================================================="
GetPossibleDESIAdducts <- function(
    ion_type, peak_list
) {
    ### Get the peak list
    lcms_mzs <- peak_list[["lcms_mzs"]]  ### This is a list
    exact_mzs <- peak_list[["exact_mzs"]] ### This is an array

    ### exact_mzs may be null
    if (!is.null(exact_mzs)) {
        ### Makesure exact_mzs is a list
        if (class(exact_mzs) != "list") {
            stop("exact_mzs must be a list")
        }

        ### Check to make sure they match
        if (any(names(lcms_mzs) != names(exact_mzs))) {
            stop("The peak IDs of lcms_mzs and exact_mzs do not match.")
        }
    }

    ### What are the possible adducts for the source?
    #source_adducts <- GetLCMSIons(ion_type)

    ### What are the possible adducts for the target?
    #target_adducts <- GetDESIIons(ion_type)

    ### Get the ion_sign character
    ion_sign <- ifelse(ion_type == "pos", "+", "-")

    ### Is exact_mzs provided?
    adduct_mat <- if (!is.null(exact_mzs)) {
        ### "------------------------------------------------------------------"
        ### Both exact and LC/MS mzs are provided.
        ### * We will always use the lcms mzs, called [X]
        ### * For those peaks with dissimilar exact and LC/MS mzs, we will also
        ###   include their exact mzs
        ### "------------------------------------------------------------------"
        ### Get peaks for the lcms mzs. We assume all of them are adducts!!
        orig_adduct_mat <- mass2mz(lcms_mzs, paste0("[M]", ion_sign))
        colnames(orig_adduct_mat) <- paste0("[X]", ion_sign, "/[X]", ion_sign)

        ### Get peaks for the lcms mzs. We assume all of them are adducts,
        ### We will also use them to calculate the original
        #orig_adduct_mat <- EstimateAdductsFromAdducts(lcms_mzs, ion_type)
        #colnames(orig_adduct_mat) <- gsub("M", "X", colnames(orig_adduct_mat))

        ### Generate the adducts for all the exact mzs
        exact_adduct_mat <- EstimateAdductsFromExactMzs(
            lcms_mzs, exact_mzs, ion_type
        )

        ### Combine them together
        cbind(orig_adduct_mat, exact_adduct_mat)

    } else {
        ### Estimate all possible adducts based on adducts
        EstimateAdductsFromAdducts(lcms_mzs, ion_type)
    }

    ### Rearrange as a long array
    ### Note: adduct_mat may contain NAs
    adduct_mzs <- NULL
    for (mz_id in rownames(adduct_mat)) {
        tmp <- adduct_mat[mz_id, ]
        names(tmp) <- paste0(mz_id, ":", names(tmp))
        tmp <- tmp[!is.na(tmp)]
        adduct_mzs <- c(adduct_mzs, tmp)
    }

    ### Make sure it is sorted
    adduct_mzs <- sort(adduct_mzs)

    ### Remove collided ion
    adduct_mzs <- FixMZCollisions(adduct_mzs)

    ### Return
    adduct_mzs
}

### "========================================================================="
### Estimate mzs of adducts from exact mzs
### * Note that different peaks may have different number of exact mzs
###   so we need to use a loop
### "========================================================================="
EstimateAdductsFromExactMzs <- function(
    lcms_mzs, exact_mzs, ion_type
) {
    ### What are the possible adducts for the target?
    target_adducts <- GetDESIIons(ion_type)

    ### Get the ion_sign character
    ion_sign <- ifelse(ion_type == "pos", "+", "-")

    ### Find the maximum number of lcms_mzs
    max_exact_num <- max(sapply(exact_mzs, length))

    ### Generate the adduct_names
    adduct_labels <- NULL
    for (idx in 1:max_exact_num) {
        label_cur <- paste0("[M]", ion_sign, "/", target_adducts)
        adduct_labels <- c(
            adduct_labels,
            gsub("M", paste0("M", idx), label_cur)
        )
    }

    ### Init the mat
    exact_adduct_mat <- matrix(
        NA,
        nrow = length(exact_mzs),
        ncol = length(adduct_labels),
        dimnames = list(names(exact_mzs), adduct_labels)
    )

    ### Loop through all peaks
    for (peak_id in names(exact_mzs)) {
        ### Generate adduct for all the exact mzs
        adduct_cur <- mass2mz(exact_mzs[[peak_id]], target_adducts)

        ### Convert to a vector row by row
        adduct_cur_mzs <- as.vector(t(adduct_cur))

        ### Find out which adducts are dissimilar from the lcms mzs
        similar_idx <- abs(adduct_cur_mzs - lcms_mzs[peak_id]) <= 1.5

        ### We don't want those similar one
        adduct_cur_mzs[similar_idx] <- NA

        ### Save the data
        exact_adduct_mat[
            peak_id, seq_along(adduct_cur_mzs)
        ] <- adduct_cur_mzs
    }

    exact_adduct_mat
}

### "========================================================================="
### Estimate mzs of adducts from adducts
### "========================================================================="
EstimateAdductsFromAdducts <- function(
    lcms_mzs, ion_type
) {
    ### What are the possible adducts for the source?
    source_adducts <- GetLCMSIons(ion_type)

    ### What are the possible adducts for the target?
    target_adducts <- GetDESIIons(ion_type)

    ### Estimate all possible neutral masses
    neutral_masses <- mz2mass(lcms_mzs, source_adducts)

    ### Estimate all possible adducts for the masses
    adducts <- mass2mz(neutral_masses, target_adducts)

    ### Collate the data
    adduct_mat <- NULL
    for (target_adduct in target_adducts) {
        ### Get the current adducts
        cur_adducts <- adducts[, , target_adduct]

        ### Remove the same ions, except for [M]+ and [M]-
        if ((target_adduct != "[M]+") && (target_adduct != "[M]-")) {
            cur_adducts <- cur_adducts[
                , colnames(cur_adducts) != target_adduct,
                drop = FALSE
            ]
        }

        colnames(cur_adducts) <- paste0(
            colnames(cur_adducts), "/", target_adduct
        )

        adduct_mat <- cbind(adduct_mat, cur_adducts)
    }

    adduct_mat
}

### "========================================================================="
### Check and fix mz collisions
### "========================================================================="
FixMZCollisions <- function(
    adduct_mzs
) {
    ### Get all the unique mzs
    unique_adduct_mzs <- unique(adduct_mzs)

    remove_names <- NULL

    ### Check any collision
    for (unique_adduct_mz in unique_adduct_mzs) {
        is_same <- adduct_mzs == unique_adduct_mz
        if (sum(is_same) > 1) {
            ### We just retained the first one
            collided_idx <- which(is_same)
            collided_names <- names(adduct_mzs)[is_same]

            message(
                "[Warning] mz collisions found for mz = ",
                unique_adduct_mz, " (",
                paste0(collided_names, collapse = ", "),
                ")"
            )

            ### Try to retain the first X
            is_X_found <- grepl("\\[X\\]", collided_names)

            if (any(is_X_found)) {
                ### Only keep the first one
                keep_idx <- which.first(is_X_found)
                remove_names <- c(remove_names, collided_names[-keep_idx])
                message(
                    "Only retained the LC/MS adducts = ",
                    names(adduct_mzs)[collided_idx[keep_idx]]
                )
            } else {
                message(
                    "Only retained the first adducts = ",
                    names(adduct_mzs)[collided_idx[1]]
                )
                remove_names <- c(remove_names, collided_names[-1])
            }
        }
    }

    if (!is.null(remove_names)) {
        adduct_mzs <- adduct_mzs[!(names(adduct_mzs) %in% remove_names)]
    }

    ### return the value
    adduct_mzs
}

### "========================================================================"
### Summarize the peak abundances
### "========================================================================"
SummarizePeakAbd <- function(
    reg_tissue_peaks, ion_type, summary_type
) {
    ### Make sure summary type is valid
    if (!(
        summary_type %in% c(
            "lcms_only",
            "lcms_or_all",
            "pred_or_all",
            "all"
        )
    )) {
        stop("Unknown summary type - ", summary_type)
    }

    ### Get the peak subset
    ### Note: peak_ids may be shorter than peak_names
    peak_idx_list <- GetPeakSubset(reg_tissue_peaks, ion_type, summary_type)
    peak_ids <- names(peak_idx_list)

    ### We need to summarize according to peak_ids
    msi_all <- matrix(
        0,
        nrow = length(peak_ids),
        ncol = ncol(reg_tissue_peaks),
        dimnames = list(
            peak_ids,
            colnames(reg_tissue_peaks)
        )
    )

    for (peak_id in peak_ids) {
        ### Get the current peak_idx
        peak_idx <- peak_idx_list[[peak_id]]

        if (sum(peak_idx) > 0) {
            msi_all[peak_id, ] <- colSums(
                reg_tissue_peaks[peak_idx, ,drop = FALSE]
            )
        }
    }

    msi_all
}

### "========================================================================"
### Get the peak subset based on the peak_summary_mode
### A list of index for each peak_id will be returned.
### "========================================================================"
GetPeakSubset <- function(
    reg_tissue_peaks, ion_type, peak_summary_mode
) {
    ### Get the current row names
    peak_names <- rownames(reg_tissue_peaks)

    ### Get the parent compound only name
    parent_ion  <- ifelse(ion_type == "pos", "[M]+/[M]+", "[M]-/[M]-")
    lcms_adduct <- ifelse(ion_type == "pos", "[X]+/[X]+", "[X]-/[X]-")

    ### Check if is adduct
    is_adduct <- any(grepl(lcms_adduct, peak_names, fixed = TRUE))

    ### Get the peak_ids
    ### Note: peak_ids may be shorter than peak_names, because peak_names
    ### may include adducts!! peak_ids are always just the parents.
    if (is_adduct) {
        ### Get the peak_ids from the lcms adducts
        lcms_adduct_names <- peak_names[
            sub("^[^:]+:(.+)$", "\\1", peak_names) == lcms_adduct
        ]
        peak_ids <- sub("^([^:]+):.+$", "\\1", lcms_adduct_names)

    } else {
        peak_ids <- peak_names

    }

    ### Get the list peak subsets
    peak_idx_list <- list()

    for (peak_id in peak_ids) {

        ### If this an adduct?
        if (is_adduct) {
            ### Generate the common logical arrays
            is_lcms_adduct <- (
                peak_names == paste0(peak_id, ":", lcms_adduct)
            )
            is_pred_adduct <- startsWith(
                peak_names, paste0(peak_id, ":[M")
            )
            is_all <- startsWith(
                peak_names, paste0(peak_id, ":")
            )
            is_parent_ion <- (
                peak_names == paste0(peak_id, ":", parent_ion)
            )

            ### Which feature to be included?
            peak_idx_list[[peak_id]] <- if (
                peak_summary_mode == "lcms_only"
            ) {
                ### LC/MS adducts only
                is_lcms_adduct
            } else if (peak_summary_mode == "pred_only") {
                ### LC/MS adducts only
                is_pred_adduct
            } else if (peak_summary_mode == "lcms_or_all") {
                ### Prefer lcms adduct if it is available
                if (sum(reg_tissue_peaks[is_lcms_adduct, ]) > 0) {
                    is_lcms_adduct
                } else {
                    is_all
                }
            } else if (peak_summary_mode == "pred_or_all") {
                ### Prefer predicted adduct if it is available
                if (sum(reg_tissue_peaks[is_pred_adduct, ]) > 0) {
                    is_pred_adduct
                } else {
                    is_all
                }
            } else {
                ### all
                is_all
            }
        } else {
            ### For other type of peaks
            peak_idx_list[[peak_id]] <- (peak_names == peak_id)
        }
    }

    ### Return
    peak_idx_list
}

### "========================================================================"
### Get the raw peaks
### "========================================================================"
GetRawPeaks <- function(
    ion_type, msi_conf, peak_type, peak_summary_mode,
    stat_type = "tissue"   ## roi or tissue
) {
    ### Get the list of current samples
    samples <- GetSamples(msi_conf)
    cur_samples <- samples[[ion_type]]

    ### Init output
    all_tissue_peaks <- list()

    ### Loop through all the samples
    for (cur_sample in cur_samples) {
        message("Processing sample ", cur_sample$sample_id)

        ### Load tissue_peaks and peak_mzs or peaks
        ### Note: In (older) common peaks, peaks is returned, but it newer
        ###       ref_peaks, peak_mzs is returned.
        load(file = here(
            "data", "msi",
            peak_type, cur_sample$ion_type,
            paste0(cur_sample$sample_id, ".Rdata")
        ))

        ### Summarize the adducts for each region
        for (reg_name in names(tissue_peaks)) {
            ### The reg_tissue_peaks
            reg_tissue_peaks <- tissue_peaks[[reg_name]]

            ### Get peak subset according to peak_summary mode
            peak_idx_list <- GetPeakSubset(
                reg_tissue_peaks,
                ion_type,
                peak_summary_mode
            )

            ### Combine all the indices
            peak_idx <- Reduce("|", peak_idx_list)
            x <- reg_tissue_peaks[peak_idx, , drop = FALSE]

            ### Cannot be negative
            x[x < 0] <- 0

            ### Transform
            x <- log2(x + 1)

            ### Save the data
            all_tissue_peaks[[cur_sample$sample_id]][[reg_name]] <- x
        }
    }

    ### return the results
    all_tissue_peaks
}

### "========================================================================"
### Get the mean abundance levels for a peak on each tissues
### * peak_type =
###      - top_common_peaks = common peaks for the whole tissue regions
###      - ref_peaks = lcms peaks for the whole tissue regions
### * stat_type =
###      - roi = roi regions
###      - tissue = tissue sections
### * peak_summary_mode (only applicable to ref_peaks)
###      - lcms_only
### "========================================================================"
GetPeakStats <- function(
    ion_type, msi_conf, peak_type, peak_summary_mode,
    min_abd = -Inf,
    stat_type = "tissue",   ## roi or tissue
    is_remove_bad_cases = TRUE
) {
    ### Get the list of current samples
    samples <- GetSamples(msi_conf)
    cur_samples <- samples[[ion_type]]

    study_stats <- NULL

    ### Loop through all the samples
    for (cur_sample in cur_samples) {
        message("Processing sample ", cur_sample$sample_id)

        ### Load tissue_peaks and peak_mzs or peaks
        load(file = here(
            "data", "msi",
            peak_type, cur_sample$ion_type,
            paste0(cur_sample$sample_id, ".Rdata")
        ))

        ### Summarize the adducts for each region
        for (reg_name in names(tissue_peaks)) {
            tissue_peaks[[reg_name]] <- SummarizePeakAbd(
                reg_tissue_peaks = tissue_peaks[[reg_name]],
                ion_type = ion_type,
                summary_type = peak_summary_mode
            )
        }

        ### Get the tissue section ids
        section_ids <- sapply(
            names(tissue_peaks),
            function(x) {
                str_split(x, "-")[[1]][[1]]
            }
        )

        ### Load the regions and loop through it
        ### There are two types of regions.
        ### - freehand are the tissue section regions
        ### - rect are the ROI regions
        regs <- read.ijzip(here(cur_sample$roi_path), names = TRUE)

        cur_stats <- list()

        for (reg_name in names(regs)) {
            ### Get the current reg and also tissue_peaks
            reg <- regs[[reg_name]]

            ### Find the tissue section name (T1, F10, N4, etc)
            reg_id <- str_split(reg_name, "-")[[1]][1]

            ### The region may be a section or ROI, thus reg_id maybe
            ### a roi_id or section_id. Nevertheless, we need to find 
            ### the section_id of the region.

            ### Get the section_idx
            section_idx <- sub('.', '', reg_id)

            ### First, we need to find out if this is from a normal section
            if (section_idx == '') {
                ### Region without a section idx comes from the normal section
                section_id <- 'N'

            } else if (section_idx == '0') {
                ### Region with section idx = 0 comes from the normal section
                section_id <- 'N'

            } else {
                ### All other regions comes from the TX section
                section_id <- sub(".", "T", reg_id)
            }

            ### section id must be valid
            if (!(section_id %in% section_ids)) {
                message("Region name = ", reg_name)
                print(section_ids)
                stop(
                    "[Error] section_id (", section_id,
                    ") cannot be found in tissue_peaks"
                )
            }

            ### Get tissue_abd
            section_name <- names(section_ids)[
                section_id == section_ids
            ]
            tissue_abd <- tissue_peaks[[section_name]]

            ### Find out what type is the region?
            if (reg$strType == "rect" && stat_type == "roi") {
                ### We will need to load the ROI
                ### Since ROI is a rect, we just need to check the bound
                ### Get the coordinates of all the pixels
                region_coord <- sapply(
                    str_split(colnames(tissue_abd), "_"),
                    function(a) {
                        c(x = as.numeric(a[2]), y = as.numeric(a[3]))
                    }
                )

                good_idx <-
                    (region_coord["x", ] >= reg$xrange[1]) &
                        (region_coord["x", ] <= reg$xrange[2]) &
                        (region_coord["y", ] >= reg$yrange[1]) &
                        (region_coord["y", ] <= reg$yrange[2])

                # Make sure the size is correct
                if (
                    sum(good_idx) !=
                        (diff(reg$xrange) + 1) * (diff(reg$yrange) + 1)
                ) {
                    stop("good_idx has invalid size")
                }

                # We subset the tissue_abd
                tissue_abd <- tissue_abd[, good_idx, drop = FALSE]
            } else if (reg$strType == "freehand" && stat_type == "tissue") {
                ### We will use all the levels

            } else {
                next
            }

            ### Define the tissue region name
            tissue_region_name <- paste0(
                cur_sample$sample_id, ".", reg_name
            )

            message(paste0("Processing ", tissue_region_name, " ..."))

            ### Calculate the stats
            stat_cur <- t(apply(
                tissue_abd, 1,
                function(x) {
                    ### Cannot be negative
                    x[x < 0] <- 0

                    ### Value must be larger than min_abd
                    if (min_abd > -Inf) {
                        x[x < min_abd] <- 0
                    }

                    coverage <- sum(x > 0, na.rm = TRUE) / length(x)

                    ### Transform
                    x <- log2(x + 1)

                    ### Calculate the stats
                    c(
                        msi_sum = sum(x, na.rm = TRUE),
                        msi_median = median(x, na.rm = TRUE),
                        msi_mean = mean(x, na.rm = TRUE),
                        msi_sd = sd(x, na.rm = TRUE),
                        msi_pct90 = quantile(
                            x, 0.9,
                            na.rm = TRUE, names = FALSE
                        ),
                        msi_coverage = coverage
                    )
                }
            ))

            peak_num <- nrow(stat_cur)

            ### Return the results as a tibble
            if (exists("peak_mzs")) {
                cur_stats[[reg_name]] <- bind_cols(
                    tibble(
                        sample_id = rep(cur_sample$sample_id, peak_num),
                        reg_name = rep(reg_name, peak_num),
                        peak_id = rownames(stat_cur)
                    ),
                    as_tibble(stat_cur)
                )

            } else {
                ### Cannot remove this, mz needed by top_common_peaks
                cur_stats[[reg_name]] <- bind_cols(
                    tibble(
                        sample_id = rep(cur_sample$sample_id, peak_num),
                        reg_name = rep(reg_name, peak_num),
                        peak_id = rownames(stat_cur),
                        mz = peaks
                    ),
                    as_tibble(stat_cur)
                )
            }
        }

        study_stats <- bind_rows(
            study_stats,
            do.call(rbind, cur_stats)
        )
    }

    ### Remove the bad samples
    if (is_remove_bad_cases) {
        study_stats <- study_stats %>% filter(
            !(sample_id %in% bad_sample_ids)
        )
    }

    study_stats
}

### "========================================================================="
### Load the LCMS raw data
### mz_names = the names of the lcms features to be loaded. Please note that
### this must include the retention time info.
### "========================================================================="
LoadLCMSRawData <- function(
    lcms_mode,
    mz_names = NULL
) {
    ### Define the modes
    if (!(lcms_mode %in% c("LPOS", "HPOS", "LNEG"))) {
        stop("Unknown LC/MS mode")
    }

    ### Load all the raw data
    readRDS(here("data", "lcms", "lcms_raw.rds")) |>
        filter(endsWith(mz_id, lcms_mode))
}

### "========================================================================="
### Load the selected LCMS peak data for the specific LCMS mode
### 
### "========================================================================="
LoadLCMSHitData <- function(
    lcms_mode,
    is_ctrl_missing_ok = FALSE,
    min_abd = 10
) {
    ### Load the data
    # load(here("data", "lcms", "lcms_peaks.Rdata"))
    ion_type <- tolower(substr(lcms_mode, 2, 4))

    ### Construct the hit_table
    #peak_table <- tibble(
    #    mz = peak_list$lcms_mzs,
    #    peak_id = names(peak_list$lcms_mzs),
    #    ion_type = !! ion_type
    #)

    ### Load the data
    load(here("data", "lcms", "lcms_peaks.Rdata"))

    ### Please note that not all mzs are uniques in LC/MS data
    peak_table <- lcms_hit_peaks %>%
        filter(`ion_type` == !!ion_type) %>%
        arrange(peak_id) %>%
        dplyr::select(peak_id, mz, mz_id)

    ### Load the raw data
    raw_data <- LoadLCMSRawData(lcms_mode)
    
    ### Only get the selected peaks
    hit_data <- raw_data %>%
        left_join(peak_table, by = c("mz_id")) %>%
        filter(
            !is.na(peak_id)
        ) %>%
        group_by(
            case_id, section_id, peak_id, mz
        ) %>%
        summarize(
            ### Take the mean of all the mz related to the peaks!
            ### because we don't know which lcms mz is the real one
            lcms_abd = mean(lcms_abd),
            .groups = "drop"
        )

    ### Check min_value
    #min_abd_cur <- min(hit_data$lcms_abd, na.rm = TRUE)
    #if (min_abd_cur < min_abd) {
    #    message("Current min value = ", min_abd_cur)
    #    message("Set min value     = ", min_abd)
    #    stop("[Error] current min values is smaller than set min value")
    #}

    if (!is_ctrl_missing_ok) {
        hit_data <- hit_data %>%
            mutate(
                ### We will not analyze LC/MS without a control data
                lcms_abd = replace(
                    lcms_abd,
                    startsWith(section_id, "N") & (lcms_abd == 0),
                    NA
                )
            )
    }

    ### Do the final normalization
    hit_data <- hit_data %>%
        ### Values that are too small are likely just noise
        mutate(
            ### Values that are too small are likely just noise
            lcms_abd = replace(
                lcms_abd, lcms_abd < min_abd, min_abd
            )
        ) %>%
        ### Normalize the data
        group_by(case_id, peak_id) %>%
        mutate(
            rel_lcms_abd =
                lcms_abd - dplyr::first(lcms_abd[startsWith(section_id, "N")])
        ) %>%
        ungroup()

    hit_data
}

### Top Prominent Peak Analysis ==============================================
### "========================================================================="
### Get refpeak filename
### "========================================================================="
GetProminentPeaksFilePath <- function(case_conf) {
    filepath <- here(
        "data", "msi", "prominent_peaks",
        paste0(case_conf$sample_id, "_promi.Rdata")
    )
    return(filepath)
}

### "========================================================================"
### Merge prominent peaks
###
### https://rformassspectrometry.github.io/MsCoreUtils/reference/group.html
### https://rformassspectrometry.github.io/Spectra/reference/combinePeaks.html
### "========================================================================"
FindCommonProminentPeaks <- function(
    case_confs, align_res
) {
    # align_res <- 5
    ### Load all the peaks
    case_promi_peaks <- NULL

    for (case_conf in case_confs) {
        ### Load the peaks
        filepath <- GetProminentPeaksFilePath(case_conf)
        load(filepath)

        ### Find prominent peaks that can be found in all samples
        cur_peaks <- promi_peaks[promi_sample_num == max(promi_sample_num)]
        names(cur_peaks) <- rep(case_conf$sample_id, length(cur_peaks))

        case_promi_peaks <- c(case_promi_peaks, cur_peaks)
    }

    ### Align the prominent peaks across all the samples
    case_promi_peaks <- sort(case_promi_peaks)
    matched_idx <- MsCoreUtils::group(case_promi_peaks, ppm = align_res)
    common_promi_peak_list <- split(case_promi_peaks, matched_idx)

    ### Align the peaks
    common_promi_peaks <- sapply(common_promi_peak_list, base::mean)
    common_promi_sample_num <- sapply(common_promi_peak_list, length)

    ### Return the results
    list(
        common_promi_peaks = common_promi_peaks,
        common_promi_sample_num = common_promi_sample_num
    )
}

### "========================================================================="
### Get top common tissue prominent peaks
### * Common peaks that can be found in all ROIs in (sample_num - 1) samples
### "========================================================================="
FindTopCommonProminentPeaks <- function(
    ion_type, align_res, msi_conf
) {
    ### Identify all the samples
    samples <- GetSamples(msi_conf)

    ### Get all the results
    res <- FindCommonProminentPeaks(
        case_confs = samples[[ion_type]],
        align_res = align_res
    )

    ### Find the top common prominent peaks in at least x - 1
    common_promi_peaks <- res$common_promi_peaks[
        res$common_promi_sample_num >= (length(samples[[ion_type]]) - 1)
    ]

    ### return
    common_promi_peaks
}

### "========================================================================"
### Plot top common peaks
### "========================================================================"
PlotTopCommonTissuePeaks <- function(
    ion_type, align_res, peak_summary_mode, msi_conf
) {
    ### Identify all the samples
    samples <- GetSamples(msi_conf)

    #### Load the tissue common peaks
    study_stats <- GetPeakStats(
        ion_type = ion_type,
        msi_conf = msi_conf,
        peak_type = "top_common_peaks",
        peak_summary_mode = peak_summary_mode,
        stat_type = "tissue"
    )

    ### Have we already defined the top common peaks?
    is_selected <- factor(
        x = rep("No", length(unique(study_stats$mz))),
        levels = c("Yes", "No")
    )

    if (
        exists("top_common_peak_idx") &&
        !is.null(top_common_peak_idx[[ion_type]])
    ) {
        ### Update the colors
        is_selected[top_common_peak_idx[[ion_type]]] <- "Yes"

        ### Remove the labels
    }

    ### Generate the summary plot
    plot_data <- study_stats %>%
        mutate(
            CV = msi_sd / msi_mean,
            region_id = paste0(sample_id, ".", reg_name)
        ) %>%
        group_by(
            mz
        ) %>%
        summarise(
            s_mean = mean(msi_mean),
            s_CV   = mean(msi_coverage) * 100,
            .groups = "drop"
        ) %>% mutate(
            label = signif(mz, 4),
            is_selected = is_selected
        )

    if (any(is_selected == "Yes")) {
        plot_data <- plot_data %>% mutate(
            label = replace(label, is_selected != "Yes", "")
        )
    }

    p <- plot_data %>%
        ggplot(
            aes(
                x = s_mean,
                y = s_CV,
                label = label,
                col = is_selected
            )
        ) +
        geom_point() +
        geom_text_repel(
            max.time = 5,
            box.padding = 0.75,
            size = 2.75,
            max.overlaps = Inf
        ) +
        scale_y_continuous(
            labels = c("40", "", "60", "", "80", "", "100"),
            breaks = c(40, 50, 60, 70, 80, 90, 100),
            limits = c(35, 110)
        ) +
        scale_color_manual(
            #guide="none",
            name = NULL,
            values = c("Yes" = "#2b8cbe", "No" = "black"),
            labels = c("Yes" = "Reference CPPs", "No" = "Other CPPs")
        ) +
        ylab("Average section coverage (%)") +
        xlab("Average section mean abundance level (log2)") +
        theme_classic() +
        theme(
            legend.position = c(0.7, 0.2),
            axis.text = element_text(colour="black")
        )

    ### save the plots
    for (img_type in c(".png", ".pdf")) {
        ggsave(
            filename = here(
                "figures", "peak_normalization",
                paste0("top_common_peak_stat_", ion_type, img_type)
            ),
            plot = p,
            width = 3.75,
            height = 2.75,
            units = "in",
            dpi = 300
        )
    }
}

### "========================================================================"
### Load the normalization coefficient
### "========================================================================"
LoadNormCoeff <- function(ion_type) {
    readRDS(
        file = here("results", paste0("norm_coeff_", ion_type, ".rds"))
    )
}

### "========================================================================"
### Find the normalization coefficient
### is_plot = plot the mean levels of top common peaks for all samples
### "========================================================================"
FindNormCoeff <- function(
    ion_type, peak_summary_mode, msi_conf, is_plot = FALSE
) {
    samples <- GetSamples(msi_conf)

    common_promi_peaks <- FindTopCommonProminentPeaks(
        ion_type = ion_type,
        align_res = top_common_align_reses[[ion_type]],
        msi_conf
    )
    common_promi_peaks <- common_promi_peaks[top_common_peak_idx[[ion_type]]]

    ### Form the image ID
    sample_id <- sapply(msi_conf[["samples"]], function(x) x$sample_id)
    image_ids <- paste0(
        sapply(msi_conf[["samples"]], function(x) x$case_id), "-",
        substring(sample_id, 9, 9)
    )
    names(image_ids) <- sample_id

    study_stats <- GetPeakStats(
            ion_type = ion_type,
            msi_conf = msi_conf,
            peak_type = "top_common_peaks",
            peak_summary_mode = peak_summary_mode,
            stat_type = "tissue",
            is_remove_bad_cases = FALSE
        ) %>%
        mutate(
            sample_id = as_factor(sample_id),
            region_id = paste0(sample_id, ".", reg_name)
        ) %>%
        filter(mz %in% common_promi_peaks) %>%
        group_by(sample_id, reg_name) %>%
        summarize(section_mean = mean(msi_mean), .groups = "drop")

    ### Do we want to plot?
    if (is_plot) {

        sample_means <- study_stats %>%
            group_by(sample_id) %>%
            summarize(sample_mean = mean(section_mean), .groups = "drop") %>%
            dplyr::pull(sample_mean, name = sample_id)

        sample_stats <- tibble(
            median = median(sample_means),
            pct25 = quantile(sample_means, 0.25),
            pct75 = quantile(sample_means, 0.75),
            thres_lo = pct25 - 1.5*(pct75 - pct25),
            thres_hi = pct75 + 1.5*(pct75 - pct25)
        )

        cols <- rep("#666666", length(sample_means))
        names(cols) <- names(sample_means)
        cols[sample_means < sample_stats$thres_lo] <- "#FF0000"
        cols[sample_means > sample_stats$thres_hi] <- "#FF0000"

        p <- study_stats %>%
            #mutate(sample_id = image_ids[sample_id]) %>%
            # group_by(sample_id) %>%
            # summarize(
            #    sample_mean = mean(section_mean)
            # ) %>%
            ggplot(
                aes(x = sample_id, y = section_mean)
            ) +
            geom_bar(
                position = "dodge", stat = "summary", fun = "mean",
                fill = "#aaaaaa", width = 0.75
            ) +
            scale_x_discrete(
                labels = image_ids
            ) +
            # geom_errorbar(stat = 'summary', position = 'dodge', width = 0.9) +
            geom_hline(data = sample_stats, aes(yintercept = median)) +
            geom_hline(
                data = sample_stats, aes(yintercept = pct25),
                linetype = "dashed"
            ) +
            geom_hline(
                data = sample_stats, aes(yintercept = pct75),
                linetype = "dashed"
            ) +
            geom_hline(
                data = sample_stats, aes(yintercept = thres_lo),
                linetype = "dotted"
            ) +
            geom_hline(
                data = sample_stats, aes(yintercept = thres_hi),
                linetype = "dotted"
            ) +
            geom_point(
                aes(x = sample_id),
                shape = 21, fill = "white", size = 2
            ) +
            #ylim(c(0, 12)) +
            xlab("DESI-MSI images ") +
            ylab("Mean CPP\nabundance (log2)") +
            theme_classic() +
            theme(
                axis.text = element_text(colour="black"),
                axis.text.x = element_text(
                    angle = 45, hjust = 1 #vjust = 0.5, 
                )
            )

        ### save the plots
        for (img_type in c(".png", ".pdf")) {
            ggsave(
                filename = here(
                    "figures", "peak_normalization",
                    paste0("mean_sample_level_", ion_type, img_type)
                ),
                plot = p,
                width = length(unique(study_stats$sample_id))/4,
                height = 2.75,
                units = "in",
                dpi = 300
            )
        }
    }

    ### Determine the normalization coefficients
    norm_coeff <- study_stats %>%
        group_by(sample_id) %>%
        summarize(sample_mean = mean(section_mean), .groups = "drop") %>%
        mutate(norm_coeff = sample_mean) %>%
        dplyr::pull(norm_coeff, name = sample_id)

    ### Save the data
    saveRDS(
        norm_coeff,
        file = here("results", paste0("norm_coeff_", ion_type, ".rds")),
        compress = "xz"
    )

    ## return
    norm_coeff
}

### Marker analysis ==========================================================
### "========================================================================="
### Compute metabolites
### is_log2 = TRUE or FALSE
### ion_type = "pos" or "neg"
### "========================================================================="
AnalyzeMetabolites <- function(
    ion_type, slide_ids, msi_conf, peak_summary_mode, is_normalized
) {
    ### Make sure slide_ids is sorted
    slide_ids <- sort(slide_ids)

    ### Load the data
    data_all <- LoadMSIMat(
        ion_type  = ion_type,
        slide_ids = slide_ids,
        peak_type = "ref_peaks",
        peak_summary_mode = peak_summary_mode,
        stat_type = "roi",
        msi_conf  = msi_conf,
        is_normalized = is_normalized
    )

    ### Find ROI types
    roi_types <- FindROITypes(data_all)

    ### Get the mz values
    peak_ids <- rownames(data_all)

    ### Subtract peaks from non-transformed regions from other ROIs
    #for (slide_id in roi_types$slide_ids) {
    #    ### Get all regions from the slide
    #    is_slide <- grepl(paste0("^", slide_id, "_.*$"), colnames(data_all))
    #
    #    ### There must be normal
    #    if (!any(roi_types$is_nontransformed & is_slide)) {
    #        stop("Cannot find control in ", slide_id)
    #    }
    #
    #    ### Find the average value from the normal roi (in log2 unit)
    #    normal_mean <- rowMeans(
    #        data_all[, roi_types$is_nontransformed & is_slide, drop = FALSE]
    #    )
    #
    #    ### Subtract the normal mean from all the log2-transformed data_all
    #    data_all[, is_slide] <- sweep(
    #        data_all[, is_slide, drop = FALSE], 1,
    #        normal_mean
    #    )
    #}

    ### Perform a two-sided one-sample t-test for each metabolite under
    ### different groups of ROIs
    ### because the value should be zero under null hypothesis

    ### Perform t-test
    TestROIDiff <- function(is_grp, is_ctrl) {
        res <- if ((sum(is_grp) > 0)  && (sum(is_ctrl) > 0)) {

            ### Make sure is_grp do not overlap with is_ctrl
            is_ctrl[is_grp] <- FALSE

            ### Return the result
            apply(
                data_all, 1,
                function(x) {
                    ### If the values are all zero but the normal is not zero,
                    ### x[is_grp] will be all the same. So, we try to detect
                    ### this percular case.
                    if (
                        all(is.na(x[is_ctrl])) ||
                        all(is.na(x[is_grp])) ||
                        (sd(x[is_grp], na.rm = TRUE) == 0 ||
                         sd(x[is_ctrl], na.rm = TRUE) == 0)
                    ) {
                        c(
                            pval = NaN,
                            log2FC = NaN
                        )
                    } else {
                        test_res <- t.test(x[is_grp], x[is_ctrl])
                        c(
                            pval = test_res$p.value,
                            log2FC = test_res$estimate[1] - test_res$estimate[2]
                        )
                    }
                }
            )

            #p.adjust(pval, method = "fdr")

        } else {
            matrix(NaN, 2, nrow(data_all))
        }

        rownames(res) <- c("pval", "log2FC")

        res
    }

    transformed_test <- TestROIDiff(
        roi_types$is_transformed,
        roi_types$is_nontransformed
    )

    fibrotic_test <- TestROIDiff(
        roi_types$is_fibrotic,
        #roi_types$is_normal
        roi_types$is_nontransformed   ### Could be S or F
    )

    necrotic_test <- TestROIDiff(
        roi_types$is_necrotic,
        roi_types$is_nontransformed
    )

    steatotic_test <- TestROIDiff(
        roi_types$is_steatotic,
        #roi_types$is_normal
        roi_types$is_nontransformed  ## Could be S or F
    )

    diff_grade2_test <- TestROIDiff(
        roi_types$is_diff_grade2,
        roi_types$is_nontransformed
    )

    diff_grade3_test <- TestROIDiff(
        roi_types$is_diff_grade3,
        roi_types$is_nontransformed
    )

    ### Do a two-sample t-test for grade 2 vs grade 3
    diff_grade_test <- TestROIDiff(
        roi_types$is_diff_grade2,
        roi_types$is_diff_grade3
    )

    ### We will not do the adjustment now, but do it at the end
    #tumor_grades_padj <- p.adjust(tumor_grades_pval, method = "fdr")

    ### Get the mean values for all the slides
    GetROI_mean <- function(is_grp) {
        list(
            all = apply(
                data_all, 1, function(x) mean(x[is_grp])
            )
        )
    }

    nontransformed_mean <- GetROI_mean(roi_types$is_nontransformed)
    transformed_mean <- GetROI_mean(roi_types$is_transformed)
    fibrotic_mean <- GetROI_mean(roi_types$is_fibrotic)
    necrotic_mean <- GetROI_mean(roi_types$is_necrotic)
    steatotic_mean <- GetROI_mean(roi_types$is_steatotic)

    diff_grade2_mean <- GetROI_mean(roi_types$is_diff_grade2)
    diff_grade3_mean <- GetROI_mean(roi_types$is_diff_grade3)

    ### Form the results
    metabolites <- tibble(
        peak_id = peak_ids,

        transformed_pval = transformed_test["pval", ],
        fibrotic_pval = fibrotic_test["pval", ],
        necrotic_pval = necrotic_test["pval", ],
        steatotic_pval = steatotic_test["pval", ],

        diff_grade_pval = diff_grade_test["pval", ],
        diff_grade2_pval = diff_grade2_test["pval", ],
        diff_grade3_pval = diff_grade3_test["pval", ],

        transformed_log2FC = transformed_test["log2FC", ],
        fibrotic_log2FC = fibrotic_test["log2FC", ],
        necrotic_log2FC = necrotic_test["log2FC", ],
        steatotic_log2FC = steatotic_test["log2FC", ],

        diff_grade_log2FC = diff_grade_test["log2FC", ],
        diff_grade2_log2FC = diff_grade2_test["log2FC", ],
        diff_grade3_log2FC = diff_grade3_test["log2FC", ]
    ) 

    metabolites$nontransformed_mean <- nontransformed_mean[["all"]]
    metabolites %>% arrange(
        transformed_pval
    )
}

### Graphing plotting ========================================================
### "========================================================================="
### Plot marker to pathological score comparison
### "========================================================================="
PlotSgMEvsPathMarkers <- function(
    lcms_subset, info_cat, score_name_in,
    x_label, y_label, y_lim = NULL,
    test_method = "wilcox.test",
    test_alt = "greater",
    test_comparisons = NULL,
    is_label = FALSE,
    is_save = FALSE,
    img_name = NULL
) {
    ### Prepare the data
    plot_data <- inner_join(
        lcms_subset %>%
            filter(mz %in% all_mkrs$mz, !is.na(.data[[info_cat]])),
        all_mkrs %>%
            filter(`score_name` == score_name_in) %>%
            dplyr::select(mz, mkr_name, norm_coeff, score_name),
        by = "mz"
    )

    ### Filter HCV data
    if (info_cat == "hbv_status") {
        plot_data <- plot_data %>% filter(case_id != "HEP261G")
    }

    plot_data <- plot_data %>% group_by(
        case_id, score_name, mkr_name
    ) %>% summarize(
        norm_T_score = T_mean / norm_coeff[which.max(abs(norm_coeff))],
        norm_N_score = N_mean / norm_coeff[which.max(abs(norm_coeff))],
        .groups = "drop"
    )

    ### Sum tumors
    plot_data <- plot_data %>% group_by(
        case_id, score_name
    ) %>% summarize(
        final_T_score = mean(norm_T_score, na.rm = TRUE),
        final_N_score = mean(norm_N_score, na.rm = TRUE),
        .groups = "drop"
    ) %>% filter(
        !is.na(final_T_score) | !is.na(final_N_score)
    )

    plot_data <- inner_join(
        plot_data,
        case_cats,
        by = c("case_id")
    )

    ### Find the control info
    final_N_scores <- plot_data %>% dplyr::pull(final_N_score) 
    ctr_mean <- mean(final_N_scores, na.rm = TRUE)

    ### Do a t-test to control
    score_data <- NULL
    cat_data <- NULL
    for (idx in 1:nrow(plot_data)) {
        ### Add the final T score
        #if (!is.na(plot_data[["final_T_score"]][idx])) {
            cat_data <- c(cat_data, as.character(plot_data[[info_cat]][idx]))
            score_data <- c(score_data, plot_data[["final_T_score"]][idx])
        #}
        ### Add the final N score
        #if (!is.na(plot_data[["final_N_score"]][idx])) {
            cat_data <- c(cat_data, "Nontransformed")
            score_data <- c(score_data, plot_data[["final_N_score"]][idx])
        #}
    }
    stat_data <- tibble(
        score = score_data,
        cat = factor(
            cat_data,
            ### This is needed to make sure the order is correct
            levels = c(levels(plot_data[[info_cat]]), "Nontransformed")
        )
    )

    comparisons <- lapply(
        levels(plot_data[[info_cat]]),
        function(x) c(x, "Nontransformed")
    )

    stat_test <- stat_data %>%
        t_test(score ~ cat, comparisons = comparisons) %>%
        adjust_pvalue(method = "fdr") %>%
        mutate(
            "{info_cat}" := group1,
            p.label = case_when(
                p.adj < 0.001 ~ '***',
                p.adj < 0.01  ~ '**',
                p.adj < 0.05  ~ '*',
                TRUE ~ as.character(signif(p.adj, 2))),
            y.position = signif(
                max(plot_data$final_T_score, na.rm = TRUE) + 0.2,
                2
            )
        ) 
        #%>% mutate(
        #    `.y.` = info_cat
        #)

    p <- plot_data %>%
        ggplot(
            aes(x = .data[[info_cat]], y = final_T_score)
        ) +
        # geom_boxplot(
        #    outlier.shape = NA
        # ) +
        # geom_jitter(width = 0.2) +
        geom_hline(
            yintercept = ctr_mean, linetype = 2, linewidth = 0.3, color = "blue"
        ) +
        geom_hline(
            yintercept = 0, linetype = 2, linewidth = 0.3
        ) +
        stat_summary(
            fun = mean,
            fun.min = mean,
            fun.max = mean,
            geom = "crossbar", width = 0.5, color = "gray"
        ) +
        geom_beeswarm(
            #method = 'pseudorandom',
            size = 1, cex = 3
        ) +
        stat_compare_means(
            comparisons = test_comparisons,
            size = 3,
            step.increase = 0.05,
            method = test_method,
            method.args = list(alternative = test_alt)
            #label.x.npc = "center"
        ) +
        ylab(y_label) +
        xlab(x_label)

    p <- p + stat_pvalue_manual(
        data = stat_test,
        label = "p.label",
        y.position = max(plot_data$final_T_score, na.rm = TRUE) + 0.5,
        xmin = "group1",
        xmax = NULL,
        size = 3,
        vjust = 1
    )

    if (!is.null(y_lim)) {
        p <- p + ylim(y_lim[1], y_lim[2])
    }

    if (is_label) {
        p <- p +
            geom_label_repel(
                aes(label = case_id),
                min.segment.length = unit(0, "lines")
                # force = 2,
                # nudge_y = 1.5
            )
    }

    p <- p +
        # facet_wrap(
        #    "mkr_name", ncol = 3
        #    #scales = "free_y"
        # ) +
        theme_classic() +
        theme(
            # axis.line=element_line(),
            strip.background = element_rect(colour = NA, fill = NA),
            panel.border = element_blank(),
            #panel.border = element_rect(fill = NA, color = "black"),
            text = element_text(color = "black"),
            axis.text = element_text(color = "black")
        )

    ### Save?
    if (is_save) {
        for (img_type in c(".pdf", ".png")) {
            if (is.null(img_name)) {
                img_name <- paste0(info_cat, "-", score_name_in)
            }

            ggsave(
                filename = here(
                    "figures", "markers", paste0(img_name, img_type)
                ),
                plot = p,
                width = 0.75 * nlevels(lcms_subset[[info_cat]]),
                height = 3.5, units = "in", dpi = 300
            )
        }

    } else{
        print(p)
    }
}

### "========================================================================="
### Plot MS overlap
### "========================================================================="
PlotMSOverlap <- function(peak_summary_mode)
{
    ### Load MSI data
    pos_msi_stats <- LoadMSIData(
        ion_type = "pos",
        peak_type = "ref_peaks",
        peak_summary_mode = peak_summary_mode,
        stat_type = "tissue",
        msi_conf
    ) %>% filter(
        !is.na(msi_mean)
    )

    neg_msi_stats <- LoadMSIData(
        ion_type = "neg",
        peak_type = "ref_peaks",
        peak_summary_mode = peak_summary_mode,
        stat_type = "tissue",
        msi_conf
    ) %>% filter(
        !is.na(msi_mean)
    )

    ### Normalize the data
    lpos_lcms_stats <- LoadLCMSHitData(lcms_mode = "LPOS", min_abd = 0) %>%
        filter(!is.na(lcms_abd), lcms_abd > 0)
    hpos_lcms_stats <- LoadLCMSHitData(lcms_mode = "HPOS", min_abd = 0) %>%
        filter(!is.na(lcms_abd), lcms_abd > 0)
    lneg_lcms_stats <- LoadLCMSHitData(lcms_mode = "LNEG", min_abd = 0) %>%
        filter(!is.na(lcms_abd), lcms_abd > 0)

    ### Get the mz lists
    lcms_lpos_mz <- unique(lpos_lcms_stats$peak_id)
    lcms_hpos_mz <- unique(hpos_lcms_stats$peak_id)
    lcms_lneg_mz <- unique(lneg_lcms_stats$peak_id)
    msi_pos_mz <- unique(pos_msi_stats$peak_id)
    msi_neg_mz <- unique(neg_msi_stats$peak_id)

    ### Construct a tible
    overlap_stat <- tibble(
        mass_num = c(
            sum(lcms_lpos_mz %in% msi_pos_mz),
            sum(!(lcms_lpos_mz %in% msi_pos_mz)),
            sum(lcms_lneg_mz %in% msi_neg_mz),
            sum(!(lcms_lneg_mz %in% msi_neg_mz)),
            sum(lcms_hpos_mz %in% msi_pos_mz),
            sum(!(lcms_hpos_mz %in% msi_pos_mz))
        ),
        lcms_type = factor(
            c(
                "Organic extract\n(pos. ESI)",
                "Organic extract\n(pos. ESI)",
                "Organic extract\n(neg. ESI)",
                "Organic extract\n(neg. ESI)",
                "Aqueous extract\n(pos. ESI)",
                "Aqueous extract\n(pos. ESI)"
            ),
            levels = c(
                "Organic extract\n(pos. ESI)",
                "Organic extract\n(neg. ESI)",
                "Aqueous extract\n(pos. ESI)"
            )),
        msi_found = factor(
            c("Yes", "No", "Yes", "No", "Yes", "No"),
            levels = c("No", "Yes")
        )
    )

    ### Calculate % found
    summary_stat <- overlap_stat %>%
        group_by(lcms_type) %>%
        summarize(
            msi_found,
            total = sum(mass_num),
            overlap_pct = mass_num / sum(mass_num) * 100,
            .groups = "drop"
        ) %>%
        filter(
            msi_found == "Yes"
        ) %>%
        mutate(
            overlap_pct = paste0(signif(overlap_pct, 3), "%")
        ) %>%
        dplyr::select(
            lcms_type,
            total,
            overlap_pct
        )

    print(overlap_stat)
    print(summary_stat)

    ### Plot the comparison
    p <- overlap_stat %>%
        ggplot(aes(x = lcms_type, y = mass_num, fill = msi_found)) +
        geom_bar(stat = "identity", color = "black", width = 0.6) +
        ylab("PM number") +
        xlab("LC-MS mode") +
        scale_fill_manual(
            name = "Also\ndetected\nusing\nDESI-MSI?",
            values = c("Yes" = "black", "No" = "#ffffff"),
            breaks = c("Yes", "No")
        ) +
        scale_y_continuous(
            breaks = c(0, 20, 40, 60, 80, 100),
            limits = c(0, 100),
            expand = c(0, 0)
        ) +
        #guides(fill = guide_legend(ncol = 2)) +
        theme_classic() +
        theme(
            legend.position = "right",
            legend.title.align = 0.5,
            legend.key.size = unit(0.75, "line"),
            axis.text = element_text(colour = "black"),
            axis.text.x = element_text(
                angle = 45, hjust = 1 #vjust = 0.5,
            )
        ) +
        geom_text(
            aes(
                x = lcms_type,
                y = total + 5,
                label = overlap_pct,
                fill = NULL
            ),
            data = summary_stat,
            size = 3
        )

    ### save the plots
    SavePlot(
        p,
        width = 3, height = 2,
        dir_name = "lcms", file_name = "lcms_overlap"
    )
}

### "========================================================================="
### Load MSI and LCMS tibblePlot MSI vs LCMS
### "========================================================================="
LoadMSILCMSTMeanTibble <- function(
    ion_type, peak_summary_mode, is_normalized, slide_ids
) {

    ### Load the clinical information
    clinical_info <- LoadClinicalInfo()

    valid_case_ids <- clinical_info %>%
        filter(slide_id %in% slide_ids) %>%
        dplyr::pull(case_id)

    ### Load MSI data
    msi_stats <- LoadMSIData(
        ion_type,
        peak_type = "ref_peaks",
        peak_summary_mode = peak_summary_mode,
        stat_type = "tissue", msi_conf
    ) %>% filter(
        case_id %in% valid_case_ids
    ) #%>% mutate(
    #    ### We retain zero
    #    msi_mean = replace_na(msi_mean, 0),
    #    msi_norm_mean = replace_na(msi_norm_mean, 0)
    #)

    ### Load the LCMS data
    lcms_stats <- if (ion_type == "pos") {
        bind_rows(
            LoadLCMSHitData(
                lcms_mode = "LPOS", min_abd = 0, is_ctrl_missing_ok = FALSE
            ),
            LoadLCMSHitData(
                lcms_mode = "HPOS", min_abd = 0, is_ctrl_missing_ok = FALSE
            )
        )

    } else {
        LoadLCMSHitData(
            lcms_mode = "LNEG", min_abd = 0, is_ctrl_missing_ok = FALSE
        )
    }

    ### Try to filter the LC/MS data
    # load(here("data", "lcms", "lcms_peaks.Rdata"))
    # new_mzs <- raw_peaks %>%
    #    filter(
    #        vip >= 0.995, ### Up to two decimal points
    #        iqr >= log2(2),
    #        med >= 10
    #    ) %>%
    #    dplyr::pull(mz)

    ### Some of the lcms may be replicated from LPOS and HPOS,
    ### thus, we need to take an average
    lcms_stats <- lcms_stats %>%
        group_by(
            case_id, section_id, peak_id
        ) %>%
        summarize(
            lcms_abd = mean(lcms_abd),
            .groups = "drop"
        )

    ### Get the validated peaks
    #validated_peaks <- lcms_peaks %>%
    #    filter(!is.na(putative_name)) %>%
    #    dplyr::pull(mz)

    ### Combine the lcms and msi stats
    all_stats <- left_join(
        msi_stats,
        lcms_stats,
        by = c("case_id", "section_id", "peak_id")
    ) %>%
    filter(
        ### Some of the mass does not have LC/MS
        !is.na(lcms_abd) #, lcms_abd > 0
    )

    ### Which columns to summarize?
    mean_msi_var <- ifelse(is_normalized, "msi_norm_mean", "msi_mean")

    ### Mean stats
    ### (all stats have a mix of tumor and normal samples)
    mean_stats <- all_stats %>%
        mutate(
            section_type = if_else(
                startsWith(section_id, "T"), "T", "N"
            ),
            .after = "section_id"
        ) %>%
        #filter(
        #    ### Don't consider normal
        #    !startsWith(section_id, "N")
        #) %>%
        #mutate(
        #    ### Don't consider metabolites that cannot be detected
        #    msi_norm_mean = replace(
        #        msi_norm_mean, msi_mean == 0, NA
        #    )
        #) %>%
        group_by(case_id, peak_id, section_type) %>%
        summarize(
            ave_msi_mean = mean(.data[[mean_msi_var]]),
            sd_msi_mean = sd(.data[[mean_msi_var]]),
            ave_lcms_mean = mean(lcms_abd),
            sd_lcms_mean = sd(lcms_abd),
            section_num = n(),
            .groups = "drop"
        ) %>%
        filter(
            ### There is no point comparing peaks that cannot be detected
            ### by MSI
            ave_msi_mean != 0,
            !is.na(ave_msi_mean),
            !is.na(ave_lcms_mean)
        )

    ### Add annotations
    mean_stats <- AddAnnotations(mean_stats)

    mean_stats

}

### "========================================================================="
### Add MSMS annotation to a data_tibble
### "========================================================================="
AddAnnotations <- function(data_tibble) {

    load(here("data", "lcms", "lcms_peaks.Rdata"))

    left_join(
        data_tibble,
        lcms_hit_peaks %>% dplyr::select(peak_id, MSMS_annotation, mz),
        by = "peak_id"
    ) %>%
    relocate(
        any_of(c("MSMS_annotation", "mz")),
        .after = "peak_id"
    )
}

### "========================================================================="
### Plot MSI vs LCMS
### * peak_id will check only for the specific mz
### "========================================================================="
PlotMSIvsLCMS <- function(
    ion_type, peak_summary_mode, is_normalized, slide_ids,
    x_lim = c(0, 30), peak_id = NULL, is_legend = TRUE
) {
    ### Get mean tibble
    mean_stats <- LoadMSILCMSTMeanTibble(
        ion_type, peak_summary_mode, is_normalized, slide_ids
    )

    ### filter MZ if needed
    if (!is.null(peak_id)) {
        mz_name <- as.character(peak_id)
        mean_stats <- mean_stats %>% filter(peak_id == !!peak_id)
    } else {
        mz_name <- "all"
    }

    ### Do correlation test
    cor_res <- cor.test(
        mean_stats$ave_msi_mean,
        mean_stats$ave_lcms_mean, method = "spearman"
    )

    ### Do we need a legend
    legend_position <- if (is_legend) {c(0.02, 1)} else {"none"}

    ### label
    if (ion_type == "pos") {
        title_txt <- "Positive PMs"
        y_lab <- "Mean normalized DESI-MSI\nabundance level (log2)"
        x_lab <- "Mean normalized LC-MS\nabundance level (log2)"

    } else {
        title_txt <- "Negative PMs"
        y_lab <- "Mean normalized DESI-MSI\nabundance level (log2)"
        x_lab <- "Mean normalized LC-MS\nabundance level (log2)"
    }

    p <- mean_stats %>%
        ggplot(
            aes(x = ave_lcms_mean, y = ave_msi_mean, col = case_id)
        ) +
        geom_point(size = 0.8) +
        scale_color_manual(
            name = "Patients",
            values = case_annot_colors
        ) +
        ylab(y_lab) +
        xlab(x_lab) +
        ggtitle(
            paste0(
                title_txt, "\n",
                "(Spearman's 𝜌 = ", signif(cor_res[["estimate"]], 3), ", ",
                "P = ", signif(cor_res[["p.value"]], 3), ")"
            )
        ) +
        coord_fixed() +
        xlim(x_lim) +
        #ylim(c(-8, 4)) +
        scale_y_continuous(limits=c(-8,4), breaks=c(-8, -6, -4, -2, 0, 2, 4)) +
        guides(col = guide_legend(ncol = 3)) +
        theme_classic(base_size = 11) +
        theme(
            plot.title = element_text(hjust = 0.5, size = 11),
            axis.text = element_text(colour="black"),
            legend.justification = c(0, 1),
            legend.position = legend_position,
            legend.title.align = 0.5
        )

    is_norm_str <- ifelse(is_normalized, "_norm", "")

    ### save the plots
    SavePlot(
        p, 
        width = 8, height = 2.5,
        dir_name = "lcms",
        file_name = paste0(
            "MSI_vs_LCMS_",
            ion_type, "_", mz_name, "_",
            peak_summary_mode, is_norm_str
        )
    )
}

### "========================================================================="
### Plot metabolites
### "========================================================================="
PlotMetabolites <- function(
    metabolites, ion_type,
    filename, padj_thres = 0.001
) {
    ### Generate the plot_data
    plot_data <- metabolites %>%
        filter(padj <= padj_thres) %>%
        mutate(
            Diff = tumor_mean - normal_mean
        ) %>%
        mutate_at(vars(ID), factor) %>%
        mutate(
            ID = fct_reorder(ID, Diff)
        ) %>%
        arrange(Diff)

    pos_plot <- plot_data %>%
        dplyr::select(
            ID, mz, normal_mean, tumor_mean, Diff
        ) %>%
        pivot_longer(
            c("normal_mean", "tumor_mean"),
            names_to = "Tissue",
            values_to = "Abundance"
        ) %>%
        ggplot(
            aes(x = ID, y = Abundance, fill = Tissue)
        ) +
        geom_bar(
            stat = "identity", position = "dodge", colour = "black"
        ) +
        scale_fill_discrete(labels = c("Normal", "Tumor")) +
        ylab("Abundance [log10]") +
        xlab("Metabolite ID") +
        theme_classic()

    ggsave(
        filename=here("figures", filename),
        plot = pos_plot,
        width = 6,
        height = 3,
        units = "in",
        dpi = 300
    )

    plot_data
}

### SgME Classifiers ========================================================
### "--------------------------------------------------------------------------"
### Common function to preprocess the training data
### * training_idx = subset to use for training
### * data_norm = full normalize data with the same rows as data_raw
### "--------------------------------------------------------------------------"
PreprocessData <- function(
    data_raw, grps, training_idx = NULL,
    is_balancing = FALSE, smote_nn = 3, smote_over_ratio = 0.5
) {
    ### If training_idx not provided, we will use all of data_raw
    if (is.null(training_idx)) {
        training_idx <- seq_len(nrow(data_raw))
    }

    ### Get the normalization coefficeints
    peak_means <- apply(
        data_raw[training_idx, ], 2, mean,
        na.rm = TRUE #, trim = 0.05 #mean_trim
    )

    peak_sds <- apply(
        data_raw[training_idx, ], 2, sd,
        na.rm = TRUE
    )

    ### Normalize the data
    data_norm <- scale(data_raw, center = peak_means, scale = peak_sds)

    ### Do we want to smoote imbalance data?
    if (is_balancing) {
        ### Balancing will only work for factor grps
        if (!is.factor(grps)) {
            stop("Balancing will only work for classifiers")
        }

        ### smote can only work with non-na data
        training_data_filled <- data_norm[training_idx, , drop = FALSE]
        for (i in 1:ncol(training_data_filled)) {
            training_data_filled[
                is.na(training_data_filled[, i]), i
            ] <- mean(training_data_filled[, i], na.rm = TRUE)
        }

        ### Resample the training data to make sure it is balanced
        sampled_data <- themis::smote(
            df = as_tibble(training_data_filled) %>%
                mutate(grps = grps[training_idx]),
            var = "grps",
            k = smote_nn,
            over_ratio = smote_over_ratio
        )

        training_data_balanced <- sampled_data %>%
            dplyr::select(where(is.numeric)) %>%
            as.matrix()

        training_grps_balanced <- sampled_data$grps

    } else {
        ### Just copy the data and grps
        training_data_balanced <- data_norm[training_idx, , drop = FALSE]
        training_grps_balanced <- grps[training_idx]
    }

    ### return
    list(
        data_norm = data_norm,
        training_data_balanced = training_data_balanced,
        training_grps_balanced = training_grps_balanced,
        peak_means = peak_means,
        peak_sds = peak_sds
    )
}

### "--------------------------------------------------------------------------"
### Function to perform and save PLS-DA
### * We need to use PLS-DA because they are missing data
### Note: for tunning just set msi_tissue_mat to zero
### "--------------------------------------------------------------------------"
RunSgMEClassifier <- function(
    msi_roi_mat, ### MSI for the ROIs
    msi_tissue_mat = NULL, ### If msi_tissue_mat is null, we will not apply
    selected_peak_ids = NULL, ### Only analyze these peaks
    pos_types, neg_types,
    model_type,
    model_opts, #orthoI, predI, permI = 20,
    is_post_processing = FALSE,
    median_filter_size = 1,
    is_impute = FALSE,
    impute_nn = 5,
    is_balancing = FALSE,
    smote_nn = 3, smote_over_ratio = 0.5,
    rnd_seed = 0,  ### Change this if stucked!
    file_name
) {
    ### "---------------------------------------------------------------------"
    ### Do we need to impute data?
    ### "---------------------------------------------------------------------"
    is_roi_na <- is.na(msi_roi_mat)
    roi_na_num <- sum(is_roi_na)
    roi_nona_num <- sum(!is_roi_na) 
    message(
        paste0(
            "msi_roi_mat has ", roi_na_num, "/", roi_na_num + roi_nona_num,
            " NA values (",
            round(roi_na_num / (roi_na_num + roi_nona_num) * 100, 2), "%)"
        )
    )

    if (!is.null(msi_tissue_mat)) {
        is_tissue_na <- is.na(msi_tissue_mat)
        tissue_na_num <- sum(is_tissue_na)
        tissue_nona_num <- sum(!is_tissue_na)
        message(
            paste0(
                "msi_tissue_mat has ", tissue_na_num, "/", tissue_na_num + tissue_nona_num,
                " NA values (",
                round(tissue_na_num / (tissue_na_num + tissue_nona_num) * 100, 2), "%)"
            )
        )
    }

    if (is_impute) {
        if (roi_na_num > 0) {
            message("Imputing msi_roi_mat ...")
            roi_impute_res <- impute::impute.knn(
                msi_roi_mat,
                k = impute_nn, rng.seed = rnd_seed
            )
            msi_roi_mat <- roi_impute_res$data
        }

        if (!is.null(msi_tissue_mat) && tissue_na_num > 0) {
            message("Imputing msi_tissue_mat ...")
            tissue_impute_res <- impute::impute.knn(
                msi_tissue_mat,
                k = impute_nn, rng.seed = rnd_seed
            )
            msi_tissue_mat <- tissue_impute_res$data
        }
    }

    ### "---------------------------------------------------------------------"
    ### Prepare roi_mat and grps
    ### "---------------------------------------------------------------------"
    ### init the groups (depends on roi_types)
    grps <- rep(NA, length(roi_types$slide_ids))

    for (pos_label in names(pos_types)) {
        for (pos_type in pos_types[[pos_label]]) {
            grps[roi_types[[pos_type]]] <- pos_label
        }
    }
    for (neg_label in names(neg_types)) {
        for (neg_type in neg_types[[neg_label]]) {
            grps[roi_types[[neg_type]]] <- neg_label
        }
    }

    ### The negative types should always be the first one, so that we will
    ### not have problem when interpolating the images. But caret assumes
    ### the positive level to be the first one, so we need to specific the
    ### positive class
    grps <- factor(grps, levels = c(names(neg_types), names(pos_types)))
    subset_idx <- !is.na(grps)

    grps <- grps[subset_idx]

    ### Get a local copy of the data matrix
    if (is.null(selected_peak_ids)) {
        selected_peak_ids <- rownames(msi_roi_mat)
    }
    data_raw <- t(msi_roi_mat[selected_peak_ids, subset_idx])

    ### "---------------------------------------------------------------------"
    ### If this is the final application, user CV to estimate its accuracy
    ### "---------------------------------------------------------------------"
    if (!is.null(msi_tissue_mat)) {
        PerformCV(
            data_raw = data_raw,
            grps = grps,
            model_type = model_type,
            model_opts = model_opts,
            file_name = "roi",
            path_score = file_name,
            cv_k = 10, cv_times = 10,
            rnd_seed = rnd_seed,
            is_balancing = is_balancing,
            smote_nn = smote_nn, smote_over_ratio = smote_over_ratio
        )
    }

    ### Determine the normalization coefficient
    ### Note: this must be done after subset selection, because different
    ### subsets may have different ranges of values
    pre_res <- PreprocessData(
        data_raw = data_raw,
        grps = grps,
        is_balancing = is_balancing,
        smote_nn = smote_nn,
        smote_over_ratio = smote_over_ratio
    )

    if (!is.null(msi_tissue_mat)) {
        ### Normalize the whole tissue data
        tissue_mat <- t(msi_tissue_mat[selected_peak_ids, ])
        tissue_mat <- scale(
            tissue_mat,
            center = pre_res[["peak_means"]],
            scale = pre_res[["peak_sds"]]
        )
    
        ### What is the data_type
        data_types <- factor(c(
            rep("balanced", nrow(pre_res[["training_data_balanced"]])),
            rep("training", nrow(pre_res[["data_norm"]])),
            rep("test", nrow(tissue_mat))
        ))

        ### Build the final model
        if (model_type == "plsda") {
            final_model <- opls(
                ### Note: If we don't include tissue_mat here, predict cannot
                ### handle NA values !!!
                rbind(
                    pre_res[["training_data_balanced"]],
                    pre_res[["data_norm"]],
                    tissue_mat
                ),
                c(
                    pre_res[["training_grps_balanced"]],
                    grps,
                    factor(
                        rep(levels(grps)[1], nrow(tissue_mat)),
                        levels = levels(grps)
                    )
                ),
                subset = 1:nrow(pre_res[["training_data_balanced"]]),
                crossvalI = 5, 
                permI = model_opts[["permI"]],
                predI = model_opts[["predI"]],
                orthoI = model_opts[["orthoI"]],
                scaleC = "none"
            )

        } else if (model_type == "svm") {
            final_model <- TrainClassifier(pre_res, model_type, model_opts)

        } else {
            stop("Unknown model type for final model buildin")
        }

    } else {
        ### "----------------------------------------------------------------"
        ### Model tuning only
        ### "----------------------------------------------------------------"
        ### Just generate empty outputs
        tissue_mat <- matrix(nrow = 0, ncol = ncol(data_raw))

        data_types <- factor(c(
            rep("balanced", nrow(pre_res[["training_data_balanced"]]))
        ))

        ### Build the model
        if (model_type == "plsda") {
            ### PLS-DA. This will plot the tuning plot, so that we can decide
            ### what predI is optimum. We don't really care about the final model
            final_model <- opls(
                pre_res[["training_data_balanced"]],
                pre_res[["training_grps_balanced"]],
                crossvalI = 5,
                permI = model_opts[["permI"]],
                predI = model_opts[["predI"]],
                orthoI = model_opts[["orthoI"]],
                scaleC = "none"
            )
        
        } else if (model_type == "svm") {
            ### SVM
            final_model <- TuneSVM(pre_res, model_type, model_opts, rnd_seed)
            ### Plot the output
            plot(final_model)
            print(final_model$bestTune)

        } else {
            stop("Unknown model type")
        }
    }

    ### Do we want to apply the model ?
    if (!is.null(msi_tissue_mat)) {
    #    ### Normalize the whole tissue data
    #    tissue_mat <- t(msi_tissue_mat[selected_peak_ids, ])
    #    tissue_mat <- scale(
    #        tissue_mat,
    #        center = pre_res[["peak_means"]],
    #        scale = pre_res[["peak_sds"]]
    #    )
    #
    #    ### Apply the model to predict the region types
    #    Note: if we call predict directly, it cannot handle missing values!!
        pred_labels <- predict(final_model, tissue_mat)
        names(pred_labels) <- rownames(tissue_mat)
                
        #pred_labels <- final_model@suppLs[["y"]][
        #    nrow(pre_res[["training_data_balanced"]]) + 1:nrow(tissue_mat)
        #]

        #pred_labels["A-03-01_2.T6A-0_2.xy_181_63"]

        ### Check to make sure pred_labels have not NA
        if (any(is.na(pred_labels))) {
            stop("pred_labels have NAs")
        }

    } else {
        pred_labels <- NULL
    }

    ### Process the pred_labels, which may be NULL
    ### This will generate the pred_img_list and update pred_labels if needed
    processed_res <- ProcessPredLabelsToImg(
        pred_labels, is_post_processing, median_filter_size 
    )

    ### Do some sanity check
    if (!is_post_processing) {
        if (any(pred_labels != processed_res[["pred_labels"]])) {
            stop("pred_labels is not the same as processed pred_labels")
        }
    }

    ### Return the results
    SgMEMap_res <- list(
        model = final_model,
        grps = grps,
        data_types = data_types,
        ### pred_labels will be with or without processing depending on the
        ### user option
        pred_labels = processed_res[["pred_labels"]],
        pred_img_list = processed_res[["pred_img_list"]],
        level_values = processed_res[["level_values"]]
    )

    ### Save the results
    fs::dir_create(path = here("results", model_type))

    saveRDS(
        SgMEMap_res,
        file = here("results", model_type, paste0(file_name, ".rds")),
        compress = "xz"
    )

    SgMEMap_res
}

### Load the trained model
LoadSgMEClassifier <- function(model_type, file_name) {
    readRDS(file = here("results", model_type, paste0(file_name, ".rds")))
}

### Get the distribution of all the groups from a traind model
GetGrpDist <- function(
    pred_labels,
    is_separate_tumor = TRUE,
    is_by_case = TRUE  ### Summarize by case
) {
    ### Get clinical info
    clinical_info <- LoadClinicalInfo()
    
    ### Get the prediction tibble
    tissue_pred_tibble <- GetPredTibble(pred_labels) %>%
        mutate(
            section_type = factor(
                if_else(
                    startsWith(section_id, "T"), "Tumor", "Normal"
                ),
                levels = c("Tumor", "Normal")
            ),
            slide_id = sub("^(.+)_[^_]$", "\\1", sample_id)
        )
    
    ### Do we want to separate by case?
    if (is_by_case) {
        ### Yes. Each case only one entry
        if (is_separate_tumor) {
            tissue_pred_tibble <- tissue_pred_tibble %>%
                group_by(slide_id, section_type)
        } else {
            tissue_pred_tibble <- tissue_pred_tibble %>%
                group_by(slide_id)
        }
    } else {
        ### No. Each section has one entry
        tissue_pred_tibble <- tissue_pred_tibble %>%
            mutate(
                section_id = sub("^([A-Za-z][0-9]*).*", "\\1", section_id)
            ) %>%
            group_by(slide_id, section_id)
    }

    tissue_pred_tibble <- tissue_pred_tibble %>%
        dplyr::count(pred_label, .drop = FALSE) %>%
        filter(!is.na(pred_label)) %>%
        mutate(
            percent = n / sum(n) * 100,
            tumor_percent = n / sum(n[pred_label != "Necrotic"]) * 100
        ) %>%
        ungroup() %>%
        left_join(
            clinical_info %>% dplyr::select(
                slide_id, edmonson_grade, metavir_fibrosis_score, tumor_stage_AJCC_V8,
                necrosis_score, steatosis_score, case_id
            ),
            by = "slide_id"
        ) %>%
        arrange(edmonson_grade, case_id) %>%
        mutate(
            tumor_stage_AJCC_V8 = case_when(
                tumor_stage_AJCC_V8 == "TNM Stage IIIA" ~ "TNM Stage III",
                tumor_stage_AJCC_V8 == "TNM Stage IIIB" ~ "TNM Stage III",
                TRUE ~ tumor_stage_AJCC_V8
            ),
            sample_id = paste0(case_id, "_", section_id)
        )
    
    tissue_pred_tibble
}

### "--------------------------------------------------------------------------"
### Common function to get region matrix
### "--------------------------------------------------------------------------"
GetRegMat <- function(
    spatial_regs, lcms_mat, feature_type, top_feature_num
) {
    ### Load the reference peaks
    load(here("data", "lcms", "lcms_peaks.Rdata"))

    ### Which types
    if (feature_type == "SgME") {
        ### Load the peak spatial distributions
        peak_spatial_dists <- readRDS(here("results", "SgME_markers.rds"))

        ### Load the peaks for the specific regions
        reg_peaks <- names(peak_spatial_dists)[
            peak_spatial_dists %in% spatial_regs
        ]

        ### Find the mz_ids
        if (is.null(top_feature_num)) {
            mz_ids <- lcms_hit_peaks %>%
                filter(peak_id %in% reg_peaks) %>%
                dplyr::pull(mz_id)
        
        } else {
            mz_ids <- lcms_hit_peaks %>%
                filter(peak_id %in% reg_peaks) %>%
                mutate(padj = p.adjust(TvsN_pval, method = "fdr")) %>%
                slice_min(padj, n = top_feature_num) %>%
                dplyr::pull(mz_id)
        }

        reg_mat <- lcms_mat[, mz_ids, drop = FALSE]

    } else if (feature_type == "HighAbd") {
        ### Find the mz_ids
        if (is.null(top_feature_num)) {
            mz_ids <- lcms_hit_peaks %>%
                dplyr::pull(mz_id)
        
        } else {
            mz_ids <- lcms_hit_peaks %>%
                mutate(padj = p.adjust(TvsN_pval, method = "fdr")) %>%
                slice_min(padj, n = top_feature_num) %>%
                dplyr::pull(mz_id)
        }

        ### Get the data
        reg_mat <- lcms_mat[, mz_ids, drop = FALSE]

    } else if (feature_type == "All") {
        ### Use all the features
        if (is.null(top_feature_num)) {
            reg_mat <- lcms_mat
        
        } else {
            mz_ids <- lcms_raw_peaks %>%
                mutate(padj = p.adjust(TvsN_pval)) %>%
                slice_min(padj, n = top_feature_num) %>%
                dplyr::pull(mz_id)
            
            reg_mat <- lcms_mat[, mz_ids, drop = FALSE]
        }

    } else {
        ### Unknown
        stop("Unknown feature type")
    }

    ### return
    reg_mat
}

### "--------------------------------------------------------------------------"
### Get subset index
### "--------------------------------------------------------------------------"
GetSubsetIdx <- function(
    lcms_grps, path_score, reg_mat, is_tumor_only
) {
    ### Sometimes, we will only analyze a subset of the samples
    sample_idx <- !is.na(lcms_grps[[path_score]])

    tissue_types <- sub("^.+_([^_]).*", "\\1", rownames(reg_mat))
    if (is_tumor_only) {
        sample_idx <- sample_idx & (tissue_types == "T")
    }

    sample_idx
}

### "--------------------------------------------------------------------------"
### Train the final Classifier
### "--------------------------------------------------------------------------"
TrainClassifier <- function (
    pre_res, model_type, model_opts, rnd_seed = -1
) {
    ### Set the seed
    if (rnd_seed > -1) {
        set.seed(rnd_seed)
    }

    if (model_type == "svm" && model_opts[["kernel"]] == "radial") {
        if (!("sigma" %in% names(model_opts))) {
            stop("sigma is not defined in model_opts")
        }
        if (!("C" %in% names(model_opts))) {
            stop("C is not defined in model_opts")
        }
        
        caret_method <- "svmRadial"
        tune_grid <- expand.grid(
            sigma = model_opts[["sigma"]],
            C = model_opts[["C"]]
        )

    } else if (model_type == "svm" && model_opts[["kernel"]] == "linear") {
        if (!("C" %in% names(model_opts))) {
            stop("C is not defined in model_opts")
        }

        caret_method <- "svmLinear"
        tune_grid <- expand.grid(
            C = model_opts[["C"]]
        )

    } else {
        stop("Unknown kernel type")
    }

    caret::train(
        x = pre_res[["training_data_balanced"]],
        y = pre_res[["training_grps_balanced"]],        
        method = caret_method,            
        trControl = trainControl(
            method = "none",
            classProbs = TRUE
        ),
        tuneGrid = tune_grid,
        preProcess = NULL,
        verbose = TRUE
    )

}

### "--------------------------------------------------------------------------"
### TuneSVM
### "--------------------------------------------------------------------------"
TuneSVM <- function (
    pre_res, model_type, model_opts, rnd_seed
) {
    ### Set the seed
    if (rnd_seed > -1) {
        set.seed(rnd_seed)
    }

    if (model_opts[["kernel"]] == "radial") {
        caret_method <- "svmRadial"
        tune_grid <- expand.grid(
            sigma = 10^(-4:2),
            C = 2^(-2:7)
        )

    } else if (model_opts[["kernel"]] == "linear") {
        caret_method <- "svmLinear"
        tune_grid <- expand.grid(
            C = 2^(-2:7)
        )

    } else {
        stop("Unknown kernel type")
    }

    caret::train(
        x = pre_res[["training_data_balanced"]],
        y = pre_res[["training_grps_balanced"]],        
        method = caret_method,
        trControl = trainControl(
            "repeatedcv", number = 10, repeats = 10, 
            classProbs = TRUE
        ),
        tuneGrid = tune_grid,
        preProcess = NULL,
        tuneLength = 10
    )
}

### "--------------------------------------------------------------------------"
### Perform CV
### * grp should be reversed ordered. meaning the negative class first, and
###   then positive class (eg. "Control", "Tumor")
### "--------------------------------------------------------------------------"
PerformCV <- function(
    data_raw, grps, model_type, model_opts, file_name, path_score,
    cv_k, cv_times,
    rnd_seed = 1234,
    is_balancing = FALSE, smote_nn = 3, smote_over_ratio = 0.5
) {
    ### Need to check and make sure data_raw and grps has NO missing values
    if (sum(is.na(data_raw)) > 0) {
        warning("PerformCV is handling data with missing values")
    }

    ### Make sure the model_type is valid
    if (!(model_type %in% c("plsda", "svm"))) {
        stop("Unknown model type for CV")
    }

    ### if grps have NA, we will need to remove the samples
    if (sum(is.na(grps)) > 0) {
        non_na_idx <- !is.na(grps)

        grps <- grps[non_na_idx]
        data_raw <- data_raw[non_na_idx, , drop = FALSE]
    }

    ### Set the seed
    set.seed(rnd_seed)
    
    ### Generate the CV fold
    cv_folds <- caret::createMultiFolds(grps, k = cv_k, times = cv_times)

    ### Init all the empty storages for performance measurements
    cv_train_class_perf <- cv_test_class_perf <- list()

    cv_test_class_pred_val <- cv_train_class_pred_val <-
        cv_test_class_true_val <- cv_train_class_true_val <- array()

    cv_train_perf <- cv_test_perf <- cv_pred <- matrix()

    cv_train_pred_val <- cv_train_true_val <-
        cv_test_pred_val <- cv_test_true_val <- matrix()

    ### Define the metric names
    classifier_metrics <- c("Sensitivity", "Specificity", "Balanced Accuracy")
    regression_metrics <- c("RMSE", "Rsquared", "MAE")

    ### Prepare the storages for multi-class classifier performance measurements
    ### Note: binary classifiers don't use this
    if (is.factor(grps) && nlevels(grps)>2) {
        for (idx in 1:nlevels(grps)) {
            cv_train_class_perf[[levels(grps)[idx]]] <- 
            cv_test_class_perf[[levels(grps)[idx]]] <-
                matrix(
                    NA, length(cv_folds), 3,
                    dimnames = list(names(cv_folds), classifier_metrics)
                )
        }

        ### Prepare the storage for multiple
        cv_test_class_pred_val <- cv_train_class_pred_val <- 
        cv_test_class_true_val <- cv_train_class_true_val <- array(
            NA, c(length(cv_folds), length(grps), nlevels(grps)),
            dimnames = list(
                names(cv_folds),
                rownames(data_raw),
                levels(grps)
            )
        )
    }

    ### Prepare the storages for binary classifers performance measurements
    if (is.factor(grps) && nlevels(grps) == 2) {
        cv_train_perf <- cv_test_perf <- matrix(
            NA, length(cv_folds), 3,
            dimnames = list(names(cv_folds), classifier_metrics)
        )
        cv_test_pred_val <- cv_train_pred_val <-
            cv_test_true_val <- cv_train_true_val <-
            matrix(NA, length(cv_folds), length(grps))
    }

    ### Prepare the storages for regression models performance measurements
    if (!is.factor(grps)) {
        cv_train_perf <- cv_test_perf <- matrix(
            NA, length(cv_folds), 3,
            dimnames = list(names(cv_folds), regression_metrics)
        )
        cv_test_pred_val <- cv_train_pred_val <-
            cv_test_true_val <- cv_train_true_val <-
            matrix(NA, length(cv_folds), length(grps))
    }

    for (fold_idx in seq_along(cv_folds)) {
    #mclapply(
    #    seq_along(cv_folds),
    #    function(fold_idx) {
        message(
            "Testing fold_idx = ", fold_idx,
            " (", fold_idx / length(cv_folds) * 100, "%)"
        )
        training_idx <- cv_folds[[fold_idx]]

        ### Preprocess the training data
        pre_res <- PreprocessData(
            data_raw, grps, training_idx,
            is_balancing = is_balancing,
            smote_nn = smote_nn,
            smote_over_ratio = smote_over_ratio
        )

        if (model_type == "plsda") {
            ### "-------------------------------------------------------------"
            ### Perform the opls model
            ### "-------------------------------------------------------------"
            ### In order to be able to retrieve the predicted values, we need
            ### to concadinate three datasets: 1) balanced, 2) train and
            ### 3) test
            final_model <- opls(
                rbind(
                    pre_res[["training_data_balanced"]],
                    pre_res[["data_norm"]][training_idx, ],
                    pre_res[["data_norm"]][-training_idx, ]
                ),
                c(
                    pre_res[["training_grps_balanced"]],
                    grps[training_idx],
                    grps[-training_idx]
                ),
                ### training_subset
                subset = 1:nrow(pre_res[["training_data_balanced"]]),
                ### This is used to estimate Q2Y. SInce we are not using it,
                ### we disable it by setting it to 1
                crossvalI = 1,
                predI = model_opts[["predI"]],
                ### THis is used to estimate the p-value of Q2Y. We also not
                ### using it
                permI = 0,
                scaleC = "none",
                fig.pdfC = "none",
                info.txtC = "none"
            )

        } else if (model_type == "svm") {
            ### "-------------------------------------------------------------"
            ### Perform SVM
            ### "-------------------------------------------------------------"
            if (is.factor(grps)) {
                ### We need to first tune the model
                ### Note: we will NOT set the random seed
                final_model <- TrainClassifier(pre_res, model_type, model_opts)

            } else {
                stop("SVM regression has not been implemented")
            }

        } else if (model_type == "glmnet") {
            ### "-------------------------------------------------------------"
            ### Use elastic net to do feature selection
            ### "-------------------------------------------------------------"
            if (is.factor(grps)) {
                final_model <- cv.glmnet(
                    pre_res[["training_data_balanced"]],
                    pre_res[["training_grps_balanced"]],
                    family = ifelse(nlevels(grps)>2, "multinomial", "binomial"),
                    type.multinomial = "grouped",
                    standardize = FALSE
                )

            } else {
                final_model <- cv.glmnet(
                    pre_res[["training_data_balanced"]],
                    pre_res[["training_grps_balanced"]],
                    standardize = FALSE
                )
            }

        } else {
            stop("Unknown model_type")

        }

        ### "-------------------------------------------------------------"
        ### Apply the model
        ### "-------------------------------------------------------------"
        if (model_type == "glmnet" && is.factor(grps)) {
            train_pred <- factor(
                as.character(
                    predict(
                        final_model,
                        pre_res[["data_norm"]][training_idx, ],
                        type = "class"
                    )
                ),
                levels = levels(grps)
            )

            test_pred <- factor(
                as.character(
                    predict(
                        final_model,
                        pre_res[["data_norm"]][-training_idx, ],
                        type = "class"
                    )
                ),
                levels = levels(grps)
            )

        } else {
            ### PLS-DA or SVM
            ### Note: ropls predict cannot handle NA values
            #balance_num <- nrow(pre_res[["training_data_balanced"]])
            #train_num <- length(training_idx)
            #test_num <- nrow(data_raw) - train_num

            #train_pred <- final_model@suppLs[["y"]][
            #    balance_num + 1:train_num
            #]

            #test_pred <- final_model@suppLs[["y"]][
            #    balance_num + train_num + 1:test_num
            #]
        
            train_pred <- predict(
                final_model,
                pre_res[["data_norm"]][training_idx, , drop = FALSE]
            )

            test_pred <- predict(
                final_model,
                pre_res[["data_norm"]][-training_idx, , drop = FALSE]
            )

            ### Sanity check
            if (any(is.na(train_pred)) || any(is.na(test_pred))) {
                stop("Either train_pred or test_pred have NA values")
            }
        }

        ### Get the labels
        train_label <- grps[training_idx]
        test_label <- grps[-training_idx]

        ### Measure the performance
        if (is.factor(grps)) {
            ### "--------------------------------------------------------------"
            ### Binary or multi-class classifiers
            ### "--------------------------------------------------------------"
            ### Compute the performance
            ### https://support.bioconductor.org/p/9143305/
            train_perf_cur <- caret::confusionMatrix(
                data = train_pred, reference = train_label,
                positive = if (nlevels(grps) == 2) {
                    levels(grps)[2]
                } else {
                    NULL
                }
            )
            test_perf_cur <- caret::confusionMatrix(
                data = test_pred, reference = test_label,
                positive = if (nlevels(grps) == 2) {
                    levels(grps)[2]
                } else {
                    NULL
                }
            )

            ### Get the predicted values for the classifier, so that we can
            ### also plot ROC and AUC
            if (model_type == "plsda") {
                ### The training predicted values are
                ### final_model@suppLs[["yPreMN"]], and the test predicted
                ### values are final_model@suppLs[["yTesMN"]]

                ### This is the predicted values for the balanced training data
                ### train_pred_val <- final_model@suppLs[["yPreMN"]]
                train_num <- length(train_label)
                test_num <- length(test_label)

                train_pred_val <- final_model@suppLs[["yTesMN"]][
                    1:train_num, , drop = FALSE
                ]
                rownames(train_pred_val) <- rownames(
                    pre_res[["data_norm"]][training_idx, ]
                )

                test_pred_val <- final_model@suppLs[["yTesMN"]][
                    (train_num + 1:test_num), , drop = FALSE
                ]

                rownames(test_pred_val) <- rownames(
                    pre_res[["data_norm"]][-training_idx, ]
                )

                if (nlevels(grps) > 2) {
                    ### There will be multiple column, each column for a class
                    colnames(train_pred_val) <- levels(grps)
                    colnames(test_pred_val) <- levels(grps)
                }

                #train_true_val <- final_model@suppLs[[".char2numF"]](
                #    matrix(grps[training_idx], ncol = 1)
                #)                
                #test_true_val <- final_model@suppLs[[".char2numF"]](
                #    matrix(grps[-training_idx], ncol = 1)
                #)

            } else if (model_type == "glmnet"){
                ### Assume is glmnet
                train_pred_val <- predict(
                    final_model,
                    pre_res[["data_norm"]][training_idx, ]
                )[, , 1]

                test_pred_val <- predict(
                    final_model,
                    pre_res[["data_norm"]][-training_idx, ]
                )[, , 1]

            } else if (model_type == "svm") {
                train_pred_val <- predict(
                    final_model,
                    pre_res[["data_norm"]][training_idx, , drop = FALSE],
                    type = "prob"
                )

                test_pred_val <- predict(
                    final_model,
                    pre_res[["data_norm"]][-training_idx, , drop = FALSE],
                    type = "prob"
                )

            } else {
                stop("Unknown model_type")
            }

            ### Save all the info
            if (nlevels(grps) > 2) {
                ###"-----------------------------------------------------------"
                ### Save the predicted values for ROC plotting
                ###"-----------------------------------------------------------"
                ### Prepare the train label
                train_true_val <- matrix(
                    0, length(grps[training_idx]), nlevels(grps),
                    dimnames = list(
                        rownames(pre_res[["data_norm"]][training_idx, ]),
                        levels(grps)
                    )
                )
                grp_numeric <- as.numeric(grps[training_idx])
                for (x in 1:length(grp_numeric)) {
                    train_true_val[x, grp_numeric[x]] <- 1
                }
                
                ### Prepare the test label
                test_true_val <- matrix(
                    0, length(grps[-training_idx]), nlevels(grps),
                    dimnames = list(
                        rownames(pre_res[["data_norm"]][-training_idx, ]),
                        levels(grps)
                    )
                )
                grp_numeric <- as.numeric(grps[-training_idx])
                for (x in 1:length(grp_numeric)) {
                    test_true_val[x, grp_numeric[x]] <- 1
                }

                for (class_name in levels(grps)) {
                    cv_train_class_true_val[
                        fold_idx, training_idx, class_name
                    ] <- train_true_val[, class_name]

                    cv_train_class_pred_val[
                        fold_idx, training_idx, class_name
                    ] <- train_pred_val[, class_name]

                    cv_test_class_true_val[
                        fold_idx, -training_idx, class_name
                    ] <- test_true_val[, class_name]

                    cv_test_class_pred_val[
                        fold_idx, -training_idx, class_name
                    ] <- test_pred_val[, class_name]
                }
            
                ### "----------------------------------------------------------"
                ### Save the classifier performance
                ### "----------------------------------------------------------"
                for (class_name in levels(grps)) {
                    cv_train_class_perf[[class_name]][fold_idx, ] <-
                        train_perf_cur$byClass[
                            paste0("Class: ", class_name), classifier_metrics
                        ]

                    cv_test_class_perf[[class_name]][fold_idx, ] <-
                        test_perf_cur$byClass[
                            paste0("Class: ", class_name), classifier_metrics
                        ]
                }

            } else {
                ### Binary classifiers
                ###"-----------------------------------------------------------"
                ### Save the predicted values for ROC plotting
                ###"-----------------------------------------------------------"
                ### Note: the original numeric values for grps is 1 and 2,
                ### thus we need to subtract one to get 0 and 1
                cv_train_true_val[fold_idx, training_idx] <-
                    as.numeric(grps[training_idx]) - 1

                cv_test_true_val[fold_idx, -training_idx] <-
                    as.numeric(grps[-training_idx]) - 1

                cv_train_pred_val[fold_idx, training_idx] <-
                    train_pred_val[, 1]
                
                cv_test_pred_val[fold_idx, -training_idx] <-
                    test_pred_val[, 1]

                ### "-----------------------------------------------------------"
                ### Save the classifier performance
                ### "-----------------------------------------------------------"
                cv_train_perf[fold_idx, ] <-
                    train_perf_cur$byClass[classifier_metrics]
                
                cv_test_perf[fold_idx, ] <-
                    test_perf_cur$byClass[classifier_metrics]
            }

        } else {
            ### "--------------------------------------------------------------"
            ### Regression models
            ### "--------------------------------------------------------------"
            ### Regression predicted values for scatter plot
            cv_train_pred_val[fold_idx, training_idx] <- train_pred
            cv_test_pred_val[fold_idx, -training_idx] <- test_pred

            cv_train_true_val[fold_idx, training_idx] <- train_label
            cv_test_true_val[fold_idx, -training_idx] <- test_label

            ### Regression performance
            cv_test_perf[fold_idx, ] <- caret::postResample(
                pred = test_pred, obs = test_label
            )

            cv_train_perf[fold_idx, ] <- caret::postResample(
                pred = train_pred, obs = train_label
            )
        }
    #)
    }

    ### Print the average performance
    if (is.factor(grps) && nlevels(grps) > 2) {
        ### Loop through all classes
        all_perf <- NULL
        for (class_name in levels(grps)) {
            all_perf <- rbind(
                all_perf,
                colMeans(cv_test_class_perf[[class_name]], na.rm = TRUE)
            )
        }
        rownames(all_perf) <- levels(grps)
        print(all_perf)

    } else {
        print(colMeans(cv_test_perf, na.rm = TRUE))
    }

    ### Get the tissue types
    tissue_types <- sub("^.+_([^_]).*", "\\1", rownames(data_raw))

    ### Save the results
    fs::dir_create(here("results", "final_models"))
    save(
        grps, tissue_types,
        
        cv_test_perf, cv_train_perf,        
        cv_train_true_val, cv_test_true_val,
        cv_train_pred_val, cv_test_pred_val,

        cv_test_class_perf, cv_train_class_perf,        
        cv_test_class_pred_val, cv_train_class_pred_val,
        cv_test_class_true_val, cv_train_class_true_val,
        
        file = here(
            "results", "final_models",
            paste0(file_name, "_", path_score, ".Rdata")
        ),
        compress = "xz"
    )
}

### "--------------------------------------------------------------------------"
### Build and test the final opls models
### "--------------------------------------------------------------------------"
TestFinalModel <- function(
    model_type, 
    file_name,
    feature_type, # "All", "HighAbd", "SgME", "SgME_Top"
    top_feature_num,
    spatial_regs, path_score, is_tumor_only,
    lcms_mat, lcms_grps,
    predI_opt, cv_k = 10, cv_times = 10,
    is_balancing = FALSE, smote_nn = 3, smote_over_ratio = 0.5
) {
    ### Get the region matrix and sample subset idx
    reg_mat <- GetRegMat(spatial_regs, lcms_mat, feature_type, top_feature_num)
    sample_idx <- GetSubsetIdx(lcms_grps, path_score, reg_mat, is_tumor_only)

    ### We need to take an average for each case!!!
    ### Continue here
    data_mean_tibble <- as_tibble(
            reg_mat[sample_idx, ],
            rownames = "case_id"
        ) %>%
        mutate(
            case_id = sub("^(.+)_.+$", "\\1", case_id),
            grps = lcms_grps[[path_score]][sample_idx]
        ) %>%
        group_by(case_id, grps) %>%
        summarize_all(mean, na.rm = TRUE) %>%
        ungroup()

    grps <- data_mean_tibble %>% dplyr::pull(grps, name = case_id)

    data_raw <- data_mean_tibble %>%
        mutate(case_id = NULL, grps = NULL) %>%
        #column_to_rownames("case_id") %>%
        as.matrix()

    ### Do the tunning, it will plot the data
    if (model_type == "pls") {
        tune_model <- opls(
            data_raw,
            grps,
            crossvalI = 5,
            predI = predI_opt
        )
    }

    ### Show the grp distrubtion
    print(table(grps))

    ### Perform CV
    PerformCV(
        data_raw, grps, model_type, file_name, path_score,
        predI_opt, cv_k, cv_times,
        is_balancing = is_balancing,
        smote_nn = smote_nn, smote_over_ratio = smote_over_ratio
    )
}

### "--------------------------------------------------------------------------"
### Load the CV performance
### "--------------------------------------------------------------------------"
LoadCVPerf <- function(
    data_type, model_name
) {
    ### Load the results
    load(
        here(
            "results", "final_models",
            paste0(data_type, "_", model_name, ".Rdata")
        )
    )

    ### If this PLS-DA?
    if (!is.factor(grps)) {
        stop("Model is not a classifier")
    }

    ### Generate the list of ROC
    perf_list <- list()

    if (nlevels(grps) > 2) {
        ### Multi-class classifier
        for (class_name in levels(grps)) {
            perf_list[["training"]] <-
                rbind(
                    perf_list[["training"]],
                    colMeans(cv_train_class_perf[[class_name]])
                )

            perf_list[["test"]] <-
                rbind(
                    perf_list[["test"]],
                    colMeans(cv_test_class_perf[[class_name]])
                )
        }

        ### Add row names
        rownames(perf_list[["training"]]) <-
            rownames(perf_list[["test"]]) <- levels(grps)
    
    } else {
        perf_list[["training"]] <- matrix(
            colMeans(cv_train_perf),
            nrow = 1,
            dimnames = list(levels(grps)[2], colnames(cv_train_perf))
        )

        perf_list[["test"]] <- matrix(
            colMeans(cv_test_perf),
            nrow = 1,
            dimnames = list(levels(grps)[2], colnames(cv_test_perf))
        )
    }

    perf_list

}

### "--------------------------------------------------------------------------"
### Load ROC results for a final model
### "--------------------------------------------------------------------------"
LoadROC <- function(
    data_type, model_name
) {
    ### Load the results
    load(
        here(
            "results", "final_models",
            paste0(data_type, "_", model_name, ".Rdata")
        )
    )

    ### If this PLS-DA?
    if (!is.factor(grps)) {
        stop("Model is not a classifier")
    }

    ### Generate the list of ROC
    roc_list <- list()

    if (nlevels(grps) > 2) {
        ### Multi-class classifier
        for (class_name in levels(grps)) {
            roc_list[["training"]][[class_name]] <- roc(
                as.numeric(cv_train_class_true_val[, , class_name]),
                as.numeric(cv_train_class_pred_val[, , class_name]),
                percent = TRUE
            )

            roc_list[["test"]][[class_name]] <- roc(
                as.numeric(cv_test_class_true_val[, , class_name]),
                as.numeric(cv_test_class_pred_val[, , class_name]),
                percent = TRUE
            )
        }

    } else {
        ### Binary classifier
        class_name <- levels(grps)[2]

        roc_list[["training"]][[class_name]] <- roc(
            as.numeric(cv_train_true_val),
            as.numeric(cv_train_pred_val),
            percent = TRUE
        )

        roc_list[["test"]][[class_name]] <- roc(
            as.numeric(cv_test_true_val),
            as.numeric(cv_test_pred_val),
            percent = TRUE
        )
    }

    roc_list
}

### "------------------------------------------------------------------------"
### Function to get the significantly correlated genes for a set of metabolites
### "------------------------------------------------------------------------"
GetCEGs <- function(
    reg_name, top_lcms_list, RNA_mat,
    corr_padj_thres = 0.1,
    min_metabolite_coverages
) {
    ### Get the metabolites
    reg_mat <- top_lcms_list[[reg_name]]

    ### Form a pairwise correlation matrix
    corr_val_mat <- matrix(
        NA,
        nrow = ncol(reg_mat),
        ncol = ncol(RNA_mat),
        dimnames = list(colnames(reg_mat), colnames(RNA_mat))
    )

    corr_pval_mat <- matrix(
        NA,
        nrow = ncol(reg_mat),
        ncol = ncol(RNA_mat),
        dimnames = list(colnames(reg_mat), colnames(RNA_mat))
    )

    for (lcms_idx in 1:ncol(reg_mat)) {
        ## message(lcms_idx)
        reg_data <- reg_mat[, lcms_idx]

        for (RNA_idx in 1:ncol(RNA_mat)) {
            cor_res <- cor.test(
                reg_data,
                RNA_mat[, RNA_idx],
                method = "spearman",
                #method = "pearson",
                alternative = "greater"
                #alternative = "two.sided"
            )
            ## Disable it if we don't need it
            corr_pval_mat[lcms_idx, RNA_idx] <- cor_res$p.value
            corr_val_mat[lcms_idx, RNA_idx] <- cor_res$estimate[[1]]
        }
    }

    corr_padj_mat <- matrix(
        p.adjust(corr_pval_mat, method = "BH"),
        nrow = ncol(reg_mat),
        ncol = ncol(RNA_mat),
        dimnames = list(colnames(reg_mat), colnames(RNA_mat))
    )

    corr_sig_mat <- corr_val_mat
    corr_sig_mat[corr_padj_mat > corr_padj_thres] <- NA

    sig_nums <- sort(colSums(!is.na(corr_sig_mat)))
    good_idx <- sig_nums > nrow(corr_sig_mat) * min_metabolite_coverages[reg_name]
    sig_ensembls <- names(sig_nums)[good_idx]
    sig_ensembls <- unique(sig_ensembls)
    message(
        reg_name, ", total metabolites = ", nrow(corr_sig_mat), 
        ", passed RNAs = ", sum(good_idx), " (", min_metabolite_coverages[reg_name], ")"
    )

    ### Convert to ids
    sig_entrezs <- unique(
        all_gene_info %>%
            filter(ENSEMBL %in% sig_ensembls) %>%
            dplyr::pull(ENTREZID)
    )

    ### Generate the export table
    export_data <- tibble(
            ENSEMBL = sig_ensembls
        ) %>% left_join(
            all_gene_info
        ) %>% group_by(
            ENSEMBL
        ) %>% summarize(
            ENTREZID = paste(ENTREZID, collapse = ", "),
            SYMBOL = paste(SYMBOL, collapse = ", ")
        ) %>% mutate(
            corr_median = apply(corr_val_mat[, sig_ensembls], 2, median),
            corr_mean = apply(corr_val_mat[, sig_ensembls], 2, mean),
            signif_num = colSums(corr_padj_mat[, sig_ensembls] <= corr_padj_thres)
        ) %>% arrange(-corr_median)

    list(
        CEG_entrezs = sig_entrezs,
        CEG_ensembls = sig_ensembls,
        CEG_corr = corr_val_mat[, sig_ensembls],
        CEG_padj = corr_padj_mat[, sig_ensembls],
        export_data = export_data
    )
}

### "------------------------------------------------------------------------"
### Function to get the significantly correlated genes for a set of metabolites
### for a SgME region
### "------------------------------------------------------------------------"
GetSgMECEGs <- function(
    spatial_regs, lcms_mat, RNA_mat,
    corr_padj_thres,
    min_metabolite_coverages
) {
    ### Get the lcms data for a specific metabolic domain
    reg_mat <- GetRegMat(
        spatial_regs = spatial_regs,
        lcms_mat = lcms_mat,
        feature_type = "SgME",
        top_feature_num = NULL
    )

    ### Get the CEGs
    GetCEGs(
        reg_mat, RNA_mat,
        corr_padj_thres,
        min_metabolite_coverages
    )
}


### "========================================================================="
### Process the pred labels and convert it into an image
### "========================================================================="
ProcessPredLabelsToImg <- function(
    pred_labels,
    is_post_processing = FALSE,
    median_filter_size = 1
) {
    ### Init an empty container
    pred_img_list <- list()
    level_values <- NULL
    processed_pred_labels <- NULL

    ### Is pred_labels NULL? Sometimes, we just want to train the model, and
    ### this, pred_labels will be NULL
    if (!is.null(pred_labels)) {
        message("Generating SgME maps ...")

        ### Convert to tibble
        pred_tibble <- GetPredTibble(pred_labels)
        pred_levels <- levels(pred_tibble$pred_label)
        sample_ids <- unique(pred_tibble$sample_id)

        ### Perform filtering
        processed_pred_labels <- factor(levels = pred_levels)

        for (sample_id in sample_ids) {
            ### Prepare the full pred img, so that we can do a median filter
            sample_pred <- ProcessPredTibbleToImg(
                pred_tibble,
                sample_id = sample_id,
                median_filter_size = median_filter_size,
                is_post_processing = is_post_processing
            )

            ### Copy the image
            pred_img_list[[sample_id]] <- sample_pred$pred_img
            level_values <- sample_pred$level_values

            ### factors cannot be concantenate by using c()
            processed_pred_labels <- unlist(
                list(processed_pred_labels, sample_pred$pred_labels)
            )
        }

        ## Re-arrange processed_pred_labels
        processed_pred_labels <- processed_pred_labels[names(pred_labels)]
    }

    ### Return
    list(
        ### Is no processing, processed_pred_labels should be exactly the
        ### same as pred_labels
        pred_labels = processed_pred_labels,
        pred_img_list = pred_img_list,
        level_values = level_values
    )
}

### "========================================================================="
### Convert a PLS-DA prediction tibble to image
### Note: the function will also do median filtering
### "========================================================================="
ProcessPredTibbleToImg <- function(
    pred_tibble, sample_id,
    median_filter_size = 1,
    is_post_processing = TRUE
) {
    ### Filter the pred_tibble
    pred_tibble <- pred_tibble %>% filter(sample_id == !!sample_id)

    ### Make sure the tibble only has one sample/slide/image
    if (n_distinct(pred_tibble$sample_id) > 1) {
        stop("pred_tibble has more than one sample")
    }

    ### Get the prediction levels of the original image
    pred_levels <- levels(pred_tibble$pred_label)
    max_val <- nlevels(pred_tibble$pred_label)

    ### Convert to an image matrix
    pred_img <- ConvertPredTibbleToMat(
        pred_tibble,
        sample_id = sample_id,
        pred_label_name = "pred_label"
    )

    ### Do we want to do post processing?
    if (is_post_processing) {
        ### EBImage must be between 0 and 1, thus we need to scale the image
        ### 0 = background
        ### 1, 2, .. = the grp levels
        pred_img <- pred_img / max_val

        ### Perform the filter
        ### Note: the returned values have rounding errors, so we need to do
        ### the rounding
        pred_img_smooth <- round(
            medianFilter(pred_img, median_filter_size) * max_val
        )

        ### Opening only work on binary images!!! but pred_img may be multiple
        ### categories
        #pred_img_smooth <- round(
        #    opening(pred_img_smooth, makeBrush(3, shape = 'disc')) * max_val
        #)

        ### At the end, we convert the results back to pred_label format
        ### Note: processed_pred_tibble may be completely empty!!
        ### when there is no pixel in the image
        processed_pred_tibble <- pred_img_smooth %>%
            as_tibble(rownames = "y") %>%
            pivot_longer(
                !y,
                names_to = "x",
                values_to = "processed_pred_label"
            ) %>%
            ### We don't need the background
            filter(processed_pred_label > 0)

        ### Convert pred_label back to factor
        processed_pred_tibble$processed_pred_label <- factor(
            pred_levels[processed_pred_tibble$processed_pred_label],
            levels = pred_levels
        )

        ### Add the processed back to pred_tibble
        pred_tibble <- left_join(
            pred_tibble,
            processed_pred_tibble,
            by = c("x", "y")
        )

        ### A small number of pixels will be NA because they are closer to the
        ### boundaries. So, we will revert them back to the original prediction
        pix_num <- nrow(pred_tibble)
        na_num <- sum(is.na(pred_tibble$processed_pred_label))

        if (na_num / pix_num > 0.05) {
            message("Total number of pixels = ", pix_num)
            message("Total number of NAs    = ", na_num)
            stop("More than 5% pixels are NAs")
        }

        if (na_num > 0) {
            na_idx <- which(is.na(pred_tibble$processed_pred_label))
            pred_tibble$processed_pred_label[na_idx] <-
                pred_tibble$pred_label[na_idx]
        }

        ### Finally, we just repliacte the pred_label
        #pred_tibble <- pred_tibble %>%
        #    mutate(pred_label = processed_pred_label)

        ### We need to regenerate the image matrix, because pred_label
        ### may have changed
        pred_img <- ConvertPredTibbleToMat(
            pred_tibble,
            sample_id = sample_id,
            pred_label_name = "processed_pred_label"
        )

    } else {
        pred_tibble <- pred_tibble %>%
            mutate(processed_pred_label = pred_label)
    }

    ### Get back the processed pred labels
    ### Convert the id to A-11-01_2.T5-0_3.xy_351_101
    pred_labels <- pred_tibble %>%
        unite(pos_id, c("x", "y"), sep = "_") %>%
        mutate(pos_id = paste0("xy_", pos_id)) %>%
        unite(id, c("sample_id", "section_id", "pos_id"), sep = ".") %>%
        dplyr::pull(processed_pred_label, name = id)

    ### Note: full_pred_img_filtered's values may be different from
    ### full_pred_img due to rounding error in EDImage, so we need to 
    ### rebuild the pred_label key tables
    level_values <- 0:max_val
    names(level_values) <- c("Background", pred_levels)

    ### return the data
    list(
        pred_img = pred_img,
        level_values = level_values,
        pred_labels = pred_labels
    )
}

### "========================================================================="
### Convert pred tibble to image matrix
### "========================================================================="
ConvertPredTibbleToMat <- function(
    pred_tibble, sample_id, pred_label_name
) {
    ### Get the dimension of the image
    img_dim <- GetMSIDim(sample_id, "pos")

    ### Get the levels
    pred_levels <- levels(pred_tibble[[pred_label_name]])
    max_val <- nlevels(pred_tibble[[pred_label_name]])

    ### Convert to a matrix
    pred_img <- pred_tibble %>%
        dplyr::select(all_of(c("x", "y", pred_label_name))) %>%
        #rename(pred_label = all_of(pred_label_name)) %>%
        mutate(
            x = as.integer(x),
            y = as.integer(y),
            "{pred_label_name}" := as.numeric(.data[[pred_label_name]])
        ) %>%
        complete(x = 1:img_dim[1], y = 1:img_dim[2]) %>%
        #complete(x = min(x):max(x), y = min(y):max(y)) %>%
        pivot_wider(
            names_from = "x",
            values_from = !!pred_label_name,
            values_fill = NA
        ) %>%
        arrange(y) %>%
        column_to_rownames("y") %>%
        as.matrix()

    ### Rearrange the column
    pred_img <- pred_img[
        , as.character(sort(as.integer(colnames(pred_img))))
    ]

    ## Leave the NA alone
    pred_img[is.na(pred_img)] <- 0

    ### return
    pred_img
}
