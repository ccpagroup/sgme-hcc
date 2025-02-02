### "========================================================================="
### Plot Fig 1
###
### Copyright (c) 2021-2025. Bioinformatics Institute, A*STAR, Singapore
### License: MIT License. See LICENSE for the full license.
### Authors: Lit-Hsin Loo (loolh@bii.a-star.edu.sg)
### "========================================================================="
library("here")
source(here("lib", "MSI_lib.R"))
source(here("lib", "Study_lib.R"))

### Get clinical info
clinical_info <- LoadClinicalInfo()

### Load the peaks
load(here("data", "lcms", "lcms_peaks.Rdata"))

### "=========================================================================="
### Fig. 1c The distribution of the LC/MS samples
### "=========================================================================="
### Load the clinical information
sample_list <- read_excel(
        here("data", "lcms", "Supp_Tables_S1 to S4_v3.xlsx"),
        sheet = "Table_S1b Sample_Usages"
    ) %>%
    dplyr::select("Case ID", "Section ID") %>%
    dplyr::rename(case_id = `Case ID`, section_id = `Section ID`) %>%
    slice_head(n=117)

sample_dist <- left_join(
        sample_list,
        clinical_info %>%
            dplyr::select(
                case_id, 
                histological_diagnosis, tumor_stage_AJCC_V8,
                edmonson_grade
            ),
        by = c("case_id")
    ) %>%
    mutate(
        sample_type = factor(
            case_when(
                startsWith(section_id, "T") & histological_diagnosis ==
                    "HCC" ~ "HCC tumor",
                startsWith(section_id, "T") & histological_diagnosis ==
                    "CCA" ~ "iCCA tumor",
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
                edmonson_grade == 1 ~ "G1",
                edmonson_grade == 2 ~ "G2",
                edmonson_grade == 3 ~ "G3",
                edmonson_grade == 4 ~ "G4",
                TRUE ~ "Unknown"
            ),
            levels = c("G1", "G2", "G3", "G4")
        )
    ) %>%
    mutate(section_id = NULL) %>%
    group_by(case_id, sample_type) %>%
    mutate(sample_num = n()) %>%
    distinct() %>%
    arrange(
        #histological_diagnosis, tumor_stage
        tumor_stage, case_id        
    )

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
            roi_annot_colors[["G3"]], "#4daf4a", roi_annot_colors[["Normal"]]
        )
    ) +
    ylab("Sector number") +
    xlab("Case ID") +
    theme_classic() +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 8)) +
    theme(
        # plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
        legend.key.size = unit(0.4, 'cm'),
        legend.title = element_text(size = 8, face = "bold"),
        legend.text = element_text(size = 8),
        legend.position = "top",
        #legend.position = c(0.8, 0.8),
        #legend.direction="horizontal",
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
    )

### save the plots
fs::dir_create(here("figures", "lcms"))
for (img_type in c(".png", ".pdf")) {
    ggsave(
        filename = here(
            "figures", "lcms",
            paste0("lcms_sample_dist", img_type)
        ),
        plot = sample_dist_plot,
        width = 4.1,
        height = 2.1,
        units = "in",
        dpi = 300
    )
}

### "=========================================================================="
### Fig. S2 Can we distinguish between normal, HCC, and CCA?
### * Tumors are different from Normal
### * HCC is different from CCA. Thus CCA can be removed.
### "=========================================================================="
load(
    file = here("results", "plsda", "lcms_tumor.Rdata")
)

SavePLSDAPlot(
    width = 5, height = 5.5,
    final_model = list(
        model = lcms_tumor,
        data_types = rep("balanced", length(lcms_grps$tumor_grps))
    ),
    grp_ids = lcms_grps$tumor_grps,
    plot_components = c(1, 2),
    grp_labels = c(
        "Normal (21)" = "N",
        "HCC Tumor (80)" = "H",
        "CCA Tumor (6)" = "C"
    ),
    grp_colors = c(
        "Normal (21)" = roi_annot_colors[["Normal"]],
        "HCC Tumor (80)" = roi_annot_colors[["G3"]],
        "CCA Tumor (6)" = roi_annot_colors[["Necrotic"]]
    ),
    title_txt = paste0(
        "Metabolomics (",
        nrow(ropls::getLoadingMN(lcms_tumor)),
        " MFs)"
    ),
    filename = "LCMS_tumor.pdf"
)

### "=========================================================================="
### Fig. 1d the lcms hcc stage metabolomics landscape
### "=========================================================================="
load(
    file = here("results", "plsda", "lcms_hcc_stage.Rdata")
)

SavePLSDAPlot(
    width = 5, height = 5.5,
    final_model = list(
        model = lcms_hcc_stage,
        data_types = rep("balanced", length(lcms_grps$hcc_stage_grps))
    ),
    grp_ids = factor(
        lcms_grps$hcc_stage_grps,
        levels = c(
            "HCC Stage IB",
            "HCC Stage II",
            "HCC Stage III",
            "Normal"
        )
    ),
    plot_components = c(1, 2),
    grp_labels = c(
        "Normal" = "0",
        "HCC Stage IB" = "1",
        "HCC Stage II" = "2",
        "HCC Stage III" = "3"
    ),
    grp_colors = c(
        "Normal" = roi_annot_colors[["Normal"]],
        "HCC Stage IB" = roi_annot_colors[["Fibrotic"]],
        "HCC Stage II" = roi_annot_colors[["G2"]],
        "HCC Stage III" = roi_annot_colors[["G3"]]
    ),
    title_txt = paste0(
        "Metabolomics (",
        nrow(ropls::getLoadingMN(lcms_hcc_stage)),
        " MFs)"
    ),
    filename = "LCMS_hcc_stage.pdf"
)

table(lcms_grps$hcc_stage_grps)

### "=========================================================================="
### Fig. 1f The RNA transcriptomics landscape
### "=========================================================================="
load(
    file = here("results", "plsda", "RNA_hcc_stage.Rdata")
)

SavePLSDAPlot(
    width = 5, height = 5.5,
    final_model = list(
        model = RNA_hcc_stage,
        data_types = rep("balanced", length(RNA_grps$hcc_stage_grps))
    ),
    grp_ids = factor(
        RNA_grps$hcc_stage_grps,
        levels = c(
            "HCC Stage IB",
            "HCC Stage II",
            "HCC Stage III",
            "Normal"
        )
    ),
    plot_components = c(1, 2),
    grp_labels = c(
        "Normal" = "0",
        "HCC Stage IB" = "1",
        "HCC Stage II" = "2",
        "HCC Stage III" = "3"
    ),
    grp_colors = c(
        "Normal" = roi_annot_colors[["Normal"]],
        "HCC Stage IB" = roi_annot_colors[["Fibrotic"]],
        "HCC Stage II" = roi_annot_colors[["G2"]],
        "HCC Stage III" = roi_annot_colors[["G3"]]
    ),
    title_txt = paste0(
        "Transcriptomics (",
        nrow(ropls::getLoadingMN(RNA_hcc_stage)),
        " RNAs)"
    ),
    "RNA_hcc_stage.pdf"
)

table(RNA_grps$hcc_stage_grps)

