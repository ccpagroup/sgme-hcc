### "========================================================================="
### Find altered MFs and DEGs
###
### Copyright (c) 2021-2025. Bioinformatics Institute, A*STAR, Singapore
### License: MIT License. See LICENSE for the full license.
### Authors: Lit-Hsin Loo (loolh@bii.a-star.edu.sg)
### "========================================================================="
library("here")
source(here("lib", "MSI_lib.R"))
source(here("lib", "Study_lib.R"))
library(ComplexHeatmap)
library(DESeq2)

### Reload the lcms and RNA data
lcms_stats <- readRDS(here("data", "lcms", "lcms_stats.rds"))
RNA_stats <- readRDS(here("data", "RNA", "RNA_stats.rds"))

### Plot the heatmap
PlotStage <- function(cur_stats, log2FC_thres, padj_thres, omics_name, col_title) {

    plot_mat <- cur_stats %>%
        filter(
            S1_vs_N_padj <= padj_thres |
            S2_vs_N_padj <= padj_thres |
            S3_vs_N_padj <= padj_thres,
            abs(S1_log2FC) > log2FC_thres |
            abs(S2_log2FC) > log2FC_thres |
            abs(S3_log2FC) > log2FC_thres 
        ) %>%
        mutate(
            S1N = case_when(
                (S1_log2FC > log2FC_thres) & (S1_vs_N_padj <= padj_thres) ~ "+",
                (S1_log2FC < -log2FC_thres) & (S1_vs_N_padj <= padj_thres) ~ "-",
                TRUE ~ "0"
            ),
            S2N = case_when(
                (S2_log2FC > log2FC_thres) & (S2_vs_N_padj <= padj_thres) ~ "+",
                (S2_log2FC < -log2FC_thres) & (S2_vs_N_padj <= padj_thres) ~ "-",
                TRUE ~ "0"
            ),
            S3N = case_when(
                (S3_log2FC > log2FC_thres) & (S3_vs_N_padj <= padj_thres) ~ "+",
                (S3_log2FC < -log2FC_thres) & (S3_vs_N_padj <= padj_thres) ~ "-",
                TRUE ~ "0"
            ),
            type = factor(
                paste0(S1N, S2N, S3N)
            ),
            #S1_log2FC = if_else(S1_vs_N_padj <= padj_thres, S1_log2FC, NA),
            #S2_log2FC = if_else(S2_vs_N_padj <= padj_thres, S2_log2FC, NA),
            #S3_log2FC = if_else(S3_vs_N_padj <= padj_thres, S3_log2FC, NA),
            S1N = case_when(S1N == "+" ~ 1, S1N == "-" ~ -1, TRUE ~ 0),
            S2N = case_when(S2N == "+" ~ 1, S2N == "-" ~ -1, TRUE ~ 0),
            S3N = case_when(S3N == "+" ~ 1, S3N == "-" ~ -1, TRUE ~ 0)
        )

    ### Level definitions
    early_dec_types <- c("-00", "--0", "0-0")
    late_dec_types <- c("---", "0--", "00-", "-0-")
    early_inc_types <- c("+00", "++0", "0+0")
    late_inc_types <- c("+++", "0++", "00+", "+0+")
    other_types   <- c(
        "0+-", "0-+", "-+0", "+-0", "-0+", "+0-",
        "--+", "-+-", "-++", "+-+", "++-"
    )

    ### Relevel for RNA
    plot_mat <- plot_mat %>% mutate(
        type = fct_relevel(type,c(
            early_dec_types, late_dec_types,
            early_inc_types, late_inc_types,
            other_types
        )),
        row_cat = factor(
            case_when(
                type %in% early_dec_types ~ "ED",
                type %in% late_dec_types ~ "LD",
                type %in% early_inc_types ~ "EI",
                type %in% late_inc_types ~ "LI",
                TRUE ~ "O"
            ),
            levels = c("ED", "LD", "EI", "LI", "O")
        )
    ) %>% arrange(type)

    ### Rename the row_cat
    cat_dist <- table(plot_mat$row_cat)
    cat_pct <- round(cat_dist / sum(cat_dist) * 100, digits = 1)

    row_titles <- c(
        "ED" = paste0("Early\ndec. ", col_title, "\n(", cat_dist[["ED"]], ",\n", cat_pct[["ED"]],"%)"),
        "LD" = paste0("Late\ndec. ", col_title, "\n(", cat_dist[["LD"]], ",\n", cat_pct[["LD"]],"%)"),
        "EI" = paste0("Early\ninc. ", col_title, "\n(", cat_dist[["EI"]], ",\n", cat_pct[["EI"]],"%)"),
        "LI" = paste0("Late\ninc. ", col_title, "\n(", cat_dist[["LI"]], ",\n", cat_pct[["LI"]],"%)"),
        "O"  = paste0("Other ", col_title, "\n(", cat_dist[["O"]], ", ", cat_pct[["O"]],"%)")
    )


    plot_mat <- plot_mat 

    col_fun <- circlize::colorRamp2(
            c(-1, 0, 1),
            c("#469CE8", "#FFE680", "#F85B5B")
        )

    column_labels <- c("Stage IB", "Stage II", "Stage III")

    plot_mat %>%
        dplyr::select(
            S1N, S2N, S3N
            # S1_log2FC, S2_log2FC, S3_log2FC, #T_log2FC
        ) %>%
        as.matrix() %>%
        ComplexHeatmap::Heatmap(
            cluster_rows = FALSE, cluster_columns = FALSE,
            row_split = plot_mat$row_cat,
            border = TRUE,
            rect_gp = gpar(lty = 0),
            col = col_fun,
            height = unit(3, "inch"),
            width = unit(0.7, "inch"),
            column_title = paste0(
                omics_name, " (", nrow(plot_mat), "/", nrow(cur_stats),
                " ", col_title, ")"
            ),
            row_title = row_titles,
            row_title_rot = 0,
            row_title_side = "right",
            column_labels = column_labels,
            column_split = c("I", "II", "III"),
            #column_title = " ",
            show_heatmap_legend = TRUE,
            heatmap_legend_param = list(
                title = "Abundance log2FC wrt. normal tissue",
                color_bar = "discrete",
                border = "black",
                at = c(-1, 0, 1),
                labels = c(
                    paste0("Negative (Padj < ", round(padj_thres, 2)," and |log2FC| > ", round(log2FC_thres, 2),")"),
                    paste0("N.S. (Padj > ", round(padj_thres, 2)," or |log2FC| < ", round(log2FC_thres, 2),")"),
                    paste0("Positive (Padj < ", round(padj_thres, 2)," and |log2FC| > ", round(log2FC_thres, 2),")")
                )            
            ),
            use_raster = TRUE,
            raster_quality = 10
        )
}

### Create the output folder
fs::dir_create(here("figures", "lcms"))
fs::dir_create(here("figures", "RNA"))

### Fig. 1e
lcms_p <- PlotStage(
    lcms_stats, LCMS_log2FC_thres, LCMS_padj_thres, "Metabolomics", "MFs"
)
pdf(here("figures", "lcms", "lcms_stages.pdf"), width = 8, height = 6)
draw(lcms_p)
dev.off()

### Fig. 2d
load(here("data", "lcms", "lcms_peaks.Rdata"))
lcms_hit_p <- PlotStage(
    lcms_stats %>% filter(mz_id %in% lcms_hit_peaks$mz_id),
    LCMS_log2FC_thres, LCMS_padj_thres, "Metabolomics", "MFs"
)
pdf(here("figures", "lcms", "lcms_hit_stages.pdf"), width = 6, height = 6)
draw(lcms_hit_p)
dev.off()

### Fig. 1g
RNA_p <- PlotStage(RNA_stats, RNA_log2FC_thres, RNA_padj_thres, "Transcriptomics", "RNAs")
pdf(here("figures", "RNA", "RNA_stages.pdf"), width = 6, height = 6)
draw(RNA_p)
dev.off()

### "--------------------------------------------------------------------------"
### Fig. 1h. Find pathway enrichment
### "--------------------------------------------------------------------------"
library(clusterProfiler)
library(enrichplot)

DEG_ensembl_list <- list(
    "Stage IB" = unique(RNA_stats %>% filter(S1_vs_N_padj <= RNA_padj_thres, abs(S1_log2FC) > RNA_log2FC_thres) %>% pull(gene_id)),
    "Stage II" = unique(RNA_stats %>% filter(S2_vs_N_padj <= RNA_padj_thres, abs(S2_log2FC) > RNA_log2FC_thres) %>% pull(gene_id)),
    "Stage III" = unique(RNA_stats %>% filter(S3_vs_N_padj <= RNA_padj_thres, abs(S3_log2FC) > RNA_log2FC_thres) %>% pull(gene_id))
)

### Rename the labels
RenameGeneSet <- function(msigdb_t2g) {
    msigdb_t2g %>%
    mutate(
        db_name = sub("([^_]+)_.+", "\\1", gs_name),
        gs_name = if_else(startsWith(gs_name, "HALLMARK_"), paste0(gs_name, "_hallmark"), gs_name),
        gs_name = sub("[^_]+_(.+)", "\\1", gs_name),
        gs_name = stringr::str_to_sentence(gs_name),
        gs_name = gsub("process", "proc.", gs_name),
        gs_name = gsub("metabolic", "metab.", gs_name),
        gs_name = gsub("catabolic", "catab.", gs_name),
        gs_name = gsub("cytochrome.p450", "P450", gs_name),
        gs_name = gsub("Tnfa", "TNFa", gs_name),
        gs_name = gsub("line_leucine", "line,_leucine", gs_name),
        gs_name = gsub("nine_aspar", "nine,_aspar", gs_name),
        gs_name = gsub("cine_serine", "cine,_serine", gs_name),
        gs_name = gsub("E2f", "E2F", gs_name),
        gs_name = gsub("g1_s", "G1/S", gs_name),
        gs_name = gsub("G1_s", "G1/S", gs_name),
        gs_name = gsub("g1", "G1", gs_name),    
        gs_name = gsub("G2m", "G2/M", gs_name),
        gs_name = gsub("nfkb", "NFkB", gs_name),
        gs_name = gsub("Wnt", "WNT", gs_name),
        gs_name = gsub("Ppar", "PPAR", gs_name),
        gs_name = gsub("hase_ii", "hase_II", gs_name),
        gs_name = gsub("hase_i", "hase_I", gs_name),
        gs_name = gsub("Fceri", "FceRI", gs_name),
        gs_name = gsub("ca_2", "Ca2+", gs_name),
        gs_name = gsub("c4_and_c2", "C4_and_C2", gs_name),
        gs_name = gsub("Cd22", "CD22", gs_name),
        gs_name = gsub("dna", "DNA", gs_name),
        gs_name = gsub("_bcr_", "_BCR_", gs_name),
        gs_name = gsub("Dna", "DNA", gs_name),
        gs_name = gsub("Ecm", "ECM", gs_name),
        gs_name = gsub("il10", "IL10", gs_name),
        gs_name = gsub("lat2_ntal_lab", "LAT2/NTAL/LAB", gs_name),
        gs_name = gsub("Il2_stat5", "IL2_STAT5", gs_name),
        gs_name = gsub("Mtorc1", "mTORC1", gs_name),
        gs_name = gsub("Kras", "KRAS", gs_name),
        gs_name = gsub("Fcgr3a", "FCGR3A", gs_name),
        gs_name = gsub("Fcgr", "FCGR", gs_name),
        gs_name = gsub("gtpase", "GTPase", gs_name)
    )
}

### Get all the potential gene sets
msigdb_t2g <- rbind(
    msigdbr::msigdbr(
        species = "Homo sapiens", category = "H"
    ) %>% dplyr::select(gs_name, ensembl_gene),
    msigdbr::msigdbr(
        species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG"
    ) %>% dplyr::select(gs_name, ensembl_gene),
    msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
        dplyr::select(gs_name, ensembl_gene),
    msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
        dplyr::select(gs_name, ensembl_gene)
) 

msigdb_t2g <- RenameGeneSet(msigdb_t2g)

### Find enriched terms
Stage_MsigDB <- compareCluster(
    geneCluster = DEG_ensembl_list,
    fun = "enricher",
    TERM2GENE = msigdb_t2g[, c("gs_name", "ensembl_gene")],
    pAdjustMethod = "BH",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05
)

Stage_MsigDB <- setReadable(
    Stage_MsigDB, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENSEMBL"
)

PlotEnrich <- function(MsigDB, pval_thres, left_margin, show_category){
    dotplot(
        MsigDB %>% filter(pvalue < pval_thres),
        showCategory = show_category,
        color = "pvalue",
        label_format = 50
    ) + ylab("") + xlab("") +
    geom_count(aes(color = -log10(pvalue))) +
    scale_size_area(max_size = 6) +
    scale_colour_gradientn(
        name = "-log10(P)",
        limits = c(-log10(pval_thres), 40),
        colors = c("blue", "red")
    ) +
    theme(
        axis.text.y = element_text(lineheight = 0.65),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.margin = margin(t=0, r=0, b=0, l=left_margin, "in")
        #legend.direction = "horizontal", 
        #legend.position = "bottom",
        #legend.box = "horizontal"
    )
}

### Fig. 1h. Plot the top pathways
SavePlot(
    PlotEnrich(Stage_MsigDB, pval_thres = 1e-12, left_margin = 2.2, show_category = 14),
    width = 8,
    height = 6.5,
    dir_name = "RNA",
    file_name = "Stage_Pathway_General_Top"
)

### Fig. S3a. Plot the all pathways
SavePlot(
    PlotEnrich(Stage_MsigDB, pval_thres = 1e-12, left_margin = 1.9, show_category = 50),
    width = 8,
    height = 15,
    dir_name = "RNA",
    file_name = "Stage_Pathway_General_All"
)

### Get all the potential gene sets
metabo_t2g <- rbind(
    msigdbr::msigdbr(
        species = "Homo sapiens", category = "H"
    ) %>% dplyr::select(gs_name, ensembl_gene),
    msigdbr::msigdbr(
        species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG"
    ) %>% dplyr::select(gs_name, ensembl_gene),
    msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
        dplyr::select(gs_name, ensembl_gene),
    msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
        dplyr::select(gs_name, ensembl_gene)
) %>% filter(
    grepl("metab|lipid|catab|tty_acid|synth|egrad|TCA|genesis", gs_name, ignore.case = TRUE)
)

metabo_t2g <- RenameGeneSet(metabo_t2g)

### Find enriched terms
Metabo_MsigDB <- compareCluster(
    geneCluster = DEG_ensembl_list,
    fun = "enricher",
    TERM2GENE = metabo_t2g[, c("gs_name", "ensembl_gene")],
    pAdjustMethod = "BH",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05
)

Metabo_MsigDB <- setReadable(
    Metabo_MsigDB, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENSEMBL"
)

### Fig. 1i. Plot the top pathways
SavePlot(
    PlotEnrich(Metabo_MsigDB, pval_thres = 1e-5, left_margin = 2.22, show_category = 16),
    width = 8,
    height = 6.5,
    dir_name = "RNA",
    file_name = "Stage_Pathway_Meta_Top"
)

### Fig. S3b. Plot all the pathways 
SavePlot(
    PlotEnrich(Metabo_MsigDB, pval_thres = 1e-5, left_margin = 2.11, show_category = 40),
    width = 8,
    height = 15,
    dir_name = "RNA",
    file_name = "Stage_Pathway_Meta_All"
)

### Find out total number of processes
length(unique(as_tibble(Stage_MsigDB) %>% pull(ID)))
length(unique(as_tibble(Metabo_MsigDB) %>% pull(ID)))

### Save as csv
readr::write_excel_csv(as_tibble(Stage_MsigDB), here("outputs", "Stage_MsigDB.csv"))
readr::write_excel_csv(as_tibble(Metabo_MsigDB), here("outputs", "Metabo_MsigDB.csv"))
