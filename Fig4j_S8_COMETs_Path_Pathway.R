### "========================================================================="
### Plot COMET's Path analysis results
###
### Copyright (c) 2021-2024. Bioinformatics Institute, A*STAR, Singapore
### License: MIT License. See LICENSE for the full license.
### Authors: Lit-Hsin Loo (loolh@bii.a-star.edu.sg)
### "========================================================================="
library(here)
source(here("lib", "MSI_lib.R"))
source(here("lib", "Study_lib.R"))
library(clusterProfiler)
library(enrichplot)

### MSIG
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
) %>%
    ### Capitalize the name
    mutate(
        gs_name = sub("[^_]+_(.+)", "\\1", gs_name),
        gs_name = stringr::str_to_sentence(gs_name),
        gs_name = gsub("cytochrome_p450", "P450", gs_name),
        gs_name = gsub("Cytochrome_p450", "P450", gs_name),
        gs_name = gsub("line_leucine", "line,_leucine", gs_name),
        gs_name = gsub("Tnfa", "TNFa", gs_name),
        gs_name = gsub("E2f", "E2F", gs_name),
        gs_name = gsub("G2m", "G2/M", gs_name),
        gs_name = gsub("nfkb", "NFkB", gs_name),
        gs_name = gsub("Ngf", "NGF", gs_name),
        gs_name = gsub("Ppar", "PPAR", gs_name),
        gs_name = gsub("Dna", "DNA", gs_name),
        gs_name = gsub("dna", "DNA", gs_name),
        gs_name = gsub("G1_s", "G1/S", gs_name),
        gs_name = gsub("pre_repli", "pre-repli", gs_name),
        gs_name = gsub("Alpha_lino", "α-lino", gs_name),
        gs_name = gsub("7alpha", "7α", gs_name),
        gs_name = gsub("of_atr_", "of_ATR_", gs_name),
        gs_name = gsub("tors_spms", "tors_(SPMs)", gs_name),
        gs_name = gsub("Ecm", "ECM", gs_name),
        gs_name = gsub("Il2_stat5", "IL2_STAT5", gs_name),
        gs_name = gsub("gtpase", "GTPase", gs_name),
        gs_name = gsub("Metabolic", "Metab.", gs_name),
        gs_name = gsub("hase_i_", "hase_I_", gs_name),
        gs_name = gsub("hase_ii_", "hase_II_", gs_name),
        gs_name = gsub("srebf_srebp", "SREBF (SREBP)", gs_name),
        gs_name = gsub("mapk", "MAPK", gs_name),
        gs_name = gsub("Fcgamma_receptor_fcgr_dependent", "FcγR-dependent", gs_name),
        gs_name = gsub("Fcgr3a mediated", "FcγR3a-mediated", gs_name),
        gs_name = gsub("Fcgr", "FcγR", gs_name),
        gs_name = gsub("fcgr", "FcγR", gs_name),
        gs_name = gsub("Tp53", "TP53", gs_name),
        gs_name = gsub("Fceri_mediated", "FcεRI-mediated", gs_name),
        gs_name = gsub("il10", "IL10", gs_name),
        gs_name = gsub("_ca_2_", "_Ca2+_", gs_name),
        gs_name = gsub("Chk1_chk2_cds1", "Chk1/Chk2(Cds1)", gs_name),
        gs_name = gsub("cyclin_b_cdk1", "Cyclin_B:Cdk1", gs_name),
        gs_name = gsub("Nadph_", "NADPH_", gs_name),
        gs_name = gsub("Nadp_", "NADP_", gs_name),        
        gs_name = gsub("lat2_ntal_lab", "LAT2/NTAL/LAB", gs_name),
        gs_name = gsub("srebp_srebf", "SREBP (SREBF)", gs_name),        
        gs_name = gsub("_c_strand_lagging_strand", "_C-strand_(lagging strand)_", gs_name)
    )

### "========================================================================="
### Plot pathway enrichment for top VIP genes
### "========================================================================="
CEG_list_vip <- readRDS(here("results", "RNA", "CEG_list_vip.rds"))

### Split into different id types
CEG_ensembl_list <- lapply(CEG_list_vip, function(x) x$CEG_ensembls)
CEG_entrez_list <- lapply(CEG_list_vip, function(x) x$CEG_entrezs)

### Capitalize the model name
#names(CEG_ensembl_list) <- tools::toTitleCase(names(CEG_ensembl_list))

### Re-arrange the models
CEG_ensembl_list <- CEG_ensembl_list[c(
    "Normal", #"Transformed", 
    "G2", "G3",  
    "Necrotic",  "Fibrotic", "Steatotic"
)]

names(CEG_ensembl_list) <- c(
    "ME-normal\nregions",
    #"ME-transformed\nregions", 
    "ME-low-grade\nregions", "ME-high-grade\nregions",
    "ME-necrotic\nregions",
    "ME-fibrotic\nregions", "ME-steatotic\nregions"
)

sapply(CEG_ensembl_list, length)

cc_vip_msigdb <- compareCluster(
    geneCluster = CEG_ensembl_list,
    fun = "enricher",
    #universe = unique(rownames(dds)),
    TERM2GENE = msigdb_t2g[, c("gs_name", "ensembl_gene")],
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.10
)
cc_vip_msigdb <- setReadable(
    cc_vip_msigdb,
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    keyType = "ENSEMBL"
)

### Save as csv
readr::write_excel_csv(
    as_tibble(cc_vip_msigdb),
    here("outputs", "SgME_Map_MsigDB.csv")
)

### Capitalize the pathway name
pval_thres <- 0.05
vip_enrich_short_p <- dotplot(
        cc_vip_msigdb %>% filter(pvalue < pval_thres),
        showCategory = 5, label_format = 40,
        color = "pvalue"
    ) +
    xlab("") + ylab("") +
    geom_count(aes(fill = -log10(pvalue))) +
    scale_size_area(max_size = 5.1) +
    scale_fill_gradientn(
        name = "-log10(P)",
        limits = c(0, 80),
        colors = c("blue", "red")
    ) +
    #coord_flip() +
    #scale_x_discrete(limits = rev, labels = ~ gsub("\n", " ", .x)) +
    #scale_y_discrete(limits = rev, labels = ~ gsub("_", " ", .x)) +
    scale_x_discrete(labels = ~ gsub("\n", " ", .x)) +
    scale_y_discrete(labels = ~ gsub("_", " ", .x)) +    
    theme(
        #axis.text.y = element_text(lineheight = 0.7),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        #plot.margin = margin(t=0, r=0, b=0, l=0, "in")
        #legend.direction = "horizontal", 
        #legend.position = "top",
        #legend.box = "horizontal"
    )

### Save the plots
SavePlot(
    vip_enrich_short_p,
    width = 7.75,
    height = 6.75,
    dir_name = "RNA",
    file_name = "vip_enriched_pathway_short"
)

vip_enrich_full_p <- dotplot(
        cc_vip_msigdb %>% filter(pvalue < pval_thres),
        showCategory = 10, label_format = 40,
        color = "pvalue"
    ) +
    xlab("") + ylab("") +
    geom_count(aes(fill = -log10(pvalue))) +
    scale_size_area(max_size = 5.5) +
    scale_fill_gradientn(
        name = "-log10(P)",
        limits = c(0, 80),
        colors = c("blue", "red")
    ) +
    #coord_flip() +
    #scale_x_discrete(limits = rev, labels = ~ gsub("\n", " ", .x)) +
    #scale_y_discrete(limits = rev, labels = ~ gsub("_", " ", .x)) +
    scale_x_discrete(labels = ~ gsub("\n", " ", .x)) +
    scale_y_discrete(labels = ~ gsub("_", " ", .x)) +    
    theme(
        #axis.text.y = element_text(lineheight = 0.7),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
        #plot.margin = margin(t=0, r=0, b=0, l=0, "in")
        #legend.direction = "horizontal", 
        #legend.position = "top",
        #legend.box = "horizontal"
    )

### Save the plots
SavePlot(
    vip_enrich_full_p,
    width = 8,
    height = 11,
    dir_name = "RNA",
    file_name = "vip_enriched_pathway_full"
)