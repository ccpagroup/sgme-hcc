### Study configurations
require(here)

### Define the location of the STx profiles for F008_CA and F011_CA
### visium_raw_dir <- ""

peak_summary_mode_study <- "lcms_only"
is_normalized_study <- TRUE

### Reference slides for training (each of them must have normal sections)
ref_slide_ids <- c(
    "A-01-01",
    #"A-02-01", ##  HEP0194 - Control has no hepatocytes
    "A-03-01",
    "A-04-01",
    # A-05-01,  ## CCA
    #"A-06-01", ## Failed QC
    "A-07-01",
    #"A-08-01", ## HEP0319 - Too much necrosis 
    "A-09-01",
    "A-10-01",
    "A-11-01"
)

test_slide_ids <- c(
    "A-01-01",
    "A-02-01", ## HEP0194 - Control has no hepatocytes
    "A-03-01",
    "A-04-01", 
    # A-05-01,  ## CCA
    #"A-06-01", ## Failed QC
    "A-07-01",
    "A-08-01", ## HEP0319 - Too much necrosis 
    "A-09-01",
    "A-10-01",
    "A-11-01"
)

### Align resolutions for top common peaks
top_common_align_reses <- c(
    "pos" = 50,
    "neg" = 30
)

### Selected CPPs
top_common_peak_idx <- list(
    pos = c(
        4, # 104.1
        30, # 380.3
        31, # 429.2
        36, # 607.5  <-
        37, # 617.5  <-
        38, # 618.5  <-
        39, # 643.5  <-
        40, # 707.5
        41, # 728.5
        42 # 781.6
    ),
    neg = c(
        21, # 187.1
        35, # 253.2
        37, # 255.2
        38, # 256.2
        46, # 281.3
        47, # 282.3
        56, # 742.5
        57, # 773.5
        60, # 857.5
        61  # 887.6
    )
)

### Bad samples
bad_sample_ids <- c()

vip_thres_study <- 0.95
coeff_thres_study <- 0.03

### Metabolite class colors
metabolite_classes_colors <- c(
    Carbonyls = "#8dd3c7",
    Cholines = "#ffffb3",
    "Fatty acyls" = "#bebada",
    Glycerolipids = "#fb8072",
    Glycosides = "#fdb462",
    Glycerophospholipids = "#80b1d3",
    Others = "#ccebc5",
    "Prenol lipids" = "#b3de69",
    Sphingolipids = "#fccde5",
    Steroids = "#bc80bd",
    Unannotated = "#d9d9d9",
    "Amino acids" = "#bc80bd",
    "Nucleotides" = "#ccebc5"
    # #ffed6f (backup color)
)

### ROI class colors
roi_annot_colors <- c(
    Normal = "#b8b8b8",
    Fibrotic = "#40ff13",
    G3 = "#ff3c31",
    G2 = "#F9A603",
    Necrotic = "#fd61f9",
    Steatotic = "#42fff6"
)

steatotic_annot_colors <- c(
    Steatotic = "#42fff6",
    Nonsteatotic = "white"
)

fibrotic_annot_colors <- c(
    Fibrotic = "#40ff13",
    Nonfibrotic = "white"
)

case_annot_colors <- c(
    B003 = "#1b9e77",
    B008 = "#d95f02",
    HEP0152 = "#7570b3",
    HEP0194 = "#377eb8",
    HEP0214 = "#e7298a",
    HEP0268 = "#66a61e",
    HEP0277 = "#e6ab02",
    HEP0321 = "#a6761d",
    HEP0319 = "#999999"
)

### Thresholds
RNA_padj_thres <- 0.05
RNA_log2FC_thres <- log2(1.2)

LCMS_padj_thres <- 0.05
LCMS_log2FC_thres <- log2(1.2)

### Load local conf if available
if (file.exists(here("conf", "study_conf_local.R"))) {
    source(here("conf", "study_conf_local.R"))
}

### Warn the users if key variables are not defined
if (!exists(visium_raw_dir) || visium_raw_dir == "") {
    stop(
        "visium_raw_dir is not defined. "
        "Please update conf/study_conf.R to the downloaded STX_profiles folder."
    )
}