library(here)
source(here("lib", "MSI_lib.R"))

# https://www.bioconductor.org/packages/release/bioc/vignettes/Cardinal/inst/doc/Cardinal-2-stats.html

### Define the dataset
msi_conf <- read_json(here("conf", "msi_conf.json"))

peak_summary_mode <- "lcms_only"

### "========================================================================"
### We need to find common peaks that are the highest and the most
### homogenous one
### "========================================================================"
### Generate the stat scatter plot
PlotTopCommonTissuePeaks(
    ion_type = "pos",
    align_res = top_common_align_reses[["pos"]],
    peak_summary_mode = peak_summary_mode,
    msi_conf
)

common_promi_peaks <- FindTopCommonProminentPeaks(
    ion_type = "pos",
    align_res = top_common_align_reses[["pos"]],
    msi_conf
)

### Generate the stat scatter plot
PlotTopCommonTissuePeaks(
    ion_type = "neg",
    align_res = top_common_align_reses[["neg"]],
    peak_summary_mode = peak_summary_mode,
    msi_conf
)

common_promi_peaks <- FindTopCommonProminentPeaks(
    ion_type = "neg",
    align_res = top_common_align_reses[["neg"]],
    msi_conf
)

### "========================================================================"
### Check the average levels of all the peaks across all the reference peaks
### across all the slides and sections to:
### 1) Flag out potential outliers and remove them from further analysis
### 2) Determine the spectra normalization coefficients. We may normalize
###    at the pixel, section, or slide levels
### 3) Save the norm_coeff
### "========================================================================"
norm_coeff <- FindNormCoeff(
    ion_type = "pos",
    peak_summary_mode = peak_summary_mode,
    msi_conf, is_plot = TRUE
)

norm_coeff <- FindNormCoeff(
    ion_type = "neg",
    peak_summary_mode = peak_summary_mode,
    msi_conf, is_plot = TRUE
)