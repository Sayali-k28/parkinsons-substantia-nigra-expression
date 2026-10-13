# =====================================================================
# Project: A reproducible transcriptomic signature of Parkinson's disease
#          in the human substantia nigra
#
# Script 1 of 3: Download data and inspect sample groups
#
# Discovery dataset : GSE7621  (substantia nigra, 16 PD vs 9 controls)
# Validation dataset: GSE49036 (substantia nigra, PD and controls)
# Both use the same platform: Affymetrix HG-U133 Plus 2.0 (GPL570)
# =====================================================================

# ---- 1. Install packages (run once) ----
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc_pkgs <- c("GEOquery", "limma", "clusterProfiler", "org.Hs.eg.db", "enrichplot")
for (p in bioc_pkgs) if (!require(p, character.only = TRUE, quietly = TRUE)) BiocManager::install(p)
cran_pkgs <- c("tidyverse", "pheatmap", "ggrepel", "pROC")
for (p in cran_pkgs) if (!require(p, character.only = TRUE, quietly = TRUE)) install.packages(p)

library(GEOquery)
library(tidyverse)

# ---- 2. Create project folders ----
dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ---- 3. Download both datasets from GEO (processed expression + metadata) ----
# This can take a few minutes the first time
gse_disc <- getGEO("GSE7621",  GSEMatrix = TRUE, destdir = "data")[[1]]
gse_val  <- getGEO("GSE49036", GSEMatrix = TRUE, destdir = "data")[[1]]

# ---- 4. Basic checks: dimensions and value range ----
# Values around 2-16 mean data are already log2-transformed;
# values in the thousands mean we need to log2-transform them in Script 2
for (nm in c("disc", "val")) {
  eset <- if (nm == "disc") gse_disc else gse_val
  cat("\n==========", nm, "==========\n")
  cat("Probes x samples:", nrow(eset), "x", ncol(eset), "\n")
  print(summary(as.vector(exprs(eset))))
}

# ---- 5. Show the sample information columns that describe groups ----
# Send me this output: I will use it to define PD vs control correctly
show_groups <- function(eset, label) {
  pd <- pData(eset)
  cols <- grep("characteristics|source_name|title", colnames(pd), value = TRUE)
  cat("\n==========", label, ": group-related columns ==========\n")
  for (cl in cols) {
    cat("\n--", cl, "--\n")
    print(table(pd[[cl]]))
  }
}
show_groups(gse_disc, "GSE7621 (discovery)")
show_groups(gse_val,  "GSE49036 (validation)")

# ---- 6. Save for the next script ----
saveRDS(gse_disc, "data/gse7621.rds")
saveRDS(gse_val,  "data/gse49036.rds")
cat("\nSaved both datasets in the data/ folder. Send the output above to continue.\n")
