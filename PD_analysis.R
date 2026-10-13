# Parkinson's disease vs control: gene expression in the substantia nigra
# Data: GSE7621 (9 controls, 16 PD), already downloaded by Script 1

library(GEOquery)
library(limma)
library(tidyverse)
library(ggrepel)

# STEP 1: load the data
gse <- readRDS("data/gse7621.rds")

# STEP 2: label each sample as Control or PD
group <- ifelse(pData(gse)$characteristics_ch1 == "Parkinson's Disease", "PD", "Control")
group <- factor(group, levels = c("Control", "PD"))

# STEP 3: clean the numbers
# remove missing values, then log2-transform (like log scaling in DESeq2's VST)
expr <- exprs(gse)
expr <- expr[complete.cases(expr), ]
expr <- log2(pmax(expr, 1))

# STEP 4: replace probe IDs with gene names
genes <- fData(gse)[rownames(expr), "Gene Symbol"]
genes <- sub(" ///.*", "", genes)            # keep the first name if there are several
expr  <- expr[genes != "" & !is.na(genes), ]
genes <- genes[genes != "" & !is.na(genes)]
expr  <- avereps(expr, ID = genes)            # average probes that measure the same gene
# make samples comparable, then remove genes with very low expression
expr <- normalizeBetweenArrays(expr, method = "quantile")
expr <- expr[rowMeans(expr) > median(rowMeans(expr)), ]

# STEP 5: compare PD vs Control (limma works like DESeq2, but for microarray data)
design <- model.matrix(~ group)
fit    <- eBayes(lmFit(expr, design))
res    <- topTable(fit, coef = "groupPD", number = Inf) %>% rownames_to_column("gene")
write.csv(res, "results_PD_vs_control.csv", row.names = FALSE)

# how many genes changed? (adjusted p < 0.05 and at least 1.5-fold change)
sig <- res %>% filter(adj.P.Val < 0.05, abs(logFC) > log2(1.5))
cat("Genes down in PD:", sum(sig$logFC < 0), "\n")
cat("Genes up in PD:  ", sum(sig$logFC > 0), "\n")

# STEP 6: sanity check - dopamine-neuron genes should be LOWER in PD,
# because these neurons die in Parkinson's disease
markers <- c("TH", "SLC6A3", "DDC", "SLC18A2")
print(res %>% filter(gene %in% markers) %>% select(gene, logFC, adj.P.Val))

# STEP 7: volcano plot
res$status <- ifelse(res$logFC < 0, "Lower in PD (n.s.)", "Higher in PD (n.s.)")
res$status[res$gene %in% sig$gene & res$logFC < 0] <- "Down in PD"
res$status[res$gene %in% sig$gene & res$logFC > 0] <- "Up in PD"

volcano <- ggplot(res, aes(logFC, -log10(P.Value), colour = status)) +
  geom_point(size = 1, alpha = 0.6) +
  scale_colour_manual(values = c("Down in PD" = "#1F4E9C",
                                    "Up in PD" = "#C0392B",
                                    "Lower in PD (n.s.)" = "#BFD3EC",
                                    "Higher in PD (n.s.)" = "#F2C4BE")) +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed", colour = "grey50") +
  geom_text_repel(data = filter(res, gene %in% c(markers, "ALDH1A1", "SV2C", "P2RX7", "PNPT1")),
                  aes(label = gene), colour = "black", size = 3.5) +
  labs(title = "Parkinson's disease vs control (GSE7621)",
       x = "log2 fold change", y = "-log10 p-value", colour = NULL) +
  theme_bw()
ggsave("volcano_PD.png", volcano, width = 7, height = 5)

# STEP 8: boxplots of key genes, one dot per person
key_genes <- c("SLC18A2", "ALDH1A1", "SV2C", "TH", "SLC6A3", "P2RX7")

box_data <- as.data.frame(t(expr[key_genes, ])) %>%
  mutate(group = group) %>%
  pivot_longer(-group, names_to = "gene", values_to = "expression") %>%
  mutate(gene = factor(gene, levels = key_genes))

boxplot_pd <- ggplot(box_data, aes(group, expression, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.5) +
  facet_wrap(~ gene, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c("Control" = "#5B9BD5", "PD" = "#E07B6A")) +
  labs(title = "Key genes in Parkinson's disease vs control",
       x = NULL, y = "log2 expression") +
  theme_bw() +
  theme(legend.position = "none")

ggsave("boxplots_key_genes.png", boxplot_pd, width = 9, height = 5.5)

# STEP 9: heatmap of the 16 significant genes
library(pheatmap)

heat_genes <- sig %>% arrange(logFC) %>% pull(gene)      # order: most down to most up
heat_data  <- expr[heat_genes, ]

# label each column (person) as Control or PD, and sort columns by group
sample_info <- data.frame(Group = group, row.names = colnames(heat_data))
heat_data   <- heat_data[, order(group)]

pheatmap(heat_data,
         scale = "row",                       # z-score each gene so genes are comparable
         cluster_rows = FALSE,                # keep genes in fold-change order
         cluster_cols = FALSE,                # keep Control and PD as separate blocks
         annotation_col = sample_info,
         annotation_colors = list(Group = c(Control = "#5B9BD5", PD = "#E07B6A")),
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         show_colnames = FALSE,
         gaps_col = sum(group == "Control"),  # white gap between the two groups
         main = "16 significant genes: Parkinson's disease vs control",
         filename = "heatmap_significant_genes.png",
         width = 8, height = 6)
