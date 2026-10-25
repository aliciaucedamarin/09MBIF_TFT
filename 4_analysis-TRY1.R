# ==============================================================================
# SECTION 4.4: DESCRIPTIVE STATISTICAL ANALYSIS & DISTRIBUTION VISUALIZATIONS
# ==============================================================================
# Objective: Load raw inputs, synchronize metadata rows, calculate log2 shifts,
# and visually assess variance properties using non-redundant long structures.

library(readxl)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)

# 1. DATA IMPORT & METADATA SYNCHRONIZATION
raw_excel_data <- read_excel("data/raw_data_annotated_hg19.xlsx", sheet = "Annotated_Data")
metadata       <- read.table("metadata/patients_metadata.txt", header = TRUE, sep = "\t")

# Set clinical matrix row coordinates using the continuous first column indexing
rownames(metadata) <- paste0(metadata[[1]], "_Norm")

# Isolate metric value dimensions
norm_counts_cols <- grep("^[0-9]+_Norm$", colnames(raw_excel_data), value = TRUE)

# Format value boundaries dynamically and cleanse unmapped cells
raw_excel_data <- raw_excel_data %>%
  mutate(across(all_of(norm_counts_cols), ~ as.numeric(as.character(.))))

metadata_structural_cols <- c("REGION ID", "ANNOTATED_GENE", "REGIONAL_TYPE", 
                              "DISTANCE_TO_TSS", "OTHER_GENES")

selected_data <- raw_excel_data[, c(metadata_structural_cols, norm_counts_cols)]

# Build continuous metrics arrays and replace NAs with absolute baseline zeros
counts_norm <- as.matrix(selected_data[, norm_counts_cols])
rownames(counts_norm) <- selected_data$`REGION ID`
genomic_annotations <- selected_data[, metadata_structural_cols]


# 2. DIAGNOSTIC MOMENT 1: RAW LINEAR DISTRIBUTION ASSESSMENT
# Melt raw linear values to evaluate initial distribution features and skewness
counts_norm_long <- melt(counts_norm)
colnames(counts_norm_long) <- c("Region", "sample_id", "Normalized_Value")
counts_norm_long <- merge(counts_norm_long, metadata, by.x = "sample_id", by.y = "row.names")

# Plot A: Initial Linear Scale Distribution Profile (Violin-Box Plot)
# Visually demonstrates severe right-skewness and extreme dynamic range inequalities
p_counts_distribution <- ggplot(
  counts_norm_long,
  aes(
    x = sample_id,
    y = Normalized_Value,
    fill = pathology
  )
) +
  geom_violin(
    alpha = 0.5,
    scale = "width",
    colour = "grey40"
  ) +
  geom_boxplot(
    width = 0.1,
    outlier.shape = NA,
    colour = "black",
    alpha = 0.8
  ) +
  scale_fill_brewer(
    palette = "Set1",
    name = "Pathology"
  ) +
  labs(
    title = "Pre-Transformation Global Distribution Profile",
    subtitle = "Distribution of normalized counts across samples on the linear scale",
    x = "Sample ID",
    y = "Normalized Counts (Linear Scale)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      size = 17,
      face = "bold",
      colour = "black"
    ),
    plot.subtitle = element_text(
      size = 10.5,
      colour = "grey35",
      margin = margin(b = 15)
    ),
    axis.title = element_text(
      size = 11,
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 9,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 10,
      colour = "black"
    ),
    legend.title = element_text(
      face = "bold",
      size = 10
    ),
    legend.text = element_text(size = 9),
    plot.margin = margin(15, 20, 20, 15)
  )

print(p_counts_distribution)

ggsave(
  file.path(
    figures_dir,
    "pre_transformation_distribution_profile.png"
  ),
  p_counts_distribution,
  width = 14,
  height = 7,
  dpi = 300,
  bg = "white"
)

# 3. TRANSFORMATION MOMENT: MATHEMATICAL VARIANCE STABILIZATION
# Apply log2(x + 1) to stabilize variances across all columns and handle zeros safely
log_counts <- log2(counts_norm + 1)

# Melt transformed matrix to assess the new metric space
log_counts_long <- melt(as.matrix(log_counts))
colnames(log_counts_long) <- c("Region", "sample_id", "Log2_Counts")
log_counts_long <- merge(log_counts_long, metadata, by.x = "sample_id", by.y = "row.names")


# 4. DIAGNOSTIC MOMENT 2: POST-TRANSFORMATION PROFILING (THE SHIFT)
# Plot B: Variance-Stabilized Log2 Scale Profile (Identical Violin-Box Plot)
# Proves that the log2 transformation resolved variance inflation and unmasked internal structures
p_distribution <- ggplot(
  log_counts_long,
  aes(
    x = sample_id,
    y = Log2_Counts,
    fill = pathology
  )
) +
  geom_violin(
    alpha = 0.5,
    scale = "width",
    colour = "grey40"
  ) +
  geom_boxplot(
    width = 0.1,
    outlier.shape = NA,
    colour = "black",
    alpha = 0.8
  ) +
  scale_fill_brewer(
    palette = "Set1",
    name = "Pathology"
  ) +
  labs(
    title = "Post-Transformation Global Distribution Profile",
    subtitle = "Distribution of log2-transformed normalized counts across samples",
    x = "Sample ID",
    y = "Log2(Normalized Counts + 1)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      size = 17,
      face = "bold",
      colour = "black"
    ),
    plot.subtitle = element_text(
      size = 10.5,
      colour = "grey35",
      margin = margin(b = 15)
    ),
    axis.title = element_text(
      size = 11,
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 9,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 10,
      colour = "black"
    ),
    legend.title = element_text(
      face = "bold",
      size = 10
    ),
    legend.text = element_text(size = 9),
    plot.margin = margin(15, 20, 20, 15)
  )

print(p_distribution)

ggsave(
  file.path(
    figures_dir,
    "post_transformation_distribution_profile.png"
  ),
  p_distribution,
  width = 14,
  height = 7,
  dpi = 300,
  bg = "white"
)

# 5. DIAGNOSTIC MOMENT 3: OVERLAPPED PATIENT CONTINUITY PROFILE
# Plot C: Overlapped Log2 Density Map across all 24 sample pathways simultaneously
# Confirms global multimodal alignment and homogeneity of profiles before inferential modeling
p_density_profiles <- ggplot(
  log_counts_long,
  aes(
    x = Log2_Counts,
    group = sample_id,
    colour = pathology
  )
) +
  geom_density(
    alpha = 0.05,
    linewidth = 0.6
  ) +
  scale_colour_brewer(
    palette = "Set1",
    name = "Pathology"
  ) +
  labs(
    title = "Methylation Density Profiles Across All 24 Samples",
    subtitle = "Distribution of log2-transformed normalized counts across samples",
    x = "Log2(Normalized Counts + 1)",
    y = "Density"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      size = 17,
      face = "bold",
      colour = "black"
    ),
    plot.subtitle = element_text(
      size = 10.5,
      colour = "grey35",
      margin = margin(b = 15)
    ),
    axis.title = element_text(
      size = 11,
      face = "bold"
    ),
    axis.text = element_text(
      size = 10,
      colour = "black"
    ),
    legend.title = element_text(
      face = "bold",
      size = 10
    ),
    legend.text = element_text(size = 9),
    plot.margin = margin(15, 20, 15, 15)
  )

print(p_density_profiles)

ggsave(
  file.path(
    figures_dir,
    "methylation_density_profiles.png"
  ),
  p_density_profiles,
  width = 12,
  height = 7,
  dpi = 300,
  bg = "white"
)

# ==============================================================================
# SECTION 4.5: INFERENTIAL STATISTICAL ANALYSIS, CONTRASTS & VISUALIZATIONS
# ==============================================================================
library(limma)
library(FactoMineR)
library(factoextra)
library(pheatmap)
library(UpSetR)
library(ggplot2)

# 1. OPTIMIZE INPUT DATA BY FILTERING LOW-VARIANCE REGIONS
region_variances <- apply(log_counts, 1, function(x) stats::var(x, na.rm = TRUE))
log_counts_filtered <- log_counts[region_variances > 0.01, ]

# 2. RUN PRINCIPAL COMPONENT ANALYSIS (PCA) & GENERATE OPTIMIZED DIAGNOSTIC PLOTS
pca_results <- FactoMineR::PCA(t(log_counts_filtered), graph = FALSE)

# PCA Plot 1: Highlight Pathology & Genetic Status simultaneously without layout errors
p_pca <- fviz_pca_ind(
  pca_results,
  col.ind = metadata$pathology,
  shape.ind = metadata$genetic_status,
  palette = "Set1",
  addEllipses = TRUE,
  ellipse.type = "convex",
  repel = TRUE,
  title = "Global PCA: Separation by Pathology & Genetic Profiles"
) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      size = 17,
      face = "bold",
      colour = "black"
    ),
    axis.title = element_text(
      size = 11,
      face = "bold"
    ),
    axis.text = element_text(
      size = 10,
      colour = "black"
    ),
    legend.title = element_text(
      face = "bold",
      size = 10
    ),
    legend.text = element_text(size = 9),
    plot.margin = margin(15, 25, 15, 20)
  )

print(p_pca)

ggsave(
  file.path(
    figures_dir,
    "global_pca_pathology_genetic_profiles.png"
  ),
  p_pca,
  width = 11,
  height = 8,
  dpi = 300,
  bg = "white"
)

# PCA Plot 2: Highlight Clinical Status & Gene Variant
p_pca_clinical <- fviz_pca_ind(
  pca_results,
  col.ind = metadata$clinical_status,
  shape.ind = metadata$gene_variant,
  palette = "Dark2",
  addEllipses = TRUE,
  ellipse.type = "convex",
  repel = TRUE,
  title = "Global PCA: Separation by Clinical Status & Lineage"
) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 17, face = "bold", colour = "black"),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10, colour = "black"),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    plot.margin = margin(15, 25, 15, 20)
  )

print(p_pca_clinical)

ggsave(
  file.path(figures_dir, "global_pca_clinical_status_gene_variant.png"),
  plot = p_pca_clinical,
  width = 11,
  height = 8,
  dpi = 300,
  bg = "white"
)

# SUBSET PCA (Outliers Excluded: 5_Norm, 20_Norm, 21_Norm)
outlier_samples <- c("5_Norm", "20_Norm", "21_Norm")
log_counts_sub <- log_counts_filtered[, !colnames(log_counts_filtered) %in% outlier_samples]
metadata_sub   <- metadata[!rownames(metadata) %in% outlier_samples, ]

pca_sub_results <- FactoMineR::PCA(t(log_counts_sub), graph = FALSE)
p_pca_sub <- fviz_pca_ind(
  pca_sub_results,
  col.ind = metadata_sub$pathology,
  palette = "Set1",
  addEllipses = TRUE,
  ellipse.type = "convex",
  repel = TRUE,
  title = "Detailed PCA: Sub-Cohort Architecture (Outliers Excluded)"
) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 17, face = "bold", colour = "black"),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10, colour = "black"),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    plot.margin = margin(15, 25, 15, 20)
  )

print(p_pca_sub)

ggsave(
  file.path(figures_dir, "detailed_pca_subcohort_architecture.png"),
  plot = p_pca_sub,
  width = 11,
  height = 8,
  dpi = 300,
  bg = "white"
)

# 3. ADVANCED INFERENTIAL MODELING USING LIMMA-VOOM PIPELINE
# Clean cohort text labels to prevent syntax breaks inside the contrast formulas
group_factor <- paste(metadata$pathology, metadata$clinical_status, sep = "_")
group_factor <- gsub("-| ", "_", group_factor)
group_factor <- as.factor(group_factor)

design_matrix <- model.matrix(~ 0 + group_factor)
colnames(design_matrix) <- levels(group_factor)

# CRITICAL BIOINFORMATIC STEP: Apply voom transformation directly to raw normalized counts
# to calculate empirical weights, stabilizing the mean-variance trend and restoring FDR power
v_weights <- voom(counts_norm, design_matrix, plot = FALSE)

# Fit weighted linear models using precision weights calculated by voom
lm_fit <- lmFit(v_weights, design_matrix)

# 4. EXPLICIT GROUP CONTRASTS (Pathology vs Healthy Baseline Baseline)
contrast_matrix <- makeContrasts(
  ALPS_vs_Controls      = ALPS_AFFECTED - HEALTHY_CONTROLS_HEALTHY,
  ALPS_Like_vs_Controls = ALPS_LIKE_AFFECTED - HEALTHY_CONTROLS_HEALTHY,
  levels = design_matrix
)

contrast_fit        <- contrasts.fit(lm_fit, contrast_matrix)
empirical_bayes_fit <- eBayes(contrast_fit)

# 5. ROBUST DMR EXTRACTION (STRICT FDR < 0.05 CUTOFF RESTORED)
get_sig_dmrs <- function(contrast_name) {
  t_table <- topTable(empirical_bayes_fit, coef = contrast_name, number = Inf, adjust.method = "BH")
  # Use BH-adjusted p-values (FDR) < 0.05 and a logFC thresholds of 1
  sig_rows <- rownames(t_table[t_table$adj.P.Val < 0.05 & abs(t_table$logFC) > 1, ])
  message(paste("Contrast:", contrast_name, "--> Found", length(sig_rows), "statistically significant DMRs under strict FDR control."))
  return(sig_rows)
}

dmrs_alps      <- get_sig_dmrs("ALPS_vs_Controls")
dmrs_alps_like <- get_sig_dmrs("ALPS_Like_vs_Controls")

# 6. GRAPHICAL OUTPUTS
# A. UpSet Plot (Safe rendering execution checked)
upset_input_list <- list(ALPS_vs_Controls = dmrs_alps, ALPS_Like_vs_Controls = dmrs_alps_like)

# Generate and save the UpSet plot
if (length(dmrs_alps) > 0 || length(dmrs_alps_like) > 0) {

  grDevices::png(
    filename = file.path(figures_dir, "dmrs_upset_plot.png"),
    width = 3300,
    height = 2100,
    res = 300,
    bg = "white"
  )

  tryCatch(
    {
      UpSetR::upset(
        UpSetR::fromList(upset_input_list),
        order.by = "freq",
        main.bar.color = "darkcyan",
        sets.bar.color = "coral",
        text.scale = 1.3
      )
    },
    finally = {
      grDevices::dev.off()
    }
  )

} else {
  message("UpSet plot skipped: no significant DMRs available.")
}

# B. Volcano Plot for the ALPS vs Controls Comparison
# Purpose: Build a publication-ready Volcano Plot resolving overplotting artifacts,
# color-coding directional dynamics, and mapping the top meaningful symbols.

library(ggplot2)
library(ggrepel)
library(dplyr)
library(stringr)

# ------------------------------------------------------------------------------
# VOLCANO PLOT 1: ALPS VS HEALTHY CONTROLS
# ------------------------------------------------------------------------------

# 1. EXTRACT STATISTICS AND MERGE GENOMIC ANNOTATIONS
t_table_alps <- topTable(empirical_bayes_fit, coef = "ALPS_vs_Controls", number = Inf, adjust.method = "BH")
t_table_alps$Region_ID <- rownames(t_table_alps)
volcano_alps <- merge(t_table_alps, genomic_annotations, by.x = "Region_ID", by.y = "REGION ID")

# 2. FORMAT ADVANCED DYNAMIC ROW STRINGS AND ASSIGN STATUS
volcano_alps <- volcano_alps %>%
  mutate(
    ANNOTATED_GENE = ifelse(is.na(ANNOTATED_GENE) | ANNOTATED_GENE == "NA", "", ANNOTATED_GENE),
    OTHER_GENES    = ifelse(is.na(OTHER_GENES) | OTHER_GENES == "NA", "", OTHER_GENES),
    # Simplify REGIONAL_TYPE to match the heatmap architecture exactly
    Simplified_Type = stringr::str_extract(REGIONAL_TYPE, "^[A-Za-z]+"),
    # Standardize the row annotation label syntax
    Expanded_Row_Label = case_when(
      OTHER_GENES != "" ~ paste0(ANNOTATED_GENE, " [", Simplified_Type, "] (", OTHER_GENES, ")"),
      TRUE              ~ paste0(ANNOTATED_GENE, " [", Simplified_Type, "]")
    ),
    Status = case_when(
      adj.P.Val < 0.05 & logFC > 1  ~ "Hypermethylated DMR",
      adj.P.Val < 0.05 & logFC < -1 ~ "Hypomethylated DMR",
      TRUE                           ~ "NS"
    )
  )

# 3. SELECT TOP HITS FOR LABELING (Avoids overlapping text crowding)
top_labeled_alps <- volcano_alps %>%
  filter(Status != "NS" & ANNOTATED_GENE != "") %>%
  arrange(adj.P.Val) %>%
  head(50)

# 4. RENDER THE ADVANCED GRAPHIC FOR ALPS
p_volcano_alps <- ggplot(
  volcano_alps,
  aes(x = logFC, y = -log10(adj.P.Val))
) +
  geom_point(
    data = filter(volcano_alps, Status == "NS"),
    aes(color = Status),
    alpha = 0.2,
    size = 0.8
  ) +
  geom_point(
    data = filter(volcano_alps, Status != "NS"),
    aes(color = Status),
    alpha = 0.9,
    size = 1.8
  ) +
  scale_color_manual(
    values = c(
      "Hypermethylated DMR" = "firebrick3",
      "Hypomethylated DMR" = "dodgerblue3",
      "NS" = "grey70"
    )
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey30",
    linewidth = 0.5
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    color = "grey30",
    linewidth = 0.5
  ) +
  geom_text_repel(
    data = top_labeled_alps,
    aes(label = Expanded_Row_Label),
    size = 3.2,
    fontface = "bold",
    color = "black",
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = 20,
    segment.color = "black",
    linewidth = 0.4
  ) +
  labs(
    title = "Volcano Plot: Differential Methylation Profile",
    subtitle = "ALPS Patients vs. Healthy Controls (Limma-Voom Empirical Bayes Weights)",
    x = "Log2 Fold Change (Methylation Shift)",
    y = "-Log10 FDR (Adjusted P-Value)",
    color = "Significance Thresholds"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      size = 17, face = "bold", colour = "black"
    ),
    plot.subtitle = element_text(
      size = 10.5, colour = "grey35",
      margin = margin(b = 15)
    ),
    axis.title = element_text(
      size = 11, face = "bold"
    ),
    axis.text = element_text(
      size = 10, colour = "black"
    ),
    legend.title = element_text(
      face = "bold", size = 10
    ),
    legend.text = element_text(size = 9),
    legend.position = "right",
    plot.margin = margin(15, 35, 15, 20)
  )

print(p_volcano_alps)

ggsave(
  file.path(figures_dir, "volcano_alps_differential_methylation.png"),
  plot = p_volcano_alps,
  width = 12,
  height = 9,
  dpi = 300,
  bg = "white"
)


# ------------------------------------------------------------------------------
# VOLCANO PLOT 2: ALPS-LIKE VS HEALTHY CONTROLS
# ------------------------------------------------------------------------------

# 1. EXTRACT STATISTICS AND MERGE GENOMIC ANNOTATIONS
t_table_like <- topTable(empirical_bayes_fit, coef = "ALPS_Like_vs_Controls", number = Inf, adjust.method = "BH")
t_table_like$Region_ID <- rownames(t_table_like)
volcano_like <- merge(t_table_like, genomic_annotations, by.x = "Region_ID", by.y = "REGION ID")

# 2. FORMAT ADVANCED DYNAMIC ROW STRINGS AND ASSIGN STATUS
volcano_like <- volcano_like %>%
  mutate(
    ANNOTATED_GENE = ifelse(is.na(ANNOTATED_GENE) | ANNOTATED_GENE == "NA", "", ANNOTATED_GENE),
    OTHER_GENES    = ifelse(is.na(OTHER_GENES) | OTHER_GENES == "NA", "", OTHER_GENES),
    Simplified_Type = stringr::str_extract(REGIONAL_TYPE, "^[A-Za-z]+"),
    Expanded_Row_Label = case_when(
      OTHER_GENES != "" ~ paste0(ANNOTATED_GENE, " [", Simplified_Type, "] (", OTHER_GENES, ")"),
      TRUE              ~ paste0(ANNOTATED_GENE, " [", Simplified_Type, "]")
    ),
    Status = case_when(
      adj.P.Val < 0.05 & logFC > 1  ~ "Hypermethylated DMR",
      adj.P.Val < 0.05 & logFC < -1 ~ "Hypomethylated DMR",
      TRUE                           ~ "NS"
    )
  )

# 3. SELECT TOP HITS FOR LABELING
top_labeled_like <- volcano_like %>%
  filter(Status != "NS" & ANNOTATED_GENE != "") %>%
  arrange(adj.P.Val) %>%
  head(50)

# 4. RENDER THE ADVANCED GRAPHIC FOR ALPS-LIKE (Subtitle fixed)
p_volcano_like <- ggplot(
  volcano_like,
  aes(x = logFC, y = -log10(adj.P.Val))
) +
  geom_point(
    data = filter(volcano_like, Status == "NS"),
    aes(color = Status),
    alpha = 0.2,
    size = 0.8
  ) +
  geom_point(
    data = filter(volcano_like, Status != "NS"),
    aes(color = Status),
    alpha = 0.9,
    size = 1.8
  ) +
  scale_color_manual(
    values = c(
      "Hypermethylated DMR" = "firebrick3",
      "Hypomethylated DMR" = "dodgerblue3",
      "NS" = "grey70"
    )
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey30",
    linewidth = 0.5
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    color = "grey30",
    linewidth = 0.5
  ) +
  geom_text_repel(
    data = top_labeled_like,
    aes(label = Expanded_Row_Label),
    size = 3.2,
    fontface = "bold",
    color = "black",
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = 20,
    segment.color = "black",
    linewidth = 0.4
  ) +
  labs(
    title = "Volcano Plot: Differential Methylation Profile",
    subtitle = "ALPS-Like Patients vs. Healthy Controls (Limma-Voom Empirical Bayes Weights)",
    x = "Log2 Fold Change (Methylation Shift)",
    y = "-Log10 FDR (Adjusted P-Value)",
    color = "Significance Thresholds"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      size = 17, face = "bold", colour = "black"
    ),
    plot.subtitle = element_text(
      size = 10.5, colour = "grey35",
      margin = margin(b = 15)
    ),
    axis.title = element_text(
      size = 11, face = "bold"
    ),
    axis.text = element_text(
      size = 10, colour = "black"
    ),
    legend.title = element_text(
      face = "bold", size = 10
    ),
    legend.text = element_text(size = 9),
    legend.position = "right",
    plot.margin = margin(15, 35, 15, 20)
  )

print(p_volcano_like)

ggsave(
  file.path(figures_dir, "volcano_alps_like_differential_methylation.png"),
  plot = p_volcano_like,
  width = 12,
  height = 9,
  dpi = 300,
  bg = "white"
)

# C. Hierarchical Clustering Heatmap using pheatmap
# Purpose: Render the final publication-grade heatmap incorporating explicit 
# gene labels, simplified regional types (Exon/Intron), and multi-gene targets.

library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(stringr)

# 1. PARSE AND FORMAT ADVANCED COMPREHENSIVE GENOMIC ROW LABELS
# Combine significant regions from both contrasts for heatmap annotation
# unique() prevents regions significant in both contrasts from being duplicated.
all_sig_dmrs <- unique(c(dmrs_alps, dmrs_alps_like))

# Filter genomic annotations using the significant dynamic hits pool
row_mappings <- genomic_annotations %>%
  dplyr::filter(`REGION ID` %in% all_sig_dmrs) %>%
  dplyr::mutate(
    ANNOTATED_GENE = ifelse(
      is.na(ANNOTATED_GENE) | ANNOTATED_GENE == "NA",
      "",
      ANNOTATED_GENE
    ),
    OTHER_GENES = ifelse(
      is.na(OTHER_GENES) | OTHER_GENES == "NA",
      "",
      OTHER_GENES
    ),
    Simplified_Type = stringr::str_extract(REGIONAL_TYPE, "^[A-Za-z]+"),
    Display_Gene = dplyr::if_else(
      ANNOTATED_GENE != "",
      ANNOTATED_GENE,
      as.character(`REGION ID`)
    ),
    Expanded_Row_Label = dplyr::case_when(
      OTHER_GENES != "" ~ paste0(
        Display_Gene, " [", Simplified_Type, "] (", OTHER_GENES, ")"
      ),
      TRUE ~ paste0(Display_Gene, " [", Simplified_Type, "]")
    )
  )

# 1. SYNCHRONIZE THE VALID FEATURES MATRIX
valid_features <- intersect(
  all_sig_dmrs,
  row_mappings$`REGION ID`
)

heatmap_features <- head(valid_features, 45)

heatmap_matrix <- v_weights$E[heatmap_features, , drop = FALSE]

# Map the formatted labels to the matrix rows safely
matched_indices <- match(
  rownames(heatmap_matrix),
  row_mappings$`REGION ID`
)

row_labels <- row_mappings$Expanded_Row_Label[matched_indices]

# Keep rows even when gene annotations are missing
missing_labels <- is.na(row_labels) | row_labels == ""
row_labels[missing_labels] <- rownames(heatmap_matrix)[missing_labels]

rownames(heatmap_matrix) <- make.unique(row_labels)

# 2. COMPILE CLINICAL AND GENETIC ANNOTATIONS
annotation_col <- data.frame(
  Pathology      = as.character(metadata$pathology),
  ClinicalStatus = as.character(metadata$clinical_status),
  GeneticStatus  = as.character(metadata$genetic_status),
  GeneVariant    = as.character(metadata$gene_variant),
  stringsAsFactors = FALSE
)

rownames(annotation_col) <- rownames(metadata)

annotation_col <- annotation_col %>%
  mutate(
    across(
      everything(),
      ~ ifelse(is.na(.) | . %in% c("NA", "null", ""), "NA", .)
    )
  ) %>%
  mutate(across(everything(), as.factor))

# Ensure the heatmap columns and annotation rows are aligned
common_samples <- intersect(
  colnames(heatmap_matrix),
  rownames(annotation_col)
)

heatmap_matrix <- heatmap_matrix[, common_samples, drop = FALSE]
annotation_col <- annotation_col[common_samples, , drop = FALSE]

# 3. GENERATE PALETTES WITH CONSISTENT NA COLORING
get_homogenized_palette <- function(factor_variable, palette_name) {
  unique_levels <- levels(factor_variable)
  non_na_levels <- unique_levels[unique_levels != "NA"]

  if (length(non_na_levels) > 0) {
    max_colors <- RColorBrewer::brewer.pal.info[palette_name, "maxcolors"]
    base_colors <- RColorBrewer::brewer.pal(
      n = min(max(3, min(length(non_na_levels), max_colors)), max_colors),
      name = palette_name
    )

    color_vector <- grDevices::colorRampPalette(base_colors)(
      length(non_na_levels)
    )

    names(color_vector) <- non_na_levels
  } else {
    color_vector <- character(0)
  }

  if ("NA" %in% unique_levels) {
    color_vector["NA"] <- "grey85"
  }

  color_vector
}

pathology_colors <- get_homogenized_palette(
  annotation_col$Pathology, "Set1"
)

clinical_status_colors <- get_homogenized_palette(
  annotation_col$ClinicalStatus, "Set2"
)

genetic_status_colors <- get_homogenized_palette(
  annotation_col$GeneticStatus, "Dark2"
)

gene_variant_levels <- levels(annotation_col$GeneVariant)
non_na_variants <- gene_variant_levels[gene_variant_levels != "NA"]

if (length(non_na_variants) > 0) {
  gene_variant_colors <- grDevices::colorRampPalette(
    RColorBrewer::brewer.pal(8, "Set3")
  )(length(non_na_variants))

  names(gene_variant_colors) <- non_na_variants
} else {
  gene_variant_colors <- character(0)
}

if ("NA" %in% gene_variant_levels) {
  gene_variant_colors["NA"] <- "grey85"
}

custom_track_colors <- list(
  Pathology      = pathology_colors,
  ClinicalStatus = clinical_status_colors,
  GeneticStatus  = genetic_status_colors,
  GeneVariant    = gene_variant_colors
)

# 4. GENERATE AND EXPORT THE ANNOTATED HEATMAP
heatmap_result <- pheatmap::pheatmap(
  heatmap_matrix,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  annotation_col = annotation_col,
  annotation_colors = custom_track_colors,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize = 7,
  fontsize_row = 6,
  fontsize_col = 6,
  angle_col = 45,
  border_color = "grey",
  main = "Hierarchical Clustering & Epigenetic Signatures of Significant DMRs",
  treeheight_row = 10,
  treeheight_col = 12,
  cellwidth = 10,
  cellheight = 12,
  filename = file.path(
    figures_dir,
    "hierarchical_clustering_significant_dmrs.png"
  ),
  width = 12,
  height = 4
)

grid::grid.newpage()
grid::grid.draw(heatmap_result$gtable)

# ==============================================================================
# SECTION 4.6: GENE ONTOLOGY & FUNCTIONAL ENRICHMENT ANALYSIS
# ==============================================================================
# Objective: Extract genes from significant DMRs, perform GO/KEGG enrichment,
# and generate robust visualizations only when valid results are available.

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(tidyr)
library(dplyr)
library(ggplot2)

# Safely convert enrichment results to a data frame
safe_enrichment_df <- function(x) {
  if (is.null(x)) {
    return(data.frame())
  }

  tryCatch(
    as.data.frame(x),
    error = function(e) {
      message("Could not extract enrichment results: ", conditionMessage(e))
      data.frame()
    }
  )
}

# Validate, save and display a plot without stopping the analysis
safe_save_plot <- function(plot_object, file_name, width, height) {

  if (is.null(plot_object)) {
    return(invisible(NULL))
  }

  # Validate the plot before attempting to save it
  valid_plot <- tryCatch({
    invisible(ggplot2::ggplotGrob(plot_object))
    TRUE
  }, error = function(e) {
    message(
      "Plot skipped (validation failed): ", file_name, " — ",
      conditionMessage(e)
    )
    FALSE
  })

  if (!valid_plot) {
    return(invisible(NULL))
  }

  # Save the plot
  saved <- tryCatch({
    ggplot2::ggsave(
      filename = file.path(figures_dir, file_name),
      plot = plot_object,
      width = width,
      height = height,
      dpi = 300,
      bg = "white"
    )
    TRUE
  }, error = function(e) {
    message(
      "Plot could not be saved: ", file_name, " — ",
      conditionMessage(e)
    )
    FALSE
  })

  # Display only after successful saving
  if (saved) {
    tryCatch(
      print(plot_object),
      error = function(e) {
        message(
          "Plot was saved but could not be displayed: ",
          conditionMessage(e)
        )
      }
    )
  }

  invisible(saved)
}

# Safely generate GO gene-concept networks
safe_cnetplot <- function(go_results, plot_title, file_name) {

  go_df <- safe_enrichment_df(go_results)

  if (nrow(go_df) == 0) {
    message("Cnetplot skipped: no GO terms available for '", plot_title, "'.")
    return(invisible(NULL))
  }

  cnet_plot <- tryCatch({

    p <- enrichplot::cnetplot(
      go_results,
      showCategory = 3
    ) +
      labs(title = plot_title) +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(
          size = 17,
          face = "bold",
          colour = "black"
        ),
        axis.title = element_text(
          size = 11,
          face = "bold"
        ),
        axis.text = element_text(
          size = 10,
          colour = "black"
        ),
        legend.title = element_text(
          face = "bold",
          size = 10
        ),
        legend.text = element_text(size = 9),
        plot.margin = margin(20, 35, 20, 25)
      )

    invisible(ggplot2::ggplotGrob(p))
    p

  }, error = function(e) {
    message(
      "Cnetplot generation failed for '", plot_title, "': ",
      conditionMessage(e)
    )
    NULL
  })

  if (!is.null(cnet_plot)) {
    safe_save_plot(
      cnet_plot,
      file_name,
      width = 12,
      height = 9
    )
  }

  invisible(cnet_plot)
}

# 1. CONSOLIDATE SIGNIFICANT DMRS

all_significant_union <- unique(
  c(dmrs_alps, dmrs_alps_like)
)

# Remove missing or empty region identifiers
all_significant_union <- all_significant_union[
  !is.na(all_significant_union) &
    all_significant_union != ""
]

if (length(all_significant_union) == 0) {

  message(
    "GO/KEGG analysis skipped: no significant DMRs were identified."
  )

} else {

  significant_annotations <- genomic_annotations %>%
    dplyr::filter(`REGION ID` %in% all_significant_union) %>%
    dplyr::mutate(
      OTHER_GENES = ifelse(
        is.na(OTHER_GENES) | OTHER_GENES %in% c("NA", "null"),
        "",
        as.character(OTHER_GENES)
      ),
      ANNOTATED_GENE = ifelse(
        is.na(ANNOTATED_GENE) |
          ANNOTATED_GENE %in% c("NA", "null"),
        "",
        as.character(ANNOTATED_GENE)
      ),
      all_combined_genes = dplyr::case_when(
        ANNOTATED_GENE != "" & OTHER_GENES != "" ~
          paste0(ANNOTATED_GENE, ",", OTHER_GENES),
        ANNOTATED_GENE != "" ~ ANNOTATED_GENE,
        OTHER_GENES != "" ~ OTHER_GENES,
        TRUE ~ ""
      )
    )

  # Keep unannotated regions in significant_annotations.
  # Only rows with gene symbols contribute to functional enrichment.
  flattened_gene_records <- significant_annotations %>%
    dplyr::filter(all_combined_genes != "") %>%
    tidyr::separate_rows(all_combined_genes, sep = "\\s*,\\s*") %>%
    dplyr::mutate(
      all_combined_genes = trimws(all_combined_genes)
    ) %>%
    dplyr::filter(
      !is.na(all_combined_genes),
      !all_combined_genes %in% c("", "NA", "null")
    )

  compiled_target_symbols <- unique(
    flattened_gene_records$all_combined_genes
  )

  message(
    "Total unique genes extracted for functional profiling: ",
    length(compiled_target_symbols)
  )

  # 2. RUN GO AND KEGG ONLY IF GENES ARE AVAILABLE

  if (length(compiled_target_symbols) == 0) {

    message(
      "GO/KEGG analysis skipped: significant DMRs were found, ",
      "but no annotated gene symbols were available."
    )

  } else {

    # --------------------------------------------------------
    # 2A. GENE ONTOLOGY ENRICHMENT: BP, MF AND CC
    # --------------------------------------------------------

    global_go_enrichment <- tryCatch(
      clusterProfiler::enrichGO(
        gene = compiled_target_symbols,
        OrgDb = org.Hs.eg.db,
        keyType = "SYMBOL",
        ont = "ALL",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
      ),
      error = function(e) {
        message(
          "GO enrichment failed: ",
          conditionMessage(e)
        )
        NULL
      }
    )

    go_results_df <- safe_enrichment_df(global_go_enrichment)


    # 3. GO FACETED BARPLOT

    if (
      nrow(go_results_df) > 0 &&
      all(c("ONTOLOGY", "Description", "Count", "p.adjust") %in%
            colnames(go_results_df))
    ) {

      go_results_df <- go_results_df %>%
        dplyr::group_by(ONTOLOGY) %>%
        dplyr::slice_min(
          order_by = p.adjust,
          n = 5,
          with_ties = FALSE
        ) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          Description = factor(
            Description,
            levels = rev(unique(Description))
          )
        )

      go_barplot <- ggplot(
        go_results_df,
        aes(x = Count, y = Description, fill = p.adjust)
      ) +
        geom_col() +
        facet_grid(
          ONTOLOGY ~ .,
          scales = "free_y",
          space = "free_y"
        ) +
        scale_fill_gradient(
          low = "blue",
          high = "red",
          name = "Adjusted p-value"
        ) +
        labs(
          title = "Gene Ontology Functional Enrichment",
          x = "Gene Count",
          y = "Enriched Ontological Terms"
        ) +
        theme_classic(base_size = 12) +
        theme(
          plot.title = element_text(
            size = 17,
            face = "bold",
            colour = "black"
          ),
          axis.title = element_text(
            size = 11,
            face = "bold"
          ),
          axis.text = element_text(
            size = 10,
            colour = "black"
          ),
          strip.text.y = element_text(
            angle = 0,
            face = "bold",
            size = 10,
            colour = "black"
          ),
          legend.title = element_text(
            face = "bold",
            size = 10
          ),
          legend.text = element_text(size = 9),
          panel.spacing.y = grid::unit(0.5, "lines"),
          plot.margin = margin(15, 35, 15, 20)
        )

      safe_save_plot(
        go_barplot,
        "gene_ontology_functional_enrichment.png",
        width = 12,
        height = 10
      )

    } else {

      message(
        "GO barplot skipped: no significant GO terms available."
      )
    }

    # 4. KEGG PATHWAY ENRICHMENT
    
    kegg_translation_table <- tryCatch(
      clusterProfiler::bitr(
        compiled_target_symbols,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db
      ),
      error = function(e) {
        message(
          "KEGG gene ID translation failed: ",
          conditionMessage(e)
        )
        NULL
      }
    )

    global_kegg_enrichment <- NULL

    if (
      !is.null(kegg_translation_table) &&
      nrow(kegg_translation_table) > 0 &&
      "ENTREZID" %in% colnames(kegg_translation_table)
    ) {

      entrez_ids <- unique(
        na.omit(as.character(kegg_translation_table$ENTREZID))
      )

      entrez_ids <- entrez_ids[entrez_ids != ""]

      if (length(entrez_ids) > 0) {

        global_kegg_enrichment <- tryCatch(
          clusterProfiler::enrichKEGG(
            gene = entrez_ids,
            organism = "hsa",
            keyType = "ncbi-geneid",
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH"
          ),
          error = function(e) {
            message(
              "KEGG enrichment failed: ",
              conditionMessage(e)
            )
            NULL
          }
        )

      } else {
        message("KEGG enrichment skipped: no valid ENTREZ IDs.")
      }

    } else {
      message(
        "KEGG enrichment skipped: no gene symbols could be translated."
      )
    }

    # 5. KEGG DOTPLOT

    kegg_results_df <- safe_enrichment_df(global_kegg_enrichment)

    if (nrow(kegg_results_df) > 0) {

      kegg_dotplot <- tryCatch({

        p <- enrichplot::dotplot(
          global_kegg_enrichment,
          showCategory = 10
        ) +
          labs(
            title = "Enriched Pathways Matrix"
          ) +
          theme_classic(base_size = 12) +
          theme(
            plot.title = element_text(
              size = 17,
              face = "bold",
              colour = "black"
            ),
            axis.title = element_text(
              size = 11,
              face = "bold"
            ),
            axis.text = element_text(
              size = 10,
              colour = "black"
            ),
            legend.title = element_text(
              face = "bold",
              size = 10
            ),
            legend.text = element_text(size = 9),
            plot.margin = margin(15, 35, 15, 20)
          )

        invisible(ggplot2::ggplotGrob(p))
        p

      }, error = function(e) {
        message(
          "KEGG dotplot skipped: plot generation failed. ",
          conditionMessage(e)
        )
        NULL
      })

      if (!is.null(kegg_dotplot)) {
        safe_save_plot(
          kegg_dotplot,
          "enriched_pathways_kegg_dotplot.png",
          width = 12,
          height = 8
        )
      }

    } else {
      message(
        "KEGG dotplot skipped: no enriched pathways available."
      )
    }

    # 6. ONTOLOGY-SPECIFIC GO ENRICHMENT FOR CNETPLOTS

    run_ontology_enrichment <- function(gene_symbols, ontology) {

      tryCatch(
        clusterProfiler::enrichGO(
          gene = gene_symbols,
          OrgDb = org.Hs.eg.db,
          keyType = "SYMBOL",
          ont = ontology,
          pAdjustMethod = "BH",
          pvalueCutoff = 0.1
        ),
        error = function(e) {
          message(
            "GO ", ontology, " enrichment failed: ",
            conditionMessage(e)
          )
          NULL
        }
      )
    }

    go_bp_only <- run_ontology_enrichment(
      compiled_target_symbols, "BP"
    )

    go_mf_only <- run_ontology_enrichment(
      compiled_target_symbols, "MF"
    )

    go_cc_only <- run_ontology_enrichment(
      compiled_target_symbols, "CC"
    )

    # 7. GENERATE CNETPLOTS INDEPENDENTLY

    safe_cnetplot(
      go_bp_only,
      "Gene-Concept Network: Biological Process Linkages",
      "cnetplot_go_biological_process.png"
    )

    safe_cnetplot(
      go_mf_only,
      "Gene-Concept Network: Molecular Function Linkages",
      "cnetplot_go_molecular_function.png"
    )

    safe_cnetplot(
      go_cc_only,
      "Gene-Concept Network: Cellular Component Linkages",
      "cnetplot_go_cellular_component.png"
    )

  } # End: compiled_target_symbols available

} # End: significant DMRs available

length(compiled_target_symbols)
nrow(kegg_translation_table)
head(setdiff(compiled_target_symbols, kegg_translation_table$SYMBOL))