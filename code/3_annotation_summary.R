# 0. Load packages

library(readxl)
library(dplyr)
library(ggplot2)
library(readr)
library(forcats)
library(stringr)
library(tidyr)

# 1. Paths and output directories

input_file <- "data/raw_data_annotated_hg19.xlsx"

figures_dir <- "results/figures"
tables_dir  <- "results/tables"

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  tables_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# 2. Read the annotated data and comparison tables

annotated_data <- read_excel(
  input_file,
  sheet = "Annotated_Data"
) %>%
  as.data.frame()

comparison_chipseeker <- read_excel(
  input_file,
  sheet = "Comparison_ChIPseeker"
) %>%
  as.data.frame()

comparison_overlapping <- read_excel(
  input_file,
  sheet = "Comparison_Overlapping"
) %>%
  as.data.frame()

# Check required columns
required_annotated <- c(
  "REGIONAL_TYPE",
  "ANNOTATED_GENE",
  "OTHER_GENES"
)

required_chipseeker <- c(
  "STATUS",
  "ANNOTATED_GENE"
)

required_overlapping <- c(
  "STATUS",
  "GENES_OVERLAPPING"
)

if (!all(required_annotated %in% names(annotated_data))) {
  stop(
    "The Annotated_Data sheet is missing required columns: ",
    paste(
      setdiff(required_annotated, names(annotated_data)),
      collapse = ", "
    )
  )
}

if (!all(required_chipseeker %in% names(comparison_chipseeker))) {
  stop(
    "The Comparison_ChIPseeker sheet is missing required columns: ",
    paste(
      setdiff(required_chipseeker, names(comparison_chipseeker)),
      collapse = ", "
    )
  )
}

if (!all(required_overlapping %in% names(comparison_overlapping))) {
  stop(
    "The Comparison_Overlapping sheet is missing required columns: ",
    paste(
      setdiff(required_overlapping, names(comparison_overlapping)),
      collapse = ", "
    )
  )
}

# 3. Common helper functions

clean_text <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == ""] <- NA_character_
  x
}

# Convert STATUS values into readable labels
status_labels <- c(
  "Exact_match" = "Exact match",
  "Partial_match" = "Partial match",
  "New_annotation" = "New annotation",
  "Discordant_annotations" = "Discordant annotations",
  "Old_annotation_no_new_symbol" =
    "Original annotation, no new symbol",
  "Unannotated_before_and_after" =
    "Unannotated before and after"
)

# Summarize annotation concordance
make_status_summary <- function(data) {

  data %>%
    mutate(
      STATUS = clean_text(STATUS),
      STATUS = ifelse(
        is.na(STATUS),
        "Unknown",
        STATUS
      )
    ) %>%
    count(STATUS, name = "n") %>%
    mutate(
      category = recode(
        STATUS,
        !!!status_labels,
        .default = STATUS
      ),
      percentage = 100 * n / sum(n),
      label = paste0(
        sprintf("%.1f", percentage),
        "% (n=",
        n,
        ")"
      )
    ) %>%
    arrange(desc(percentage)) %>%
    mutate(
      category = factor(
        category,
        levels = rev(category)
      )
    )
}

# Plot annotation concordance
plot_status_summary <- function(
    summary_data,
    plot_title,
    plot_subtitle,
    output_file
) {

  max_percentage <- max(
    summary_data$percentage,
    na.rm = TRUE
  )

  p <- ggplot(
    summary_data,
    aes(
      x = percentage,
      y = category,
      fill = category
    )
  ) +
    geom_col(
      width = 0.68,
      colour = "black",
      linewidth = 0.35
    ) +
    geom_text(
      aes(label = label),
      hjust = -0.08,
      size = 3.5,
      fontface = "bold",
      colour = "black"
    ) +
    scale_x_continuous(
      limits = c(0, max_percentage * 1.35),
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_fill_brewer(
      palette = "Dark2",
      guide = "none"
    ) +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = "Percentage of regions",
      y = NULL
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
      axis.title.x = element_text(
        size = 11,
        face = "bold",
        margin = margin(t = 10)
      ),
      axis.text = element_text(
        size = 10,
        colour = "black"
      ),
      plot.margin = margin(15, 45, 15, 15)
    )

  ggsave(
    filename = file.path(figures_dir, output_file),
    plot = p,
    width = 11,
    height = 7,
    dpi = 300,
    bg = "white"
  )

  p
}

# 4. Genomic region annotation distribution
# This plot uses REGIONAL_TYPE, which represents the primary
# ChIPseeker annotation assigned to each region.

region_annotation_summary <- annotated_data %>%
  mutate(
    REGIONAL_TYPE = clean_text(REGIONAL_TYPE),
    category = case_when(

      is.na(REGIONAL_TYPE) ~ "Not annotated",

      grepl(
        "^Promoter \\(<=1kb\\)",
        REGIONAL_TYPE,
        ignore.case = TRUE
      ) ~ "Promoter (<=1 kb)",

      grepl(
        "^Promoter \\(1-2kb\\)",
        REGIONAL_TYPE,
        ignore.case = TRUE
      ) ~ "Promoter (1-2 kb)",

      grepl(
        "^Promoter",
        REGIONAL_TYPE,
        ignore.case = TRUE
      ) ~ "Promoter (other)",

      grepl(
        "^Distal Intergenic",
        REGIONAL_TYPE,
        ignore.case = TRUE
      ) ~ "Distal intergenic",

      grepl(
        "^Downstream",
        REGIONAL_TYPE,
        ignore.case = TRUE
      ) ~ "Downstream",

      grepl(
        "^5['’] UTR",
        REGIONAL_TYPE,
        ignore.case = TRUE
      ) ~ "5' UTR",

      grepl(
        "^3['’] UTR",
        REGIONAL_TYPE,
        ignore.case = TRUE
      ) ~ "3' UTR",

      grepl(
        "^Exon",
        REGIONAL_TYPE,
        ignore.case = TRUE
      ) ~ "Exon",

      grepl(
        "^Intron",
        REGIONAL_TYPE,
        ignore.case = TRUE
      ) ~ "Intron",

      TRUE ~ "Other"
    )
  ) %>%
  count(category, name = "n") %>%
  mutate(
    percentage = 100 * n / sum(n),
    label = paste0(
      sprintf("%.1f", percentage),
      "%\n(n=",
      n,
      ")"
    )
  ) %>%
  arrange(percentage) %>%
  mutate(
    category = factor(
      category,
      levels = category
    )
  )

print(
  as.data.frame(region_annotation_summary),
  row.names = FALSE
)

p_region_annotation <- ggplot(
  region_annotation_summary,
  aes(
    x = category,
    y = percentage,
    fill = category
  )
) +
  geom_col(
    width = 0.68,
    colour = "black",
    linewidth = 0.35
  ) +
  geom_text(
    aes(label = label),
    hjust = -0.08,
    size = 3.2,
    fontface = "bold"
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(
      0,
      max(region_annotation_summary$percentage) * 1.30
    ),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_brewer(palette = "Dark2") +
  labs(
    title = "Genomic region annotation",
    subtitle = paste(
      "Distribution of regions across broad genomic",
      "categories assigned by ChIPseeker"
    ),
    x = NULL,
    y = "Percentage of regions",
    fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      colour = "black"
    ),
    plot.subtitle = element_text(
      size = 11,
      colour = "grey35",
      margin = margin(b = 15)
    ),
    axis.title.x = element_text(
      size = 11,
      face = "bold",
      margin = margin(t = 10)
    ),
    axis.text = element_text(
      size = 10,
      colour = "black"
    ),
    legend.position = "none",
    plot.margin = margin(15, 40, 15, 15)
  )

ggsave(
  file.path(
    figures_dir,
    "genomic_region_annotation.png"
  ),
  p_region_annotation,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)

write_tsv(
  region_annotation_summary,
  file.path(
    tables_dir,
    "genomic_region_annotation_summary.tsv"
  )
)

# 5. Original genes vs primary ChIPseeker gene
# This comparison evaluates the original annotation against
# ANNOTATED_GENE only. OTHER_GENES is not used here.

annotation_summary_chipseeker <- make_status_summary(
  comparison_chipseeker
)

print(
  as.data.frame(annotation_summary_chipseeker),
  row.names = FALSE
)

p_annotation_chipseeker <- plot_status_summary(
  summary_data = annotation_summary_chipseeker,
  plot_title = "Gene annotation concordance: ChIPseeker",
  plot_subtitle = paste(
    "Comparison between original gene annotations",
    "and the primary ChIPseeker gene"
  ),
  output_file = "gene_annotation_chipseeker.png"
)

write_tsv(
  annotation_summary_chipseeker,
  file.path(
    tables_dir,
    "gene_annotation_chipseeker_summary.tsv"
  )
)

# 6. Original genes vs physically overlapping genes
# This comparison evaluates the original annotation against
# the genes overlapping the genomic region.
annotation_summary_overlapping <- make_status_summary(
  comparison_overlapping
)

print(
  as.data.frame(annotation_summary_overlapping),
  row.names = FALSE
)

p_annotation_overlapping <- plot_status_summary(
  summary_data = annotation_summary_overlapping,
  plot_title = "Gene annotation concordance: overlapping genes",
  plot_subtitle = paste(
    "Comparison between original gene annotations",
    "and genes physically overlapping each region"
  ),
  output_file = "gene_annotation_overlapping.png"
)

write_tsv(
  annotation_summary_overlapping,
  file.path(
    tables_dir,
    "gene_annotation_overlapping_summary.tsv"
  )
)


# 7. Number of additional genes per region
# OTHER_GENES contains the union of original genes and
# overlapping genes, excluding ANNOTATED_GENE.
# This section counts the distinct genes in OTHER_GENES.

split_gene_list <- function(x) {

  x <- clean_text(x)

  if (is.na(x)) {
    return(character(0))
  }

  genes <- unlist(
    strsplit(x, "\\s*[;,|]\\s*")
  )

  genes <- trimws(genes)
  genes <- genes[!is.na(genes) & genes != ""]

  unique(genes)
}

additional_gene_counts <- annotated_data %>%
  mutate(
    N_OTHER_GENES_PLOT = vapply(
      OTHER_GENES,
      function(x) length(split_gene_list(x)),
      integer(1)
    ),
    additional_gene_category = case_when(
      N_OTHER_GENES_PLOT == 0 ~ "0",
      N_OTHER_GENES_PLOT == 1 ~ "1",
      N_OTHER_GENES_PLOT == 2 ~ "2",
      N_OTHER_GENES_PLOT == 3 ~ "3",
      N_OTHER_GENES_PLOT == 4 ~ "4",
      N_OTHER_GENES_PLOT >= 5 ~ "5+"
    )
  )

additional_gene_summary <- additional_gene_counts %>%
  count(
    additional_gene_category,
    name = "n"
  ) %>%
  mutate(
    percentage = 100 * n / sum(n),
    label = paste0(
      sprintf("%.1f", percentage),
      "% (n=",
      n,
      ")"
    ),
    additional_gene_category = factor(
      additional_gene_category,
      levels = c("0", "1", "2", "3", "4", "5+")
    )
  ) %>%
  arrange(additional_gene_category)

print(
  as.data.frame(additional_gene_summary),
  row.names = FALSE
)

p_additional_genes <- ggplot(
  additional_gene_summary,
  aes(
    x = additional_gene_category,
    y = percentage,
    fill = additional_gene_category
  )
) +
  geom_col(
    width = 0.68,
    colour = "black",
    linewidth = 0.35
  ) +
  geom_text(
    aes(label = label),
    vjust = -0.35,
    size = 3.5,
    fontface = "bold"
  ) +
  scale_y_continuous(
    limits = c(
      0,
      max(additional_gene_summary$percentage) * 1.25
    ),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_brewer(
    palette = "Dark2",
    guide = "none"
  ) +
  labs(
    title = "Additional genes per genomic region",
    subtitle = paste(
      "Number of genes in OTHER_GENES after excluding",
      "the primary ChIPseeker gene"
    ),
    x = "Number of additional genes",
    y = "Percentage of regions"
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
    plot.margin = margin(15, 20, 15, 15)
  )

ggsave(
  file.path(
    figures_dir,
    "additional_genes_distribution.png"
  ),
  p_additional_genes,
  width = 9,
  height = 6,
  dpi = 300,
  bg = "white"
)

write_tsv(
  additional_gene_summary,
  file.path(
    tables_dir,
    "additional_genes_distribution_summary.tsv"
  )
)

# 8. Compare the two annotation methods
# Comparison of STATUS distributions between:
#   1. Original genes vs primary ChIPseeker gene
#   2. Original genes vs physically overlapping genes

status_comparison <- bind_rows(
  annotation_summary_chipseeker %>%
    mutate(
      method = "Primary ChIPseeker gene"
    ),

  annotation_summary_overlapping %>%
    mutate(
      method = "Physical gene overlap"
    )
) %>%
  mutate(
    category = as.character(category),
    method = factor(
      method,
      levels = c(
        "Primary ChIPseeker gene",
        "Physical gene overlap"
      )
    )
  )

p_status_comparison <- ggplot(
  status_comparison,
  aes(
    x = percentage,
    y = fct_rev(factor(category)),
    fill = method
  )
) +
  geom_col(
    position = position_dodge(width = 0.78),
    width = 0.68,
    colour = "black",
    linewidth = 0.25
  ) +
  scale_x_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Comparison of gene annotation methods",
    subtitle = paste(
      "Concordance of original annotations with the",
      "primary ChIPseeker gene or physically overlapping genes"
    ),
    x = "Percentage of regions",
    y = NULL,
    fill = "Annotation method"
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
    axis.title.x = element_text(
      size = 11,
      face = "bold"
    ),
    axis.text = element_text(
      size = 9.5,
      colour = "black"
    ),
    legend.position = "bottom",
    plot.margin = margin(15, 20, 15, 15)
  )

ggsave(
  file.path(
    figures_dir,
    "gene_annotation_methods_comparison.png"
  ),
  p_status_comparison,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)

write_tsv(
  status_comparison,
  file.path(
    tables_dir,
    "gene_annotation_methods_comparison.tsv"
  )
)

# 9. Print all figures

print(p_region_annotation)
print(p_annotation_chipseeker)
print(p_annotation_overlapping)
print(p_additional_genes)
print(p_status_comparison)

# 10. Final checks and summary

cat("\n============================================\n")
cat("ANNOTATION VISUALIZATION SUMMARY\n")
cat("============================================\n")

cat(
  "Total annotated regions:",
  nrow(annotated_data),
  "\n"
)

cat(
  "Regions with a primary ChIPseeker gene:",
  sum(!is.na(clean_text(annotated_data$ANNOTATED_GENE))),
  "\n"
)

cat(
  "Regions with additional genes:",
  sum(additional_gene_counts$N_OTHER_GENES_PLOT > 0),
  "\n"
)

cat(
  "Regions without additional genes:",
  sum(additional_gene_counts$N_OTHER_GENES_PLOT == 0),
  "\n"
)

cat(
  "\nFigures saved to:",
  figures_dir,
  "\n"
)

cat(
  "Tables saved to:",
  tables_dir,
  "\n"
)

message(
  "All annotation figures and summary tables have been generated."
)