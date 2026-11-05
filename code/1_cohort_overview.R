# 0. Load required libraries
library(readr)   # Provides fast and friendly tools to read flat and tabular data (like CSV, TSV) into R
library(dplyr)   # Grammar of data manipulation providing a consistent set of verbs to solve data hurdles
library(tidyr)   # Tools to reshape and tidy data layout to achieve consistent and rectangular structures
library(ggplot2) # System for declaratively creating elegant and publication-ready data visualizations

# 1. Read metadata

metadata <- read_tsv(
  "metadata/patients_metadata.txt",
  na = c("NA"),
  show_col_types = FALSE
)

# 2. Variables represented within pathology

variables_to_plot <- c(
  "clinical_status",
  "gene_variant",
  "relationship",
  "genetic_status"
)

# 3. Figure titles

variable_labels <- c(
  clinical_status = "Clinical status",
  gene_variant = "Gene variant",
  relationship = "Relationship",
  genetic_status = "Genetic status"
)

# 4. Prepare data

pathology_characterization <- metadata %>%
  dplyr::select(
    pathology,
    all_of(variables_to_plot)
  ) %>%
  mutate(
    pathology = ifelse(
      is.na(pathology),
      "NA",
      as.character(pathology)
    )
  ) %>%
  pivot_longer(
    cols = all_of(variables_to_plot),
    names_to = "variable",
    values_to = "category"
  ) %>%
  mutate(
    category = ifelse(
      is.na(category),
      "NA (Not Available)",
      as.character(category)
    )
  ) %>%
  count(
    pathology,
    variable,
    category,
    name = "n"
  ) %>%
  group_by(
    pathology,
    variable
  ) %>%
  mutate(
    percentage = 100 * n / sum(n)
  ) %>%
  ungroup()

# 5. Generate one PNG figure per variable
for (var in variables_to_plot) {
  
  plot_data <- pathology_characterization %>%
    filter(variable == var)
  
  # Order pathology groups as ALPS, ALPS-LIKE and HEALTHY CONTROLS
  pathology_order <- c(
    "ALPS",
    "ALPS-LIKE",
    "HEALTHY CONTROLS"
  )
  
  plot_data <- plot_data %>%
    mutate(
      pathology = factor(
        pathology,
        levels = pathology_order
      )
    )
  
  # Labels
  plot_data <- plot_data %>%
    mutate(
      label = ifelse(
        percentage > 0,
        paste0(round(percentage, 1), "%"),
        ""
      )
    )
  
  # Plot
  p <- ggplot(
    plot_data,
    aes(x = pathology, y = percentage, fill = category)
  ) +
    geom_col(
      width = 0.68,
      colour = "black",
      linewidth = 0.35
    ) +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      size = 3.5,
      fontface = "bold",
      colour = "black" 
    ) +
    scale_y_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, 20),
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0.01))
    ) +
    scale_fill_brewer(palette = "Dark2") + 
    labs(
      title = variable_labels[var],
      subtitle = "Distribution within each pathology group",
      x = "Pathology",
      y = "Percentage of patients",
      fill = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(size = 18, face = "bold", colour = "black", hjust = 0),
      plot.subtitle = element_text(size = 11, colour = "grey35", hjust = 0, margin = margin(b = 15)),
      axis.title.x = element_text(size = 11, face = "bold", colour = "black", margin = margin(t = 10)),
      axis.title.y = element_text(size = 11, face = "bold", colour = "black", margin = margin(r = 10)),
      axis.text.x = element_text(size = 10, colour = "black"),
      axis.text.y = element_text(size = 10, colour = "black"),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 10, colour = "black"),
      legend.key.height = unit(0.45, "cm"),
      legend.key.width = unit(0.45, "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      plot.margin = margin(15, 20, 15, 15)
    )
  
  output_png <- file.path("results/figures", paste0(var, ".png"))
  
  ggsave(
    filename = output_png,
    plot = p,
    width = 10,
    height = 7,
    dpi = 300
  )
  
  message("Chart created: ", output_png)
}

# 7. Final confirmation
message("All cohort overview have been successfully saved in results/figures/")

cohort_charts <- file.path("results/figures", paste0(variables_to_plot, ".png"))
knitr::include_graphics(cohort_charts)