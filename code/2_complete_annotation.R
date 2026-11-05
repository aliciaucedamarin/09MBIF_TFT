# 0. Load required libraries
library(readxl)
library(openxlsx)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)
library(GenomicFeatures)
library(AnnotationDbi)
library(IRanges)
library(S4Vectors)
library(tidyverse)

# 1. Load original data and validate REGION IDs
# Loads the Excel and BED files, validates required columns, and checks region IDs.

input_file <- "data/original/all_24_content_clean_annot_nodup.bed.xlsx"
input_sheet <- "all_24_content_clean_annot_nodu"
bed_path <- "data/original/Panel-Regions.bed"

raw_data_input <- readxl::read_excel(input_file, sheet = input_sheet)

bed <- read.table(
  bed_path, header = FALSE, sep = "\t", skip = 2,
  stringsAsFactors = FALSE
)
colnames(bed)[1:4] <- c("#CHROM", "START", "END", "REGION ID")

required_columns <- c("GEN", "REGION ID", "#CHROM", "START", "END")
stopifnot(all(required_columns %in% names(raw_data_input)))

raw_data_input <- raw_data_input %>%
  dplyr::mutate(
    dplyr::across(where(is.character), ~ {
      x <- trimws(.x)
      x[!is.na(x) & x == ""] <- NA_character_
      x
    }),
    GEN = as.character(GEN),
    `REGION ID` = as.character(`REGION ID`)
  )

norm_columns <- grep("_Norm$", names(raw_data_input), value = TRUE)
meth_columns <- grep("_Met$", names(raw_data_input), value = TRUE)
stopifnot(length(norm_columns) == 24)
signal_columns <- unique(c(meth_columns, norm_columns))

cat("\n===== ORIGINAL DATA =====\n")
cat("Unique REGION IDs in BED file:", dplyr::n_distinct(bed$`REGION ID`), "\n")
cat("Unique REGION IDs in raw data:", dplyr::n_distinct(raw_data_input$`REGION ID`), "\n")
cat("Normalized Methylation Reads columns:", length(norm_columns), "\n")
cat("Methylation Reads columns:", length(meth_columns), "\n")

regions_bed <- unique(bed$`REGION ID`)
regions_raw <- unique(raw_data_input$`REGION ID`)
only_in_bed <- setdiff(regions_bed, regions_raw)
only_in_raw <- setdiff(regions_raw, regions_bed)

cat("Regions only in BED:", length(only_in_bed), "\n")
if (length(only_in_bed) > 0) print(only_in_bed)
cat("Regions only in raw_data_input:", length(only_in_raw), "\n")
if (length(only_in_raw) > 0) print(only_in_raw)
same_regions <- setequal(regions_bed, regions_raw)
cat("Same REGION IDs in both files:", same_regions, "\n")

# 2. Assess duplicated REGION IDs before removing them
# Checks whether duplicate IDs have different signal profiles.

valid_region_id <- !is.na(raw_data_input$`REGION ID`) &
  trimws(raw_data_input$`REGION ID`) != ""

duplicate_id_assessment <- raw_data_input %>%
  dplyr::filter(valid_region_id) %>%
  dplyr::group_by(`REGION ID`) %>%
  dplyr::summarise(
    N_ROWS = dplyr::n(),
    N_DISTINCT_PROFILES = nrow(unique(as.data.frame(
      dplyr::pick(dplyr::all_of(signal_columns))
    ))),
    DIFFERENT_SIGNAL_PROFILES = N_DISTINCT_PROFILES > 1,
    .groups = "drop"
  ) %>%
  dplyr::filter(N_ROWS > 1)

cat("\n===== DUPLICATED REGION IDs =====\n")
cat("Duplicated REGION IDs:", nrow(duplicate_id_assessment), "\n")
cat("Rows associated with duplicated IDs:", sum(duplicate_id_assessment$N_ROWS), "\n")
cat(
  "Duplicated IDs with different signal profiles:",
  sum(duplicate_id_assessment$DIFFERENT_SIGNAL_PROFILES), "\n"
)
if (nrow(duplicate_id_assessment) > 0) str(duplicate_id_assessment)

# 3. Remove duplicated REGION IDs and retain the first occurrence
# Keeps the first row for each valid ID; rows with missing IDs are not collapsed.

duplicate_mask <- rep(FALSE, nrow(raw_data_input))
duplicate_mask[valid_region_id] <- duplicated(
  raw_data_input$`REGION ID`[valid_region_id]
)
duplicate_rows_removed <- raw_data_input[duplicate_mask, , drop = FALSE]
raw_data <- raw_data_input[!duplicate_mask, , drop = FALSE]

cat("\n===== REGION ID DEDUPLICATION =====\n")
cat("Rows before removing duplicates:", nrow(raw_data_input), "\n")
cat("Rows removed:", nrow(duplicate_rows_removed), "\n")
cat("Rows remaining:", nrow(raw_data), "\n")

cat(
  "Unique non-empty REGION IDs:",
  dplyr::n_distinct(raw_data$`REGION ID`[
    !is.na(raw_data$`REGION ID`) & trimws(raw_data$`REGION ID`) != ""
  ]),
  "\n"
)

if (
  nrow(duplicate_id_assessment) > 0 &&
  any(duplicate_id_assessment$DIFFERENT_SIGNAL_PROFILES)
) {
  warning(
    "Some duplicated REGION IDs have different signal profiles. ",
    "The first row was retained for each ID. ",
    "Review Duplicate_REGION_IDs_Removed before using the final matrix."
  )
}
if (nrow(raw_data) != 3481) {
  warning(
    "The dataset does not contain the expected rows after deduplication. ",
    "Review the duplicate assessment."
  )
}
raw_data_original <- raw_data

# 4. Clean gene lists and preserve original annotations
# Removes invalid and duplicated symbols, returning NA when no valid symbol remains.

clean_gene_list <- function(x) {
  if (length(x) == 0 || is.na(x) || !nzchar(trimws(as.character(x)))) {
    return(NA_character_)
  }
  genes <- trimws(unlist(strsplit(as.character(x), ",", fixed = TRUE)))
  genes <- genes[!is.na(genes) & genes != "" & genes != "." & genes != "-1"]
  genes <- unique(genes)
  if (length(genes) == 0) return(NA_character_)
  paste(genes, collapse = ", ")
}

has_annotation <- function(x) {
  if (length(x) == 0 || is.na(x) || !nzchar(trimws(as.character(x)))) {
    return(FALSE)
  }
  genes <- trimws(unlist(strsplit(as.character(x), ",", fixed = TRUE)))
  any(!is.na(genes) & genes != "" & genes != "." & genes != "-1")
}

raw_data$GEN_ORIGINAL <- as.character(raw_data$GEN)
raw_data$GEN <- vapply(raw_data$GEN, clean_gene_list, character(1))

cat("\n===== ORIGINAL GENE ANNOTATION =====\n")
cat("Rows with an original gene:", sum(vapply(raw_data$GEN, has_annotation, logical(1))), "\n")
cat("Rows without an original gene:", sum(!vapply(raw_data$GEN, has_annotation, logical(1))), "\n")

# 5. Report regions without a gene before annotation
# Summarizes regions that lack valid original gene symbols.

report_unannotated <- function(data, stage) {
  gene_values <- as.character(data$GEN)
  missing_gene <- is.na(gene_values) | trimws(gene_values) %in% c("", ".", "-1")
  missing_gene[is.na(missing_gene)] <- TRUE
  result <- data[missing_gene, , drop = FALSE]

  cat("\n===== UNANNOTATED REGIONS:", stage, "=====\n")
  cat("Total rows checked:", nrow(data), "\n")
  cat("Regions without an assigned gene:", nrow(result), "\n")
  cat("GEN values that are NA:", sum(is.na(result$GEN)), "\n")

  columns_to_show <- intersect(
    c("REGION ID", "#CHROM", "START", "END", "GEN_ORIGINAL", "GEN",
      "gene_start", "gene_end", "REGIONAL_TYPE", "DISTANCE_TO_TSS"),
    names(result)
  )
  report <- result %>% dplyr::select(dplyr::all_of(columns_to_show))
  if (nrow(report) > 0) str(report) else cat("No unannotated regions found.\n")
  report
}
unannotated_before <- report_unannotated(raw_data, "BEFORE ANNOTATION")

# 6. Check genomic coordinates before annotation
# Ensures coordinates are present, integer-valued, and ordered correctly.

stopifnot(
  all(!is.na(raw_data$`#CHROM`)),
  all(!is.na(raw_data$START)),
  all(!is.na(raw_data$END)),
  all(as.numeric(raw_data$START) == floor(as.numeric(raw_data$START))),
  all(as.numeric(raw_data$END) == floor(as.numeric(raw_data$END))),
  all(as.numeric(raw_data$START) <= as.numeric(raw_data$END))
)

# 7. Annotate regions using hg19 and identify physically overlapping genes
# ChIPseeker assigns the primary gene; findOverlaps identifies genes intersecting each region.

txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene
all_regions_gr <- GenomicRanges::GRanges(
  seqnames = as.character(raw_data$`#CHROM`),
  ranges = IRanges::IRanges(
    start = as.integer(raw_data$START),
    end = as.integer(raw_data$END)
  )
)

annotated_peaks_all <- ChIPseeker::annotatePeak(
  all_regions_gr, tssRegion = c(-2000, 2000), TxDb = txdb_hg19,
  level = "gene", annoDb = "org.Hs.eg.db", verbose = TRUE
)
annotation_temp_df <- as.data.frame(annotated_peaks_all)

# Verify that annotation preserves the original region count, coordinates, and order.
stopifnot(
  nrow(annotation_temp_df) == nrow(raw_data),
  identical(as.character(annotation_temp_df$seqnames), as.character(raw_data$`#CHROM`)),
  identical(as.integer(annotation_temp_df$start), as.integer(raw_data$START)),
  identical(as.integer(annotation_temp_df$end), as.integer(raw_data$END))
)

genes_hg19 <- GenomicFeatures::genes(txdb_hg19)
GenomeInfoDb::seqlevelsStyle(genes_hg19) <-
  GenomeInfoDb::seqlevelsStyle(all_regions_gr)[1]

overlaps <- GenomicRanges::findOverlaps(
  query = all_regions_gr, subject = genes_hg19, ignore.strand = TRUE
)

gene_id_map <- AnnotationDbi::select(
  org.Hs.eg.db, keys = unique(names(genes_hg19)),
  keytype = "ENTREZID", columns = "SYMBOL"
) %>%
  dplyr::mutate(
    ENTREZID = as.character(ENTREZID),
    SYMBOL = trimws(as.character(SYMBOL))
  ) %>%
  dplyr::filter(
    !is.na(SYMBOL), SYMBOL != "", SYMBOL != ".", SYMBOL != "-1"
  ) %>%
  dplyr::distinct(ENTREZID, SYMBOL)

overlap_gene_df <- data.frame(
  REGION_ROW = S4Vectors::queryHits(overlaps),
  GENE_ID = names(genes_hg19)[S4Vectors::subjectHits(overlaps)],
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(GENE_ID = as.character(GENE_ID)) %>%
  dplyr::left_join(gene_id_map, by = c("GENE_ID" = "ENTREZID")) %>%
  dplyr::mutate(SYMBOL = trimws(as.character(SYMBOL))) %>%
  dplyr::filter(
    !is.na(SYMBOL), SYMBOL != "", SYMBOL != ".", SYMBOL != "-1"
  ) %>%
  dplyr::distinct(REGION_ROW, SYMBOL)

overlap_gene_list <- overlap_gene_df %>%
  dplyr::group_by(REGION_ROW) %>%
  dplyr::summarise(
    GENES_OVERLAPPING = paste(sort(unique(SYMBOL)), collapse = ", "),
    N_OVERLAPPING_GENES = dplyr::n_distinct(SYMBOL),
    .groups = "drop"
  )

raw_data$GENES_OVERLAPPING <- NA_character_
raw_data$N_OVERLAPPING_GENES <- NA_integer_
matched_rows <- overlap_gene_list$REGION_ROW
raw_data$GENES_OVERLAPPING[matched_rows] <- overlap_gene_list$GENES_OVERLAPPING
raw_data$N_OVERLAPPING_GENES[matched_rows] <- overlap_gene_list$N_OVERLAPPING_GENES

invalid_overlap <- is.na(raw_data$GENES_OVERLAPPING) |
  trimws(raw_data$GENES_OVERLAPPING) %in% c("", ".", "-1")
raw_data$GENES_OVERLAPPING[invalid_overlap] <- NA_character_
raw_data$N_OVERLAPPING_GENES[invalid_overlap] <- NA_integer_

cat("\n===== OVERLAPPING GENE SYMBOLS =====\n")
cat("Total regions:", nrow(raw_data), "\n")
cat("Regions with overlapping genes:", sum(!is.na(raw_data$GENES_OVERLAPPING)), "\n")
cat("Regions without overlapping genes:", sum(is.na(raw_data$GENES_OVERLAPPING)), "\n")
cat(
  "Maximum overlapping genes per region:",
  if (all(is.na(raw_data$N_OVERLAPPING_GENES))) NA_integer_
  else max(raw_data$N_OVERLAPPING_GENES, na.rm = TRUE),
  "\n"
)
str(raw_data[c("GENES_OVERLAPPING", "N_OVERLAPPING_GENES")])

# 8. Add the primary ChIPseeker gene annotation
# ANNOTATED_GENE is the sole primary gene column linked to REGIONAL_TYPE and DISTANCE_TO_TSS.

annotation_temp_df <- as.data.frame(annotated_peaks_all)
if ("SYMBOL" %in% names(annotation_temp_df)) {
  annotated_gene <- as.character(annotation_temp_df$SYMBOL)
} else if ("geneId" %in% names(annotation_temp_df)) {
  gene_ids <- as.character(annotation_temp_df$geneId)
  gene_id_map_chipseeker <- AnnotationDbi::select(
    org.Hs.eg.db, keys = unique(gene_ids[!is.na(gene_ids)]),
    keytype = "ENTREZID", columns = "SYMBOL"
  ) %>%
    dplyr::mutate(
      ENTREZID = as.character(ENTREZID),
      SYMBOL = trimws(as.character(SYMBOL))
    ) %>%
    dplyr::filter(
      !is.na(SYMBOL), SYMBOL != "", SYMBOL != ".", SYMBOL != "-1"
    ) %>%
    dplyr::distinct(ENTREZID, .keep_all = TRUE)
  annotated_gene <- gene_id_map_chipseeker$SYMBOL[
    match(gene_ids, gene_id_map_chipseeker$ENTREZID)
  ]
} else {
  stop("No SYMBOL or geneId column was found in the ChIPseeker output.")
}

raw_data$ANNOTATED_GENE <- vapply(annotated_gene, clean_gene_list, character(1))
raw_data$REGIONAL_TYPE <- trimws(as.character(annotation_temp_df$annotation))
raw_data$REGIONAL_TYPE[
  is.na(raw_data$REGIONAL_TYPE) | raw_data$REGIONAL_TYPE == ""
] <- NA_character_

raw_data$DISTANCE_TO_TSS <- suppressWarnings(
  as.numeric(as.character(annotation_temp_df$distanceToTSS))
)
raw_data$DISTANCE_TO_TSS[!is.finite(raw_data$DISTANCE_TO_TSS)] <- NA_real_

stopifnot(
  nrow(raw_data) == nrow(annotation_temp_df),
  length(raw_data$ANNOTATED_GENE) == nrow(raw_data)
)

cat("\n===== CHIPSEEKER ANNOTATION =====\n")
cat(
  "Regions with a ChIPseeker gene symbol:",
  sum(vapply(raw_data$ANNOTATED_GENE, has_annotation, logical(1))), "\n"
)
cat(
  "Regions without a ChIPseeker gene symbol:",
  sum(!vapply(raw_data$ANNOTATED_GENE, has_annotation, logical(1))), "\n"
)

# 9. Identify duplicated genomic coordinates
# Reports repeated chromosome/start/end combinations without collapsing regions.

duplicate_coordinates <- raw_data %>%
  dplyr::count(`#CHROM`, START, END, name = "N_ROWS") %>%
  dplyr::filter(N_ROWS > 1)
cat("\n===== DUPLICATED COORDINATES =====\n")
cat("Coordinate combinations occurring more than once:", nrow(duplicate_coordinates), "\n")

# 10. Define gene-set and comparison functions
# Converts comma-separated symbols to unique sets and compares two annotations.

get_gene_set <- function(x) {
  if (length(x) == 0 || is.na(x) || !nzchar(trimws(as.character(x)))) {
    return(character(0))
  }
  genes <- trimws(unlist(strsplit(as.character(x), ",", fixed = TRUE)))
  genes <- genes[!is.na(genes) & genes != "" & genes != "." & genes != "-1"]
  unique(genes)
}

compare_gene_sets <- function(old_genes, new_genes) {
  old_set <- get_gene_set(old_genes)
  new_set <- get_gene_set(new_genes)
  shared_genes <- intersect(old_set, new_set)
  missing_old_genes <- setdiff(old_set, new_set)
  additional_genes <- setdiff(new_set, old_set)
  all_genes <- union(old_set, new_set)

  status <- if (length(old_set) == 0 && length(new_set) == 0) {
    "Unannotated_before_and_after"
  } else if (length(old_set) == 0) {
    "New_annotation"
  } else if (length(new_set) == 0) {
    "Old_annotation_no_new_symbol"
  } else if (setequal(old_set, new_set)) {
    "Exact_match"
  } else if (length(shared_genes) > 0) {
    "Partial_match"
  } else {
    "Discordant_annotations"
  }

  data.frame(
    STATUS = status,
    N_OLD_GENES = length(old_set),
    N_NEW_GENES = length(new_set),
    N_SHARED_GENES = length(shared_genes),
    SHARED_GENES = if (length(shared_genes) > 0) paste(sort(shared_genes), collapse = ", ") else NA_character_,
    MISSING_OLD_GENES = if (length(missing_old_genes) > 0) paste(sort(missing_old_genes), collapse = ", ") else NA_character_,
    ADDITIONAL_GENES = if (length(additional_genes) > 0) paste(sort(additional_genes), collapse = ", ") else NA_character_,
    ORIGINAL_GENE_RETENTION_PERCENT = if (length(old_set) > 0) 100 * length(shared_genes) / length(old_set) else NA_real_,
    JACCARD_PERCENT = if (length(all_genes) > 0) 100 * length(shared_genes) / length(all_genes) else NA_real_,
    stringsAsFactors = FALSE
  )
}

# 11. Build OTHER_GENES by combining original and overlapping genes
# ANNOTATED_GENE remains the primary ChIPseeker gene.
# OTHER_GENES combines original and physically overlapping gene symbols.
# Any gene already present in ANNOTATED_GENE is excluded from OTHER_GENES.

raw_data$OTHER_GENES <- vapply(
  seq_len(nrow(raw_data)),
  function(i) {

    # Primary gene assigned by ChIPseeker
    primary_genes <- get_gene_set(
      raw_data$ANNOTATED_GENE[i]
    )

    # Original gene annotations
    original_genes <- get_gene_set(
      raw_data$GEN_ORIGINAL[i]
    )

    # Genes physically overlapping the genomic region
    overlapping_genes <- get_gene_set(
      raw_data$GENES_OVERLAPPING[i]
    )

    # Combine original and overlapping genes, removing duplicates
    combined_genes <- unique(
      c(original_genes, overlapping_genes)
    )

    # Exclude genes already represented in the primary annotation
    other_genes <- setdiff(
      combined_genes,
      primary_genes
    )

    # Return NA if no additional genes remain
    if (length(other_genes) == 0) {
      return(NA_character_)
    }

    # Store additional genes in alphabetical order
    paste(
      sort(other_genes),
      collapse = ", "
    )
  },
  character(1)
)

# Count additional genes independently of the overlap count
raw_data$N_OTHER_GENES <- vapply(
  raw_data$OTHER_GENES,
  function(x) {
    length(get_gene_set(x))
  },
  integer(1)
)

raw_data$N_OTHER_GENES[
  is.na(raw_data$OTHER_GENES)
] <- 0L

cat("\n===== FINAL GENE COLUMNS =====\n")

cat(
  "Regions with a primary ChIPseeker gene:",
  sum(vapply(
    raw_data$ANNOTATED_GENE,
    has_annotation,
    logical(1)
  )),
  "\n"
)

cat(
  "Regions with additional genes:",
  sum(!is.na(raw_data$OTHER_GENES)),
  "\n"
)

cat(
  "Regions without additional genes:",
  sum(is.na(raw_data$OTHER_GENES)),
  "\n"
)

cat(
  "Maximum number of additional genes:",
  max(raw_data$N_OTHER_GENES, na.rm = TRUE),
  "\n"
)

# Verify that primary genes do not appear in OTHER_GENES
for (i in seq_len(nrow(raw_data))) {

  primary_genes <- get_gene_set(
    raw_data$ANNOTATED_GENE[i]
  )

  other_genes <- get_gene_set(
    raw_data$OTHER_GENES[i]
  )

  stopifnot(
    length(intersect(primary_genes, other_genes)) == 0
  )
}

# 12. Compare original annotations with the ChIPseeker primary annotation
# Compares original GEN symbols against ANNOTATED_GENE for each region.

comparison_chipseeker_metrics <- dplyr::bind_rows(lapply(seq_len(nrow(raw_data)), function(i) {
  compare_gene_sets(raw_data$GEN[i], raw_data$ANNOTATED_GENE[i])
}))

annotation_comparison_chipseeker <- data.frame(
  ROW_ID = seq_len(nrow(raw_data)),
  `REGION ID` = raw_data$`REGION ID`,
  `#CHROM` = raw_data$`#CHROM`,
  START = raw_data$START,
  END = raw_data$END,
  OLD_GENES = raw_data$GEN,
  ANNOTATED_GENE = raw_data$ANNOTATED_GENE,
  REGIONAL_TYPE = raw_data$REGIONAL_TYPE,
  DISTANCE_TO_TSS = raw_data$DISTANCE_TO_TSS,
  comparison_chipseeker_metrics,
  check.names = FALSE
)

# 13. Compare original annotations with physically overlapping genes
# Calculates concordance between GEN and GENES_OVERLAPPING.

comparison_overlap_metrics <- dplyr::bind_rows(lapply(seq_len(nrow(raw_data)), function(i) {
  compare_gene_sets(raw_data$GEN[i], raw_data$GENES_OVERLAPPING[i])
}))

annotation_comparison_overlapping <- data.frame(
  ROW_ID = seq_len(nrow(raw_data)),
  `REGION ID` = raw_data$`REGION ID`,
  `#CHROM` = raw_data$`#CHROM`,
  START = raw_data$START,
  END = raw_data$END,
  OLD_GENES = raw_data$GEN,
  GENES_OVERLAPPING = raw_data$GENES_OVERLAPPING,
  N_OVERLAPPING_GENES = raw_data$N_OVERLAPPING_GENES,
  comparison_overlap_metrics,
  check.names = FALSE
)

# 14. Summarize both annotation comparisons
# Counts concordance categories and calculates gene-set similarity statistics.

annotation_summary_chipseeker <- annotation_comparison_chipseeker %>%
  dplyr::count(STATUS, name = "N_REGIONS") %>%
  dplyr::mutate(PERCENT = 100 * N_REGIONS / nrow(raw_data)) %>%
  dplyr::arrange(dplyr::desc(N_REGIONS))

annotation_summary_overlapping <- annotation_comparison_overlapping %>%
  dplyr::count(STATUS, name = "N_REGIONS") %>%
  dplyr::mutate(PERCENT = 100 * N_REGIONS / nrow(raw_data)) %>%
  dplyr::arrange(dplyr::desc(N_REGIONS))

mean_or_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
median_or_na <- function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)

similarity_summary <- data.frame(
  METHOD = c("ChIPseeker", "Physical_gene_overlap"),
  MEAN_GENE_RETENTION_PERCENT = c(
    mean_or_na(annotation_comparison_chipseeker$ORIGINAL_GENE_RETENTION_PERCENT),
    mean_or_na(annotation_comparison_overlapping$ORIGINAL_GENE_RETENTION_PERCENT)
  ),
  MEDIAN_GENE_RETENTION_PERCENT = c(
    median_or_na(annotation_comparison_chipseeker$ORIGINAL_GENE_RETENTION_PERCENT),
    median_or_na(annotation_comparison_overlapping$ORIGINAL_GENE_RETENTION_PERCENT)
  ),
  MEAN_JACCARD_PERCENT = c(
    mean_or_na(annotation_comparison_chipseeker$JACCARD_PERCENT),
    mean_or_na(annotation_comparison_overlapping$JACCARD_PERCENT)
  ),
  MEDIAN_JACCARD_PERCENT = c(
    median_or_na(annotation_comparison_chipseeker$JACCARD_PERCENT),
    median_or_na(annotation_comparison_overlapping$JACCARD_PERCENT)
  )
)

cat("\n===== CHIPSEEKER COMPARISON SUMMARY =====\n")
print(annotation_summary_chipseeker)
cat("\n===== OVERLAPPING GENES COMPARISON SUMMARY =====\n")
print(annotation_summary_overlapping)
cat("\n===== GENE-SET SIMILARITY SUMMARY =====\n")
print(similarity_summary)

# 15. Identify regions that remain unannotated after both methods
# Reports regions without original, ChIPseeker, and physical-overlap gene symbols.

original_missing <- !vapply(raw_data$GEN, has_annotation, logical(1))
chipseeker_missing <- !vapply(raw_data$ANNOTATED_GENE, has_annotation, logical(1))
overlap_missing <- !vapply(raw_data$GENES_OVERLAPPING, has_annotation, logical(1))

unannotated_before <- raw_data[original_missing, , drop = FALSE]
unannotated_after <- raw_data[
  original_missing & chipseeker_missing & overlap_missing, , drop = FALSE
]

cat("\n===== ANNOTATION COVERAGE =====\n")
cat("Regions without original annotation:", nrow(unannotated_before), "\n")
cat("Regions without annotation from either new method:", nrow(unannotated_after), "\n")

# 16. Verify row counts, original columns, and gene-list construction
# Checks that annotation did not alter region order or coordinates and that cleaning preserves valid symbols.

stopifnot(
  nrow(raw_data) == nrow(raw_data_original),
  nrow(annotation_comparison_chipseeker) == nrow(raw_data),
  nrow(annotation_comparison_overlapping) == nrow(raw_data),
  identical(raw_data$`REGION ID`, raw_data_original$`REGION ID`),
  identical(raw_data$`#CHROM`, raw_data_original$`#CHROM`),
  identical(raw_data$START, raw_data_original$START),
  identical(raw_data$END, raw_data_original$END)
)

original_columns_to_check <- setdiff(names(raw_data_original), "GEN")
stopifnot(identical(
  raw_data[, original_columns_to_check, drop = FALSE],
  raw_data_original[, original_columns_to_check, drop = FALSE]
))

for (i in seq_len(nrow(raw_data))) {
  original_genes <- get_gene_set(raw_data$GEN_ORIGINAL[i])
  cleaned_genes <- get_gene_set(raw_data$GEN[i])
  if (length(original_genes) > 0) stopifnot(setequal(original_genes, cleaned_genes))

  primary_genes <- get_gene_set(raw_data$ANNOTATED_GENE[i])
  other_genes <- get_gene_set(raw_data$OTHER_GENES[i])
  stopifnot(length(intersect(primary_genes, other_genes)) == 0)
}

raw_data <- raw_data %>%
  dplyr::mutate(dplyr::across(where(is.character), ~ {
    x <- trimws(.x)
    x[!is.na(x) & x == ""] <- NA_character_
    x
  }))

# 17. Prepare and export the final dataset and quality-control reports
# The final Annotated_Data sheet keeps ANNOTATED_GENE and OTHER_GENES as the only gene columns.
# Intermediate original and overlap columns remain available in QC comparison tables.

output_dir <- "results/tables/gene_annotation"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write.table(
  annotation_comparison_chipseeker,
  file = file.path(output_dir, "Gene_annotation_comparison_ChIPseeker.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)
write.table(
  annotation_summary_chipseeker,
  file = file.path(output_dir, "Gene_annotation_summary_ChIPseeker.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)
write.table(
  annotation_comparison_overlapping,
  file = file.path(output_dir, "Gene_annotation_comparison_Overlapping.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)
write.table(
  annotation_summary_overlapping,
  file = file.path(output_dir, "Gene_annotation_summary_Overlapping.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)
write.table(
  similarity_summary,
  file = file.path(output_dir, "Gene_annotation_similarity_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)

# Remove intermediate gene columns from the final user-facing dataset.
# Keep only ANNOTATED_GENE and OTHER_GENES as gene annotation columns.

intermediate_gene_columns <- c(
  "GEN",
  "GEN_ORIGINAL",
  "GENES_OVERLAPPING",
  "N_OVERLAPPING_GENES"
)

final_annotated_data <- raw_data %>%
  dplyr::select(
    -dplyr::any_of(intermediate_gene_columns)
  )

# Verify the final structure
stopifnot(
  "ANNOTATED_GENE" %in% names(final_annotated_data),
  "OTHER_GENES" %in% names(final_annotated_data),
  !any(
    intermediate_gene_columns %in%
      names(final_annotated_data)
  ),
  nrow(final_annotated_data) == nrow(raw_data)
)

# Verify that the primary gene is not duplicated among additional genes
for (i in seq_len(nrow(final_annotated_data))) {

  primary_genes <- get_gene_set(
    final_annotated_data$ANNOTATED_GENE[i]
  )

  other_genes <- get_gene_set(
    final_annotated_data$OTHER_GENES[i]
  )

  stopifnot(
    length(intersect(primary_genes, other_genes)) == 0
  )
}

wb <- openxlsx::createWorkbook()
sheets_to_export <- list(
  "Annotated_Data" = final_annotated_data,
  "Unannotated_Before" = unannotated_before,
  "Unannotated_After" = unannotated_after,
  "Duplicate_ID_Assessment" = duplicate_id_assessment,
  "Duplicate_REGION_IDs_Removed" = duplicate_rows_removed,
  "Comparison_ChIPseeker" = annotation_comparison_chipseeker,
  "Summary_ChIPseeker" = annotation_summary_chipseeker,
  "Comparison_Overlapping" = annotation_comparison_overlapping,
  "Summary_Overlapping" = annotation_summary_overlapping,
  "Similarity_Summary" = similarity_summary,
  "Duplicate_Coordinates" = duplicate_coordinates
)

for (sheet_name in names(sheets_to_export)) {
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(
    wb, sheet = sheet_name, x = sheets_to_export[[sheet_name]],
    withFilter = TRUE, keepNA = TRUE, na.string = "NA"
  )
  openxlsx::freezePane(wb, sheet = sheet_name, firstRow = TRUE)
}

output_excel <- file.path("data", "raw_data_annotated_hg19.xlsx")
openxlsx::saveWorkbook(wb, output_excel, overwrite = TRUE)

# 18. Reload the workbook and verify the exported results
# Confirms that the expected worksheets exist and that the final table has the expected structure.

check_wb <- openxlsx::loadWorkbook(output_excel)
expected_sheets <- names(sheets_to_export)
stopifnot(
  all(expected_sheets %in% names(check_wb)),
  nrow(final_annotated_data) == nrow(raw_data_original),
  "ANNOTATED_GENE" %in% names(final_annotated_data),
  "OTHER_GENES" %in% names(final_annotated_data)
)

cat("\n===== FINAL RESULTS =====\n")
cat("Final number of regions:", nrow(final_annotated_data), "\n")
cat("Excel output:", output_excel, "\n")
cat("TSV files saved in:", output_dir, "\n")
cat("Worksheets exported:", length(expected_sheets), "\n")
cat("\n===== FINAL DATA STRUCTURE =====\n")
str(final_annotated_data)
cat("\n===== EXPORTED WORKSHEETS =====\n")
print(names(check_wb))