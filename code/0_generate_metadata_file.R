# 0. Load required libraries
library(readxl)   # Imports and reads data directly from Excel files (.xls/.xlsx) into R data frames
library(dplyr)    # Grammar of data manipulation providing a consistent set of verbs to transform datasets
library(stringr)  # Provides a cohesive and user-friendly set of tools for regular expressions and string manipulation

# 1. Define input / output files location
input_file <- file.path("metadata/original/INFO-PATIENTS_Methylseq.xlsx")
output_file <- file.path("metadata/patients_metadata.txt")

# 2. Read original patient information
patients <- read_excel(input_file, sheet = "Hoja 1", skip = 1)

dim(patients)
print(as.data.frame(head(patients, n = 24)))

# 3. Build metadata table
patients_metadata <- patients %>%
  dplyr::select(
    pathology = PATHOLOGY,
    gene_variant = `GENE VARIANT`,
    sample_id = `SAMPLE ID`,
    relationship = RELATIONSHIP,
    genetic_status = `GENETIC STATUS`,
    clinical_status = `CLINICAL STATUS`
  ) %>%
  mutate(
    patient_id = as.integer(str_remove(sample_id, "^P")),
    across(
      c(pathology, gene_variant, relationship, genetic_status, clinical_status), 
      str_squish
    )
  ) %>%
  dplyr::select(
    patient_id,
    pathology,
    gene_variant,
    relationship,
    genetic_status,
    clinical_status
  ) %>%
  arrange(patient_id)

# 4. Write tab-delimited TXT
write.table(
  patients_metadata,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = "NA" 
)

# 5. Final Confirmation
message("Metadata file created: ", output_file)
message("Number of patients: ", nrow(patients_metadata))
message("Variables:\n", paste("-", colnames(patients_metadata), collapse = "\n"))

head(patients_metadata, n = 24)