
library(PharmacoGx)
library(CoreGx)
library(SummarizedExperiment)
library(MultiAssayExperiment)

# ============================================
# STEP 3: ASSEMBLE THE PHARMACOSET
# ============================================

output_dir <- "../outputs"

# Load outputs from Steps 1 and 2
sample_metadata <- readRDS(
  file.path(output_dir, "sample_metadata.rds")
)

se <- readRDS(
  file.path(output_dir, "mutation_experiment.rds")
)

# Load the PharmacoSet created in Task 1
task1_pset <- readRDS(
  "../../task_1/outputs/AnnotationDB_pset.rds"
)

# Extract its existing treatmentResponse object
treatment_response <- slot(task1_pset, "treatmentResponse")

# Verify its class
print(class(treatment_response))

stopifnot(
  inherits(treatment_response, "LongTable") ||
    inherits(treatment_response, "LongTableDataMapper")
)

# ============================================
# 3.1 PREPARE THE SAMPLE TABLE
# ============================================

sample_table <- as.data.frame(sample_metadata)

# Use raw cell line identifiers as row names
rownames(sample_table) <- sample_table$sampleid

stopifnot(
  "sampleid" %in% names(sample_table),
  !anyNA(sample_table$sampleid),
  !anyDuplicated(sample_table$sampleid),
  identical(rownames(sample_table), sample_table$sampleid)
)

cat("Sample table rows:", nrow(sample_table), "\n")


# ============================================
# 3.2 CURATE IDENTIFIERS
# ============================================

# Raw identifiers are already the canonical sample IDs.
# Both columns therefore contain the same values.

curation_sample <- data.frame(
  sampleid = sample_table$sampleid,
  sampleid_dataset = sample_table$sampleid,
  row.names = sample_table$sampleid,
  stringsAsFactors = FALSE
)

# No raw tissue labels require curation.
curation_tissue <- data.frame()

# No drug treatments in this task.
curation_treatment <- data.frame()

curation_list <- list(
  sample = curation_sample,
  tissue = curation_tissue,
  treatment = curation_treatment
)

stopifnot(
  nrow(curation_sample) == nrow(sample_table),
  identical(
    rownames(curation_sample),
    rownames(sample_table)
  ),
  identical(
    curation_sample$sampleid,
    curation_sample$sampleid_dataset
  )
)

cat("Sample curation completed.\n")


# ============================================
# 3.3 BUILD MULTIASSAYEXPERIMENT
# ============================================

# Verify identifiers match before construction.
# This prevents silent loss of cell lines.

stopifnot(
  identical(
    colnames(se),
    rownames(sample_table)
  )
)

# One named experiment containing mutation profiles
experiment_list <- list(
  mutation.exome = se
)

mae <- MultiAssayExperiment(
  experiments = experiment_list,
  colData = sample_table
)

# Validate that harmonization preserved every sample
stopifnot(
  validObject(mae),
  length(experiments(mae)) == 1L,
  identical(
    rownames(colData(mae)),
    rownames(sample_table)
  ),
  identical(
    colnames(mae[["mutation.exome"]]),
    colnames(se)
  ),
  ncol(mae[["mutation.exome"]]) == nrow(sample_table),
  nrow(mae[["mutation.exome"]]) == nrow(se),
  nrow(sampleMap(mae)) == ncol(se)
)

cat("\n========== MULTIASSAYEXPERIMENT ==========\n")
print(mae)

cat("Samples before:", nrow(sample_table), "\n")
cat("Samples after:", nrow(colData(mae)), "\n")
cat("Mutation experiment columns:",
    ncol(mae[["mutation.exome"]]), "\n")


# ============================================
# 3.4 CREATE PHARMACOSET
# ============================================

pset <- PharmacoGx::PharmacoSet2(
  name = "Mutation Exome Cell Line Dataset",
  sample = sample_table,
  treatment = data.frame(),
  molecularProfiles = mae,
  treatmentResponse = treatment_response,
  curation = curation_list,
  datasetType = "sensitivity"
)

cat("\nPharmacoSet constructed.\n")
print(pset)


# ============================================
# 3.5 VALIDATE AND SAVE
# ============================================

# 3.5.1 Validate S4 object
stopifnot(validObject(pset))


# 3.5.2 Retrieve the full MultiAssayExperiment
pset_samples <- PharmacoGx::sampleInfo(pset)

pset_mae <- PharmacoGx::molecularProfilesSlot(pset)

# Extract the complete SummarizedExperiment
pset_mutations <- pset_mae[["mutation.exome"]]

pset_names <- PharmacoGx::mDataNames(pset)

pset_curation <- PharmacoGx::curation(pset)

# Confirm the retrieved object types
stopifnot(
  inherits(pset_mae, "MultiAssayExperiment"),
  inherits(pset_mutations, "SummarizedExperiment")
)


# 3.5.3 Validate sample preservation
stopifnot(
  nrow(pset_samples) == nrow(sample_table),
  setequal(
    rownames(pset_samples),
    sample_table$sampleid
  ),
  ncol(pset_mutations) == ncol(se),
  identical(
    colnames(pset_mutations),
    colnames(se)
  )
)


# 3.5.4 Validate mutation experiment
stopifnot(
  "mutation.exome" %in% pset_names,
  
  # Compare mutation assay values
  identical(
    SummarizedExperiment::assay(
      pset_mutations,
      "value"
    ),
    SummarizedExperiment::assay(
      se,
      "value"
    )
  ),
  
  # Verify gene identifiers
  identical(
    rownames(pset_mutations),
    rownames(se)
  ),
  
  # Verify gene metadata
  identical(
    as.data.frame(rowData(pset_mutations)),
    as.data.frame(rowData(se))
  ),
  
  # Verify sample metadata
  identical(
    as.data.frame(colData(pset_mutations)),
    as.data.frame(colData(se))
  )
)

cat("Mutation experiment validated.\n")


# 3.5.5 Validate identifier curation
stopifnot(
  setequal(
    pset_curation$sample$sampleid,
    sample_table$sampleid
  ),
  all(
    pset_curation$sample$sampleid ==
      pset_curation$sample$sampleid_dataset
  ),
  nrow(pset_curation$tissue) == 0L,
  nrow(pset_curation$treatment) == 0L
)

cat("Identifier curation validated.\n")


# 3.5.6 Save output
saveRDS(
  pset,
  file.path(output_dir, "pharmacoset.rds")
)

saveRDS(
  sample_table,
  file.path(output_dir, "sample_table.rds")
)


# 3.5.7 Final summary
cat("\n========== FINAL VALIDATION ==========\n")

cat("Valid PharmacoSet:", validObject(pset), "\n")
cat("Samples:", nrow(pset_samples), "\n")
cat("Genes:", nrow(pset_mutations), "\n")
cat("Mutation experiment:",
    paste(pset_names, collapse = ", "), "\n")

cat("Mutation presence count:",
    sum(
      SummarizedExperiment::assay(
        pset_mutations,
        "value"
      )
    ), "\n")

cat("Sample curation rows:",
    nrow(pset_curation$sample), "\n")

cat("\nAll Step 3 checks passed.\n")
cat("Saved to:", file.path(output_dir, "pharmacoset.rds"), "\n")
