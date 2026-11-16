library(CoreGx)
library(PharmacoGx)
library(SummarizedExperiment)
library(MultiAssayExperiment)
library(S4Vectors)

# ============================================
# 1. CREATE PHARMACOSET
# ============================================

# Load FDA-approved compounds
fda_compounds <- readRDS("../outputs/fda_compounds.rds")

# Load cell lines
cell_lines <- readRDS("../outputs/cell_lines.rds")

#Check that the identifiers are unique and nonmissing
stopifnot(
  !anyNA(cell_lines$sampleid),
  !anyDuplicated(cell_lines$sampleid),
  !anyNA(fda_compounds$treatmentid),
  !anyDuplicated(fda_compounds$treatmentid)
)

cat("All sample and treatment identifiers are unique and nonmissing.\n")

rownames(fda_compounds) <- fda_compounds$treatmentid
rownames(cell_lines) <- cell_lines$sampleid

placeholder_ids <- data.frame(treatmentid = "placeholder", sampleid = "placeholder")

one_row_response <- CoreGx::TreatmentResponseExperiment(
  rowData = placeholder_ids["treatmentid"], rowIDs = "treatmentid",
  colData = placeholder_ids["sampleid"],    colIDs = "sampleid",
  assays = list(raw = placeholder_ids),     assayIDs = list(raw = names(placeholder_ids))
)
empty_response <- one_row_response[0, 0]

placeholder_experiment <- SummarizedExperiment::SummarizedExperiment(
  colData = S4Vectors::DataFrame(sampleid = character(), batchid = character())
)
placeholder_profiles <- MultiAssayExperiment::MultiAssayExperiment(list(placeholder = placeholder_experiment))

pset <- PharmacoSet2(
  name = "AnnotationDB",
  sample = cell_lines,
  treatment = fda_compounds,
  molecularProfiles = placeholder_profiles,
  treatmentResponse = empty_response,
  datasetType = "sensitivity"
)

# ============================================
# 2. VALIDATE PHARMACOSET
# ============================================

cat("\n========== PHARMACOSET VALIDATION ==========\n")

# 1. Check S4 object validity
stopifnot(isTRUE(validObject(pset)))
cat("S4 object validity: PASSED\n")

# 2. Retrieve metadata stored in the PharmacoSet
pset_samples <- CoreGx::sampleInfo(pset)
pset_treatments <- CoreGx::treatmentInfo(pset)

# 3. Verify row counts
stopifnot(
  nrow(pset_samples) == nrow(cell_lines),
  nrow(pset_treatments) == nrow(fda_compounds)
)

cat("Sample count:", nrow(pset_samples), "\n")
cat("Treatment count:", nrow(pset_treatments), "\n")

# 4. Verify required metadata columns
stopifnot(
  all(c(
    "accession",
    "cell_line_name",
    "age_at_sampling",
    "sex_of_cell"
  ) %in% names(pset_samples)),
  
  all(c(
    "name",
    "cid",
    "fda_approval"
  ) %in% names(pset_treatments))
)

cat("Required metadata columns: PASSED\n")

# 5. Verify unique, nonmissing identifiers
stopifnot(
  !anyNA(pset_samples$accession),
  !anyDuplicated(pset_samples$accession),
  
  !anyNA(pset_treatments$cid),
  !anyDuplicated(pset_treatments$cid)
)

cat("Identifier uniqueness: PASSED\n")

# 6. Verify data was preserved exactly
stopifnot(
  isTRUE(all.equal(
    pset_samples,
    cell_lines,
    check.attributes = FALSE
  )),
  
  isTRUE(all.equal(
    pset_treatments,
    fda_compounds,
    check.attributes = FALSE
  ))
)

cat("Metadata preservation: PASSED\n")

# 7. Verify all treatments are FDA-approved
stopifnot(
  all(pset_treatments$fda_approval %in% TRUE)
)

cat("FDA approval: PASSED\n")

# 8. Inspect sample and treatment identifier conventions
cat("\nSample identifier column present:",
    "sampleid" %in% names(pset_samples), "\n")

cat("Treatment identifier column present:",
    "treatmentid" %in% names(pset_treatments), "\n")

# 9. Verify placeholder experiments contain no observations
stopifnot(
  nrow(CoreGx::treatmentResponse(pset)) == 0,
  ncol(CoreGx::treatmentResponse(pset)) == 0
)

cat("Treatment response is empty: PASSED\n")

cat("\nAll structural and metadata tests passed.\n")

# ============================================
# 3. SAVE PHARMACOSET
# ============================================

saveRDS(
  pset,
  "../outputs/AnnotationDB_pset.rds"
)
