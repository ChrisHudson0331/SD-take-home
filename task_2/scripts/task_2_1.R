
library(httr)
library(jsonlite)

# ============================================
# 1. LOAD DATA
# ============================================

mutations_df <- read.delim("../raw_data/mutations.tsv")
clines_df <- readRDS("../../task_1/outputs/cell_lines.rds")

# ============================================
# 2. CONSTRUCT SAMPLE IDENTIFIERS
# ============================================

# Preserve original sample IDs from mutations.tsv
sample_metadata <- data.frame(
  sampleid = unique(mutations_df$cell_line),
  stringsAsFactors = FALSE
)

# Remove supplier suffix to obtain cell line names
sample_metadata$name <- sub(
  "_[^_]+$",
  "",
  sample_metadata$sampleid
)

# ============================================
# 3. MATCH CELL LINES TO ANNOTATIONDB CATALOG
# ============================================

# Normalize names by removing punctuation and capitalization
normalize_name <- function(x) {
  toupper(gsub("[^A-Za-z0-9]", "", x))
}

normalized_samples <- normalize_name(sample_metadata$name)
normalized_catalog <- normalize_name(clines_df$cell_line_name)

# First attempt exact-name matching
catalog_idx <- match(
  sample_metadata$name,
  clines_df$cell_line_name
)

# Use normalized matching for names without exact matches
unmatched <- which(is.na(catalog_idx))

for (i in unmatched) {
  
  candidates <- which(
    normalized_catalog == normalized_samples[i]
  )
  
  if (length(candidates) == 1) {
    catalog_idx[i] <- candidates
    
  } else if (length(candidates) > 1) {
    warning(
      "Ambiguous catalog match for ",
      sample_metadata$name[i],
      ": ",
      paste(clines_df$cell_line_name[candidates], collapse = ", ")
    )
  }
}

# Retrieve Cellosaurus accessions
sample_metadata$accession <- clines_df$accession[catalog_idx]

# Known alias documented in assignment example
sample_metadata$accession[
  sample_metadata$name == "U266"
] <- "CVCL_0566"

# ============================================
# 4. QUERY ANNOTATIONDB FOR CELL LINE DETAILS
# ============================================

# Use names recognized by AnnotationDB
query_names <- unique(na.omit(
  clines_df$cell_line_name[catalog_idx]
))

# Include U266 alias
if ("U266" %in% sample_metadata$name) {
  query_names <- unique(c(query_names, "U266B1"))
}

url <- "https://annotationdb.bhklab.ca/cell_line/many"

response <- GET(
  url,
  query = list(cell_lines = paste(query_names, collapse = ",")),
  timeout(60)
)

stop_for_status(response)

json_clines_details <- fromJSON(
  content(response, "text", encoding = "UTF-8")
)

if (length(json_clines_details) == 0) {
  stop("AnnotationDB returned no cell line details.")
}

# ============================================
# 5. MATCH API RESULTS BY CELLOSAURUS ACCESSION
# ============================================

match_idx <- match(
  sample_metadata$accession,
  json_clines_details$accession
)

# ============================================
# 6. EXTRACT CELL LINE METADATA
# ============================================

sample_metadata$age <- json_clines_details$age_at_sampling[match_idx]
sample_metadata$sex <- json_clines_details$sex_of_cell[match_idx]
sample_metadata$population <- json_clines_details$population[match_idx]

# Extract first available disease description
sample_metadata$disease <- vapply(
  match_idx,
  function(i) {
    
    if (is.na(i)) return(NA_character_)
    
    diseases <- json_clines_details$diseases[[i]]
    
    if (is.null(diseases) || nrow(diseases) == 0) {
      return(NA_character_)
    }
    
    descriptions <- diseases$description
    descriptions <- descriptions[
      !is.na(descriptions) & nzchar(trimws(descriptions))
    ]
    
    if (length(descriptions) == 0) {
      return(NA_character_)
    }
    
    descriptions[1]
  },
  character(1)
)

# ============================================
# 7. ADD REQUIRED FIELDS
# ============================================

# TRUE only when detailed API metadata was found
sample_metadata$In.AnnotationDB <- !is.na(match_idx)

# Preserve the exact title returned by AnnotationDB
sample_metadata$AnnotationDB.Title <-
  json_clines_details$cell_line_name[match_idx]

# Tissue identifier specified by assignment
sample_metadata$tissueid <- "haematopoietic_and_lymphoid"

# Set row names to sample IDs
rownames(sample_metadata) <- sample_metadata$sampleid

# ============================================
# 8. CLEAN AND VALIDATE SAMPLE METADATA
# ============================================

# --------------------------------------------
# 8.1 Standardize missing metadata values
# --------------------------------------------

metadata_cols <- c("age", "sex", "population", "disease")

for (col in metadata_cols) {
  x <- sample_metadata[[col]]
  x[!is.na(x) & trimws(x) == ""] <- NA_character_
  sample_metadata[[col]] <- x
}

# --------------------------------------------
# 8.2 Validate sample identifiers
# --------------------------------------------

stopifnot(
  nrow(sample_metadata) == length(unique(mutations_df$cell_line)),
  !anyNA(sample_metadata$sampleid),
  !anyDuplicated(sample_metadata$sampleid),
  identical(rownames(sample_metadata), sample_metadata$sampleid),
  setequal(sample_metadata$sampleid, unique(mutations_df$cell_line))
)

cat("Sample identifier checks passed.\n")

# --------------------------------------------
# 8.3 Validate required columns
# --------------------------------------------

required_cols <- c(
  "sampleid",
  "name",
  "accession",
  "age",
  "sex",
  "population",
  "disease",
  "In.AnnotationDB",
  "AnnotationDB.Title",
  "tissueid"
)

stopifnot(
  all(required_cols %in% names(sample_metadata)),
  !anyNA(sample_metadata$name),
  all(sample_metadata$tissueid == "haematopoietic_and_lymphoid")
)

cat("Required column checks passed.\n")

# --------------------------------------------
# 8.4 Validate AnnotationDB matches
# --------------------------------------------

matched <- sample_metadata$In.AnnotationDB

stopifnot(
  !anyNA(matched),
  all(!is.na(sample_metadata$accession[matched])),
  all(!is.na(sample_metadata$AnnotationDB.Title[matched])),
  all(sample_metadata$accession[matched] %in%
        json_clines_details$accession)
)

# Verify returned titles match the corresponding accessions
api_idx <- match(
  sample_metadata$accession[matched],
  json_clines_details$accession
)

stopifnot(
  identical(
    sample_metadata$AnnotationDB.Title[matched],
    json_clines_details$cell_line_name[api_idx]
  )
)

cat("AnnotationDB match checks passed.\n")

# --------------------------------------------
# 8.5 Check duplicate accessions
# --------------------------------------------

# Different sample IDs can legitimately share an accession.
# Report duplicates without rejecting them.

duplicate_accessions <- sample_metadata[
  !is.na(sample_metadata$accession) &
    (duplicated(sample_metadata$accession) |
       duplicated(sample_metadata$accession, fromLast = TRUE)),
  c("sampleid", "name", "accession"),
  drop = FALSE
]

cat("\nDuplicate accession records:",
    nrow(duplicate_accessions), "\n")

if (nrow(duplicate_accessions) > 0) {
  print(duplicate_accessions)
}

# --------------------------------------------
# 8.6 Check ambiguous API accessions
# --------------------------------------------

api_accessions <- json_clines_details$accession
api_accessions <- api_accessions[!is.na(api_accessions)]

ambiguous_api <- unique(api_accessions[
  duplicated(api_accessions)
])

ambiguous_api <- intersect(
  ambiguous_api,
  sample_metadata$accession[matched]
)

cat("Matched accessions with multiple API records:",
    length(ambiguous_api), "\n")

if (length(ambiguous_api) > 0) {
  print(ambiguous_api)
}

# --------------------------------------------
# 8.7 Report unmatched cell lines
# --------------------------------------------

unmatched_samples <- sample_metadata[
  !matched,
  c("sampleid", "name", "accession"),
  drop = FALSE
]

cat("\nUnmatched cell lines:",
    nrow(unmatched_samples), "\n")

if (nrow(unmatched_samples) > 0) {
  print(unmatched_samples)
}

# --------------------------------------------
# 8.8 Report missing metadata
# --------------------------------------------

missing_counts <- colSums(
  is.na(sample_metadata[metadata_cols])
)

cat("\nMissing metadata values:\n")
print(missing_counts)

incomplete_samples <- sample_metadata[
  matched &
    !complete.cases(sample_metadata[metadata_cols]),
  c("sampleid", metadata_cols),
  drop = FALSE
]

cat("\nMatched cell lines with incomplete metadata:",
    nrow(incomplete_samples), "\n")

if (nrow(incomplete_samples) > 0) {
  print(incomplete_samples)
}

# --------------------------------------------
# 8.9 Final validation summary
# --------------------------------------------

cat("\n========== VALIDATION SUMMARY ==========\n")

cat("Total samples:", nrow(sample_metadata), "\n")
cat("Matched samples:", sum(matched), "\n")
cat("Unmatched samples:", sum(!matched), "\n")
cat("Match rate:",
    round(100 * mean(matched), 2), "%\n")

cat("Unique sample IDs:",
    length(unique(sample_metadata$sampleid)), "\n")

cat("\nAll mandatory validation checks passed.\n")

# ============================================
# 9. SAVE RESULTS
# ============================================

saveRDS(
  sample_metadata,
  "../outputs/sample_metadata.rds"
)

cat("\nSaved sample metadata to ../outputs/sample_metadata.rds\n")
