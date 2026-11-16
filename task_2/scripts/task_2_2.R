
library(rtracklayer)
library(S4Vectors)
library(SummarizedExperiment)

# ============================================
# 2.1 INSPECT RAW MUTATION DATA
# ============================================

mutations_df <- read.delim("../raw_data/mutations.tsv")
sample_metadata <- readRDS("../outputs/sample_metadata.rds")

cat("Unique cell lines:",
    length(unique(mutations_df$cell_line)), "\n")

cat("Unique gene symbols:",
    length(unique(mutations_df$SNPEFF_GENE_NAME)), "\n")

cat("Variant records:", nrow(mutations_df), "\n")

# Inspect predicted effects
cat("\nEffect counts:\n")
print(sort(
  table(mutations_df$SNPEFF_EFFECT, useNA = "ifany"),
  decreasing = TRUE
))

# Inspect predicted impacts
cat("\nImpact counts:\n")
print(sort(
  table(mutations_df$SNPEFF_IMPACT, useNA = "ifany"),
  decreasing = TRUE
))

# Depending on the nature of the experiment, only variants having a sufficiently high impact
# or type of effect should be kept. Since no context is given for the variant types to be studied
# I do not think there is any reason to filter the dataset on the basis of variant type. Likewise
# the impact is almost all moderate and 9 between low and modifier with the rest high. I therefore
# don't think there is much reason to motivate filtering on variant impact either.

# The resulting matrix will have 46 columns (cell lines) and 7466 rows (genes) based on the number
# of unique cell lines and genes in the data. I think this is a reasonably sized matrix and 
# therefore there is no need to filter out any of the variants at this stage.

# Validate identifiers before constructing the matrix
stopifnot(
  !anyNA(mutations_df$cell_line),
  all(nzchar(trimws(mutations_df$cell_line))),
  !anyNA(mutations_df$SNPEFF_GENE_NAME),
  all(nzchar(trimws(mutations_df$SNPEFF_GENE_NAME))),
  !anyDuplicated(sample_metadata$sampleid),
  setequal(unique(mutations_df$cell_line),
           sample_metadata$sampleid)
)


# ============================================
# 2.2 BUILD GENE BY CELL LINE MATRIX
# ============================================

# Unique gene symbols from mutation records
genes <- unique(mutations_df$SNPEFF_GENE_NAME)

# Use sample identifiers from Step 1
samples <- sample_metadata$sampleid

# Initialize binary matrix with zeros
mutation_matrix <- matrix(
  0L,
  nrow = length(genes),
  ncol = length(samples),
  dimnames = list(genes, samples)
)

# Identify matrix coordinates for each variant
gene_idx <- match(
  mutations_df$SNPEFF_GENE_NAME,
  rownames(mutation_matrix)
)

sample_idx <- match(
  mutations_df$cell_line,
  colnames(mutation_matrix)
)

stopifnot(
  !anyNA(gene_idx),
  !anyNA(sample_idx)
)

# Set gene-cell line combinations containing variants to 1
mutation_matrix[cbind(gene_idx, sample_idx)] <- 1L

# Validate binary matrix
stopifnot(
  all(mutation_matrix %in% c(0L, 1L)),
  !anyNA(mutation_matrix),
  !anyDuplicated(rownames(mutation_matrix)),
  identical(colnames(mutation_matrix), sample_metadata$sampleid)
)

cat("\n========== INITIAL MATRIX ==========\n")
cat("Genes:", nrow(mutation_matrix), "\n")
cat("Cell lines:", ncol(mutation_matrix), "\n")
cat("Gene-cell line combinations with variants:",
    sum(mutation_matrix), "\n")

print(mutation_matrix[1:10, 1:5, drop = FALSE])


# ============================================
# 2.3 MAP GENE SYMBOLS TO ENSEMBL GENE IDS
# ============================================

# 2.3.1 Load GENCODE v50 annotation

gtf <- rtracklayer::import(
  "../raw_data/gencode.v50.annotation.gtf"
)

# Keep only gene features
gtf_genes <- gtf[gtf$type == "gene"]

# Extract required gene annotations
gene_reference <- data.frame(
  gene_id = sub("\\..*$", "", gtf_genes$gene_id),
  gene_name = gtf_genes$gene_name,
  chromosome = as.character(
    GenomicRanges::seqnames(gtf_genes)
  ),
  start = BiocGenerics::start(gtf_genes),
  end = BiocGenerics::end(gtf_genes),
  strand = as.character(
    BiocGenerics::strand(gtf_genes)
  ),
  stringsAsFactors = FALSE
)

# Exclude reference records with invalid identifiers
gene_reference <- gene_reference[
  !is.na(gene_reference$gene_id) &
    nzchar(trimws(gene_reference$gene_id)) &
    !is.na(gene_reference$gene_name) &
    nzchar(trimws(gene_reference$gene_name)),
  ,
  drop = FALSE
]

cat("\nGENCODE gene records:",
    nrow(gene_reference), "\n")


# --------------------------------------------
# 2.3.2 Identify uniquely mapped symbols
# --------------------------------------------

# Count distinct Ensembl IDs associated with each symbol
symbol_id_counts <- tapply(
  gene_reference$gene_id,
  gene_reference$gene_name,
  function(x) length(unique(x))
)

# Retain symbols associated with exactly one Ensembl ID
unique_symbols <- names(symbol_id_counts)[
  symbol_id_counts == 1
]

# Identify ambiguous symbols for exclusion logging
ambiguous_symbols <- names(symbol_id_counts)[
  symbol_id_counts > 1
]

gene_lookup <- gene_reference[
  gene_reference$gene_name %in% unique_symbols,
  ,
  drop = FALSE
]

# Ensure one reference row per uniquely mapped symbol
gene_lookup <- gene_lookup[
  !duplicated(gene_lookup$gene_name),
  ,
  drop = FALSE
]


# --------------------------------------------
# 2.3.3 Match mutation genes to GENCODE
# --------------------------------------------

original_symbols <- rownames(mutation_matrix)

match_idx <- match(
  original_symbols,
  gene_lookup$gene_name
)

# Record symbols excluded from mapping
excluded_symbols <- original_symbols[is.na(match_idx)]

exclusion_log <- data.frame(
  gene_symbol = excluded_symbols,
  reason = ifelse(
    excluded_symbols %in% ambiguous_symbols,
    "Ambiguous symbol: maps to multiple Ensembl IDs",
    "Gene symbol not found in GENCODE v50"
  ),
  stringsAsFactors = FALSE
)

# Retain successfully mapped mutation rows
mapped_matrix <- mutation_matrix[
  !is.na(match_idx),
  ,
  drop = FALSE
]

mapped_genes <- gene_lookup[
  match_idx[!is.na(match_idx)],
  ,
  drop = FALSE
]


# --------------------------------------------
# 2.3.4 Combine duplicate Ensembl mappings
# --------------------------------------------

# If multiple symbols map to the same Ensembl ID,
# combine their binary mutation rows using logical OR.

ensembl_ids <- unique(mapped_genes$gene_id)

ensembl_matrix <- matrix(
  0L,
  nrow = length(ensembl_ids),
  ncol = ncol(mapped_matrix),
  dimnames = list(
    ensembl_ids,
    colnames(mapped_matrix)
  )
)

for (i in seq_len(nrow(mapped_matrix))) {
  
  id <- mapped_genes$gene_id[i]
  
  ensembl_matrix[id, ] <- as.integer(
    ensembl_matrix[id, ] == 1L |
      mapped_matrix[i, ] == 1L
  )
}

# Identify Ensembl IDs associated with multiple symbols
duplicate_ids <- unique(
  mapped_genes$gene_id[
    duplicated(mapped_genes$gene_id)
  ]
)

cat("\nEnsembl IDs receiving multiple symbols:",
    length(duplicate_ids), "\n")


# --------------------------------------------
# 2.3.5 Validate gene mapping
# --------------------------------------------

stopifnot(
  !anyNA(rownames(ensembl_matrix)),
  !anyDuplicated(rownames(ensembl_matrix)),
  all(grepl("^ENSG[0-9]+$", rownames(ensembl_matrix))),
  all(ensembl_matrix %in% c(0L, 1L)),
  !anyNA(ensembl_matrix),
  identical(
    colnames(ensembl_matrix),
    colnames(mutation_matrix)
  ),
  nrow(mapped_matrix) + nrow(exclusion_log) ==
    nrow(mutation_matrix)
)

# Independently check expected mutation presence after mapping
expected_pairs <- unique(data.frame(
  gene_id = mapped_genes$gene_id[
    row(mapped_matrix)[mapped_matrix == 1L]
  ],
  sampleid = colnames(mapped_matrix)[
    col(mapped_matrix)[mapped_matrix == 1L]
  ]
))

stopifnot(
  sum(ensembl_matrix) == nrow(expected_pairs)
)

cat("\n========== GENE MAPPING SUMMARY ==========\n")
cat("Original gene symbols:", nrow(mutation_matrix), "\n")
cat("Mapped gene symbols:", nrow(mapped_matrix), "\n")
cat("Unique Ensembl genes:", nrow(ensembl_matrix), "\n")
cat("Excluded symbols:", nrow(exclusion_log), "\n")
cat("Original matrix 1s:", sum(mutation_matrix), "\n")
cat("Mapped matrix 1s:", sum(ensembl_matrix), "\n")
cat("Mapping rate:",
    round(100 * nrow(mapped_matrix) /
            nrow(mutation_matrix), 2),
    "%\n")


# --------------------------------------------
# 2.3.6 Save mapping outputs
# --------------------------------------------

saveRDS(
  ensembl_matrix,
  "../outputs/ensembl_matrix.rds"
)

saveRDS(
  exclusion_log,
  "../outputs/gene_exclusions.rds"
)


# ============================================
# 2.4 BUILD ROWDATA AND COLDATA
# ============================================

# --------------------------------------------
# 2.4.1 Build gene annotation table (rowData)
# --------------------------------------------

# One annotation row per Ensembl ID
gene_table <- mapped_genes[
  !duplicated(mapped_genes$gene_id),
  ,
  drop = FALSE
]

# Preserve original symbols if multiple map to one ID
symbols_by_id <- split(
  mapped_genes$gene_name,
  mapped_genes$gene_id
)

gene_table$original_symbol <- vapply(
  gene_table$gene_id,
  function(id) {
    paste(
      unique(symbols_by_id[[id]]),
      collapse = ";"
    )
  },
  character(1)
)

# Reorder gene annotations to match matrix rows
gene_table <- gene_table[
  match(
    rownames(ensembl_matrix),
    gene_table$gene_id
  ),
  ,
  drop = FALSE
]

rownames(gene_table) <- gene_table$gene_id


# --------------------------------------------
# 2.4.2 Build cell line annotation table (colData)
# --------------------------------------------

# Reorder metadata to match matrix columns
col_table <- sample_metadata[
  match(
    colnames(ensembl_matrix),
    sample_metadata$sampleid
  ),
  ,
  drop = FALSE
]

# Add required batch identifier
col_table$batchid <- NA_character_

# Set row names to sample identifiers
rownames(col_table) <- col_table$sampleid


# --------------------------------------------
# 2.4.3 Validate annotation tables
# --------------------------------------------

stopifnot(
  # Gene metadata
  identical(
    rownames(gene_table),
    rownames(ensembl_matrix)
  ),
  identical(
    gene_table$gene_id,
    rownames(ensembl_matrix)
  ),
  !anyDuplicated(gene_table$gene_id),
  
  # Cell line metadata
  identical(
    rownames(col_table),
    colnames(ensembl_matrix)
  ),
  identical(
    col_table$sampleid,
    colnames(ensembl_matrix)
  ),
  !anyDuplicated(col_table$sampleid),
  
  # Required columns
  all(
    c(
      "gene_id", "original_symbol", "gene_name",
      "chromosome", "start", "end", "strand"
    ) %in% names(gene_table)
  ),
  all(
    c("sampleid", "batchid") %in% names(col_table)
  )
)

cat("\nRowData and colData validation passed.\n")


# --------------------------------------------
# 2.4.4 Convert to Bioconductor DataFrames
# --------------------------------------------

row_data <- S4Vectors::DataFrame(gene_table)
col_data <- S4Vectors::DataFrame(col_table)


# --------------------------------------------
# 2.4.5 Save annotation tables
# --------------------------------------------

saveRDS(
  gene_table,
  "../outputs/row_data.rds"
)

saveRDS(
  col_table,
  "../outputs/col_data.rds"
)


# --------------------------------------------
# 2.4.6 Inspect results
# --------------------------------------------

cat("\nGene annotations:\n")
print(head(row_data))

cat("\nCell line annotations:\n")
print(head(col_data))

cat("\n========== FINAL DIMENSIONS ==========\n")
cat("Mutation matrix:", dim(ensembl_matrix), "\n")
cat("RowData:", dim(row_data), "\n")
cat("ColData:", dim(col_data), "\n")

cat("\nSections 2.1-2.4 completed successfully.\n")


# ============================================
# 2.5 CREATE SUMMARIZEDEXPERIMENT
# ============================================

# 2.5.1 Construct SummarizedExperiment
se <- SummarizedExperiment(
  assays = list(value = ensembl_matrix),
  rowData = row_data,
  colData = col_data
)

# 2.5.2 Add experiment metadata
metadata(se)$annotation <- "mutation.exome"
metadata(se)$value_type <- "binary source-reported variant presence"

# 2.5.3 Validate object
stopifnot(
  validObject(se),
  identical(rownames(se), rownames(ensembl_matrix)),
  identical(colnames(se), colnames(ensembl_matrix)),
  identical(assay(se, "value"), ensembl_matrix)
)

# 2.5.4 Save SummarizedExperiment
saveRDS(
  se,
  "../outputs/mutation_experiment.rds"
)

# 2.5.5 Inspect result
print(se)
print(metadata(se))

cat("\nSummarizedExperiment created successfully.\n")


# ============================================
# 2.6 VALIDATE SUMMARIZEDEXPERIMENT
# ============================================

# 2.6.1 Validate object structure
stopifnot(
  validObject(se),
  identical(rownames(se), rownames(rowData(se))),
  identical(colnames(se), rownames(colData(se))),
  identical(colnames(se), sample_metadata$sampleid)
)

cat("Experiment structure validated.\n")


# 2.6.2 Validate mutation values
values <- assay(se, "value")

stopifnot(
  !anyNA(values),
  all(values %in% c(0L, 1L))
)

cat("Binary mutation values validated.\n")


# 2.6.3 Validate Ensembl gene identifiers
gene_ids <- rownames(se)

stopifnot(
  !anyNA(gene_ids),
  !anyDuplicated(gene_ids),
  all(grepl("^ENSG[0-9]+$", gene_ids)),
  identical(gene_ids, rowData(se)$gene_id)
)

cat("Ensembl gene identifiers validated.\n")


# 2.6.4 Validate cell line identifiers
sample_ids <- colnames(se)

stopifnot(
  !anyNA(sample_ids),
  !anyDuplicated(sample_ids),
  all(sample_ids %in% sample_metadata$sampleid),
  identical(sample_ids, colData(se)$sampleid),
  "batchid" %in% names(colData(se))
)

cat("Cell line identifiers validated.\n")


# 2.6.5 Validate mutation preservation

# Number of gene-cell line combinations before mapping
original_ones <- sum(mutation_matrix)

# Number after excluding unmapped symbols, before combining IDs
mapped_ones <- sum(mapped_matrix)

# Number in final SummarizedExperiment
final_ones <- sum(values)

# Number lost because symbols could not be mapped
excluded_ones <- original_ones - mapped_ones

# Number merged because different symbols mapped to the same ID
merged_ones <- mapped_ones - final_ones

stopifnot(
  excluded_ones >= 0,
  merged_ones >= 0,
  final_ones == sum(ensembl_matrix)
)

cat("\n========== MUTATION PRESERVATION ==========\n")
cat("Original gene-cell line combinations:", original_ones, "\n")
cat("Excluded due to gene mapping:", excluded_ones, "\n")
cat("Merged due to duplicate Ensembl IDs:", merged_ones, "\n")
cat("Final gene-cell line combinations:", final_ones, "\n")

cat("\nAccounting check:",
    original_ones == excluded_ones + merged_ones + final_ones,
    "\n")


# 2.6.6 Validate annotation metadata
stopifnot(
  identical(metadata(se)$annotation, "mutation.exome"),
  identical(
    metadata(se)$value_type,
    "binary source-reported variant presence"
  ),
  all(c(
    "gene_id", "original_symbol", "gene_name",
    "chromosome", "start", "end", "strand"
  ) %in% names(rowData(se)))
)

cat("\nExperiment metadata validated.\n")


# 2.6.7 Final validation summary
cat("\n========== FINAL VALIDATION SUMMARY ==========\n")

cat("Genes:", nrow(se), "\n")
cat("Cell lines:", ncol(se), "\n")
cat("Reported mutation combinations:", sum(values), "\n")
cat("Excluded gene symbols:", nrow(exclusion_log), "\n")
cat("Cell lines without AnnotationDB matches:",
    sum(!colData(se)$In.AnnotationDB), "\n")

cat("\nAll validation checks passed.\n")
