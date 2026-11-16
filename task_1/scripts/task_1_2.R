
library(httr)
library(jsonlite)

# ============================================
# 1. FETCH ALL CELL LINES
# ============================================

url_clines <- "https://annotationdb.bhklab.ca/cell_line/all"

response_clines <- GET(url_clines, timeout(60))
stop_for_status(response_clines)

raw_clines <- content(
  response_clines,
  "text",
  encoding = "UTF-8"
)

json_clines <- fromJSON(raw_clines)

cat("Total cell lines:", nrow(json_clines), "\n")


# ============================================
# 2. PREPARE CELL LINE BATCHES
# ============================================

url_many_base <- "https://annotationdb.bhklab.ca/cell_line/many"

cell_lines <- unique(json_clines$name[
  !is.na(json_clines$name) &
    trimws(json_clines$name) != ""
])

batch_size <- 128

cell_line_batches <- split(
  cell_lines,
  ceiling(seq_along(cell_lines) / batch_size)
)


# ============================================
# 3. FETCH CELL LINE METADATA
# ============================================

fetch_batch <- function(batch) {
  
  query_params <- list(
    cell_lines = paste(batch, collapse = ",")
  )
  
  response <- GET(
    url_many_base,
    query = query_params,
    timeout(60)
  )
  
  # Split batch if URL is too long
  if (status_code(response) == 414) {
    
    cat(
      "HTTP 414: Splitting batch of",
      length(batch),
      "cell lines\n"
    )
    
    if (length(batch) == 1) {
      stop("HTTP 414 for single cell line: ", batch)
    }
    
    midpoint <- floor(length(batch) / 2)
    
    return(rbind(
      fetch_batch(batch[seq_len(midpoint)]),
      fetch_batch(batch[seq.int(midpoint + 1, length(batch))])
    ))
  }
  
  stop_for_status(response)
  
  raw <- content(
    response,
    "text",
    encoding = "UTF-8"
  )
  
  fromJSON(raw)
}

# Fetch all batches
json_clines_details <- vector(
  "list",
  length(cell_line_batches)
)

for (i in seq_along(cell_line_batches)) {
  
  json_clines_details[[i]] <- fetch_batch(
    cell_line_batches[[i]]
  )
  
  cat(
    "Completed batch",
    i,
    "of",
    length(cell_line_batches),
    "\n"
  )
}

json_clines_details <- do.call(
  rbind,
  json_clines_details
)


# ============================================
# 4. MATCH METADATA TO CELL LINES
# ============================================

# Match using accession
match_idx <- match(
  json_clines$accession,
  json_clines_details$accession
)

cat(
  "Matched cell lines:",
  sum(!is.na(match_idx)),
  "\n"
)

cat(
  "Unmatched cell lines:",
  sum(is.na(match_idx)),
  "\n"
)

# Select requested metadata
metadata <- json_clines_details[
  match_idx,
  c(
    "cell_line_name",
    "age_at_sampling",
    "sex_of_cell"
  ),
  drop = FALSE
]

# Combine with original accessions
final_clines <- data.frame(
  accession = json_clines$accession,
  metadata,
  row.names = NULL
)

# Standardize sample identifiers
final_clines$sampleid <- final_clines$accession

# ============================================
# 5. BASIC VERIFICATION
# ============================================

cat("\n========== VERIFICATION ==========\n")

cat(
  "Original cell lines:",
  nrow(json_clines),
  "\n"
)

cat(
  "Detailed records returned:",
  nrow(json_clines_details),
  "\n"
)

cat(
  "Final dataframe rows:",
  nrow(final_clines),
  "\n"
)

cat(
  "Missing ages:",
  sum(is.na(final_clines$age_at_sampling) |
        final_clines$age_at_sampling == ""),
  "\n"
)

cat(
  "Missing sexes:",
  sum(is.na(final_clines$sex_of_cell) |
        final_clines$sex_of_cell == ""),
  "\n"
)

stopifnot(
  nrow(final_clines) == nrow(json_clines),
  !anyDuplicated(final_clines$accession),
  identical(
    final_clines$accession,
    json_clines$accession
  )
)

cat("\nBasic verification passed.\n")


# ============================================
# 6. SAVE FINAL DATASET
# ============================================

saveRDS(
  final_clines,
  "../outputs/cell_lines.rds"
)

cat("\nCell line datasets saved successfully.\n")
