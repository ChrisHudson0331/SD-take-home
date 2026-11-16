library(httr)
library(jsonlite)

# ============================================
# 1. FETCH FDA-APPROVED COMPOUNDS
# ============================================

url_fda_comp <- "https://annotationdb.bhklab.ca/compound/all"

response_fda_comp <- GET(url_fda_comp, timeout(60))
stop_for_status(response_fda_comp)

raw_fda_comp <- content(
  response_fda_comp,
  "text",
  encoding = "UTF-8"
)

json_fda_comp <- fromJSON(raw_fda_comp)

# Keep only FDA-approved compounds
json_fda_comp <- json_fda_comp[
  !is.na(json_fda_comp$fda_approval) &
    json_fda_comp$fda_approval == TRUE,
  ,
  drop = FALSE
]

json_fda_comp_original <- json_fda_comp

cat("FDA-approved compounds:", nrow(json_fda_comp), "\n")


# ============================================
# 2. PREPARE COMPOUND BATCHES
# ============================================

url_many_base <- "https://annotationdb.bhklab.ca/compound/many"

compounds <- json_fda_comp$name

# Remove missing, empty, and duplicate names
compounds <- unique(compounds[
  !is.na(compounds) & trimws(compounds) != ""
])

batch_size <- 128

compound_batches <- split(
  compounds,
  ceiling(seq_along(compounds) / batch_size)
)


# ============================================
# 3. FETCH MECHANISMS
# ============================================

fetch_batch <- function(batch) {
  
  url_many_end <- c(
    setNames(
      as.list(batch),
      rep("compound", length(batch))
    ),
    list(
      bioassay = "false",
      mechanism = "true",
      toxicity = "false",
      golden_bioassay = "false"
    )
  )
  
  response_fda_mech <- GET(
    url_many_base,
    query = url_many_end,
    timeout(60)
  )
  
  # Split batches automatically if URL is too long
  if (status_code(response_fda_mech) == 414) {
    
    cat(
      "HTTP 414: Splitting batch of",
      length(batch),
      "compounds\n"
    )
    
    if (length(batch) == 1) {
      stop("HTTP 414 for single compound: ", batch)
    }
    
    midpoint <- floor(length(batch) / 2)
    
    return(rbind(
      fetch_batch(batch[seq_len(midpoint)]),
      fetch_batch(batch[seq.int(midpoint + 1, length(batch))])
    ))
  }
  
  stop_for_status(response_fda_mech)
  
  raw_fda_mech <- content(
    response_fda_mech,
    "text",
    encoding = "UTF-8"
  )
  
  fromJSON(raw_fda_mech)
}

# Fetch all batches
json_fda_mech <- vector("list", length(compound_batches))

for (i in seq_along(compound_batches)) {
  
  json_fda_mech[[i]] <- fetch_batch(compound_batches[[i]])
  
  cat(
    "Completed batch",
    i,
    "of",
    length(compound_batches),
    "\n"
  )
}

json_fda_mech <- do.call(rbind, json_fda_mech)


# ============================================
# 4. EXTRACT FIRST MECHANISM
# ============================================

if (!all(c("query_field", "cid", "mechanisms") %in% names(json_fda_mech))) {
  stop("API response is missing required fields.")
}

# Collect all mechanism column names
mech_cols <- unique(unlist(lapply(
  json_fda_mech$mechanisms,
  names
)))

if (length(mech_cols) == 0) {
  stop("No mechanism fields found.")
}

# Extract the first mechanism for each API record
mech_df <- do.call(rbind, lapply(
  json_fda_mech$mechanisms,
  function(df) {
    
    if (is.null(df) || nrow(df) == 0) {
      return(as.data.frame(
        setNames(rep(list(NA), length(mech_cols)), mech_cols)
      ))
    }
    
    for (col in setdiff(mech_cols, names(df))) {
      df[[col]] <- NA
    }
    
    df[1, mech_cols, drop = FALSE]
  }
))

rownames(mech_df) <- NULL


# ============================================
# 5. MATCH MECHANISMS TO COMPOUNDS
# ============================================

# First match using both compound name and CID
match_idx <- match(
  paste(json_fda_comp$name, json_fda_comp$cid),
  paste(json_fda_mech$query_field, json_fda_mech$cid)
)

# Fallback: match by CID only where query_field is missing
unmatched <- which(is.na(match_idx))
missing_query <- which(is.na(json_fda_mech$query_field))

match_idx[unmatched] <- missing_query[
  match(
    json_fda_comp$cid[unmatched],
    json_fda_mech$cid[missing_query]
  )
]

cat("Matched compounds:", sum(!is.na(match_idx)), "\n")
cat("Unmatched compounds:", sum(is.na(match_idx)), "\n")

# Rename overlapping mechanism columns
overlapping_cols <- intersect(
  names(json_fda_comp),
  names(mech_df)
)

names(mech_df)[names(mech_df) %in% overlapping_cols] <-
  paste0("mechanism_", overlapping_cols)

final_mech_cols <- names(mech_df)

# Attach first mechanism to each compound
json_fda_comp <- cbind(
  json_fda_comp,
  mech_df[match_idx, , drop = FALSE]
)

rownames(json_fda_comp) <- NULL

# Standardize treatment identifiers
json_fda_comp$treatmentid <- as.character(json_fda_comp$cid)

# ============================================
# 7. BASIC VERIFICATION
# ============================================

cat("\n========== VERIFICATION ==========\n")

cat(
  "FDA-approved compounds:",
  nrow(json_fda_comp_original),
  "\n"
)

cat(
  "Unique compounds queried:",
  length(compounds),
  "\n"
)

cat(
  "Mechanism records returned:",
  nrow(json_fda_mech),
  "\n"
)

cat(
  "Compounds with mechanisms:",
  length(nonempty),
  "\n"
)

cat(
  "Final dataframe rows:",
  nrow(json_fda_comp),
  "\n"
)

cat(
  "Compounds without matching API records:",
  sum(is.na(match_idx)),
  "\n"
)

cat(
  "Duplicate returned compound names:",
  sum(duplicated(json_fda_mech$query_field)),
  "\n"
)

stopifnot(
  nrow(json_fda_comp) == nrow(json_fda_comp_original),
  !anyDuplicated(names(json_fda_comp)),
  all(json_fda_comp$fda_approval == TRUE),
  identical(
    json_fda_comp$name,
    json_fda_comp_original$name
  )
)

cat("\nBasic verification passed.\n")


# ============================================
# 8. TEST FIRST MECHANISM MATCHING
# ============================================

mismatches <- integer(0)
tested <- 0

for (i in seq_len(nrow(json_fda_comp))) {
  
  api_idx <- match_idx[i]
  
  # Skip compounds with no matching API record
  if (is.na(api_idx)) {
    next
  }
  
  mechanisms <- json_fda_mech$mechanisms[[api_idx]]
  
  # Compounds without mechanisms must have NA values
  if (is.null(mechanisms) || nrow(mechanisms) == 0) {
    
    if (!all(is.na(
      json_fda_comp[i, final_mech_cols, drop = FALSE]
    ))) {
      mismatches <- c(mismatches, i)
    }
    
    next
  }
  
  tested <- tested + 1
  
  # Compare the first mechanism entry
  expected <- mechanisms[1, , drop = FALSE]
  
  for (col in names(expected)) {
    
    final_col <- if (col %in% overlapping_cols) {
      paste0("mechanism_", col)
    } else {
      col
    }
    
    actual_value <- json_fda_comp[[final_col]][i]
    expected_value <- expected[[col]][1]
    
    if (!(
      (is.na(actual_value) && is.na(expected_value)) ||
      isTRUE(all.equal(
        actual_value,
        expected_value,
        check.attributes = FALSE
      ))
    )) {
      
      mismatches <- c(mismatches, i)
      break
    }
  }
}

mismatches <- unique(mismatches)

test_first_mechanism <- length(mismatches) == 0

cat("\n========== FIRST MECHANISM TEST ==========\n")

cat("Compounds tested:", tested, "\n")
cat("Mismatches:", length(mismatches), "\n")
cat("Test passed:", test_first_mechanism, "\n")

if (!test_first_mechanism) {
  
  print(
    head(
      json_fda_comp[mismatches, c("name", "cid"), drop = FALSE],
      20
    )
  )
  
  stop("First mechanism validation failed.")
}


# ============================================
# 9. SAVE FINAL DATASET
# ============================================

saveRDS(
  json_fda_comp,
  "../outputs/fda_compounds.rds"
)

cat("\nAll tests passed. Dataset saved successfully.\n")