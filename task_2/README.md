## Task 2

*You are not expected to finish Task 2. Get as far as you can, then use the section at the end of this README to note any difficulties and how you would continue with more time.*

This task will build off some of the AnnotationDB API workflows from **Task 1** while expanding to biological curation work.

You'll notice there's a file in `raw_data/mutations.tsv` in `task_2`. You will use this file to guide the construction of your PharmacoSet sample slot and mutation experiment.

The diagram below shows how the objects you build in this task fit together. The mutation experiment is a `SummarizedExperiment`, which is placed inside a `MultiAssayExperiment`, which becomes the `molecularProfiles` slot of the PharmacoSet:

```
PharmacoSet
├── sample              data.frame, one row per cell line            ← Step 1
├── treatment           data.frame of drugs (empty in this task)
├── molecularProfiles   MultiAssayExperiment                         ← Step 3.3
│   ├── colData         the sample table from Step 1
│   └── ExperimentList  named list of experiments
│       └── mutation.exome   SummarizedExperiment                    ← Step 2
│           ├── assay        gene × cell line matrix of 0/1 values
│           ├── rowData      one row per gene (gene information)
│           ├── colData      one row per cell line
│           └── metadata     list describing the experiment
├── treatmentResponse   drug response data (empty in this task)
└── curation            identifier mapping tables                    ← Step 3.2
```

### Step 1: Build the sample slot

Extract **cell line** identifiers from `mutations.tsv` and query AnnotationDB to get the information needed to construct your sample slot. Specifically you want to extract and store name, accession, age, sex, population, disease (select any ONE description). You also want to include a boolean column labeled `In.AnnotationDB` which would contain a `TRUE` entry if the cell line metadata is found in AnnotationDB and `FALSE` otherwise. To accompany the `In.AnnotationDB` there should be an `AnnotationDB.Title` which would strictly be the title received in the AnnotationDB query. Also, have a field for `tissueid = "haematopoietic_and_lymphoid"`.

| Field | Value | Example |
|---|---|---|
| `sampleid` | Cell line identifier in the raw data | `U266_ATCC` |
| `name` | Cell line name used to query AnnotationDB (`sampleid` without the supplier suffix) | `U266` |
| `accession` | Cellosaurus accession | `CVCL_0566` |
| `age` | Age at sampling | `53Y` |
| `sex` | Sex of the cell line | `Male` |
| `population` | Donor population | `Caucasian` |
| `disease` | One disease description (chosen arbitrarily) | `Plasma cell myeloma` |
| `In.AnnotationDB` | `TRUE` if the cell line metadata is found in AnnotationDB, `FALSE` otherwise | `TRUE` |
| `AnnotationDB.Title` | Title exactly as received in the AnnotationDB query (`cell_line_name`) | `U266B1` |
| `tissueid` | Always `"haematopoietic_and_lymphoid"` | `haematopoietic_and_lymphoid` |

*If a cell line isn't found in AnnotationDB, every field returned by AnnotationDB should be `NA` and `In.AnnotationDB` is `FALSE`.*

The sample table should be its own output file (TSV or RDS).

### Step 2: Create mutation molecular profile experiment

The goal of this step is to create a `SummarizedExperiment` with variants recorded in `mutations.tsv` and attach correct cell line and gene annotations. Each row of `mutations.tsv` is one variant found by exome sequencing in one cell line, annotated with the gene it falls in (`SNPEFF_GENE_NAME`) and its predicted effect on the protein (`SNPEFF_EFFECT`, `SNPEFF_IMPACT`).

If you are used to working with tabular data: a `SummarizedExperiment` is a two-dimensional table (a matrix) with an extra table describing its rows and another describing its columns. Its row names and column names act as keys. `rowData` holds one row of information per matrix row, and `colData` holds one row per matrix column, in the same order.

**2.1 Inspect the raw data before writing any processing code**

- Identify the mapping columns you need: the cell line identifier (`cell_line`) and the gene symbol (`SNPEFF_GENE_NAME`).
- Identify fields you may need for variant consideration: `SNPEFF_EFFECT` and `SNPEFF_IMPACT`.
	- Make a choice to omit or keep variants based on these criteria. Explain your decision.

**2.2 Build the gene by cell line matrix**

- Collapse the variant table into a matrix with one row per gene symbol and one column per cell line (rows start as gene symbols and are replaced with Ensembl gene IDs in 2.3).
- Each intersection cell is `1` if there is at least one variant for that gene in that cell line, and `0` otherwise.
- Column names must be exactly the `sampleid` values from your step 1 sample table.

**2.3 Map gene symbols to Ensembl gene IDs**

Gene symbols ARE NOT stable identifiers. Gene symbols can change over time in new GENCODE iterations, so it is important to map to a more stable identifier like Ensembl gene IDs. All PharmacoSets produced in the lab will use Ensembl gene IDs instead of gene symbols for row name identifiers in tables to ensure clear and consistent mapping.

- Use the [GENCODE v50 annotation](https://storage.googleapis.com/annotationdb_pipeline_outputs/annotation_assets/gencode.v50.annotation.gtf.gz) as your gene reference. ONLY keep rows where the feature type is `gene`, and extract `gene_id`, `gene_name`, `chromosome`, `start`, `end` and `strand` from it.
- Build a lookup from gene symbol to Ensembl ID. **Keep only symbols that match exactly one ID and exclude the rest.**
- Check whether two different symbols end up mapping to the same Ensembl ID. Remember, row names must be unique, so if this happens you need a rule for which row to keep (or how to combine them).
- Write every excluded symbol to a file together with the reason it was excluded. This exclusion log is part of your submission.

**2.4 Build `rowData` and `colData`**

In 2.2 and 2.3 you were tasked with creating a matrix where the columns were the unique `sampleid`'s and the rows were the unique Ensembl gene IDs. Now we want to provide metadata for these unique identifiers for the end user. 

- `rowData` (data.frame or DataFrame): one row per Ensembl ID, in the same order as the matrix rows. Include at least the Ensembl ID, the original symbol from the raw data, the GENCODE gene name, chromosome, start, end and strand.
	- Save a copy of the gene table (TSV or RDS) as an output to submit.
- `colData` (data.frame or DataFrame): one row per cell line, in the same order as the matrix columns, with at least a `sampleid` column. Row names should equal `sampleid`.
	- *Note: the colData constructor requires a batchid. You can use NA for that in this case.*

**2.5 Create the `SummarizedExperiment`**

Now combine the matrix, `rowData` and `colData` into one object. PharmacoGx reads the experiment's metadata to know what kind of data it holds.

- Name the assay (e.g. `assays = list(value = matrix)`).
- Store at least two entries in `metadata(se)`:
  - `annotation`: the profile type (e.g. `"mutation.exome"`).
  - `value_type`: what the values mean (e.g. `"binary source-reported variant presence"`).

**2.6 Check the experiment before building the PharmacoSet**

Errors here are much easier to find now than after the PharmacoSet is built. Confirm that:

- Every value is `0` or `1`, with no `NA`.
- Row names are unique Ensembl IDs, with no version suffix and no `NA`.
- Every column name exists in your sample table's `sampleid`.
- The number of `1`s after mapping is close to the number before mapping, minus the genes you excluded.

Things to watch out for:

- A `0` means "no variant reported", not "confirmed wild type". A gene can have no rows because there is no variant, or because the region wasn't covered by the capture kit, or because the variant was filtered out upstream. Mention this in your `task_2` README.
- Cell lines not found in AnnotationDB in step 1 still have mutation data. You must still keep them in the experiment.

### Step 3: Assemble the PharmacoSet

The goal of this step is one PharmacoSet that holds the sample slot from step 1 and the mutation experiment from step 2, created with `PharmacoGx::PharmacoSet2()`.

**3.1 Prepare the sample table**

PharmacoGx uses `sampleid` as row names to link the sample slot to everything else in the PharmacoSet.

- If not already complete, convert the step 1 table to a plain `data.frame` and set its row names to `sampleid`.

**3.2 Curate the identifiers**

In PharmacoGx, the `curation` slot records how the identifiers used in this PharmacoSet relate to the identifiers in the original data. Reviewers use it to trace any sample back to the raw file.

- Use the raw identifier from `mutations.tsv` (e.g. `U266_ATCC`) as the `sampleid` throughout the PharmacoSet. It is available for every cell line, including those not found in AnnotationDB, and it is always unique. Standard names such as the Cellosaurus accession are already stored in the sample table (`accession`, `AnnotationDB.Title`).
- Build `curation$sample` as a `data.frame` with one row per sample and two columns: `sampleid` and `sampleid_dataset` (the identifier exactly as written in `mutations.tsv`). Set its row names to `sampleid`. Because you are using the raw identifier, both columns hold the same values.
- Build `curation$tissue` as an empty `data.frame`. The raw data has no tissue labels to curate; `tissueid` is set directly in the sample table.
- Build `curation$treatment` as an empty `data.frame`, since this task has no drugs.

**3.3 Build the `MultiAssayExperiment`**

The `MultiAssayExperiment` links your mutation experiment to the sample table. **It silently drops any cell line whose name doesn't match the sample table**, printing only a short "harmonizing input" message, so check dimensions before and after.

- Put the experiment in a named list (e.g. `list(mutation.exome = se)`). Users retrieve it later with `molecularProfiles(pset, "mutation.exome")`.
- Use the sample table as the `colData`.

**3.4 Create the PharmacoSet**

- Call `PharmacoGx::PharmacoSet2()` with a descriptive `name`, `sample`, `molecularProfiles`, `curation` and `datasetType = "sensitivity"`.
  - `datasetType` describes the study the data came from (a drug sensitivity screen), not what your object contains.
- `PharmacoSet2()` expects a `treatmentResponse` object even though this task has no drug data. Reuse the empty `treatmentResponse` from Task 1.

**3.5 Validate and save**

- Run `validObject(pset)` and make the script stop if it fails.
- Check the object with accessor functions, not just by printing it: `sampleInfo(pset)`, `molecularProfiles(pset, "mutation.exome")`, `mDataNames(pset)` and `curation(pset)`.
- Save the PharmacoSet with `saveRDS()`.

Expected outputs: README.md, sample table (TSV or RDS), gene table (TSV or RDS), exclusion log, PharmacoSet.

Submit script(s), envs, and add your notes to the README within the `task_2` folder.

## YOUR COMMENTS AND EXPLANATIONS HERE

### Task 2.1

This file outputs sample_metadata.rds to the outputs folder. It takes task_1/outputs/cell_lines.rds and task_2/raw_data/mutations.tsv as inputs. The script identifies all unique cell lines in mutations.tsv and then queries the AnnotationDB for the required metadata for the cell lines: age, sex, disease, etc. cell_lines.rds is used to change the cell line names used in mutations.tsv to the format expected for the AnnotationDB API call.
The script must be run from the task_2/scripts directory to correctly read the input files and output the output file.
The script is divided into 9 parts.

1. task_1/outputs/cell_lines.rds and task_2/raw_data/mutations.tsv are loaded as dataframes into the program.

2. The sample_metadata dataframe is initialized containing the unique cell lines in mutations.tsv as sampleid. The name of each cell line is also identified by removing the suffix in the sampleid. For instance, U266_ATCC has a sampleid of U266_ATCC and a name of U266. The name is what is used to query the AnnotationsDB, not the sampleid.

3. The name field in sample_metadata is matched with the cell_line_name field in clines_df. When a match is found the index of the clines_df with the appropriate cell_line_name is saved so that it can be used to query the AnnotationDB for cell_line details. First exact matching is tried, then normalized matching wherein punctuation is removed from the matching fields and they are converted to uppercase/numeric strings. The accession field for sample_metadata is also recorded from the clines_df, this is used later to match the sample_metadata rows to the data returned from the AnnotationDB (see part 5).

4. The AnnotationDB is queried at the many endpoint for cell line details containing age, sex, population and disease information for each cell line. Of the 46 unique cell lines in mutations.tsv, 44 can be found in the AnnotationDB and their information can be fetched in a single batch call.

5. The results of AnnotationDB API call are matched to the sample_metadata table via the accession field which was used as the sampleid in task 1.

6. age, sex, population and the first listed disease are added to the sample_metadata table.

7. The final fields for the sample metadata table are added, notably the fields In.AnnotationDB, AnnotationDB.Title and tissueid. The row names are also set equal to the sampleid which is meant to be the unique identifier of each row.

8. Some standardization of missing values and validation checks are applied. Empty strings are replaced with NA. "Age unspecified" and "Sex unspecified" should have also been replaced with NA and would be if given more time. The following validation checks and statistics are derived in this section:

- sample_metadata has the same length as the mutations_df derived directly from unique cell lines in mutations.tsv
- sample_metadata's sampleid is unique, complete and matches the original mutations_df cell line values
- sample_metadata has all required columns
- sample_metadata doesn't have any NA in the In.AnnotationDB field
- There are 0 duplicate accession records in sample_metadata
- No accessions are matched with multiple API records
- 2 cell lines in mutations.tsv were not found in AnnotationDB
- An inaccurate count of data missing in AnnotationDB is reported based on the untreated "Age unspecified" and "Sex unspecified" values.

9. sample_metadata is saved as sample_metadata.rds in the output folder


























