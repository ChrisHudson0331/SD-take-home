## Task 1

1. Extract and store all FDA drugs via `Compounds` routes in [AnnotationDB](https://annotationdb.bhklab.ca/).
	- Extract **THE FIRST** mechanism entry and attach its information to the compound entry.
2. Extract all cell lines via `Cell Lines` routes in [AnnotationDB](https://annotationdb.bhklab.ca/) and the metadata for `accession`, `cell_line_name`, `age_at_sampling`, `sex_of_cell`.
3. Create an empty PharmacoSet (PSet) with [sample](https://bhklab.github.io/PharmacoGx/articles/CreatingPharmacoSet.html#sample) and [treatment](https://bhklab.github.io/PharmacoGx/articles/CreatingPharmacoSet.html#treatment) data frames using the data above and PharmacoGx.

Note: `PharmacoSet2()` will give you an error if `molecularProfiles` OR `treatmentResponse` is empty. Please use these placeholder slots upon PSet creation.

```R
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
  sample = data.frame(),     # replace with your cell line table
  treatment = data.frame(),  # replace with your FDA compound table
  molecularProfiles = placeholder_profiles,
  treatmentResponse = empty_response,
  datasetType = "sensitivity"
)
```

Expected outputs: README.md, TSV or RDS of FDA compounds, TSV or RDS of cell lines, PharmacoSet.

Submit script(s), envs, and add your notes to the README within the `task_1` folder.

## YOUR COMMENTS AND EXPLANATIONS HERE

### Task 1.1

The script task_1_1.R outputs the file fda_compounds.rds containing the first mechanism entry attached information attached to each FDA approved drug in the AnnotationDB.
The script must be run from the task_1/scripts/ folder in order to place the output file in the correct location.
The script logic is divided into 8 parts.

1. An API request is made to the AnnotationDB requesting all compound data. The resulting data is converted to a JSON object and all non FDA approved drugs are filtered out. Due to the tabular structure of the compound data, the JSON object is actually simplified to a dataframe automatically, this however should be made more clear as the name of the dataframe is json_fda_comp instead of df_fda_comp. A clarity change that was not made due to time constraints.

2. The resulting 2790 compounds need to be queried to the many endpoint to fetch their mechanism data. Each query takes about 1 second, to avoid making 2790 queries, batches of size 128 are prepared so the many endpoint returns mechanism data for 128 compounds in a single hit. This number was chosen from an experiment wherein the endpoint was hit with batch sizes in powers of 2 until the resulting url was too long, resulting in a fail. 128 was found to be the batch size that didn't fail.

3. The mechanism data is fetched for each batch and then used to construct a single JSON object containing the mechanism data. The code is more complicated than in part 1 for two reasons, first 2790/128=22 queries must be executed and the results combined into a single JSON. Second, there is some recursive logic in this section which breaks up batch sizes in half in the event the URL becomes too long for 128 compounds with unusually large names. This was found to be necessary when the script was originally written for the whole AnnotationDB dataset rather than just the FDA approved drugs. It is now extra logic which should be simplified but hasn't been due to the time constraints of this assessment.

4. The first mechanism data from the JSON object constructed in part 3 is extracted and stored in a dataframe. The rows correspond to the compounds and the columns to the first mechanism data.

5. The first mechanism data from part 4 is joined to the compound dataframe in part 1. Most of the logic is dedicated to matching the the rows of the two dataframes on the basis of compound name and CID. The CID, being a unique value for the dataset was also set to the treatmentid field which should be present and unique when constructing the PharmacoSet() object in task 1.3.

6. Some basic checks are applied to the resulting dataframe before saving it to catch any simple mistakes. The original compound table was 2790 rows, 2790 compounds were queried for mechanisms, 2811 mechanism entries were returned indicating some compounds had multiple mechanism entries, 2790-2596=194 compounds had no mechanism data, the final dataframe had 2790 rows, no compounds were missing from the mechanisms API, and 24 duplicate compounds were found in the mechanisms API going off compound name alone. These stats are largely promising but suggest some handling of duplicate compound names from the mechanism API need to be handled, the appropriate logic is applied in part 5. Lastly a few basic checks on the final dataframe will stop the program if they fail, they all pass.

7. Another verification check is applied to make sure the final dataframe is in good shape. The mechanisms entries in the final dataframe are compared to the first mechanisms found in the original JSON response from the mechanisms API. 2584 exact matches were found and 0 mismatches which is great.

8. The constructed and validated final dataframe is saved as fda_compounds.rds in the outputs folder. The script must be run from the scripts directory for this save to execute successfully.

### Task 1.2

The script task_1_2.R outputs the file cell_lines.rds containing the cell lines data and their metadata from AnnotationDB.
The script must be run from the task_1/scripts/ folder in order to place the output file in the correct location.
The script logic is divided into 6 parts and closely follows in logic to what was described in Task 1.1

1. The base cell line data is retrieved via API call from the AnnotationDB and the resulting data converted to a JSON object. Just like in Task 1.1 part 1, the resulting data is tabular in nature and so the object json_clines is actually a poorly named dataframe. Given more time this would be made more clear. 1743 cell lines are retrieved in this part.

2. As with Task 1.1 part 2, the cell lines are divided into batches of size 128. This way, cell line metadata can be fetched from the "many" AnnotationDB API endpoint for 128 cell lines at a time.

3. The metadata for each cell line is fetched in batches and the resulting data is merged together in a single dataframe. The name of the dataframe is poorly named json_clines_details. The redundant complex recursive logic in part 3 of task 1.1 is not present in this part. A total of 1743/128=14 batch api calls are made.

4. The metadata table from part 3 is combined with the cell line table from part 1 where the data is matched on the accession field which is unique. That field was also used as the sampleid for building the PharmacoSet()

5. Some basic checks on the final dataframe are done to validate the process. There were originally 1743 cell lines, 1774 cell line metadata rows, a final dataframe of 1743 rows, 72 missing ages entries, and 39 missing sex entries. The pipeline stops running if the final dataframe fails certain basic properties like having duplicate accession entries.

6. The final dataframe is saved to cell_lines.rds in the output folder. The script must be run from the scripts directory for this save to execute successfully.

### Task 1.3

The script task_1_3.R outputs the file AnnotationDB_pset.rds containing the PharmacoSet() with the tables from tasks 1.1 and 1.2.
The script must be run from the task_1/scripts/ folder in order to place the output file in the correct location and read the dataframes from tasks 1.1 and 1.2.
The script logic is divided into 3 parts

1. The dataframes containing the fda compounds and cell lines data (the outputs of task 1.1 and 1.2) are read and used to create a PharmacoSet() following the structure laid out in this README file. The treatmentid and sampleid fields from those fields are also checked to make sure they are non empty and non duplicated and then they are used to label the rows of the two dataframes.

2. Some checks are made to ensure the structure of the PharmacoSet() is sensible. The following checks are made:

   - That the PharmacoSet() is a valid S4 object
   - That the counts of the sample and treatment dataframes match those of the cell lines and fda compounds dataframes
   - That all the cell line metadata fields are present in the sample dataframe
   - That the identifiers of the two dataframes are present and unique
   - That the sample and treatment dataframes exactly match the cell lines and fda compounds dataframes
   - That all treatment dataframe entries are FDA approved
   - That the sampleid and treatmentid columns are present in the sample and treatment dataframes
   - That an empty treatment response is present

3. The resulting PharmacoSet() object is saved to the output folder as AnnotationDB_pset.rds. 
















