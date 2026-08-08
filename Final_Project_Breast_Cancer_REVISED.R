# Introduction:
# The major difference in this revision of the project is that although I ask the same question as before, namely if the clinical data on cBioPortal can predict the overall survival of the patient, the way the data is being obtained is completely different from the way I initially asked the question. 
# There were three things that turned out to be true, but that I had not expected before starting this project. 
# First, the endpoint listing of studies does not actually include all studies by default. 
# Second, the clinical endpoint does not output a patient-by-variable table, but a long table instead. 
# Third, the two encodings for missing values are not the same in different studies.
# In this report, I first gather all the public studies about breast cancer and describe their completeness when the response has been expanded. I then create a relational database in third normal form using SQL, where the primary and foreign keys are included in the file itself, as well as a diagram of my schema. Finally, I compile the data for analysis using a SQL view based on joins and use this data to train a Random Forest.

# Part 0: Exploration
library(httr)
library(jsonlite)
library(cbioportalR)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(stringr)
library(ggplot2)
library(knitr)
library(DBI)
library(RSQLite)
library(tidymodels)
library(vip)
library(DiagrammeR)

# jsonlite is attached before purrr on purpose, because both export flatten and I want purrr's version to win.

set_cbioportal_db("public")

# available_studies() only gives back one page of results, which in my case was 94 records and did not include METABRIC, although the study is public and accessible. Because this call needs all publicly available breast cancer studies, I have to page through the studies endpoint by myself. The number of records per page is specified as an integer literal because R represents large integers as 1e+05, and the server does not accept it.
fetch_all_studies <- function(page_size = 10000L) {
  out <- list()
  page <- 0L
  repeat {
    r <- GET("https://www.cbioportal.org/api/studies",
             query = list(projection = "SUMMARY",
                          pageSize  = page_size,
                          pageNumber = page))
    if (status_code(r) != 200) {
      stop("studies endpoint returned HTTP ", status_code(r), ": ",
           content(r, "text", encoding = "UTF-8"))
    }
    d <- jsonlite::fromJSON(content(r, "text", encoding = "UTF-8"), flatten = TRUE)
    if (length(d) == 0 || NROW(d) == 0) break
    out[[length(out) + 1L]] <- as_tibble(d)
    if (NROW(d) < page_size) break
    page <- page + 1L
  }
  bind_rows(out)
}

all_studies <- fetch_all_studies()
stopifnot(nrow(all_studies) > 400)

# Two aspects of this filter have been chosen deliberately. First, matching on the description brings in pan-cancer cohorts like MSK MetTropism where breast appears just as one cancer type among several and where no breast-specific clinical variables are recorded at all; with twenty five thousand such cases, the resulting analytical dataset ended up being dominated by cohorts that had genomic information only. Matching on the cancer type tag and study name makes sure that the pool includes only breast cancer studies. Secondly, since str_detect produces NA and not FALSE for NA inputs and filter deletes the row with NA, it has been decided to enclose each text variable within coalesce.
brca_studies <- all_studies %>%
  filter(
    tolower(coalesce(cancerTypeId, "")) %in% c("brca", "breast") |
      str_detect(tolower(coalesce(name, "")), "breast")
  ) %>%
  distinct(studyId, .keep_all = TRUE)

# The two largest public breast cohorts must make the cut. If not, then the discovery phase was unsuccessful, and all numbers after that will paint an incomplete picture; therefore, I must ensure they appear on the list.
# This check is made at the level of the cohort and not the study id because TCGA breast cohort is available via five separate ids and METABRIC cohort is available via one id. Deduplication later on must be allowed the freedom of selecting the version of the cohort that has the most clinical information about the patients in it; hence, demanding the survival of one single id is counterproductive.
cohort_of <- function(sid) {
  dplyr::case_when(
    str_detect(sid, "^brca_tcga") ~ "TCGA",
    str_detect(sid, "metabric")   ~ "METABRIC",
    TRUE                          ~ sid
  )
}

must_have_cohorts <- c("TCGA", "METABRIC")
stopifnot(all(must_have_cohorts %in% cohort_of(brca_studies$studyId)))

# I retrieve the full list of participants of each study once and keep it, since it is used twice: the number of the list is the denominator in all completeness percentages, and the list is the scaffold onto which I expand my clinical data table. A study of 2000 participants where there are 1200 entries for age is 60 percent complete for the age field, but that is only possible when those 800 participants without an age entry are still around as rows. A failed retrieval is logged instead of ignored, since it is just how I lost METABRIC once.
patient_roster <- map_dfr(brca_studies$studyId, function(sid) {
  pts <- tryCatch(
    available_patients(study_id = sid),
    error = function(e) {
      warning("available_patients failed for ", sid, ": ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )
  if (is.null(pts) || nrow(pts) == 0) return(NULL)
  tibble(study_id = sid, patient_id = as.character(pts$patientId))
})

patient_counts <- patient_roster %>% count(study_id, name = "n_patients")

# Cohorts smaller than fifty patients are discarded. One of the cohorts in the breast cancer studies contains only one patient that can neither contribute to a stratified split nor a cross-validation fold but will only cause noise in the combined model.
study_table <- brca_studies %>%
  select(study_id = studyId, study_name = name) %>%
  inner_join(patient_counts, by = "study_id") %>%
  filter(n_patients >= 50) %>%
  arrange(desc(n_patients))

patient_roster <- patient_roster %>% semi_join(study_table, by = "study_id")

stopifnot(all(must_have_cohorts %in% cohort_of(study_table$study_id)))

kable(study_table, caption = "Table 1. Candidate breast cancer studies, before duplicate patients are resolved")

# Here comes the actual clinical data. The form of the response is what matters: get_clinical_by_study yields one row for each entity-attribute combination, where the variable name is stored within clinicalAttributeId, and the measurement within value, with both patient-level and sample-level rows combined in a single object, marked by dataLevel. In this respect, AGE_AT_DIAGNOSIS and OS_STATUS are values of this table rather than its columns.
# That is why the completeness computation becomes irrelevant. summarise(across(everything(), ~ sum(!is.na(.)))) run for this data will assess the completeness of its own seven columns in the long form, and the values for them are constructed, so each one will get 100 percent.
# There is also a second and slightly more subtle way to fall into the same trap. The publications do not have a consensus on how a missing observation should be presented. METABRIC just omits the row, while TCGA gives a row with a placeholder value like [Not Available]. The fact that a placeholder ends up in the wide format means that a variable will look complete even though it is not, so I filter out these rows upfront.
MISSING_CODES <- c("", "NA", "N/A", "NaN", "Unknown", "unknown", "UNKNOWN",
                   "[Not Available]", "[Unknown]", "[Not Applicable]",
                   "[Not Evaluated]", "[Discrepancy]", "[Pending]",
                   "[Completed]", "[Redacted]", "[Missing]")

fetch_clinical <- function(sid) {
  dat <- tryCatch(
    get_clinical_by_study(study_id = sid),
    error = function(e) {
      warning("get_clinical_by_study failed for ", sid, ": ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )
  if (is.null(dat) || nrow(dat) == 0) return(NULL)
  # cbioportalR has used both names for the attribute column across versions, so I accept either.
  if (!"clinicalAttributeId" %in% names(dat) && "attrId" %in% names(dat)) {
    dat <- rename(dat, clinicalAttributeId = attrId)
  }
  if (!"studyId"   %in% names(dat)) dat$studyId   <- sid
  if (!"sampleId"  %in% names(dat)) dat$sampleId  <- NA_character_
  if (!"dataLevel" %in% names(dat)) dat$dataLevel <- ifelse(is.na(dat$sampleId), "PATIENT", "SAMPLE")
  dat %>%
    transmute(
      study_id   = as.character(studyId),
      patient_id = as.character(patientId),
      sample_id  = as.character(sampleId),
      attr_id    = as.character(clinicalAttributeId),
      value      = str_trim(as.character(value)),
      data_level = toupper(as.character(dataLevel))
    ) %>%
    filter(!is.na(value), !value %in% MISSING_CODES)
}

# The download takes a while across this many studies, so I cache it. Delete the rds file to force a refresh.
cache_file <- "clinical_long.rds"
if (file.exists(cache_file)) {
  clinical_long <- readRDS(cache_file)
  clinical_long <- clinical_long %>% semi_join(study_table, by = "study_id")
} else {
  clinical_long <- map_dfr(study_table$study_id, fetch_clinical)
  saveRDS(clinical_long, cache_file)
}

# If my download attempt is unsuccessful, I would prefer to quit right there rather than proceed using a dataset which is not the one that I am supposed to be working with, hence I verify the download rather than substitute it.
stopifnot(nrow(clinical_long) > 0)
stopifnot(all(must_have_cohorts %in% cohort_of(unique(clinical_long$study_id))))


# There is considerable overlap among the public studies, and this cannot be determined from the study IDs. The breast cancer data from TCGA is included in five separate studies, while the MSK cohorts include each other in a hierarchical fashion, meaning that the same patient ID may appear under multiple study IDs.  
# If left untreated, it creates two kinds of problems. Firstly, it adds to the total patient count for each cohort, and since the surrogate key is constructed by concatenating the study with the patient ID, it results in multiple keys being assigned to one physical patient, who may then end up in both training and testing sets, making the AUC in-sample statistic disguised as out-of-sample.
# I eliminate any duplication by including only one copy of each patient, using the study containing the maximal amount of clinical information about this patient.
duplicate_audit <- clinical_long %>%
  distinct(study_id, patient_id) %>%
  count(patient_id, name = "n_studies") %>%
  filter(n_studies > 1)

message("Patients appearing in more than one study: ", nrow(duplicate_audit),
        " (up to ", max(c(0, duplicate_audit$n_studies)), " copies each)")


# The ranking is done based on the number of clinical characteristics recorded for the study, which means that if a patient occurs in multiple studies, the survivor will be the one with the most characteristics. Nothing is nailed down: the TCGA breast data is actually five studies, with the attribute sets being different from each other by more than an order of magnitude, and picking any specific one will result in the one I mentioned winning over the one with the most information about the patients.
study_rank <- clinical_long %>%
  filter(data_level == "PATIENT") %>%
  distinct(study_id, attr_id) %>%
  count(study_id, name = "n_attrs") %>%
  mutate(cohort = cohort_of(study_id)) %>%
  arrange(desc(n_attrs)) %>%
  mutate(rank = row_number())

print(study_rank, n = 30)

primary_record <- clinical_long %>%
  distinct(study_id, patient_id) %>%
  left_join(study_rank, by = "study_id") %>%
  group_by(patient_id) %>%
  slice_min(rank, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(study_id, patient_id)

clinical_long  <- clinical_long  %>% semi_join(primary_record, by = c("study_id", "patient_id"))
patient_roster <- patient_roster %>% semi_join(primary_record, by = c("study_id", "patient_id"))

# The counts have to be rebuilt from the deduplicated roster, otherwise every completeness percentage would still be divided by the pre-deduplication cohort size.
patient_counts <- patient_roster %>% count(study_id, name = "n_patients")

study_table <- study_table %>%
  select(study_id, study_name) %>%
  inner_join(patient_counts, by = "study_id") %>%
  filter(n_patients >= 50) %>%
  arrange(desc(n_patients))

patient_roster <- patient_roster %>% semi_join(study_table, by = "study_id")
clinical_long  <- clinical_long  %>% semi_join(study_table, by = "study_id")

stopifnot(all(must_have_cohorts %in% cohort_of(study_table$study_id)))
stopifnot(n_distinct(clinical_long$patient_id) == nrow(distinct(clinical_long, study_id, patient_id)))

# Which version of each multiply published cohort survived, so the choice is recorded in the report rather than left implicit in the ranking.
kable(study_table %>% mutate(cohort = cohort_of(study_id)) %>%
        filter(cohort %in% must_have_cohorts) %>%
        select(cohort, study_id, study_name, n_patients),
      caption = "Table 2b. Version of each reference cohort retained after deduplication")

kable(study_table, caption = "Table 2. Studies used after removing cross-study duplicate patients")

# Variable completeness, done two ways that should agree.
# The first is the direct one the question asks for: for each study I pivot the long table into a wide patient by attribute table, join it back onto the full patient roster so that a patient holding no record for an attribute really appears as NA, and count from there.
wide_by_study <- clinical_long %>%
  filter(data_level == "PATIENT") %>%
  group_by(study_id) %>%
  group_split() %>%
  set_names(map_chr(., ~ .x$study_id[1])) %>%
  map(function(d) {
    sid <- d$study_id[1]
    roster <- patient_roster %>% filter(study_id == sid) %>% distinct(patient_id)
    d %>%
      distinct(patient_id, attr_id, .keep_all = TRUE) %>%
      pivot_wider(id_cols = patient_id, names_from = attr_id, values_from = value) %>%
      right_join(roster, by = "patient_id") %>%
      mutate(study_id = sid, .before = 1)
  })

completeness_wide <- imap_dfr(wide_by_study, function(w, sid) {
  n_pat <- nrow(w)
  w %>%
    summarise(across(-c(study_id, patient_id), ~ sum(!is.na(.x)))) %>%
    pivot_longer(everything(), names_to = "attr_id", values_to = "n_nonmissing") %>%
    mutate(study_id = sid,
           n_patients = n_pat,
           completeness_pct = 100 * n_nonmissing / n_patients)
})

# The second measure, for each characteristic, tallies the number of unique patients within the study who have an entry for that characteristic, in proportion to the size of the study. This second measure never constructs the wide table but will be used to verify the first measure using SQL.
completeness_long <- clinical_long %>%
  filter(data_level == "PATIENT") %>%
  group_by(study_id, attr_id) %>%
  summarise(n_nonmissing = n_distinct(patient_id), .groups = "drop") %>%
  left_join(study_table %>% select(study_id, n_patients), by = "study_id") %>%
  mutate(completeness_pct = 100 * n_nonmissing / n_patients)

# A guard against the problem that this entire section is trying to prevent. If all the percentages turn out to be 100, then either the pivot function failed to do what I think it should have or there’s still some dummy string being treated as an actual number.
summary(completeness_wide$completeness_pct)
mean(completeness_wide$completeness_pct >= 99.999)
stopifnot(mean(completeness_wide$completeness_pct >= 99.999) < 0.9)

completeness_table <- completeness_wide %>%
  select(study_id, attr_id, n_nonmissing, n_patients, completeness_pct) %>%
  arrange(study_id, desc(completeness_pct))

kable(head(completeness_table, 40), digits = 1,
      caption = "Table 3. Completeness of patient level clinical attributes, by study")

# This is what determines whether or not an attribute is eligible for use as a predictor based on its commonality to either study.
n_studies <- n_distinct(completeness_wide$study_id)

attr_coverage <- completeness_wide %>%
  group_by(attr_id) %>%
  summarise(n_studies_present = n_distinct(study_id),
            mean_completeness = mean(completeness_pct),
            .groups = "drop") %>%
  mutate(scope = case_when(
    n_studies_present == n_studies ~ "shared by all studies",
    n_studies_present > 1          ~ "partially shared",
    TRUE                           ~ "study-specific"
  )) %>%
  arrange(desc(n_studies_present), desc(mean_completeness))

kable(head(attr_coverage, 30), digits = 1,
      caption = "Table 4. Attribute availability landscape across studies")

# The outcomes also have their own table, since having an outcome attribute within a study does not necessarily mean that the attribute is filled out.
outcome_pattern <- "^(OS_STATUS|OS_MONTHS|DFS_STATUS|DFS_MONTHS|DSS_STATUS|DSS_MONTHS|PFS_STATUS|PFS_MONTHS|RFS_STATUS|RFS_MONTHS|VITAL_STATUS)$"

outcome_table <- completeness_wide %>%
  filter(str_detect(attr_id, outcome_pattern)) %>%
  select(study_id, attr_id, completeness_pct) %>%
  pivot_wider(names_from = attr_id, values_from = completeness_pct)

kable(outcome_table, digits = 1,
      caption = "Table 5. Outcome attributes by study, as percentage of patients populated. A blank cell means the study does not record that outcome at all")

plot_attrs <- attr_coverage %>%
  filter(n_studies_present >= 2) %>%
  slice_max(mean_completeness, n = 25) %>%
  pull(attr_id)

completeness_wide %>%
  filter(attr_id %in% plot_attrs) %>%
  ggplot(aes(x = study_id, y = reorder(attr_id, completeness_pct), fill = completeness_pct)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(completeness_pct)), size = 2.2) +
  scale_fill_gradient(low = "#ff9999", high = "#66cc66", limits = c(0, 100), name = "Completeness (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Clinical variable completeness",
       subtitle = "Measured after pivoting the long API response into a wide patient by attribute table",
       x = "Study ID", y = "Clinical variable")

# Data completeness usually has a bimodal distribution. While the median of percentage completeness is 99.5 percent, the smallest percentage completeness is 0.09 percent, owing to the fact that whereas the general fields such as age and survival status are almost completely covered, the specialized fields are captured only once.

# This section is diagnostic rather than analytical. Before formulating any harmonization rules, it is important to see how each study refers to the particular concepts I am interested in, for the naming scheme differs across each database – TCGA calls PATH_T_STAGE, AURORA calls PRIM_T, and METABRIC yet again uses its own vocabulary.
concept_probe <- clinical_long %>%
  filter(data_level == "PATIENT") %>%
  distinct(study_id, attr_id) %>%
  filter(str_detect(attr_id, "AGE|STAGE|_T$|_N$|_M$|GRADE|SIZE|ER|PR|HER2|RECEPT|NODE")) %>%
  arrange(attr_id, study_id)

print(concept_probe, n = 200)

# Harmonization plan.
# It is not the case that a clinical concept will have the same identifier across all of the studies. Instead of silently remapping identifiers in R, I use the exact identifier returned by the API as my identifier and store the harmonization mapping as a column in the same table in the database.
# The following table includes all the harmonization mappings suggested by the above diagnostic code for the studies I have worked with. If an identifier is not in the table, it remains its own harmonization identifier and thus remains as a study-specific attribute in the database, never used as a predictor.
# Extend this table after reading the concept_probe function, then recreate the database, since harmonized_id is a stored column.
harmonization_map <- tribble(
  ~raw_attr_id,                   ~harmonized_id,
  "AGE",                          "AGE",
  "AGE_AT_DIAGNOSIS",             "AGE",
  "SEX",                          "SEX",
  "AJCC_PATHOLOGIC_TUMOR_STAGE",  "TUMOR_STAGE",
  "TUMOR_STAGE",                  "TUMOR_STAGE",
  "PRIM_STAGE_DX",                "TUMOR_STAGE",
  "PATH_T_STAGE",                 "T_STAGE",
  "PRIM_T",                       "T_STAGE",
  "PATH_N_STAGE",                 "N_STAGE",
  "PRIM_N",                       "N_STAGE",
  "PATH_M_STAGE",                 "M_STAGE",
  "PRIM_M",                       "M_STAGE",
  "GRADE",                        "GRADE",
  "NEOPLASM_HISTOLOGIC_GRADE",    "GRADE",
  "TUMOR_SIZE",                   "TUMOR_SIZE",
  "ER_STATUS",                    "ER_STATUS",
  "ER_IHC",                       "ER_STATUS",
  "ER_STATUS_BY_IHC",             "ER_STATUS",
  "PRIM_ER",                      "ER_STATUS",
  "PR_STATUS",                    "PR_STATUS",
  "PR_STATUS_BY_IHC",             "PR_STATUS",
  "PRIM_PR",                      "PR_STATUS",
  "HER2_STATUS",                  "HER2_STATUS",
  "IHC_HER2",                     "HER2_STATUS",
  "PRIM_HER2_INTERPRETATION",     "HER2_STATUS",
  "MUTATION_COUNT",               "MUTATION_COUNT",
  "FRACTION_GENOME_ALTERED",      "FRACTION_GENOME_ALTERED",
  "OS_STATUS",                    "OS_STATUS",
  "OS_MONTHS",                    "OS_MONTHS"
)

# Naming the concepts identically is only half of harmonization. The studies also encode the values differently, so a stage recorded as STAGE IIA in one place and as 2 in another would join into a single column that is not comparable across studies. I normalise the value vocabulary for the staging variables, and leave a record of what the raw values were.
stage_probe <- clinical_long %>%
  inner_join(harmonization_map, by = c("attr_id" = "raw_attr_id")) %>%
  filter(harmonized_id %in% c("TUMOR_STAGE", "T_STAGE", "N_STAGE", "M_STAGE")) %>%
  count(study_id, harmonized_id, value, sort = TRUE)

print(head(stage_probe, 60), n = 60)

# Since the biggest cohort of all in the set needs its own check, and an idea that cannot be mapped on it will affect more patients than any other group, this gives the list of all staging, grading, and receptor attributes that METABRIC has in addition to whether my map includes the attribute or not.
metabric_probe <- clinical_long %>%
  filter(study_id == "brca_metabric", data_level == "PATIENT",
         str_detect(attr_id, "STAGE|GRADE|SIZE|^ER|^PR|HER2|NODE|AGE")) %>%
  count(attr_id, name = "n_patients") %>%
  left_join(harmonization_map, by = c("attr_id" = "raw_attr_id")) %>%
  mutate(mapped = if_else(is.na(harmonized_id), "NOT MAPPED", harmonized_id)) %>%
  arrange(desc(mapped == "NOT MAPPED"), attr_id)

print(metabric_probe, n = 40)

# And the raw vocabularies of the concepts that map successfully, so the regular expressions below can be tested against actual values instead of presumed ones.
clinical_long %>%
  filter(study_id == "brca_metabric",
         attr_id %in% c("TUMOR_STAGE", "GRADE", "NEOPLASM_HISTOLOGIC_GRADE",
                        "ER_STATUS", "PR_STATUS", "HER2_STATUS", "ER_IHC")) %>%
  count(attr_id, value) %>%
  print(n = 40)

normalise_value <- function(harmonized_id, value) {
  v <- toupper(str_trim(value))
  dplyr::case_when(
    harmonized_id == "TUMOR_STAGE" ~ str_extract(str_replace(v, "^STAGE\\s*", ""), "^(IV|III|II|I|[0-4])"),
    harmonized_id == "T_STAGE"     ~ str_extract(v, "T[0-4X]"),
    harmonized_id == "N_STAGE"     ~ str_extract(v, "N[0-3X]"),
    harmonized_id == "M_STAGE"     ~ str_extract(v, "M[01X]"),
    harmonized_id %in% c("ER_STATUS", "PR_STATUS", "HER2_STATUS") ~
      dplyr::case_when(str_detect(v, "^POS|^\\+|^YES") ~ "POSITIVE",
                       str_detect(v, "^NEG|^-|^NO")    ~ "NEGATIVE",
                       TRUE                            ~ NA_character_),
    TRUE ~ value
  )
}


# Normalisation may quietly delete a variable if the pattern fails to match the vocabulary that a study really used, and the result appears exactly the same way as a completely missing variable. This will report, for each study and concept, the number of values that remain, so an expression that doesn't match anything will show up here instead of three hundred lines later as an empty column.
normalisation_check <- clinical_long %>%
  inner_join(harmonization_map, by = c("attr_id" = "raw_attr_id")) %>%
  filter(harmonized_id %in% c("TUMOR_STAGE", "T_STAGE", "N_STAGE", "M_STAGE",
                              "ER_STATUS", "PR_STATUS", "HER2_STATUS")) %>%
  mutate(kept = !is.na(normalise_value(harmonized_id, value))) %>%
  group_by(study_id, harmonized_id) %>%
  summarise(n_raw = n(), n_kept = sum(kept), pct_kept = round(100 * mean(kept), 1),
            .groups = "drop") %>%
  arrange(pct_kept)

print(normalisation_check, n = 40)

# Part 1: DB Setup
# One thing to note about keys is that patientId and sampleId are unique within each study but not between studies, so using patientId as a one column primary key would be erroneous when we have multiple cohorts together. I create a surrogate key by appending the study name to the key and maintain the natural key as a UNIQUE key.
OUTCOME_ATTRS <- c("OS_STATUS", "OS_MONTHS")

clinical_norm <- clinical_long %>%
  left_join(harmonization_map, by = c("attr_id" = "raw_attr_id")) %>%
  mutate(harmonized_id = coalesce(harmonized_id, attr_id),
         value = normalise_value(harmonized_id, value)) %>%
  filter(!is.na(value), !value %in% MISSING_CODES)

tbl_studies <- study_table %>%
  transmute(study_id, study_name, cancer_type = "Breast Cancer", n_patients)

tbl_patients <- clinical_norm %>%
  distinct(study_id, patient_id) %>%
  mutate(patient_key = paste(study_id, patient_id, sep = ":")) %>%
  select(patient_key, patient_id, study_id)

tbl_samples <- clinical_norm %>%
  filter(data_level == "SAMPLE", !is.na(sample_id)) %>%
  distinct(study_id, patient_id, sample_id) %>%
  mutate(sample_key  = paste(study_id, sample_id,  sep = ":"),
         patient_key = paste(study_id, patient_id, sep = ":")) %>%
  select(sample_key, sample_id, patient_key)

NUMERIC_ATTRS <- c("AGE", "AGE_AT_DIAGNOSIS", "OS_MONTHS", "TUMOR_SIZE",
                   "MUTATION_COUNT", "FRACTION_GENOME_ALTERED", "TMB_NONSYNONYMOUS")

tbl_attributes <- clinical_norm %>%
  group_by(attr_id) %>%
  summarise(data_level = first(data_level), .groups = "drop") %>%
  mutate(datatype = if_else(attr_id %in% NUMERIC_ATTRS, "NUMBER", "STRING"))

# The bridge table is where the study and the attribute intersect. It solves the problem of availability of the study and also holds the harmonized name, which is correct because the raw identifier could be complete in one study but incomplete in the other.
tbl_study_attributes <- clinical_norm %>%
  distinct(study_id, attr_id, harmonized_id) %>%
  left_join(completeness_long %>% select(study_id, attr_id, completeness_pct),
            by = c("study_id", "attr_id")) %>%
  transmute(study_id, attr_id, harmonized_id,
            patient_completeness_pct = round(completeness_pct, 2))

# The two fact tables contain the values. The outcome attributes are intentionally excluded from patient_clinical because they have their own table, and including the same facts twice would simply bring back the redundancy.
tbl_patient_clinical <- clinical_norm %>%
  filter(data_level == "PATIENT", !attr_id %in% OUTCOME_ATTRS) %>%
  mutate(patient_key = paste(study_id, patient_id, sep = ":")) %>%
  distinct(patient_key, attr_id, .keep_all = TRUE) %>%
  select(patient_key, attr_id, value)

tbl_sample_clinical <- clinical_norm %>%
  filter(data_level == "SAMPLE", !is.na(sample_id)) %>%
  mutate(sample_key = paste(study_id, sample_id, sep = ":")) %>%
  distinct(sample_key, attr_id, .keep_all = TRUE) %>%
  select(sample_key, attr_id, value)

tbl_outcomes <- clinical_norm %>%
  filter(attr_id %in% OUTCOME_ATTRS) %>%
  mutate(patient_key = paste(study_id, patient_id, sep = ":")) %>%
  distinct(patient_key, attr_id, .keep_all = TRUE) %>%
  select(patient_key, attr_id, value) %>%
  pivot_wider(names_from = attr_id, values_from = value) %>%
  rename_with(tolower)

if (!"os_status" %in% names(tbl_outcomes)) tbl_outcomes$os_status <- NA_character_
if (!"os_months" %in% names(tbl_outcomes)) tbl_outcomes$os_months <- NA_character_
tbl_outcomes <- tbl_outcomes %>% select(patient_key, os_status, os_months)

# After foreign key constraints are enforced, any attempt to add a child record without its parent results in failure by SQLite, thus I do this step in order to reconcile children records with their parents.
tbl_patient_clinical <- tbl_patient_clinical %>%
  semi_join(tbl_patients, by = "patient_key") %>%
  semi_join(tbl_attributes, by = "attr_id")
tbl_sample_clinical <- tbl_sample_clinical %>%
  semi_join(tbl_samples, by = "sample_key") %>%
  semi_join(tbl_attributes, by = "attr_id")
tbl_outcomes         <- tbl_outcomes %>% semi_join(tbl_patients, by = "patient_key")
tbl_study_attributes <- tbl_study_attributes %>%
  semi_join(tbl_attributes, by = "attr_id") %>%
  semi_join(tbl_studies, by = "study_id")

# Writing the schema.
# Although the function dbWriteTable is useful, it creates an unadorned “CREATE TABLE” without the primary key, foreign key and NOT NULL, hence there are no keys in the database regardless of how much the report says about the design. I create the DDL myself, and then attach the frames to the tables that it created.
db_file <- "brca_clinical.sqlite"
if (exists("con") && inherits(con, "SQLiteConnection") && dbIsValid(con)) dbDisconnect(con)
if (file.exists(db_file)) file.remove(db_file)
con <- dbConnect(RSQLite::SQLite(), db_file)
dbExecute(con, "PRAGMA foreign_keys = ON;")

ddl <- c(
"CREATE TABLE studies (
   study_id    TEXT    NOT NULL PRIMARY KEY,
   study_name  TEXT    NOT NULL,
   cancer_type TEXT    NOT NULL,
   n_patients  INTEGER NOT NULL CHECK (n_patients >= 0)
 );",

"CREATE TABLE clinical_attributes (
   attr_id    TEXT NOT NULL PRIMARY KEY,
   data_level TEXT NOT NULL CHECK (data_level IN ('PATIENT','SAMPLE')),
   datatype   TEXT NOT NULL CHECK (datatype  IN ('NUMBER','STRING'))
 );",

"CREATE TABLE study_attributes (
   study_id                 TEXT NOT NULL,
   attr_id                  TEXT NOT NULL,
   harmonized_id            TEXT NOT NULL,
   patient_completeness_pct REAL,
   PRIMARY KEY (study_id, attr_id),
   FOREIGN KEY (study_id) REFERENCES studies(study_id)            ON DELETE CASCADE,
   FOREIGN KEY (attr_id)  REFERENCES clinical_attributes(attr_id) ON DELETE CASCADE
 );",

"CREATE TABLE patients (
   patient_key TEXT NOT NULL PRIMARY KEY,
   patient_id  TEXT NOT NULL,
   study_id    TEXT NOT NULL,
   UNIQUE (study_id, patient_id),
   FOREIGN KEY (study_id) REFERENCES studies(study_id) ON DELETE CASCADE
 );",

"CREATE TABLE samples (
   sample_key  TEXT NOT NULL PRIMARY KEY,
   sample_id   TEXT NOT NULL,
   patient_key TEXT NOT NULL,
   FOREIGN KEY (patient_key) REFERENCES patients(patient_key) ON DELETE CASCADE
 );",

"CREATE TABLE patient_clinical (
   patient_key TEXT NOT NULL,
   attr_id     TEXT NOT NULL,
   value       TEXT,
   PRIMARY KEY (patient_key, attr_id),
   FOREIGN KEY (patient_key) REFERENCES patients(patient_key)        ON DELETE CASCADE,
   FOREIGN KEY (attr_id)     REFERENCES clinical_attributes(attr_id) ON DELETE CASCADE
 );",

"CREATE TABLE sample_clinical (
   sample_key TEXT NOT NULL,
   attr_id    TEXT NOT NULL,
   value      TEXT,
   PRIMARY KEY (sample_key, attr_id),
   FOREIGN KEY (sample_key) REFERENCES samples(sample_key)          ON DELETE CASCADE,
   FOREIGN KEY (attr_id)    REFERENCES clinical_attributes(attr_id) ON DELETE CASCADE
 );",

"CREATE TABLE outcomes (
   patient_key TEXT NOT NULL PRIMARY KEY,
   os_status   TEXT,
   os_months   TEXT,
   FOREIGN KEY (patient_key) REFERENCES patients(patient_key) ON DELETE CASCADE
 );",

"CREATE INDEX idx_pc_attr   ON patient_clinical(attr_id);",
"CREATE INDEX idx_sc_attr   ON sample_clinical(attr_id);",
"CREATE INDEX idx_pat_study ON patients(study_id);"
)

invisible(lapply(ddl, function(s) dbExecute(con, s)))

# Parents go in before children, otherwise the foreign key check rejects the insert.
dbAppendTable(con, "studies",             tbl_studies)
dbAppendTable(con, "clinical_attributes", tbl_attributes %>% select(attr_id, data_level, datatype))
dbAppendTable(con, "study_attributes",    tbl_study_attributes)
dbAppendTable(con, "patients",            tbl_patients)
dbAppendTable(con, "samples",             tbl_samples)
dbAppendTable(con, "patient_clinical",    tbl_patient_clinical)
dbAppendTable(con, "sample_clinical",     tbl_sample_clinical)
dbAppendTable(con, "outcomes",            tbl_outcomes)

# I just print the schema directly from sqlite_master so that the keys are visible in the file itself, and I perform the integrity check which shows no results if all the foreign keys are resolved.
cat("Schema as stored in brca_clinical.sqlite\n\n")
cat(dbGetQuery(con, "SELECT sql FROM sqlite_master WHERE type = 'table' ORDER BY name;")$sql,
    sep = "\n\n")
key_summary <- bind_rows(
  dbGetQuery(con, "SELECT 'patients' AS tbl UNION ALL SELECT 'samples' UNION ALL
                   SELECT 'studies' UNION ALL SELECT 'clinical_attributes' UNION ALL
                   SELECT 'study_attributes' UNION ALL SELECT 'patient_clinical' UNION ALL
                   SELECT 'sample_clinical' UNION ALL SELECT 'outcomes'") %>%
    rowwise() %>%
    mutate(
      primary_key  = paste(dbGetQuery(con, paste0("PRAGMA table_info('", tbl, "')")) %>%
                             filter(pk > 0) %>% arrange(pk) %>% pull(name), collapse = ", "),
      foreign_keys = paste(dbGetQuery(con, paste0("PRAGMA foreign_key_list('", tbl, "')")) %>%
                             mutate(fk = paste0(from, " -> ", table, "(", to, ")")) %>% pull(fk),
                           collapse = "; ")
    ) %>%
    ungroup()
)

kable(key_summary, caption = "Table 6b. Primary and foreign keys as declared in the database file")

stopifnot(nrow(dbGetQuery(con, "PRAGMA foreign_key_check;")) == 0)

# Database diagram.
grViz("
digraph ER {
  graph [rankdir = LR, fontname = Helvetica, splines = ortho]
  node  [shape = plaintext, fontname = Helvetica, fontsize = 10]

  studies [label=<
    <TABLE BORDER='0' CELLBORDER='1' CELLSPACING='0'>
    <TR><TD BGCOLOR='#dbe9f6'><B>studies</B></TD></TR>
    <TR><TD ALIGN='LEFT'>PK study_id</TD></TR>
    <TR><TD ALIGN='LEFT'>study_name</TD></TR>
    <TR><TD ALIGN='LEFT'>cancer_type</TD></TR>
    <TR><TD ALIGN='LEFT'>n_patients</TD></TR>
    </TABLE>>]

  clinical_attributes [label=<
    <TABLE BORDER='0' CELLBORDER='1' CELLSPACING='0'>
    <TR><TD BGCOLOR='#dbe9f6'><B>clinical_attributes</B></TD></TR>
    <TR><TD ALIGN='LEFT'>PK attr_id</TD></TR>
    <TR><TD ALIGN='LEFT'>data_level</TD></TR>
    <TR><TD ALIGN='LEFT'>datatype</TD></TR>
    </TABLE>>]

  study_attributes [label=<
    <TABLE BORDER='0' CELLBORDER='1' CELLSPACING='0'>
    <TR><TD BGCOLOR='#e8f6db'><B>study_attributes</B></TD></TR>
    <TR><TD ALIGN='LEFT'>PK, FK study_id</TD></TR>
    <TR><TD ALIGN='LEFT'>PK, FK attr_id</TD></TR>
    <TR><TD ALIGN='LEFT'>harmonized_id</TD></TR>
    <TR><TD ALIGN='LEFT'>patient_completeness_pct</TD></TR>
    </TABLE>>]

  patients [label=<
    <TABLE BORDER='0' CELLBORDER='1' CELLSPACING='0'>
    <TR><TD BGCOLOR='#dbe9f6'><B>patients</B></TD></TR>
    <TR><TD ALIGN='LEFT'>PK patient_key</TD></TR>
    <TR><TD ALIGN='LEFT'>patient_id</TD></TR>
    <TR><TD ALIGN='LEFT'>FK study_id</TD></TR>
    </TABLE>>]

  samples [label=<
    <TABLE BORDER='0' CELLBORDER='1' CELLSPACING='0'>
    <TR><TD BGCOLOR='#dbe9f6'><B>samples</B></TD></TR>
    <TR><TD ALIGN='LEFT'>PK sample_key</TD></TR>
    <TR><TD ALIGN='LEFT'>sample_id</TD></TR>
    <TR><TD ALIGN='LEFT'>FK patient_key</TD></TR>
    </TABLE>>]

  patient_clinical [label=<
    <TABLE BORDER='0' CELLBORDER='1' CELLSPACING='0'>
    <TR><TD BGCOLOR='#e8f6db'><B>patient_clinical</B></TD></TR>
    <TR><TD ALIGN='LEFT'>PK, FK patient_key</TD></TR>
    <TR><TD ALIGN='LEFT'>PK, FK attr_id</TD></TR>
    <TR><TD ALIGN='LEFT'>value</TD></TR>
    </TABLE>>]

  sample_clinical [label=<
    <TABLE BORDER='0' CELLBORDER='1' CELLSPACING='0'>
    <TR><TD BGCOLOR='#e8f6db'><B>sample_clinical</B></TD></TR>
    <TR><TD ALIGN='LEFT'>PK, FK sample_key</TD></TR>
    <TR><TD ALIGN='LEFT'>PK, FK attr_id</TD></TR>
    <TR><TD ALIGN='LEFT'>value</TD></TR>
    </TABLE>>]

  outcomes [label=<
    <TABLE BORDER='0' CELLBORDER='1' CELLSPACING='0'>
    <TR><TD BGCOLOR='#f6e3db'><B>outcomes</B></TD></TR>
    <TR><TD ALIGN='LEFT'>PK, FK patient_key</TD></TR>
    <TR><TD ALIGN='LEFT'>os_status</TD></TR>
    <TR><TD ALIGN='LEFT'>os_months</TD></TR>
    </TABLE>>]

  studies             -> patients         [label='1 : N']
  studies             -> study_attributes [label='1 : N']
  clinical_attributes -> study_attributes [label='1 : N']
  patients            -> samples          [label='1 : N']
  patients            -> patient_clinical [label='1 : N']
  patients            -> outcomes         [label='1 : 1']
  samples             -> sample_clinical  [label='1 : N']
  clinical_attributes -> patient_clinical [label='1 : N']
  clinical_attributes -> sample_clinical  [label='1 : N']
}
")

# Why this schema & why it is in 3NF.
# For the first normal form, all columns consist of one atomic unit. Using a format that represents the measurements as entity, attribute, and value rows will also solve the issue of repeating groups for a single wide table containing the clinical data because the number of different attributes for the various studies is unknown and would require every possible column in the table.
# The only composite keys needed for the second normal form are study_id and attr_id for study_attributes, patient_key and attr_id for patient_clinical, and sample_key and attr_id for sample_clinical. In all cases, the other columns depend on the entire composite key. A value is a characteristic of a specific entity and a specific attribute, while a standardized name or completeness score is a characteristic of a specific study and a specific attribute because the raw identifiers may be standardized or incomplete in another cohort.
# Third normal form requires that no non-key column be functionally determined by another non-key column. Study name, type of cancer, and number of patients are dependent upon the study identifier and thus are stored in studies and not duplicated in each patient row, while data level and datatype are dependent upon the attr_id and thus are stored in clinical_attributes. This is where my prior design was incorrect, as the type of cancer was included in patients when in fact it was determined by the study and not by the patient, and was thus a transitive dependency from patient_key to study_id to cancer_type.
# Joins are necessary since no table alone contains a row of the modeling dataset. Age and stage are rows in patient_clinical, mutation count is a row in sample_clinical and can only be accessed through samples, survival is in outcomes, and the study label is in studies. Putting together one row of the modeling dataset requires joins across five tables. Below is a view that does just that.
PATIENT_PREDICTOR_IDS <- c("AGE", "SEX", "TUMOR_STAGE", "T_STAGE", "N_STAGE",
                           "M_STAGE", "GRADE", "TUMOR_SIZE",
                           "ER_STATUS", "PR_STATUS", "HER2_STATUS")

# The CASE WHEN block is generated rather than typed out, so that adding a concept to the harmonization map adds a column here automatically instead of requiring the SQL to be edited in a second place.
case_lines <- paste0(
  "    MAX(CASE WHEN sa.harmonized_id = '", PATIENT_PREDICTOR_IDS,
  "' THEN pc.value END) AS ", tolower(PATIENT_PREDICTOR_IDS),
  collapse = ",\n"
)

sql_analytic <- paste0(
"CREATE VIEW v_analytic AS
WITH pat_wide AS (
  SELECT
    pc.patient_key,\n", case_lines, "
  FROM patient_clinical pc
  JOIN patients p          ON p.patient_key = pc.patient_key
  JOIN study_attributes sa ON sa.attr_id = pc.attr_id AND sa.study_id = p.study_id
  GROUP BY pc.patient_key
),
samp_agg AS (
  SELECT
    s.patient_key,
    AVG(CASE WHEN sc.attr_id = 'MUTATION_COUNT'          THEN CAST(sc.value AS REAL) END) AS mutation_count,
    AVG(CASE WHEN sc.attr_id = 'FRACTION_GENOME_ALTERED' THEN CAST(sc.value AS REAL) END) AS fga
  FROM samples s
  JOIN sample_clinical sc ON sc.sample_key = s.sample_key
  GROUP BY s.patient_key
)
SELECT
  p.patient_key,
  p.study_id,
  st.cancer_type,
  pw.*,
  sg.mutation_count,
  sg.fga,
  o.os_status,
  CAST(o.os_months AS REAL) AS os_months
FROM patients p
JOIN studies  st ON st.study_id   = p.study_id
JOIN outcomes o  ON o.patient_key = p.patient_key
LEFT JOIN pat_wide pw ON pw.patient_key = p.patient_key
LEFT JOIN samp_agg sg ON sg.patient_key = p.patient_key
WHERE o.os_status IS NOT NULL AND o.os_status <> '';")

dbExecute(con, sql_analytic)

# The number of mutations is averaged across the patient’s samples, not from one sample, as patients may have multiple samples sequenced, and selecting any particular one will cause the result to depend on row order. That’s also the reason why sample level values are aggregated within the query itself, and not appended to the patient table as an afterthought, which is the problem that led to the error in my previous attempt.
# As the final step in Part 0 I calculate the percentages of completeness again, but this time using SQLite. When the numbers coincide with those in R, that means that the data in the database matches the exploration.
sql_completeness <- "
SELECT sa.study_id,
       sa.attr_id,
       sa.harmonized_id,
       COUNT(DISTINCT pc.patient_key) AS n_nonmissing,
       s.n_patients,
       ROUND(100.0 * COUNT(DISTINCT pc.patient_key) / s.n_patients, 1) AS completeness_pct
FROM study_attributes sa
JOIN studies  s ON s.study_id = sa.study_id
JOIN patients p ON p.study_id = sa.study_id
LEFT JOIN patient_clinical pc
       ON pc.patient_key = p.patient_key AND pc.attr_id = sa.attr_id
GROUP BY sa.study_id, sa.attr_id
ORDER BY sa.study_id, completeness_pct DESC;
"
kable(head(dbGetQuery(con, sql_completeness), 30),
      caption = "Table 7. Completeness recomputed inside SQLite from the stored tables")

# Part 2: SQL Viz
# Everything from here comes out of the database through the view, so the plotting dataset, the training data and the test data all share the same SQL provenance.
raw_analytic <- dbGetQuery(con, "SELECT * FROM v_analytic;")

# All the value columns were output as character data types since this is the way they are stored in the fact tables. For a column to be output as numeric, the conversion must not create any new missing values, hence avoiding conversion of T2 to NA.
to_numeric_if_safe <- function(x) {
  if (!is.character(x)) return(x)
  num <- suppressWarnings(as.numeric(x))
  if (mean(is.na(num)) - mean(is.na(x)) < 0.02) num else x
}

df_viz <- raw_analytic %>%
  select(-any_of("patient_key.1")) %>%
  mutate(across(everything(), to_numeric_if_safe)) %>%
  mutate(
    os_event = if_else(str_detect(toupper(os_status), "^1|DECEASED"), "1", "0"),
    os_event = factor(os_event, levels = c("0", "1")),
    study_id = factor(study_id)
  ) %>%
  filter(!is.na(os_event))

kable(df_viz %>% count(study_id, os_event) %>% pivot_wider(names_from = os_event, values_from = n),
      caption = "Table 8. Analytic cohort size and event counts by study")

if (sum(!is.na(df_viz$age)) > 0) {
  print(
    ggplot(df_viz %>% filter(!is.na(age)), aes(x = os_event, y = age, fill = os_event)) +
      geom_violin(alpha = 0.5) +
      geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA) +
      facet_wrap(~ study_id, scales = "free_y") +
      theme_minimal() +
      labs(title = "Age by survival status", x = "OS event, 1 means deceased", y = "Age")
  )
}

if (sum(!is.na(df_viz$mutation_count)) > 0) {
  print(
    ggplot(df_viz %>% filter(!is.na(mutation_count), mutation_count > 0),
           aes(x = os_event, y = mutation_count, fill = os_event)) +
      geom_violin(alpha = 0.5) +
      scale_y_log10() +
      theme_minimal() +
      labs(title = "Mutation count by survival status", x = "OS event", y = "Mutation count, log10 scale")
  )
}

if (sum(!is.na(df_viz$tumor_stage)) > 0) {
  print(
    df_viz %>%
      filter(!is.na(tumor_stage)) %>%
      count(tumor_stage, os_event) %>%
      group_by(tumor_stage) %>%
      mutate(prop = n / sum(n)) %>%
      ggplot(aes(x = tumor_stage, y = prop, fill = os_event)) +
      geom_col() +
      theme_minimal() +
      labs(title = "Proportion deceased by tumor stage", x = "Tumor stage", y = "Proportion")
  )
}

# A missingness summary for the analytic dataset, which is the one relevant for modeling, as opposed to the missingness in the long table, which is zero by definition.
df_viz %>%
  select(any_of(c(tolower(PATIENT_PREDICTOR_IDS), "mutation_count", "fga"))) %>%
  summarise(across(everything(), ~ mean(is.na(.x)) * 100)) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "pct_missing") %>%
  ggplot(aes(x = reorder(variable, pct_missing), y = pct_missing)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Missingness in the analytic dataset", x = NULL, y = "Percent missing")


# In terms of these graphs, the one that clearly differentiates between the two survivor groups is the age variable, with the dead patients being the older ones, and this holds true for the data in all the studies as well as the pooled data set.
# The receptor status is weakly associated as can be predicted by the clinical literature, whereas sex does not differentiate between the two groups as the sample comprises mostly of women.
# All the above represent unadjusted relationships with no corrections made for the studies or the follow-up period.

# Part 3 & 4: ML Modeling
# The decision on what to pick as predictors is where the availability landscape becomes relevant not just descriptively but practically as well. A variable available in one cohort but unavailable in another becomes after merging of the two not a clinical observation but an indication of the cohort origin. Imputing this variable as "Unknown" allows achieving high AUC score for the Random Forest algorithm through the cohort identification rather than by recognizing the biological signal, with the subsequent variable importance plot showing a batch effect.
# Consequently, the variables available and the cohorts kept become exactly the same choice made from two different perspectives. It is nothing other than the choice between keeping the shared variables and keeping the study-specific ones, as described in the assignment statement. Instead of setting a cut-off value, I present the trade-off in the form of a table and make a decision based on it.
candidate_vars <- intersect(tolower(PATIENT_PREDICTOR_IDS), names(df_viz))

avail <- df_viz %>%
  group_by(study_id) %>%
  summarise(across(all_of(candidate_vars), ~ mean(!is.na(.x))),
            n = n(), .groups = "drop")

kable(avail, digits = 2,
      caption = "Table 9. Within-study availability of each candidate predictor, with cohort size")

# The rankings of the variables are based on how many members of the pooled cohort have them, and then I go down that list. Requiring k number of variables means I must keep only those studies that have all k number of them, thus adding one more variable adds cost to the patients. The following is just an attempt to quantify that cost.
var_rank <- avail %>%
  summarise(across(all_of(candidate_vars), ~ sum(.x * avail$n) / sum(avail$n))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "patient_coverage") %>%
  arrange(desc(patient_coverage))

tradeoff <- map_dfr(seq_len(nrow(var_rank)), function(k) {
  vars <- var_rank$variable[seq_len(k)]
  keep <- avail %>% filter(if_all(all_of(vars), ~ .x >= 0.5))
  tibble(n_vars      = k,
         variables   = paste(vars, collapse = ", "),
         n_studies   = nrow(keep),
         n_patients  = sum(keep$n),
         information = k * sum(keep$n))
})

kable(tradeoff, caption = "Table 10. Trade-off between number of predictors and pooled cohort size")

best         <- tradeoff %>% slice_max(information, n = 1, with_ties = FALSE)
predictors   <- var_rank$variable[seq_len(best$n_vars)]
keep_studies <- avail %>% filter(if_all(all_of(predictors), ~ .x >= 0.5)) %>% pull(study_id)

message("Predictors: ", paste(predictors, collapse = ", "))
message("Studies retained: ", length(keep_studies), " covering ", best$n_patients, " patients")

df_viz <- df_viz %>% filter(study_id %in% keep_studies) %>% droplevels()

stopifnot(length(predictors) >= 2, nrow(df_viz) >= 200)

# The score optimizes variables times patients and picks the three variable row, which is the same one that I would pick for other reasons, so in this particular case the arithmetic matches the decision.
# Not shown in this table, and what in reality shaped the candidate set is that the genomic summary variables were excluded from this calculation at this point in time. 
# If they had been included, the scoring method would have picked the two variable combination of MUTATION_COUNT and sex from a larger pool of patients, with a decent score but learning only about the provenance of the patient and not anything else, as described in Part 4. 
# More rows than three cannot be used anyway: addition of stage reduces the dataset to a single study of 55 patients, because stage is recorded at patient level in TCGA but sample level in METABRIC.

# As character conversion occurs prior to the NA replacement since columns in SQLite where all the values are NULL end up being imported as logical in R, hence putting the string Unknown into a logical variable is not possible. The last select statement removes any constant variables after applying the cohort filter.
df_model <- df_viz %>%
  select(patient_key, os_event, all_of(predictors)) %>%
  mutate(across(where(is.character) | where(is.logical),
                ~ factor(replace_na(as.character(.x), "Unknown")))) %>%
  select(any_of(c("patient_key", "os_event")) | where(~ n_distinct(.x, na.rm = TRUE) > 1))

set.seed(42)
data_split <- initial_split(df_model, prop = 0.80, strata = os_event)
df_train <- training(data_split)
df_test  <- testing(data_split)
cv_folds <- vfold_cv(df_train, v = 5, strata = os_event)

# The algorithm uses the median to provide numeric predictors with numeric values and provides categorical predictors with a level called "Unknown" rather than deleting patients with no information for predictors, combines the rare levels so that a predictor level which occurs three times is not treated as another dummy variable.
rf_recipe <- recipe(os_event ~ ., data = df_train) %>%
  update_role(patient_key, new_role = "ID") %>%
  step_impute_median(all_numeric_predictors()) %>%
  step_unknown(all_nominal_predictors()) %>%
  step_other(all_nominal_predictors(), threshold = 0.02) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_zv(all_predictors())

# As for the scaling issue the assignment brings up: I do not need to scale since Random Forest splits on the ordering of values in a predictor variable, and any monotonic scaling does not change this ordering.
rf_spec <- rand_forest(mtry = tune(), min_n = tune(), trees = 500) %>%
  set_engine("ranger", importance = "permutation") %>%
  set_mode("classification")

rf_workflow <- workflow() %>% add_recipe(rf_recipe) %>% add_model(rf_spec)

n_pred <- ncol(bake(prep(rf_recipe), new_data = NULL)) - 2
rf_grid <- grid_regular(mtry(range = c(2, max(3, floor(sqrt(n_pred)) + 3))),
                        min_n(range = c(5, 40)),
                        levels = 4)

rf_tune_res <- tune_grid(rf_workflow, resamples = cv_folds, grid = rf_grid,
                         metrics = metric_set(roc_auc))

kable(collect_metrics(rf_tune_res) %>% arrange(desc(mean)) %>% head(10), digits = 3,
      caption = "Table 11. Cross-validated AUC for the best hyperparameter settings")

autoplot(rf_tune_res) + theme_minimal() +
  labs(title = "Five fold cross-validated AUC across the tuning grid")

best_auc    <- select_best(rf_tune_res, metric = "roc_auc")
final_rf_wf <- finalize_workflow(rf_workflow, best_auc)

# The test set is touched only here, once, through last_fit.
final_fit <- last_fit(final_rf_wf, split = data_split, metrics = metric_set(roc_auc, accuracy))
kable(collect_metrics(final_fit), digits = 3, caption = "Table 12. Held out test set performance")

test_preds <- collect_predictions(final_fit)

test_preds %>% roc_curve(truth = os_event, .pred_0) %>% autoplot() +
  labs(title = "ROC curve on the 20 percent held out test set")

test_preds %>% conf_mat(truth = os_event, estimate = .pred_class) %>%
  autoplot(type = "heatmap") + labs(title = "Confusion matrix on the test set")

vip(extract_fit_parsnip(final_fit), num_features = 12, fill = "steelblue") +
  theme_minimal() +
  labs(title = "Random Forest variable importance")

# There are two findings from this study which are methodological and not clinical and they are related to how things happen after pooling the cohorts.
# The first finding is that there are inconsistent levels of attribute assignment among different studies. 
# In the METABRIC study, the attributes such as TUMOR_STAGE, PR_STATUS and HER2_STATUS have been assigned as sample level attributes whereas in TCGA they have been classified as patient level attributes. 
# Since the dataset has been generated on the basis of patient level attributes, all the above three variables will be totally missing from 1981 patients of METABRIC.

# The second reason is that genomic summary variables cannot be combined across different studies. 
# The median MUTATION_COUNT is 341 in AURORA, 30 in TCGA and 5 in METABRIC, with two orders of magnitude between them, owing to differences in sequencing (targeted panel, whole exome and array based sequencing), and not due to tumour biology. 
# In the initial version of the model, this was the most important variable, but since the death rate in AURORA is six times that in TCGA, the order is more a consequence of the model's understanding of the data origin than biological information. 
# Genomic summary variables were excluded from the list of predictor variables for the final model.

# Similar caution is also warranted for the final model, although less so. 
# Death rate is 13.9 percent in TCGA and 57.7 percent in METABRIC, primarily because METABRIC patients are followed for over two decades, while TCGA's median follow-up is less than three years. 
# Thus, a binary OS_STATUS carries information on both survival and length of follow-up, and the test AUC of 0.742 is not a true estimate of the prediction power of the clinical covariates alone.
# Permutation importance is based on the actual decrease in prediction accuracy when the variable is randomly shuffled, while impurity importance, on the other hand, tends to prefer variables that take many different values.
# It should also be mentioned what variables, which could be considered clinically relevant, are not included into the model and why: stage and HER2 status are excluded due to level mismatch.


dbDisconnect(con)
