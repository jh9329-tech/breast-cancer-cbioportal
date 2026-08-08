# Breast Cancer Clinical Data Integration and Survival Modelling

In this project, I am collecting all the publicly available studies of breast cancer from 
the cBioPortal API, putting them in a normalized relational database,  and asking whether 
routinely collected covariates are predictive of overall survival when combined in a cohort.

The answer is yes, but just barely (AUC of test set 0.742). But the real answer is that 
the effort of figuring out how to pool the variables is the bulk of the work done, and 
it turned out that three of the obstacles were more informative than the model itself.

---

## Data

| | |
|---|---|
| Source | [cBioPortal](https://www.cbioportal.org) public API |
| Studies retrieved | 539, filtered to breast cancer cohorts |
| Candidate cohorts after deduplication and a ≥50 patient threshold | 20+ |
| Final analytic cohort | 3,133 patients across METABRIC, TCGA and AURORA |
| Outcome | `OS_STATUS`, binarised to deceased vs alive |

---

## Pipeline

```
cBioPortal API
      ↓  paged retrieval, placeholder-value filtering
  long-format clinical records
      ↓  cross-study patient deduplication
      ↓  attribute-name and value-vocabulary harmonisation
  SQLite database (8 tables, 3NF, declared PK/FK)
      ↓  SQL view: 5-table join + in-SQL EAV pivot
  analytic dataset
      ↓  predictor audit, 80/20 split, 5-fold CV
  tuned Random Forest
```

### Database schema

Eight tables in third normal form, with primary and foreign keys declared in the DDL rather than
left to `dbWriteTable()`, which emits no constraints at all:

`studies` · `clinical_attributes` · `study_attributes` · `patients` · `samples` ·
`patient_clinical` · `sample_clinical` · `outcomes`


Measurements are captured in entity-attribute-value rows since the studies capture different sets 
of attributes; one big table capturing all the attributes for a particular study would require 
one column for each attribute used in any given study. There is no table that holds the model row, 
and therefore the analytic dataset is put together using five tables and 
pivoting the fact rows into columns inside SQL.

---

## Three problems that shaped the analysis

None of these was anticipated at the outset, and each one changed the modelling set.

### 1. The same patient appears under several study IDs

TCGA's breast cohort is republished as five separate studies and the MSK cohorts nest inside one
another. **3,928 patients appeared in more than one study, up to six times each.**

Since the surrogate key appends the study ID to the patient ID, one single patient gets divided 
into many separate keys, and a randomly selected train/test split would end up placing duplicate patients 
on both sides of the split. Unless corrected, the AUC that has been held out becomes an in-sample value 
with an out-of-sample name tag. One copy of each patient is retained from whichever study records the most clinical attributes.

### 2. Attribute level assignment is not consistent across studies

METABRIC records `TUMOR_STAGE`, `PR_STATUS` and `HER2_STATUS` as **sample-level** attributes.
TCGA places the corresponding fields — `AJCC_PATHOLOGIC_TUMOR_STAGE` and its receptor equivalents
— at **patient level**.

In this case, since the analytic data set is constructed per-patient, all of these three variables 
seem completely absent from the 1,981 METABRIC patients in the analytic data set, despite the fact that 
this data exists in the 'sample_clinical' table. This is a case of level definition differences, 
not a case of data absence, and it does not become apparent until you look for it; the symptoms are exactly 
the same as if the variable itself didn’t exist.

METABRIC is expression array data, and has just one sample per patient. Therefore, in this particular case, 
it is a reasonable decision to collect clinical variables per-sample. This is an indicator that the problem 
has a solution: in cases when there is a one-to-one correspondence between patients and samples, 
sample-level attributes can be moved up to patient-level. This requires an explicit aggregation function in general case.

### 3. Genomic summary variables cannot be pooled across sequencing platforms

Median `MUTATION_COUNT` by cohort:

| Cohort | Median mutation count | Death rate |
|---|---:|---:|
| AURORA | 341 | 83.6% |
| TCGA | 30 | 13.9% |
| METABRIC | 5 | 57.7% |

The factor of two orders of magnitude can be attributed to targeted panel, whole-exome 
and array-based sequencing but not the biology of tumours. Together with the six-fold mortality difference 
between AURORA and TCGA, it implies that a classifier can estimate the survival by identifying the patient's group.

In the previous version of the model, the `MUTATION_COUNT` was the most important variable. 
This importance is more likely to result from patient classification than any biological effect, 
so none of the genomic summary variables were included in the prediction set.

---

## Predictor selection

Which variables can be used and which cohorts can be kept are the same decision seen from two
sides: adding a predictor forces the exclusion of every study that does not record it. The
trade-off is therefore tabulated explicitly rather than resolved by a fixed missingness threshold.

| Predictors | Studies retained | Patients |
|---|---:|---:|
| sex | 6 | 5,058 |
| sex, age | 4 | 3,203 |
| **sex, age, er_status** | **3** | **3,133** |
| + stage variables | 1 | 55 |

Rows beyond the third collapse to a single 55-patient study, for the level-mismatch reason described above.

---

## Results

| Metric | Value |
|---|---|
| Test AUC | 0.742 |
| Test accuracy | 0.693 |
| Majority-class baseline | 0.572 |

Importance is measured using permutation and not impurity, such that the order will take into account 
the actual loss in predictive accuracy upon shuffling a predictor variable, and not just its number of unique values. 
Age and ER Status contain practically all the information.

Two caveats belong with these numbers.

**Out of the three predictors, only two of the three predictors carry information.** Variable audit proves `sex` to be a predictor 
since a few males prevent it from being constant, but in a breast cancer study population,
it is almost constant and adds nothing. The model is effectively `age + er_status.`

**AUC cannot be interpreted as predictive value of covariates per se.** Mortality is 13.9% in TCGA 
and 57.7% in METABRIC, mainly due to the fact that METABRIC patients are followed for more than 
two decades while the median follow-up in TCGA is less than three years. Binary `OS_STATUS` 
thus contains information on follow-up period along with survival status, 
and any predictor correlated with belonging to a particular cohort does too.

---

## Running it

```r
# install.packages(c("cbioportalR", "tidyverse", "tidymodels", "DBI", "RSQLite",
#                    "ranger", "vip", "DiagrammeR", "httr", "jsonlite", "knitr"))

rmarkdown::render("Final_Project_Breast_Cancer_REVISED.R", output_format = "html_document")
```

First pass pulls data related to all breast cancer studies and takes about 15 minutes to complete. 
Other passes will use `clinical_long.rds` and won’t pull data anew.

Output files include `brca_clinical.sqlite` and HTML file which includes all tables, 
completeness heatmap, ER diagram, and model diagnostics.

---

## Repository contents

```
Final_Project_Breast_Cancer_REVISED.R      full pipeline
Final_Project_Breast_Cancer_REVISED.html   rendered report
README.md
```

The SQLite database is not committed, since the script rebuilds it from scratch on every run.

---

## What I would do differently

Binary encoding of `OS_STATUS` represents the most significant concession in terms of the design. 
The Cox regression model based on `OS_MONTHS` would account for all the information on follow-up 
lost due to the dichotomous variable, while remaining unconfounded by the discrepancy in the length of 
follow-up among the cohorts – a potential threat to validity of the reported finding.

It is better to solve the level mismatch issue described in section 2 than document it. 
Elevating the patient attributes available at sample level to the patient level, where possible, 
will make the stage and HER2 status information available for METABRIC’s 1,981 patients 
and increase the number of predictors significantly more than hyperparameter tuning of the current three.

---

Written in R using `cbioportalR`, `tidymodels`, `RSQLite` and `ranger`.

Mike Han
