## CDCDB Triplets Preprocessing Workflow

This readme provides a three-step pipeline for preprocessing CDCDB drug–condition–drug triplets for downstream modeling.

### 1. Download the Latest CDCDB Data
Use **`cdcdb_file_download.ipynb`** to download the latest CDCDB SQLite file.

From the unzipped database, extract and save:
- `drug_trial.csv` — contains NCT IDs and drug combinations  
- `mesh_terms.csv` — contains NCT IDs and MeSH terms

---

### 2. Map MeSH Terms to MeSH IDs (R)
Use **R** with the UMLS **MRCONSO.RRF** mapping tool via **`cdcdb_table_merge.R`** to convert MeSH terms to their MeSH IDs.

Save the merged output as:
- `drugcombo_max2_drugs_conditions_MeSH_IDs.txt`  
  (includes DrugBank codes + `nct_id` + MeSH term + MeSH ID)

---

### 3. Clean and Filter Final CDCDB Triplets (Python)
Use **`cdcdb_final_triplets_cleaning.ipynb`** to clean the merged CDCDB table and filter out triplets in which **both drugs** and the **condition** are already present in the clinical trial dataset.

