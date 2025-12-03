# libraries
library(dplyr)

# directory to the final model results folder
result_dir = "/project/pi_rachel_melamed_uml_edu/Jianfeng/Drug_combinations/06122025/final_model/results/"
# save directory
save_dir = "/project/pi_rachel_melamed_uml_edu/Jianfeng/Drug_combinations/06122025/final_result_tables_mapping/"
# this directory is where you download your umls MRCONSO.RRF reference file
mrconso_dir = "/project/pi_rachel_melamed_uml_edu/Panos/drug_combo_jianfeng/CT_20250605/"
# drug bank vocabulary unzip csv file
db_file = "/project/pi_rachel_melamed_uml_edu/Jianfeng/Drug_combinations/06122025/final_result_tables_mapping/drugbank vocabulary.csv"
# read DrugBank csv zip file
db_csv <- read.csv(db_file, sep = ",", header = TRUE, stringsAsFactors = FALSE)

# read UMLS zip file
rrf = fread(paste0(mrconso_dir, "umls-2025AA-mrconso.zip"), sep = "|", quote = "")
rrf = rrf[, -19]
# Column names: https://www.ncbi.nlm.nih.gov/books/NBK9685/table/ch03.T.concept_names_and_sources_file_mr/?report=objectonly
colnames(rrf) = c("CUI", "LAT", "TS", "LUI", "STT", "SUI", "ISPREF", "AUI", "SAUI", "SCUI", "SDUI", "SAB", "TTY", "CODE" , "STR", "SRL", "SUPPRESS", "CVF")
# Filter for english terms only
rrf = rrf %>% 
  filter(LAT == "ENG" & SAB == "MSH")
rrf$STR_lowercase = tolower(rrf$STR)
# Keep only necessary columns in rrf
rrf_sub <- rrf[, .(CODE, STR_lowercase)]


#####-----second drug dataframe-----#####
# read the result table
second_drug_df <- fread(paste0(result_dir, "clinical_trials_with_proposed_second_drugs.csv"))
second_drug_df <- second_drug_df[, -1, with = FALSE]
second_drug_df <- as.data.table(second_drug_df)

### Map mesh code to lower case mesh term
# Join and pick only the first match if multiple exist
second_drug_df <- rrf_sub[second_drug_df, on = c("CODE" = "condition"), mult = "first"]
# Add condition_name column
second_drug_df[, condition_name := STR_lowercase]
# Rename original MeSH code column
setnames(second_drug_df, "CODE", "condition_id")
# Remove temporary column
second_drug_df[, STR_lowercase := NULL]
# Check mapping
sum(is.na(second_drug_df$condition_name))  # should be 0 if all codes matched

### Map drug bank code to drug bank name
# mapping
second_drug_df <- second_drug_df %>%
  left_join(db_csv %>% dplyr::select(DrugBank.ID, Common.name),
            by = c("drug1" = "DrugBank.ID")) %>%
  rename(drug1_name = Common.name) %>%
  left_join(db_csv %>% dplyr::select(DrugBank.ID, Common.name),
            by = c("drug2" = "DrugBank.ID")) %>%
  rename(drug2_name = Common.name)
# reorder the columns
second_drug_df <- second_drug_df %>%
  rename(drug1_id = drug1, drug2_id = drug2) %>%
  dplyr::select(drug1_id, drug1_name, condition_id, condition_name, drug2_id, drug2_name, scores)

### save the table
# write.csv(second_drug_df, file = paste0(save_dir, "clinical_trials_with_proposed_second_drugs_w_names.csv"), row.names = FALSE)


#####-----top 3 predicted triplets per condition-----#####
# read the result table
top_3_df <- fread(paste0(result_dir, "top_3_triplets_per_condition.csv"))
top_3_df <- top_3_df[, -1, with = FALSE]
top_3_df <- as.data.table(top_3_df)

### Map mesh code to lower case mesh term
# Join and pick only the first match if multiple exist
top_3_df <- rrf_sub[top_3_df, on = c("CODE" = "condition"), mult = "first"]
# Add condition_name column
top_3_df[, condition_name := STR_lowercase]
# Rename original MeSH code column
setnames(top_3_df, "CODE", "condition_id")
# Remove temporary column
top_3_df[, STR_lowercase := NULL]
# Check mapping
sum(is.na(top_3_df$condition_name))  # should be 0 if all codes matched

### Map drug bank code to drug bank name
# mapping
top_3_df <- top_3_df %>%
  left_join(db_csv %>% dplyr::select(DrugBank.ID, Common.name),
            by = c("drug1" = "DrugBank.ID")) %>%
  rename(drug1_name = Common.name) %>%
  left_join(db_csv %>% dplyr::select(DrugBank.ID, Common.name),
            by = c("drug2" = "DrugBank.ID")) %>%
  rename(drug2_name = Common.name)
# reorder the columns
top_3_df <- top_3_df %>%
  rename(drug1_id = drug1, drug2_id = drug2) %>%
  dplyr::select(condition_id, condition_name, drug1_id, drug1_name, drug2_id, drug2_name, scores)

### save the table
# write.csv(top_3_df, file = paste0(save_dir, "top_3_triplets_per_condition_w_names.csv"), row.names = FALSE)


#####-----all predicted triplets passing threshold 0.1-----#####
# read the result table
all_triplets_df <- fread(paste0(result_dir, "triplets_pass_threshold.csv"))
all_triplets_df <- all_triplets_df[, -1, with = FALSE]
all_triplets_df <- as.data.table(all_triplets_df)

### Map mesh code to lower case mesh term
# Join and pick only the first match if multiple exist
all_triplets_df <- rrf_sub[all_triplets_df, on = c("CODE" = "condition"), mult = "first"]
# Add condition_name column
all_triplets_df[, condition_name := STR_lowercase]
# Rename original MeSH code column
setnames(all_triplets_df, "CODE", "condition_id")
# Remove temporary column
all_triplets_df[, STR_lowercase := NULL]
# Check mapping
sum(is.na(all_triplets_df$condition_name))  # should be 0 if all codes matched

### Map drug bank code to drug bank name
# mapping
all_triplets_df <- all_triplets_df %>%
  left_join(db_csv %>% dplyr::select(DrugBank.ID, Common.name),
            by = c("drug1" = "DrugBank.ID")) %>%
  rename(drug1_name = Common.name) %>%
  left_join(db_csv %>% dplyr::select(DrugBank.ID, Common.name),
            by = c("drug2" = "DrugBank.ID")) %>%
  rename(drug2_name = Common.name)
# reorder the columns
all_triplets_df <- all_triplets_df %>%
  rename(drug1_id = drug1, drug2_id = drug2) %>%
  dplyr::select(drug1_id, drug1_name, drug2_id, drug2_name, condition_id, condition_name, scores)

### save the table
# write.csv(all_triplets_df, file = paste0(save_dir, "triplets_pass_threshold_w_names.csv"), row.names = FALSE)





