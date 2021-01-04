setwd("~/Research 2019")

# Load dataset
pd_final_all <- readRDS("06-11-19 pd_final_all_samps.rds")
pd.qc <- readRDS("pd.qc.rds")

#select the columns to keep
covariate_file_unQC <- pd.qc %>%
  select(meth_id,
         Sample.ID,
         id,
         celltype,
         month,
         year,
         age.r,
         Male_label,
         Hispanic,
         White,
         Black,
         AmerI.AN,
         Asian,
         Other,
         race_birac,
         Sick,
         cell_count,
         viability,
         volume,
         overnight,
         hours,
         batch,
         batch_date
         )

covariate_file_QC <- pd_final_all %>%
  select(meth_id,
         Sample.ID,
         id,
         celltype,
         month,
         year,
         age.r,
         Male_label,
         Hispanic,
         White,
         Black,
         AmerI.AN,
         Asian,
         Other,
         race_birac,
         Sick,
         cell_count,
         viability,
         volume,
         overnight,
         hours,
         batch,
         batch_date
  )

#save as csv
setwd("~/Research 2020")
write.csv(covariate_file_unQC, "saliva_unQC_covariates.csv")
write.csv(covariate_file_QC, "saliva_QC_covariates.csv")