#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##############################  Determine Which Chemicals To Include In Analysis  #############################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose: This code file calculates the paired t-tests between paired epithelial and immune samples
#          
# Inputs:   pd_final   - dataframe containing information about the samples and demographics 
# 
#           beta_final - dataframe containing 
#
# Outputs:  results from a paired t-test of cell types in the beta matrix, number of significant probes

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#############################################################################################################
######################################## Load Libraries And Datasets ########################################
#############################################################################################################

library(tidyverse)
library(dplyr)
library(matrixTests)
library(genefilter)

setwd("C:/Users/HP/Documents/Research 2019")

#load in data files
#these don't include whole or oragene or c45neg
pd_final <- readRDS("06-11-19 pd_final.rds")
beta_final <- readRDS("06-11-19 beta_final.rds")

#############################################################################################################
################################### Determine Which Fractions Have a Pair ###################################
#############################################################################################################

#drop the fractions that don't have an epith/immune pair
paired_df <- pd_final %>%
  group_by(Sample_Group) %>%
  filter(n() > 1) %>%
  ungroup()

table(paired_df$celltype)
# CD45pos   large 
# 17        17

#############################################################################################################
################################### Subset Beta Matrix to These Fractions ###################################
#############################################################################################################

#38 original samples in beta matrix, need to drop that to 34
dim(beta_final)

#convert beta matrix to dataframe so tidyverse::select works
beta_final <- as.data.frame(beta_final)

#identify the meth_ids for the paired samples
meth_id_pairs <- paired_df$meth_id

#select only the paired samples
beta_paired <- beta_final %>%
  select(all_of(meth_id_pairs))

#############################################################################################################
################################ Split Beta Matrix Into Epithelial And Immune ###############################
#############################################################################################################

#split paired_df by cell type
immune_pairs <- paired_df %>%
  filter(celltype == "CD45pos")
epith_pairs <- paired_df %>%
  filter(celltype == "large")

#grab the meth_ids from both cell types
immune_meth_id <- immune_pairs$meth_id
epith_meth_id <- epith_pairs$meth_id

#create immune cells beta matrix
beta_immune <- beta_final %>%
  select(all_of(immune_meth_id))

#create epithelial cells beta matrix
beta_epith <- beta_final %>%
  select(all_of(epith_meth_id))

#############################################################################################################
########################################### Calculate Row T-tests ###########################################
#############################################################################################################

#do a paired row ttest of immune vs epithelial beta values for each probe
ttest_results <- row_t_paired(x = beta_immune,
                              y = beta_epith,
                              alternative = "two.sided",
                              mu = 0, #null hypothesis
                              conf.level = 0.95)

#############################################################################################################
################################### Calculate Number of Significant Probes ##################################
#############################################################################################################

#what is the bonferroni cutoff
num_probes <- length(rownames(beta_immune))
bonferroni_long <- 0.05 / num_probes
bonferroni <- round(bonferroni_long, digits = 10)
#6.28e-08

#how many sig probes at bonferroni level
sig_probes <- length(which(ttest_results$pvalue < bonferroni))
#111,922

#percent of total probes
(sig_probes / num_probes)*100
#14.1%

#which probes?
probes_paired <- ttest_results[which(ttest_results$pvalue < bonferroni), ]
names_probes_paired <- rownames(probes_paired)

#############################################################################################################
############################### Calculate Unpaired Row T-tests For All Samples ##############################
#############################################################################################################

pd_final$cd45pos <-ifelse(pd_final$celltype =="CD45pos", 1, 0)

#check the coding;
table(pd_final$celltype, pd_final$cd45pos)

# Run row ttests
#do a row ttest of one cell types against all other cell types within one probe
ttest_cd45pos <- rowttests(as.matrix(beta_final), as.factor(pd_final$cd45pos), tstatOnly = FALSE)

#which probes?
probes_unpaired <- ttest_cd45pos[which(ttest_cd45pos$p.value < bonferroni), ]
names_probes_unpaired <- rownames(probes_unpaired)


#overlap between paired and unpaired
length(intersect(names_probes_unpaired, names_probes_paired))
#111819 - all probes overlap from paired to unpaired