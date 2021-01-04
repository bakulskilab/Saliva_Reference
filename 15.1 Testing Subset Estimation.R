#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
###############################  ESTIMATE CELL PROPORTIONS IN TESTING SUBSET  #################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose:  This code file estimates the proportion of cells in the testing subset of saliva samples
#          
# Inputs:   beta_final   - beta matrix of sorted saliva samples
#           pd_final     - demographics file for sorted samples
#           training_pd  - demographics file for training subset of samples
#           meth_id      - vector of meth_ids from the training_pd
#
# Outputs:  est_encode_clean - Dataframe of estimates of epithelial and immune cells per sample
#           long_encode      - Long dataframe of est_encode_clean
#           est_saliva_clean - Dataframe of estimates of epithelial and immune cells per sample
#           long_saliva      - Long dataframe of est_saliva_clean

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#############################################################################################################
###################################### Load Libraries And Set Directory #####################################
#############################################################################################################

# Load global environment: w - training_testing_global_environment.RData

setwd(current_directory)

library(tidyverse)

#for the developer branch of ewastools:
# .libPaths("C:/Users/HP/Documents/R-dev")
library(ewastools)

#############################################################################################################
######################################### Subset The Testing Dataset ########################################
#############################################################################################################

#select the testing fractions based on the meth_ids from the training subset
testing_pd <- anti_join(pd_final, training_pd, by = "meth_id")

beta_final_df <- as.data.frame(beta_final)

testing_beta <- beta_final_df %>%
  select(!all_of(meth_id))

rm(beta_final_df)

#############################################################################################################
######################################## Estimate The Cell Proportions ######################################
#############################################################################################################

#use ewastools package to estimate cell proportions

#saliva
unclean_est_saliva <- estimateLC(meth = testing_beta,
                                 ref = "saliva_training_constrained",
                                 constrained = TRUE)

#encode/reinius
unclean_est_encode <- estimateLC(meth = testing_beta,
                                 ref = "lauren_encode_reinius_ref",
                                 constrained = TRUE)

#############################################################################################################
########################################### Clean Up The Estimates ##########################################
#############################################################################################################

#make it a dataframe
unclean_est_saliva <- as.data.frame(unclean_est_saliva)
unclean_est_encode <- as.data.frame(unclean_est_encode)

#make the cell type names nice
est_saliva <- unclean_est_saliva %>%
  rename(Epithelial = Epi,
         Immune = IC)
#the estimates got flipped so i'm putting it back
est_encode <- unclean_est_encode %>%
  rename(Immune = Epi,
         Epithelial = IC)

#add the meth_ids to the dataset
identical(colnames(testing_beta), testing_pd$meth_id)
#TRUE
est_saliva$meth_id <- testing_pd$meth_id
est_encode$meth_id <- testing_pd$meth_id

#add cell type to the dataset
est_saliva$celltype <- testing_pd$celltype
est_encode$celltype <- testing_pd$celltype

#fix the cell types
est_saliva_clean <- est_saliva %>%
  mutate(celltype = ifelse(celltype == "large",
                           yes = "Epithelial Fractions",
                           no = "Immune Fractions")
         )
est_encode_clean <- est_encode %>%
  mutate(celltype = ifelse(celltype == "large",
                           yes = "Epithelial Fractions",
                           no = "Immune Fractions")
         )

#############################################################################################################
############################################ Make A Long Dataset ############################################
#############################################################################################################

#make the dataset long
long_saliva <- gather(data = est_saliva_clean,
                      key = est_celltype,
                      value = percent,
                      Epithelial:Immune)
long_encode <- gather(data = est_encode_clean,
                      key = est_celltype,
                      value = percent,
                      Epithelial:Immune)

#change proportions into percentages
long_saliva$percent <- long_saliva$percent*100
long_encode$percent <- long_encode$percent*100

rm(unclean_est_saliva,
   unclean_est_encode,
   est_encode,
   est_saliva)
