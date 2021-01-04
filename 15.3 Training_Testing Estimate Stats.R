#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####################################  CALCULATE STATS FOR CELL ESTIMATES  #####################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose:  This code file calculates the stats for the cell estimates from training/testing split
#          
# Inputs:   long_encode - Long dataframe of cell proportion estimates from encode ref
#           long_saliva - Long dataframe of cell proportion estimates from saliva ref
#
# Outputs:  Statistics (range, median, iqr)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#############################################################################################################
###################################### Load Libraries And Set Directory #####################################
#############################################################################################################

#NOTE: load w - training_testing_global_environment to run this code

library(tidyverse)

setwd(current_directory)

#############################################################################################################
################################## Combine ENCODE and Saliva Long Estimates #################################
#############################################################################################################

#add label to the adult and peds datasets to specify reference panel type
est_saliva_clean$reference <- "saliva"
est_encode_clean$reference <- "encode_reinius"

#combine all datasets into one long dataset
combined_est_long <- bind_rows(est_encode_clean,
                               est_saliva_clean)

#clean up the long dataset
combined_est_long$Epithelial <- combined_est_long$Epithelial*100
combined_est_long$Immune <- combined_est_long$Immune*100

#############################################################################################################
######################################### Calculate The Statistics ##########################################
#############################################################################################################

stats <- combined_est_long %>%
  group_by(reference, celltype) %>%
  dplyr::summarise(
    n = n(),
    min_epithelial = min(Epithelial),
    max_epithelial = max(Epithelial),
    min_immune = min(Immune),
    max_immune = max(Immune),
    median_epithelial = median(Epithelial),
    median_immune = median(Immune),
    IQR_epithelial = IQR(Epithelial),
    IQR_immune = IQR(Immune)
  ) %>%
  arrange(n) %>%
  ungroup()
