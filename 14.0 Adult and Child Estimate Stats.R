#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####################################  CALCULATE STATS FOR CELL ESTIMATES  #####################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose:  This code file calculates the stats for the cell estimates of adult and peds datasets
#          
# Inputs:   adult_long - Long dataframe of cell proportion estimates from adult saliva
#           child_long - Long dataframe of cell proportion estimates from child saliva
#
# Outputs:  Statistics (range, median, iqr)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#############################################################################################################
###################################### Load Libraries And Set Directory #####################################
#############################################################################################################

#NOTE: load w - child_saliva_global_environment.RData and w - child_saliva_global_environment.RData
#      to run this code

library(tidyverse)

setwd(current_directory)

#############################################################################################################
################################### Combine Adult And Peds Long Estimates ###################################
#############################################################################################################

#add label to the adult and peds datasets to specify reference panel type
adult_est_saliva$reference <- "saliva"
adult_est_encode$reference <- "encode_reinius"
child_est_saliva$reference <- "saliva"
child_est_encode$reference <- "encode_reinius"

#add label to the adult and peds datasets to specify dataset
adult_est_saliva$dataset <- "adult"
adult_est_encode$dataset <- "adult"
child_est_saliva$dataset <- "child"
child_est_encode$dataset <- "child"

#combine all four datasets into one long dataset
combined_est_long <- bind_rows(adult_est_encode,
                               adult_est_saliva,
                               child_est_encode,
                               child_est_saliva)

#clean up the long dataset
combined_est_long$Epithelial <- combined_est_long$Epithelial*100
combined_est_long$Immune <- combined_est_long$Immune*100

#############################################################################################################
######################################### Calculate The Statistics ##########################################
#############################################################################################################

stats <- combined_est_long %>%
  group_by(dataset, reference) %>%
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