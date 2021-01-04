#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
###############################  ESTIMATE THE CELL PROPORTIONS IN CHILD SALIVA  ###############################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose:  This code file estimates the proportion of cells in the child saliva
#          
# Inputs:   beta_child - This matrix contains the cleaned DNAm data for the child saliva samples 
# 
#           pd_child   - This dataframe contains the information about the saliva samples
#
# Outputs:  child_est  - Dataframe of estimates of epithelial and immune cells per sample
#           child_long - Long dataframe of child_est

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#############################################################################################################
###################################### Load Libraries And Set Directory #####################################
#############################################################################################################

#NOTE: load w - child_saliva_global_environment.RData to run this code

library(tidyverse)
library(ewastools)

current_directory <- "C:/Users/HP/Google Drive/Colacino Lab/Saliva/Peds Saliva"
setwd(current_directory)

#############################################################################################################
######################################## Estimate The Cell Proportions ######################################
#############################################################################################################

#use ewastools package to estimate cell proportions

#saliva
child_unclean_est_saliva <- estimateLC(meth = beta_child,
                                       ref = "saliva",
                                       constrained = TRUE)

#encode/reinius
child_unclean_est_encode <- estimateLC(meth = beta_child,
                                       ref = "encode_reinius",
                                       constrained = TRUE)

#############################################################################################################
########################################### Clean Up The Estimates ##########################################
#############################################################################################################

#make it a dataframe
child_est_saliva <- as.data.frame(child_unclean_est_saliva)
child_est_encode <- as.data.frame(child_unclean_est_encode)

#make the cell type names nice - flip encode because the reference panel is backwards
child_est_saliva$Epithelial <- child_unclean_est_saliva$Epithelial.cells
child_est_saliva$Immune <- child_unclean_est_saliva$Leukocytes
child_est_encode$Epithelial <- child_unclean_est_encode$IC
child_est_encode$Immune <- child_unclean_est_encode$Epi

#add the gsm ids to the dataset
identical(colnames(beta_child), rownames(pd_child))
#TRUE
child_est_saliva$GSM <- pd_child$geo_accession
child_est_encode$GSM <- pd_child$geo_accession

#add the sex to the dataset
child_est_saliva$sex <- pd_child$sex
child_est_encode$sex <- pd_child$sex

#make the dataset long
child_long_saliva <- gather(data = child_est_saliva,
                            key = celltype,
                            value = percent,
                            Epithelial:Immune)
child_long_encode <- gather(data = child_est_encode,
                            key = celltype,
                            value = percent,
                            Epithelial:Immune)

#change proportions into percentages
child_long_saliva$percent <- child_long_saliva$percent*100
child_long_encode$percent <- child_long_encode$percent*100

rm(child_unclean_est_saliva,
   child_unclean_est_encode)
