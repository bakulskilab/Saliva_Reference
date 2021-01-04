#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
###############################  ESTIMATE THE CELL PROPORTIONS IN ADULT SALIVA  ###############################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose:  This code file estimates the proportion of cells in the adult saliva
#          
# Inputs:   beta_adult - This matrix contains the cleaned DNAm data for the adult saliva samples 
# 
#           pd_adult   - This dataframe contains the information about the saliva samples
#
# Outputs:  adult_est_saliva/encode  - Dataframe of estimates of epithelial and immune cells per sample
#           adult_long_saliva/encode - Long dataframe of adult_est_saliva/encode

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#############################################################################################################
###################################### Load Libraries And Set Directory #####################################
#############################################################################################################

#NOTE: load w - adult_saliva_global_environment.RData to run this code

library(tidyverse)
library(ewastools)

current_directory <- "C:/Users/HP/Google Drive/Colacino Lab/Saliva/Adult Saliva"
setwd(current_directory)

#############################################################################################################
######################################## Estimate The Cell Proportions ######################################
#############################################################################################################

#use ewastools package to estimate cell proportions - saliva reference
adult_unclean_est_saliva <- estimateLC(meth = beta_adult,
                                       ref = "saliva",
                                       constrained = TRUE)

#use ewastools package to estimate cell proportions - encode reference
adult_unclean_est_encode <- estimateLC(meth = beta_adult,
                                       ref = "encode_reinius",
                                       constrained = TRUE)

#############################################################################################################
########################################### Clean Up The Estimates ##########################################
#############################################################################################################

#make it a dataframe
adult_est_saliva <- as.data.frame(adult_unclean_est_saliva)
adult_est_encode <- as.data.frame(adult_unclean_est_encode)

#make the cell type names nice - and flip encode cell types because ref panel is backwards
adult_est_saliva$Epithelial <- adult_est_saliva$Epi
adult_est_saliva$Immune <- adult_est_saliva$Leukocytes
adult_est_encode$Epithelial <- adult_est_encode$IC
adult_est_encode$Immune <- adult_est_encode$Epi


#add the gsm ids to the dataset
identical(colnames(beta_adult), pd_adult$geo_accession)
#TRUE
adult_est_saliva$GSM <- pd_adult$geo_accession
adult_est_encode$GSM <- pd_adult$geo_accession

#add the sex to the dataset
adult_est_saliva$sex <- pd_adult$sex
adult_est_encode$sex <- pd_adult$sex

#make the dataset long
adult_long_saliva <- gather(data = adult_est_saliva,
                            key = celltype,
                            value = percent,
                            Epithelial:Immune)
adult_long_encode <- gather(data = adult_est_encode,
                            key = celltype,
                            value = percent,
                            Epithelial:Immune)

#change proportions into percentages
adult_long_saliva$percent <- adult_long_saliva$percent*100
adult_long_encode$percent <- adult_long_encode$percent*100

rm(adult_unclean_est_saliva,
   adult_unclean_est_encode)