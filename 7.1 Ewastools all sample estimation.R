---
title: "07-16-20 Ewastools estimation graphs"
author: "Lauren Middleton"
date: "7/16/2020"
output: html_document
---
# Purpose:  generate basic stats for all four cell-type samples
#
# Inputs:   "01-30-20 06-11-19 pd_final_all_samps.rds" - descriptives of all saliva samples
#           "06-11-19 beta_final_all_samps.rds" - beta matrix of all saliva samples
#
# Outputs:  stats
  
###############################################################################################
###################################### Load the datasets ######################################
###############################################################################################

setwd("~/Research 2019")

# Load datasets
pd_final_all <- readRDS("06-11-19 pd_final_all_samps.rds")
beta_final_all <- readRDS("06-11-19 beta_final_all_samps.rds")
  
###############################################################################################
################################## Run ewastools estimates ####################################
###############################################################################################

#run ewastools package with encode-reinius reference panel
library(ewastools)
out.1 <- estimateLC(meth = beta_final_all, ref = "lauren_saliva_ref")
head(out.1)


#make it a dataframe
est <- as.data.frame(out.1)

###############################################################################################
###################################### Fix the dataset ########################################
###############################################################################################

library(dplyr)
library(tidyverse)


#add sample names and cell types
est$sample <- pd_final_all$Sample.ID
est$celltype <- pd_final_all$celltype
table(est$celltype)
# CD45pos   large oragene   whole 
#      20      18       4      18


#make the cell type names nice
est$celltype <- gsub("large", "Epithelial Fractions", est$celltype)
est$celltype <- gsub("CD45pos", "Immune Fractions", est$celltype)
est$celltype <- gsub("whole", "Whole Samples", est$celltype)
est$celltype <- gsub("oragene", "Oragene Kits", est$celltype)
table(est$celltype)


#make the dataset into the long form
enc_rein_long <- gather(est, cell_type_est, estimate, Epi:IC, factor_key=TRUE)
# View(enc_rein_long)

#change proportions into percentages
enc_rein_long$estimate <- enc_rein_long$estimate*100

#Fix IC, Epi to real names
table(enc_rein_long$cell_type_est) #60 of each
enc_rein_long$cell_type_est <- gsub("IC", "Immune", enc_rein_long$cell_type_est)
enc_rein_long$cell_type_est <- gsub("Epi", "Epithelial", enc_rein_long$cell_type_est)
table(enc_rein_long$cell_type_est) #60 of each

###############################################################################################
######################## Ranges of cell types and make cell-type datasets #####################
###############################################################################################

############################### separate out EPITHELIAL fractions #############################
est_long_epith <- enc_rein_long[enc_rein_long$celltype == "Epithelial Fractions", ]
#View(est_long_epith)
#choose only the epithelial estimates
epi_epith_est <- est_long_epith[est_long_epith$cell_type_est == "Epithelial", ]
#View(epi_epith_est)
range(epi_epith_est$estimate)
#76.63 105.29
median(epi_epith_est$estimate)
#102.2%
IQR(epi_epith_est$estimate)
#6.1%

################################# separate out IMMUNE fractions ###############################
est_long_immune <- enc_rein_long[enc_rein_long$celltype == "Immune Fractions", ]
# View(est_long_immune)
#choose only the immune estimates
immune_est <- est_long_immune[est_long_immune$cell_type_est == "Immune", ]
# View(immune_est)
range(immune_est$estimate)
#48.88 113.49%
median(immune_est$estimate)
#105.4%
IQR(immune_est$estimate)
#28.6%

################################## separate out WHOLE fractions ###############################
est_long_whole <- enc_rein_long[enc_rein_long$celltype == "Whole Samples", ]
# View(est_long_whole)
#choose only the immune estimates
immune_est_whole <- est_long_whole[est_long_whole$cell_type_est == "Immune", ]
# View(immune_est_whole)
range(immune_est_whole$estimate)
#13.21 113.19% immune
median(immune_est_whole$estimate)
#48.2%
IQR(immune_est_whole$estimate)
#49.5%

#find the range of the epithelial cells within whole samples
epith_est_whole <- est_long_whole[est_long_whole$cell_type_est == "Epithelial", ]
range(epith_est_whole$estimate)
#0.0 90.33%
median(epith_est_whole$estimate)
#56.1%
IQR(epith_est_whole$estimate)
#48.8

################################# separate out ORAGENE fractions ##############################
est_long_ora <- enc_rein_long[enc_rein_long$celltype == "Oragene Kits",]
# View(est_long_ora)
#choose only the immune estimates
immune_est_ora <- est_long_ora[est_long_ora$cell_type_est == "Immune", ]
#View(immune_est_ora)
range(immune_est_ora$estimate)
#72.8 110.0% immune
median(immune_est_ora$estimate)
#103.8%
IQR(immune_est_ora$estimate)

#epithelial cells:
epith_est_ora <- est_long_ora[est_long_ora$cell_type_est == "Epithelial", ]
# View(immune_est_ora)
range(epith_est_ora$estimate)
#0.0 28.9%
median(epith_est_ora$estimate)
#0