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
out.1 <- estimateLC(meth = beta_final_all,
                    ref = "saliva",
                    constrained = TRUE)
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

#make the ewastools cell names nice
est_clean <- est %>%
  rename(Immune = Leukocytes,
         Epithelial = Epithelial.cells)
table(est_clean$celltype)

#make the dataset into the long form
saliva_long <- gather(est_clean, cell_type_est, estimate, Immune:Epithelial, factor_key=TRUE)
# View(saliva_long)

#change proportions into percentages
saliva_long$estimate <- saliva_long$estimate*100

###############################################################################################
######################## Ranges of cell types and make cell-type datasets #####################
###############################################################################################

stats <- saliva_long %>%
  group_by(celltype, cell_type_est) %>%
  dplyr::summarise(
    n = n(),
    min = min(estimate),
    max = max(estimate),
    median = median(estimate),
    IQR = IQR(estimate),
    mean = mean(estimate)
  ) %>%
  arrange(n) %>%
  ungroup()

# ############################### separate out EPITHELIAL fractions #############################
# est_long_epith <- enc_rein_long[enc_rein_long$celltype == "Epithelial Fractions", ]
# #View(est_long_epith)
# #choose only the epithelial estimates
# epi_epith_est <- est_long_epith[est_long_epith$cell_type_est == "Epithelial", ]
# #View(epi_epith_est)
# range(epi_epith_est$estimate)
# #76.63 105.29
# median(epi_epith_est$estimate)
# #102.2%
# IQR(epi_epith_est$estimate)
# #6.1%
# 
# ################################# separate out IMMUNE fractions ###############################
# est_long_immune <- enc_rein_long[enc_rein_long$celltype == "Immune Fractions", ]
# # View(est_long_immune)
# #choose only the immune estimates
# immune_est <- est_long_immune[est_long_immune$cell_type_est == "Immune", ]
# # View(immune_est)
# range(immune_est$estimate)
# #48.88 113.49%
# median(immune_est$estimate)
# #105.4%
# IQR(immune_est$estimate)
# #28.6%
# 
# ################################## separate out WHOLE fractions ###############################
# est_long_whole <- enc_rein_long[enc_rein_long$celltype == "Whole Samples", ]
# # View(est_long_whole)
# #choose only the immune estimates
# immune_est_whole <- est_long_whole[est_long_whole$cell_type_est == "Immune", ]
# # View(immune_est_whole)
# range(immune_est_whole$estimate)
# #13.21 113.19% immune
# median(immune_est_whole$estimate)
# #48.2%
# IQR(immune_est_whole$estimate)
# #49.5%
# 
# #find the range of the epithelial cells within whole samples
# epith_est_whole <- est_long_whole[est_long_whole$cell_type_est == "Epithelial", ]
# range(epith_est_whole$estimate)
# #0.0 90.33%
# median(epith_est_whole$estimate)
# #56.1%
# IQR(epith_est_whole$estimate)
# #48.8
# 
# ################################# separate out ORAGENE fractions ##############################
# est_long_ora <- enc_rein_long[enc_rein_long$celltype == "Oragene Kits",]
# # View(est_long_ora)
# #choose only the immune estimates
# immune_est_ora <- est_long_ora[est_long_ora$cell_type_est == "Immune", ]
# #View(immune_est_ora)
# range(immune_est_ora$estimate)
# #72.8 110.0% immune
# median(immune_est_ora$estimate)
# #103.8%
# IQR(immune_est_ora$estimate)
# 
# #epithelial cells:
# epith_est_ora <- est_long_ora[est_long_ora$cell_type_est == "Epithelial", ]
# # View(immune_est_ora)
# range(epith_est_ora$estimate)
# #0.0 28.9%
# median(epith_est_ora$estimate)
# #0

###############################################################################################
############################ Graph the estimates as a violin plot #############################
###############################################################################################

setwd("C:/Users/HP/Documents/Research 2020")
library(ggplot2)

saliva_long$cell_type_est <- factor(saliva_long$cell_type_est,
                                    levels = c("Immune", "Epithelial"))


png("07-29-20 saliva prediction panel.png", width=1000, height=1000)
ggplot(data = saliva_long, aes(x = cell_type_est, y = estimate, fill = cell_type_est))+
  geom_violin()+
  geom_jitter(width = 0.05, size = 3)+
  xlab("")+
  ylab("Percent Cell Type (%)")+
  scale_y_continuous(breaks = seq(0,100,10),
                     limits = c(-1,101))+ #set the y axis from 0-100 labeled every 10 units
  theme(text = element_text(size=35))+ #makes text bigger
  labs(fill = "Estimated Cell Types")+ #change the legend title
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))+ #change axis text color
  theme(panel.background = element_rect(fill = "white", #this is the background
                                        colour = "black", #this is the border
                                        size = 0.1, linetype = "solid"))+
  scale_fill_manual(values=c("#F8766D", "#00BFC4", "yellow"))+
  theme(legend.position="bottom")+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  facet_wrap(vars(celltype))
dev.off()