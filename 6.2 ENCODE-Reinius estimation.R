---
title: "02-05-20 epiDISH estimation"
author: "Lauren Middleton"
date: "7/16/2020"
output: html_document
---
# Purpose:  make a plot showing epidish estimates for all 4 saliva sample types
#
# Inputs:   "01-30-20 06-11-19 pd_final_all_samps.rds" - descriptives of all saliva samples
#           "06-11-19 beta_final_all_samps.rds" - beta matrix of all saliva samples
#
# Outputs:  plot - figure 5a
  
###############################################################################################
###################################### Load the datasets ######################################
###############################################################################################

setwd("C:/Users/HP/Documents/Research 2019")

# Load datasets
pd_final_all <- readRDS("06-11-19 pd_final_all_samps.rds")
beta_final_all <- readRDS("06-11-19 beta_final_all_samps.rds")
  
###############################################################################################
################################## Run ewastools estimates ####################################
###############################################################################################

#run ewastools package with encode-reinius reference panel
library(ewastools)
out.1 <- estimateLC(meth = beta_final_all, ref = "lauren_encode_reinius_ref")
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
#81.75432 97.38068
median(epi_epith_est$estimate)
#93.41753%
IQR(epi_epith_est$estimate)
#4.7%

################################# separate out IMMUNE fractions ###############################
est_long_immune <- enc_rein_long[enc_rein_long$celltype == "Immune Fractions", ]
# View(est_long_immune)
#choose only the immune estimates
immune_est <- est_long_immune[est_long_immune$cell_type_est == "Immune", ]
# View(immune_est)
range(immune_est$estimate)
#12.99822 88.94242%
median(immune_est$estimate)
#58.9%
IQR(immune_est$estimate)
#26.4%

################################## separate out WHOLE fractions ###############################
est_long_whole <- enc_rein_long[enc_rein_long$celltype == "Whole Samples", ]
# View(est_long_whole)
#choose only the immune estimates
immune_est_whole <- est_long_whole[est_long_whole$cell_type_est == "Immune", ]
# View(immune_est_whole)
range(immune_est_whole$estimate)
#10.82925 88.06632% immune
median(immune_est_whole$estimate)
#29.9%
IQR(immune_est_whole$estimate)
#24.6%

#find the range of the epithelial cells within whole samples
epith_est_whole <- est_long_whole[est_long_whole$cell_type_est == "Epithelial", ]
range(epith_est_whole$estimate)
#8.368852 88.291419%
median(epith_est_whole$estimate)
#71.75%
IQR(epith_est_whole$estimate)
#25.6%

################################# separate out ORAGENE fractions ##############################
est_long_ora <- enc_rein_long[enc_rein_long$celltype == "Oragene Kits",]
# View(est_long_ora)
#choose only the immune estimates
immune_est_ora <- est_long_ora[est_long_ora$cell_type_est == "Immune", ]
#View(immune_est_ora)
range(immune_est_ora$estimate)
#30.60612 69.61919% immune
median(immune_est_ora$estimate)
#55.4%

#epithelial cells:
epith_est_ora <- est_long_ora[est_long_ora$cell_type_est == "Epithelial", ]
# View(immune_est_ora)
range(epith_est_ora$estimate)
#19.79448 59.84410%
median(epith_est_ora$estimate)
#28.5


###############################################################################################
############################ Graph the estimates as a violin plot #############################
###############################################################################################

setwd("C:/Users/HP/Documents/Research 2020")
library(ggplot2)

enc_rein_long$cell_type_est <- factor(enc_rein_long$cell_type_est,
                                      levels = c("Immune", "Epithelial"))


png("07-29-20 encode-reinius prediction panel.png", width=1000, height=1000)
ggplot(data = enc_rein_long, aes(x = cell_type_est, y = estimate, fill = cell_type_est))+
  # geom_violin(scale="count", width=0.5, position=position_dodge(width=1.5))+
  #geom_point(size = 2)+
  geom_boxplot()+
  geom_jitter(width = 0.05, size = 3)+
  xlab("")+
  ylab("Percent Cell Type (%)")+
  scale_y_continuous(breaks = seq(0,120,10),
                     limits = c(-1,120))+ #set the y axis from 0-100 labeled every 10 units
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