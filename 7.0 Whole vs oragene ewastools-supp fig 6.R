---
title: "07-15-20 Whole vs oragene"
author: "Lauren Middleton"
date: "1/23/2020"
output: html_document
---
# Purpose:  make two spaghetti plot figures
#             1) whole vs oragene predicted immune cell %
#             2) within whole, predicted immune vs epith % by sample to check that high % immune means low % epith
#
# Inputs:   "01-30-20 06-11-19 pd_final_all_samps.rds" - descriptives of all saliva samples
#           "06-11-19 beta_final_all_samps.rds" - beta matrix of all saliva samples
#
# Outputs:  plots

###############################################################################################
###################################### Load the datasets ######################################

# Load datasets
setwd("C:/Users/HP/Documents/Research 2019")

pd_final_all <- readRDS("06-11-19 pd_final_all_samps.rds")
beta_final_all <- readRDS("06-11-19 beta_final_all_samps.rds")

###############################################################################################
##################################### Set up the datasets #####################################
library(ewastools)
library(tidyverse)
setwd("~/Research 2019")

#check that the order of pd_final_all$meth is the same as beta_final_all columns
beta_meth <- colnames(beta_final_all)
pd_meth <- pd_final_all$meth_id
beta_meth == pd_meth #all TRUE

#add the names of samples (Sample.ID) to colnames of beta matrix
if(pd_final_all$meth_id == colnames(beta_final_all))
{
   colnames(beta_final_all) <- pd_final_all$Sample.ID
}

#grab the whole and oragene samples from beta
beta_whole <- beta_final_all[, endsWith(colnames(beta_final_all), "_whole") == "TRUE"]
dim(beta_whole) #795694     18
beta_oragene <- beta_final_all[, endsWith(colnames(beta_final_all), "_oragene") == "TRUE"]
dim(beta_oragene) #795694     4

#make the cell proportion estimates
w <- estimateLC(beta_whole,
                ref = "saliva",
                constrained = TRUE)
o <- estimateLC(beta_oragene,
                ref = "saliva",
                constrained = TRUE)

#attach sample ids and cell types and ids to whole and oragene estimates
sample_id <- pd_final_all[pd_final_all$celltype == "whole",]
w$sample_id <- sample_id$Sample.ID
w$celltype <- sample_id$celltype
w$id <- sample_id$id
sample_id <- pd_final_all[pd_final_all$celltype == "oragene",]
o$sample_id <- sample_id$Sample.ID
o$celltype <- sample_id$celltype
o$id <- sample_id$id

#merge whole and oragene estimates
w_o_est <- rbind(w, o)
dim(w_o_est) #22(18+4) 5
colnames(w_o_est) #"Leukocytes"       "Epithelial.cells"        "sample_id"  "id"


#saliva percent conversion
w_o_est$Epi <- w_o_est$Epithelial.cells*100
w_o_est$IC <- w_o_est$Leukocytes*100
w_o_est$Epithelial <- w_o_est$Epi
w_o_est$Immune <- w_o_est$IC

estimates <- w_o_est %>%
              select(-c(Epi,
                        IC,
                        Leukocytes,
                        Epithelial.cells))


# Check that all whole and oragene samples passed QC
#are all whole samples in oragene
# Oragene samples: 5, 6, 21, 23
#   table(estimates$sample_id, estimates$celltype)
# Whole samples:  #5, 6, 23 - 21_whole doesn't have a good sample and was dropped due to probe fail rate >= 3%
# there are 3 whole/oragene pairs


# Drop 21 oragene
#keep all columns, drop 2 rows
estimates <- estimates[!estimates$sample_id == "SAL21_oragene", ]
dim(estimates) #21 5


###############################################################################################
######################## Make plot 1: oragene vs whole without sample 21 ######################
###############################################################################################
library(ggplot2)

# Add labels to the factor levels
estimates$celltype <- factor(estimates$celltype,
                                     levels = c("whole","oragene"),
                                     labels = c("Whole",
                                                "Oragene"))
# Change to factor levels
estimates$celltype <- as.factor(estimates$celltype)
estimates$id <- as.factor(estimates$id)

#grab only the three whole and oragene pairs
est_5 <- estimates[estimates$id == "5", ]
est_6 <- estimates[estimates$id == "6", ]
est_23 <-estimates[estimates$id == "23", ]
est_ora_whole <- rbind(est_5, est_6, est_23)
rm(est_23, est_5, est_6)

est_ora_whole$`as.factor(id)` <- est_ora_whole$id

setwd("~/Research 2020")
png("07-15-20 oragene vs whole percent immune subset.png", width = 500, height = 500)
ggplot(data = est_ora_whole, aes(x = factor(celltype, level = c("Whole", "Oragene")),
       y = Immune, group = `as.factor(id)`))+
  geom_point(aes(color = `as.factor(id)`), size = 3)+
  geom_line(aes(color = `as.factor(id)`), size = 1)+
  xlab("")+
  ylab("Percent Immune (%)")+
  scale_y_continuous(breaks = seq(0,100,10),
                     limits = c(0,100))+ #set the y axis from 0-100 labeled every 10 units
  theme(text = element_text(size=25))+ #makes text bigger
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))+ #change axis text color
  theme(panel.background = element_rect(fill = "white", #this is the background
                                colour = "black", #this is the border
                                size = 0.1, linetype = "solid"))+
  theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
dev.off()

stats_oragene_saliva <- est_ora_whole %>%
  group_by(celltype) %>%
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
  ungroup()


###############################################################################################
########################## Make a Whole saliva dataset and get stats ##########################
###############################################################################################

# Starting with "estimates" dataset
str(estimates)

# Select out only the Whole samples - should be 18
#keep all columns, keep only whole rows
whole_estimates <- estimates[estimates$celltype == "Whole", ]
str(whole_estimates)

stats_whole_saliva <- whole_estimates %>%
  group_by(celltype) %>%
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

long_whole_est <- gather(data = whole_estimates,
                         key = celltype,
                         value = percent,
                         Epithelial:Immune)

#get individual whole sample stats
stats_indiv_whole_saliva <- long_whole_est %>%
  group_by(id, celltype) %>%
  dplyr::summarise(
    min = percent
  ) %>%
  ungroup()

# ## Calculate average epithelial cells in whole samples based on ewastools estimates
# mean_epith <- mean(whole_estimates$Epithelial)
# median_epith <- median(whole_estimates$Epithelial)
# sd_epith <- sd(whole_estimates$Epithelial)
# iqr_epith <- IQR(whole_estimates$Epithelial)
# mean_epith #46.62202
# median_epith #56.10145
# sd_epith #28.93889
# iqr_epith #48.84833
# 
# ## Calculate average immune cells in whole samples based on ewastools estimates
# mean_ic <- mean(whole_estimates$Immune)
# median_ic <- median(whole_estimates$Immune)
# sd_ic <- sd(whole_estimates$Immune)
# iqr_ic <- IQR(whole_estimates$Immune)
# mean_ic #58.36656
# median_ic #48.22453
# sd_ic #30.1962
# iqr_ic #49.53824

#see below for encode-reinius estimates average (line ~307) and ttest comparison

###############################################################################################
####################### Make plot 2: estimated immune vs epith % by sample ####################
###############################################################################################

## One sample number per prediction - make the dataset long
library(tidyverse)

#make the dataset long instead of wide
whole_est_long <- gather(whole_estimates,
                         est_celltype,
                         estimates,
                         Epithelial:Immune,
                         factor_key=TRUE)
str(whole_est_long)

## Graph the whole samples - ewastools
setwd("~/Research 2020")
png("07-15-20 whole estimate comparisons_update.png", width = 500, height = 500)
ggplot(data = whole_est_long, aes(x = est_celltype, y = estimates, group = id))+
  geom_point(aes(color = as.factor(id)), size = 3)+
  geom_line(aes(color = as.factor(id)), size = 1)+
  xlab("")+
  ylab("Cell Type (%)")+
  scale_y_continuous(breaks = seq(0,100,10),
                     limits = c(0,100))+ #set the y axis from 0-100 labeled every 10 units
  theme(text = element_text(size=20))+ #makes text bigger
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))+ #change axis text color
  theme(panel.background = element_rect(fill = "white", #this is the background
                                        colour = "black", #this is the border
                                        size = 0.1, linetype = "solid"))+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

###############################################################################################
##################### Make plot 3: same as plot 2, but circle sick children ###################
###############################################################################################

# sick_id <- c(4, 6, 18)

## Graph the whole samples - ewastools
setwd("~/Research 2020")
png("11-25-20 spaghetti_plot_whole_saliva_sick.png", width = 500, height = 500)
ggplot(data = whole_est_long, aes(x = est_celltype, y = estimates, group = id))+
  geom_point(aes(color = as.factor(id)), size = 3)+
  geom_line(aes(color = as.factor(id)), size = 1)+
  geom_point(data = whole_est_long[whole_est_long$id == 4,],
             pch = 21, fill = NA, size = 5, colour = "black", stroke = 2)+
  geom_point(data = whole_est_long[whole_est_long$id == 6,],
             pch = 21, fill = NA, size = 5, colour = "black", stroke = 2)+
  geom_point(data = whole_est_long[whole_est_long$id == 18,],
             pch = 21, fill = NA, size = 5, colour = "black", stroke = 2)+
  xlab("")+
  ylab("Cell Type (%)")+
  scale_y_continuous(breaks = seq(0,100,10),
                     limits = c(0,100))+ #set the y axis from 0-100 labeled every 10 units
  theme(text = element_text(size=20))+ #makes text bigger
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))+ #change axis text color
  theme(panel.background = element_rect(fill = "white", #this is the background
                                colour = "black", #this is the border
                                size = 0.1, linetype = "solid"))+
  theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
dev.off()

###############################################################################################
#################### ENCODE-Reinius estimates and comparison to ewastools stats ###############
###############################################################################################

# Estimate cell proportions in whole samples using ENCODE-Reinius reference panel in ewastools
library(tidyverse)
library(ewastools)

#load datasets
setwd("~/Research 2019")
beta_final <- readRDS("06-11-19 beta_final_all_samps.rds")
pd_final <- readRDS("06-11-19 pd_final_all_samps.rds")


#subset pd_final to only include the whole samples
pd_samps <- pd_final[pd_final$celltype == "whole", ]
dim(pd_samps) #18 64
#add sample_id to pd_samps
sample_id <- sub("*SAL", "", pd_samps$Sample_Group)
pd_samps$sample_id <- sample_id

meth_id_whole <- as.vector(pd_samps$meth_id)

#subset beta_final to only include the whole samples using the meth_ids
beta_samps <- beta_final[ , meth_id_whole]
dim(beta_samps)

est <- as.data.frame(estimateLC(meth = beta_samps,
                                ref = "encode_reinius",
                                constrained = TRUE))

#change the proportions to percents - and flip estimates because my ref is backwards
est$Immune <- est$Epi * 100
est$Epithelial <- est$IC * 100

est <- est %>%
  select(-Epi, -IC)

#add sample_id based on the pd_final
est$id <- pd_samps$id

stats_whole_encode <- est %>%
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
  )

est_long <- gather(est,
                   est_celltype,
                   estimates,
                   Epithelial:Immune,
                   factor_key=TRUE)


# ## Calculate average epithelial/immune cells in whole samples based on estimates
# ## Calculate epithelial cells
# mean_epith_epi <- mean(est$Epi)
# median_epith_epi <- median(est$Epi)
# sd_epith_epi <- sd(est$Epi)
# iqr_epith_epi <- IQR(est$Epi)
# mean_epith_epi #65.4
# median_epith_epi #71.8
# sd_epith_epi #19.5
# iqr_epith_epi #25.5
# 
# ## Calculate immune cells
# mean_ic_epi <- mean(est$IC)
# median_ic_epi <- median(est$IC)
# sd_ic_epi <- sd(est$IC)
# iqr_ic_epi <- IQR(est$IC)
# mean_ic_epi #35.8
# median_ic_epi #29.9
# sd_ic_epi #20.0
# iqr_ic_epi #24.6

# Compare the estimates of each epithelial sample for ewastools and ENCODE/Reinius using a paired ttest
ttest <- t.test(whole_estimates$Epithelial, est$Epithelial,
                paired = TRUE, alternative = "two.sided")
ttest
# Paired t-test
# 
# data:  whole_estimates$Epithelial and est$Epithelial
# t = -7.4725, df = 17, p-value = 9.122e-07
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -33.81513 -18.92440
# sample estimates:
#   mean of the differences 
# -26.36976




# # Compare the estimates of each immune sample for ewastools and ENCODE/Reinius using a paired ttest
# ttest <- t.test(whole_estimates$Immune, est$IC,
#                 paired = TRUE, alternative = "two.sided")
# ttest
# # Paired t-test
# # 
# # data:  whole_estimates$Immune and est$IC
# # t = 7.8259, df = 17, p-value = 4.918e-07
# # alternative hypothesis: true difference in means is not equal to 0
# # 95 percent confidence interval:
# #   16.50806 28.69438
# # sample estimates:
# #   mean of the differences 
# # 22.60122


###############################################################################################
################### Make plot 3: immune vs epith % by sample in encode-reinius ################
###############################################################################################

library(tidyverse)

#make the dataset long instead of wide
enc_rein_long <- gather(est, est_celltype, estimates,
                         Epithelial:Immune, factor_key=TRUE)
#est_celltype: column with Immune or Epithelial label
#estimates: the numbers


## Graph the whole samples estimated by encode-reinius
setwd("C:/Users/HP/Documents/Research 2020")
png("07-30-20 whole estimates encode-reinius and saliva.png", width=500, height=500)
ggplot(data = enc_rein_long, aes(x = est_celltype, y = estimates, group = id))+
  geom_point(aes(color = as.factor(id)), size = 3)+
  geom_line(aes(color = as.factor(id)), size = 1)+
  xlab("")+
  ylab("Cell Type (%)")+
  scale_y_continuous(breaks = seq(0,100,10),
                     limits = c(0,100))+ #set the y axis from 0-100 labeled every 10 units
  theme(text = element_text(size=20))+ #makes text bigger
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))+ #change axis text color
  theme(panel.background = element_rect(fill = "white", #this is the background
                                colour = "black", #this is the border
                                size = 0.1, linetype = "solid"))+
  theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
dev.off()