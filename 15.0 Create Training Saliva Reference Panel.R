#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
###############################  CREATE REFERENCE PANEL FOR TRAINING/TESTING  #################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose:  This code file creates a reference panel based on the training saliva fractions
#          
# Inputs:   "06-11-19 beta_final.rds" - beta matrix of sorted saliva samples
#           "06-11-19 pd_final.rds"   - demographics file for samples
#
# Outputs:  "saliva_training.txt"     - list of probes for subsetted saliva reference panel

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#############################################################################################################
############################################### Load Datasets ###############################################
#############################################################################################################

# setwd("C:/Users/HP/Documents/Research 2019")
# beta_final <- readRDS("06-11-19 beta_final.rds")
# pd_final <- readRDS("06-11-19 pd_final.rds")
# beta_all_samps <- readRDS("06-11-19 beta_final_all_samps.rds")
# pd_all_samps <- readRDS("06-11-19 pd_final_all_samps.rds")

# Load global environment: w - training_testing_global_environment.RData

#############################################################################################################
###################################### Load Libraries And Set Directory #####################################
#############################################################################################################

# current_directory <- "C:/Users/HP/Google Drive/Colacino Lab/Saliva/Training-Testing Analysis"

setwd(current_directory)

library(minfi)
library(ewastools)
library(tidyverse)
library(forcats)
library(stringi)
library(data.table)

#############################################################################################################
########################################### Select Random Samples ###########################################
#############################################################################################################

# 20 immune samples and 18 epithelial samples, select 70% of each for training

0.7*20
# 14 immune

0.7*18
# 12 epithelial

#split the pd into immune and epithelial dataframes
pd_epith <- pd_final %>%
  filter(celltype == "large")
pd_immune <- pd_final %>%
  filter(celltype == "CD45pos")

#randomly sample from each dataframe - initial sampling
# train_pd_epith <- pd_epith %>%
#   sample_n(size = 12,
#            replace = FALSE)
# train_pd_immune <- pd_immune %>%
#   sample_n(size = 14,
#            replace = FALSE)

#identify which samples were chosen from each cell type for future replication
# train_epith_samps <- train_pd_epith$id
# # 3  9  7 10  6 12  5  1 18  2  4 11
# train_immune_samps <- train_pd_immune$id
# # 19 24 23  4 10  7  5 20  6 18  3 21 15 12

#replication sampling - epithelial
epith_ids <- c(3, 9, 7, 10, 6, 12, 5, 1, 18, 2, 4, 11)
train_pd_epith <- pd_epith %>%
  filter(id %in% epith_ids)

#replication sampling - immune
immune_ids <- c(19, 24, 23, 4, 10, 7, 5, 20, 6, 18, 3, 21, 15, 12)
train_pd_immune <- pd_immune %>%
  filter(id %in% immune_ids)

#bind the training datasets back together
training_pd <- bind_rows(train_pd_epith,
                         train_pd_immune)

#############################################################################################################
############################################## Subset The Beta ##############################################
#############################################################################################################

#identify samples based on meth_id
meth_id <- training_pd$meth_id

#change the beta matrix to df
beta_final_df <- as.data.frame(beta_final)

#subset the beta matrix based on the meth_ids
beta_training <- beta_final_df %>%
  select(all_of(meth_id))

#clean up
rm(pd_epith, pd_immune, train_pd_epith, train_pd_immune, beta_final_df)

#############################################################################################################
########################################## Clean Up The Pd Dataset ##########################################
#############################################################################################################

#rename cd45pos and large labels for cell types
training_pd <- training_pd %>%
  mutate(celltype = ifelse(celltype == "large",
                           yes = "Epi",
                           no = "IC")
         )
table(training_pd$celltype)
# Epi  IC 
# 12   14

#create a cell type vector
cell_types = as.character(training_pd$celltype)

#############################################################################################################
########################################### Create Reference Panel ##########################################
#############################################################################################################

# Code from algorithm
train_model = function(input, output)
{
  train_beta = input
  markers = list() 
  
  for(ct in c("Epi", "IC"))
  { 
    cat(ct,'\n') 
    
    j = cell_types == ct
    
    tmp = apply(train_beta,
                1,
                function(x)
                  { 
                    if(!any(is.na(x[j]))) 
                    {  
                      tmp = t.test(x[j],
                                   x[!j],
                                   var.equal=T) 
                      return(c(tmp$p.value,
                               tmp$estimate[1] - tmp$estimate[2])
                             ) 
                    }else
                      { 
                      return(c(NA,NA)) 
                      } 
                  }
                )
    
    i = which(p.adjust(tmp[1,]) <0.05) 
    o = order(tmp[2,i], na.last=NA)
    markers[[ct]] = c(i[head(o,50)], i[tail(o,50)]) 
  } 
  
  markers %<>% unlist %>% unique 
  
  coefs = sapply(c("Epi", "IC"),
                 function(ct)
                   { 
                    j = cell_types == ct
                    rowMeans(train_beta[markers,j],na.rm=TRUE) 
                   }
                 )
  rownames(coefs) = rownames(train_beta)[markers] 
  colnames(coefs) = c("Epi", "IC")
  
  write.table(coefs, file = output, row.names=TRUE) 
}

#only change this line - output is sent to current_directory
train_model(input = beta_training,
            output = "saliva_training.txt")
# train_model(input = beta_training,
#             output = "saliva_training_constrained.txt")

#put saliva_training manually into ewastools package ("C:/Users/HP/Documents/R/win-library/3.6/ewastools/data")
#putting saliva_training_constrained.txt into ("C:\Users\HP\Documents\R-dev\ewastools\data")