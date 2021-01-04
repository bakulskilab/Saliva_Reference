#######################################################################################################
#######################################################################################################
# title: "Adding saliva to estimateLC"
# author: "Lauren Middleton"
# date: "07-03-20"
#######################################################################################################
#######################################################################################################

#######################################################################################################
# Purpose:  generate the list of probes to use as the reference panel for saliva cell types
#
# Inputs:   "06-11-19 beta_final.rds" - beta matrix of sorted saliva samples
#           "06-11-19 pd_final.rds"   - demographics file for samples
#
# Outputs:  "saliva.txt" - list of probes for saliva reference panel
#######################################################################################################


setwd("C:/Users/HP/Documents/Research 2020")


################################### Load the libraries ################################################

library(minfi)
library(ewastools)
library(tidyverse)
library(forcats)
library(stringi)
library(data.table)

#######################################################################################################

# We want to create list of 100 saliva probes to add to estimateLC() function in ewastools

################################## Load the beta matrix ###############################################

setwd("C:/Users/HP/Documents/Research 2019")
beta_final <- readRDS("06-11-19 beta_final.rds")
pd_final <- readRDS("06-11-19 pd_final.rds")

#######################################################################################################

pd_final$celltype_1 <- gsub(pattern = "CD45pos", replacement = "IC", pd_final$celltype)
pd_final$celltype_recode <- gsub(pattern = "large", replacement = "Epi", pd_final$celltype_1)
table(pd_final$celltype_recode)
  
cell_types = c(as.character(pd_final$celltype_recode))
#######################################################################################################

## Code from algorithm
train_model = function(input, output){
  
  train_beta = input
  markers = list() 
  
  for(ct in c("Epi", "IC")){ 
    cat(ct,'\n') 
    
    j = cell_types == ct
    
    tmp = apply(train_beta,1,function(x){ 
      if(!any(is.na(x[j]))) 
      {  
        tmp = t.test(x[j],x[!j],var.equal=T) 
        return(c(tmp$p.value,tmp$estimate[1]-tmp$estimate[2])) 
      }else{ 
        return(c(NA,NA)) 
      } 
    }) 
    
    i = which(p.adjust(tmp[1,])<0.05) 
    o = order(tmp[2,i],na.last=NA)
    markers[[ct]] = c(i[head(o,50)],i[tail(o,50)]) 
  } 
  
  markers %<>% unlist %>% unique 
  
  coefs = sapply(c("Epi", "IC"),function(ct){ 
    j = cell_types == ct
    rowMeans(train_beta[markers,j],na.rm=TRUE) 
  }) 
  rownames(coefs) = rownames(train_beta)[markers] 
  colnames(coefs) = c("Epi", "IC")
  
  write.table(coefs,file=output,row.names=TRUE) 
}

setwd("C:/Users/HP/Documents/Research 2020")
train_model(input = beta_final, output = "saliva.txt") #only change this line



#test practice
# train_model = function(input, output){
  
  train_beta = beta_final[1,]
  markers = list() 
  
  for(ct in c("Epi", "IC")){ 
    cat(ct,'\n') 
    
    j = cell_types == ct
    
    tmp = apply(train_beta,1,function(x){ 
      if(!any(is.na(x[j]))) 
      {  
        tmp = t.test(x[j],x[!j],var.equal=T) #this x is the first row of beta_final
        return(c(tmp$p.value,tmp$estimate[1]-tmp$estimate[2])) 
      }else{ 
        return(c(NA,NA)) 
      } 
    }) 
    
    i = which(p.adjust(tmp[1,])<0.05) 
    #orders temp, takes second row: what is the second row
    o = order(tmp[2,i],na.last=NA)
    markers[[ct]] = c(i[head(o,50)],i[tail(o,50)]) 
  } 
  
  markers %<>% unlist %>% unique 
  
  coefs = sapply(c("Epi", "IC"),function(ct){ 
    j = cell_types == ct
    rowMeans(train_beta[markers,j],na.rm=TRUE) 
  }) 
  rownames(coefs) = rownames(train_beta)[markers] 
  colnames(coefs) = c("Epi", "IC")
  
  write.table(coefs,file=output,row.names=TRUE) 
}


#######################################################################################################


intersect(rownames(saliva), rownames(ewastools_saliva_probes))

#manually add saliva.txt data library of ewastools