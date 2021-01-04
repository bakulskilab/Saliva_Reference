#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
###############################  LOAD AND PREPARE THE DNAm DATA FOR PROCESSING  ###############################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose: This code file loads the GSE111631 data (adult saliva samples)
#          
# Inputs:   GSE111631 data
#
# Outputs:  ewastools meth list and minfi RGset, pd_adult_clean for sample information

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#############################################################################################################
###################################### Load Libraries And Set Directory #####################################
#############################################################################################################

library(GEOquery)
library(tidyverse)
library(minfi)
library(ewastools)

current_directory <- "C:/Users/T7920/Desktop/Lauren/Saliva - Adult dataset"
idat_directory <- "C:/Users/T7920/Desktop/Lauren/Saliva - Adult dataset/GSE111631_IDATS"
setwd(current_directory)

#############################################################################################################
######################################## Download DNAm Data And Pd Data #####################################
#############################################################################################################

#download dataset
adult_saliva_raw <- getGEO('GSE111631',
                           GSEMatrix = TRUE)

pd_adult_raw <- pData(phenoData(adult_saliva_raw[[1]]))
# saveRDS(pd_adult_raw, file='pd_adult_raw.rds')

#############################################################################################################
############################################## Clean The Pd Data ############################################
#############################################################################################################

#keep the useful columns
pd_adult_simplified <- pd_adult_raw[,c('geo_accession',
                                       'source_name_ch1',
                                       'characteristics_ch1.2',
                                       'characteristics_ch1.4',
                                       "supplementary_file",
                                       "supplementary_file.1")]
head(pd_adult_simplified)

#turn the supplementary file column into something resembling the IDAT file names
pd_adult_simplified$supplementary_file[1]
# ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3035nnn/GSM3035989/suppl/GSM3035989_201530450003_R01C01_Grn.idat.gz

#remove the ftp://... at start
pd_adult_simplified$meth_id <- gsub('ftp://.*suppl/',
                         '',
                         pd_adult_simplified$supplementary_file)
head(pd_adult_simplified$meth_id)

#Remove _Grn.idat.gz from each file
pd_adult_simplified$meth_id <- gsub('_Grn.idat.gz',
                                    '',
                                    pd_adult_simplified$meth_id)
head(pd_adult_simplified$meth_id)

#drop supplementary_file column
pd_adult_clean <- pd_adult_simplified %>%
  select(-supplementary_file)

#############################################################################################################
############################################ Download the DNAm Data #########################################
#############################################################################################################

# Download raw data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111631

#idats are stored in .tar file, so open it to extract the idats
untar("GSE111631_RAW.tar", exdir = "GSE111631_IDATS")

#manually removed the csv file - not sure if this is necessary

#############################################################################################################
############################################### Read In The IDATS ###########################################
#############################################################################################################

setwd(idat_directory)

#reading into ewastools
meth <- read_idats(pd_adult_clean$meth_id)
#read in 1051815 probes (duplicates for type 1 probes red and green)

class(meth)
#list

class(pd_adult_clean)
#data.frame

#reading into minfi
RGset <- read.metharray(basenames = pd_adult_clean$meth_id,
                        verbose = T,
                        force = TRUE)
#have to include force = TRUE line for some reason related to the EPIC array:
# https://support.bioconductor.org/p/97773/