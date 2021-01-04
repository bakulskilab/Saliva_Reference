#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
###############################  LOAD AND PREPARE THE DNAm DATA FOR PROCESSING  ###############################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose: This code file loads the GSE138279 data (child saliva samples)
#          
# Inputs:   GSE138279 data
#
# Outputs:  ewastools meth list and minfi RGset, pd_child_clean for sample information

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#############################################################################################################
###################################### Load Libraries And Set Directory #####################################
#############################################################################################################

library(GEOquery)
library(tidyverse)
library(minfi)
library(ewastools)

current_directory <- "C:/Users/T7920/Desktop/Lauren/Saliva - Peds dataset"
idat_directory <- "C:/Users/T7920/Desktop/Lauren/Saliva - Peds dataset/GSE138279_IDATS"
setwd(current_directory)

#############################################################################################################
######################################## Download DNAm Data And Pd Data #####################################
#############################################################################################################

#download dataset
child_saliva_raw <- getGEO('GSE138279',
                           GSEMatrix = TRUE)
#parsing errors, but the pd looks fine

pd_child_raw <- pData(phenoData(child_saliva_raw[[1]]))

#############################################################################################################
############################################## Clean The Pd Data ############################################
#############################################################################################################

#keep the useful columns
pd_child_simplified <- pd_child_raw[,c('geo_accession',
                                       "title",
                                       'gender:ch1',
                                       'characteristics_ch1.2',
                                       "supplementary_file",
                                       "supplementary_file.1")]
head(pd_child_simplified)

#turn the supplementary file column into something resembling the IDAT file names
pd_child_simplified$supplementary_file[1]
# ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3035nnn/GSM3035989/suppl/GSM3035989_201530450003_R01C01_Grn.idat.gz

#remove the ftp://... at start
pd_child_simplified$meth_id <- gsub('ftp://.*suppl/',
                         '',
                         pd_child_simplified$supplementary_file)
head(pd_child_simplified$meth_id)

#Remove _Grn.idat.gz from each file
pd_child_simplified$meth_id <- gsub('_Grn.idat.gz',
                                    '',
                                    pd_child_simplified$meth_id)
head(pd_child_simplified$meth_id)

#drop supplementary_file column
pd_child_clean <- pd_child_simplified %>%
  select(-supplementary_file)

#############################################################################################################
############################################ Download the DNAm Data #########################################
#############################################################################################################

# Download raw data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138279

#idats are stored in .tar file, so open it to extract the idats
untar("GSE138279_RAW.tar", exdir = "GSE138279_IDATS")

#manually removed the csv file - not sure if this is necessary

#############################################################################################################
############################################### Read In The IDATS ###########################################
#############################################################################################################

setwd(idat_directory)

#reading into ewastools
meth <- read_idats(pd_child_clean$meth_id)
#read in 622399 probes (duplicates for type 1 probes red and green)

class(meth)
#list

class(pd_child_clean)
#data.frame

#reading into minfi
RGset <- read.metharray(basenames = pd_child_clean$meth_id,
                        verbose = T,
                        force = TRUE)
#have to include force = TRUE line for some reason related to the EPIC array:
# https://support.bioconductor.org/p/97773/


