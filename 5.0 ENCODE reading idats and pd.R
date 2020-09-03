---
title: "06-25-19 ENCODE Data Pre-process"
author: "Lauren Middleton"
date: "6/25/2019"
output: html_document
---
# Purpose:  Preprocess the data from this paper in the same way the saliva 850k data was done:
#           https://www.futuremedicine.com/doi/10.2217/epi-2018-0037?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub=www.ncbi.nlm.nih.gov&
#
# Inputs:   download ENCODE datasets from the internet
#
# Outputs:  methylation dataset and RGset for QC

# Set up
library(GEOquery)
setwd("C:/Users/HP/Documents/Research 2019")

# Download dataset and descriptive dataset
# Make the pd file
#download dataset
encode <- getGEO('GSE40699',GSEMatrix = TRUE)
#downloaded sample descriptive data

pd.encode <- pData(phenoData(encode[[1]]))
saveRDS(pd.encode,file='pd.encode.rds')

#only keep a few columns
pd <- pd.encode[,c('geo_accession','source_name_ch1','characteristics_ch1.6','supplementary_file')]
head(pd)


# Remove file path from ftp: to suppl/
#we can turn the supplementary file column into something resembling the file names, which look like
# GSM999335_hg19_wgEncodeHaibMethyl450SkmcSitesRep1_Grn.idat

pd$supplementary_file[1]
#ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM999nnn/GSM999335/suppl/GSM999335_hg19_wgEncodeHaibMethyl450SkmcSitesRep1_Grn.idat.gz

#remove the ftp://... stuff at start
pd$meth_id <- gsub('ftp://.*suppl/','',pd$supplementary_file)
head(pd$meth_id)
# [1] "GSM999335_hg19_wgEncodeHaibMethyl450SkmcSitesRep1_Grn.idat.gz"  
# [2] "GSM999336_hg19_wgEncodeHaibMethyl450HipeSitesRep1_Grn.idat.gz"  
# [3] "GSM999337_hg19_wgEncodeHaibMethyl450Helas3SitesRep1_Grn.idat.gz"
# [4] "GSM999338_hg19_wgEncodeHaibMethyl450Hepg2SitesRep1_Grn.idat.gz" 
# [5] "GSM999339_hg19_wgEncodeHaibMethyl450HepatoSitesRep1_Grn.idat.gz"
# [6] "GSM999340_hg19_wgEncodeHaibMethyl450Imr90SitesRep1_Grn.idat.gz"


# Remove _Grn.idat.gz from each file
#now remove the _Grn.idat.gz
pd$meth_id <- gsub('_Grn.idat.gz','',pd$meth_id)
head(pd$meth_id)
# [1] "GSM999335_hg19_wgEncodeHaibMethyl450SkmcSitesRep1"   "GSM999336_hg19_wgEncodeHaibMethyl450HipeSitesRep1"  
# [3] "GSM999337_hg19_wgEncodeHaibMethyl450Helas3SitesRep1" "GSM999338_hg19_wgEncodeHaibMethyl450Hepg2SitesRep1" 
# [5] "GSM999339_hg19_wgEncodeHaibMethyl450HepatoSitesRep1" "GSM999340_hg19_wgEncodeHaibMethyl450Imr90SitesRep1"


# Remove supplementary_file column from dataset
#get rid of old supplementary_file column because its too long
pd <- pd[,-4]
setwd("C:/Users/HP/Documents/Research 2019")
saveRDS(pd,file='pd.reduced.rds')


# Download raw data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40699
#GSE40699_RAW.tar download link at bottom of page
#extract the contents into a folder
#make sure code matches folder you have files in

#idats are zipped, but packages don't like that,
#here I list out all the downloaded idats, then unzip the idats
idatFiles <- list.files(path="C:/Users/HP/Documents/Research 2019/ENCODE idats",
                        pattern = "idat.gz",
                        full = TRUE)

sapply(idatFiles, gunzip, overwrite = TRUE)


#read in idats in ewastools and minfi
library(minfi)
library(ewastools)


#reading into ewastools
setwd("C:/Users/HP/Documents/Research 2019/ENCODE idats")
meth <- read_idats(pd$meth_id) #read in 622399 probes (duplicates for type 1 probes red and green)
# head(meth)
class(meth) #list
class(pd) #data.frame

#reading into minfi
setwd("C:/Users/HP/Documents/Research 2019/ENCODE idats")
RGset <- read.metharray(basenames=pd$meth_id,verbose=T)



setwd("C:/Users/HP/Documents/Research 2019")
saveRDS(meth,file="07-02-19 ENCODE ewastools-meth.rds")
saveRDS(RGset,file="07-02-19 ENCODE RGset.rds")