#######################################################################################################
#######################################################################################################
# title: "Adding Reinius and ENCODE to estimateLC"
# author: "Jonah Fisher, Lauren Middleton"
# date: "07-02-20"
#######################################################################################################
#######################################################################################################

#######################################################################################################
# Purpose:  generate the list of probes to use as the reference panel for epidish/encode cell types
#
# Inputs:   "Reinius_sample_sheet_IDAT.csv" - descriptives of Reinius samples
#           IDATs
# and encode files
#
# Outputs:  "encode_reinius.txt" - list of probes for epidish/encode reference panel
#######################################################################################################


setwd("C:/Users/HP/Documents/Research 2020")

#######################################################################################################
########################################## Load the libraries #########################################
#######################################################################################################

library(minfi)
library(ewastools)
library(tidyverse)
library(forcats)
library(stringi)
library(data.table)
library(GEOquery)

#######################################################################################################

# We want to add Reinius and ENCODE data sets used in epidish to estimateLC() function in ewastools
#   as a reference in order to compare "apples to apples" style with our saliva reference panel.


# Creating the beta matrices to plug into algorithm
## Create beta1 using the Reinius data
### Download idats from GSE

#######################################################################################################

#This code is straight from Jonathans github for reinius.
# It is based on data downloaded from reinius dropbox folder
reinius = fread("Reinius_sample_sheet_IDAT.csv")
dim(reinius)
#60 samples 5 columns
head(reinius)
#      Sample        Type    Chip#ID Chip CO position Chip Row Pos
# 1:   WB 105 Whole blood 5684819001              CO1          RO1
# 2:   WB 218 Whole blood 5684819001              CO1          RO2
# 3:   WB 261 Whole blood 5684819001              CO1          RO3
# 4: PBMC 105        PBMC 5684819001              CO1          RO4

#add a column for just the donor ID, grab the last 3 characters of Sample
reinius[,donor:=stri_sub(Sample,-3,-1)]

#the zeroes in chip row pos are actually letters so make the letters into numbers
reinius[,`Chip Row Pos`:=stri_replace_all_fixed(`Chip Row Pos`,"O","0")]
reinius[,`Chip CO position`:=stri_replace_all_fixed(`Chip CO position`,"O","0")]

#combine the chipID, row, and column into the meth_id, add idat/ in front and call the column "file"
reinius[,file:=paste0(`Chip#ID`,"_",`Chip Row Pos`,`Chip CO position`)] 
# reinius[,file:=paste0(`Chip#ID`,"_",`Chip Row Pos`,`Chip CO position`)] 

#add a column for cell_type based on the Type column
reinius$cell_type = fct_recode(reinius$Type, MO="CD14+ Monocytes",
                                             NK="CD56+ NK-cells",
                                             CD8="CD8+ T-cells",
                                             CD4="CD4+ T-cells",
                                             GR="Granulocytes",
                                             B="CD19+ B-cells",
                                             WB="Whole blood",
                                             EO="Eosinophils",
                                             NE="Neutrophils")

#add a column that says that all sex is male
reinius$sex = "m"

#kept only the columns cell_type, donor, sex, and file. added study: reinius
reinius = reinius[,list(study="Reinius",cell_type,donor,sex,file)]

#keep only the sorted cell types, changed to include eosinophils
reinius = reinius[cell_type %in% c("MO","NK","CD8","GR","CD4","B", "EO")] 
dim(reinius)
#42 (6 donors, 7 cell types) 5 columns

## JONATHAN EXCLUDES THESE SAMPLES AND I DONT KNOW WHY BUT THE PAPER LOOKS LIKE IT USES
##    ALL 7 CELLTYPES FROM EACH OF THE 6 DONORS
# reinius = reinius[donor != "105"]
# reinius = reinius[! (donor=="160" & cell_type %in% c("CD8","CD4"))]
# reinius = reinius[! (donor=="261" & cell_type %in% c("NK","B"))]
# reinius = reinius[! (donor=="218" & cell_type=="NK")]

#group cell type into just immune cells
reinius$cell_type <- factor("IC", levels = "IC")

setwd("C:/Users/HP/Documents/Research 2019/reinius_idat_files")

#create a beta matrix: read the idats, detection p mask -2, correct for dye bias, don't normalize
beta1 = dont_normalize(correct_dye_bias(ewastools::mask(detectionP.neg(read_idats(reinius$file)),-2)))
dim(beta1)
#485577     42

#######################################################################################################
#######################################################################################################
#######################################################################################################

#######################################################################################################
######################################## Create ENCODE Dataset ########################################
#######################################################################################################

#######################################################################################################
######################################## Create ENCODE pd File ########################################
#######################################################################################################

#load descriptives for encode data
setwd("C:/Users/HP/Documents/Research 2019")
encode <- readRDS("pd.reduced.rds")
dim(encode) #62 4

#select only the cell types included in epidish and 
encode$cell_type <- fct_recode(encode$source_name_ch1,
                               Epi = "HIPEpiC", 
                               Epi = "SAEC", 
                               Epi = "HRE", 
                               Epi = "HAEpiC", 
                               Epi = "HRPEpiC", 
                               Epi = "PrEC", 
                               Epi = "HEEpiC", 
                               Epi = "HCPEpiC", 
                               Epi = "HNPCEpiC", 
                               Epi = "HMEC", 
                               Epi = "HRCEpiC"
                               # Fib = "IMR90",
                               # Fib = "ProgFib", 
                               # Fib = "AG04449", 
                               # Fib = "AG09319",  
                               # Fib = "AG04450", 
                               # Fib = "BJ",
                               # Fib = "NHDF-neo"
                                )

#keep only the epithelial cells
encode <- filter(encode, cell_type %in% c("Epi"))

#######################################################################################################
###################################### Create ENCODE Beta Matrix ######################################
#######################################################################################################

#download ENCODE IDAT files
setwd("C:/Users/HP/Documents/Research 2019/ENCODE idats")

#make beta matrix from encode idats
beta2 = dont_normalize(correct_dye_bias(ewastools::mask(detectionP.neg(read_idats(encode$meth_id)),-2)))

#######################################################################################################
#######################################################################################################
#######################################################################################################

#######################################################################################################
###################################### Create ENCODE/Reinius Beta #####################################
#######################################################################################################

# Combine beta matrices from reinius and encode

# setwd("C:/Users/HP/Documents/Research 2020")
# common = paste0("beta", 1:2) %>%
#   map(get) %>%
#   map(rownames) %>%
#   reduce(intersect)
# 
# common = intersect(common,
#                    ewastools:::manifest_450K[!chr%in%c("X","Y") & probe_type=="cg"]$probe_id)
# 
# beta = cbind(
#   beta1[ match(common,rownames(beta1)) ,]
#   ,beta2[ match(common,rownames(beta2)) ,]
# )
# 
# beta = beta[common,]

#load beta matrix from code file 5.5
setwd("C:/Users/HP/Documents/Research 2019")
beta <- readRDS("07-09-20 beta_reinius_epi.rds")

cell_types = c(
  as.character(reinius$cell_type),
  as.character(encode$cell_type)
)


## Code from algorithm
train_model = function(output){
  
  train_beta = beta
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
train_model("encode_reinius.txt")

#manually add encode_reinius.txt data library of ewastools