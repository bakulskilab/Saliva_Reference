setwd("~/Research 2019")

# Load files

#encode beta
beta_encode <- readRDS("07-03-19 beta_final.rds")

#encode pd
pd_encode <- readRDS("07-03-19 pd.paper.rds")

#reinius pd
setwd("~/Research 2020")
pd_reinius <- readRDS("07-09-20 pd_reinius.rds")

#reinius beta
setwd("~/Research 2019")
beta_reinius <- readRDS("07-31-19 beta_final_wbc.rds")


#BiocManager::install("FlowSorted.Blood.450k")



# Find which probes overlap between 450k and 850k
overlapping_probes <- Reduce(intersect, list(rownames(beta_encode),
                                             rownames(beta_reinius)))
length(overlapping_probes) #443574


# Cut betamatrix to only the ones that intersect - subsetting to overlapping probes
#keeps the order of the probes the same
beta_encode_overlap <- beta_encode[overlapping_probes, ]
beta_reinius_overlap <- beta_reinius[overlapping_probes, ]
dim(beta_encode_overlap) #443574     11 (11 epithelial)
dim(beta_reinius_overlap) #443574     42
identical(rownames(beta_reinius_overlap), rownames(beta_encode_overlap)) #TRUE
sum(is.na(beta_reinius_overlap)) #0 NA
sum(is.na(beta_encode_overlap)) #0 NA

# Merge the two datasets
beta_overlap <- cbind(beta_encode_overlap, beta_reinius_overlap)
beta_reinius_overlap <- as.data.frame(beta_reinius_overlap)
beta_encode_overlap <- as.data.frame(beta_encode_overlap)
dim(beta_overlap) #443574    53 (11+42)
colnames(beta_overlap)


setwd("~/Research 2020/")
saveRDS(beta_overlap, file = "07-09-20 beta_encode_reinius.rds")





############################## Combine reinius and encode pds ###############################
library(data.table)
library(tidyverse)

setwd("~/Research 2020")
reinius <- fread("Reinius_sample_sheet_IDAT.csv") #downloaded
setwd("~/Research 2019")
encode <- readRDS("pd.reduced.rds") #generated in code file 5.0
setwd("~/Research 2020")
beta_wbc_encode <- readRDS("07-09-20 beta_encode_reinius.rds") #generated in this code file

#combine reinius and encode pds
View(reinius)
View(encode)
colnames(encode)
#ENCODE
encode1 <- encode %>%
  select(c(geo_accession, source_name_ch1))
encode1 <- rename(encode1, meth_id = geo_accession)
encode1$source_name_ch1 <- as.character(encode1$source_name_ch1)
encode2 <- encode1 %>%
  filter(source_name_ch1 %in% c("HIPEpiC",
                                "SAEC",
                                "HRE",
                                "HAEpiC",
                                "HRPEpiC",
                                "PrEC",
                                "HEEpiC",
                                "HCPEpiC",
                                "HNPCEpiC",
                                "HMEC",
                                "HRCEpiC")) %>%
  add_column(cell_type = "Epi") %>%
  select(c(meth_id, cell_type))
dim(encode2) #11 2
colnames(encode2)
pd_encode <- encode2

#################################################################################

#Reinius
reinius1 <- reinius %>%
  select(c(Sample, Type)) %>%
  add_column(meth_id = "")

reinius1$meth_id <- gsub(pattern = " ",
                         replacement = "_",
                         x = reinius1$Sample)
reinius1$cell_type <- "IC"


#keep only the sorted cell types
reinius2 <- reinius1 %>%
  select(c(meth_id, Type, cell_type)) %>%
  filter(Type %in% c("CD14+ Monocytes",
                     "CD56+ NK-cells",
                     "CD8+ T-cells",
                     "CD4+ T-cells",
                     "CD19+ B-cells",
                     "Eosinophils",
                     "Neutrophils")) %>%
  select(c(meth_id, cell_type))
dim(reinius2)
pd_reinius <- reinius2
colnames(pd_reinius)

#bind the two datasets together
pd_merge <- rbind(pd_encode, pd_reinius)
dim(pd_merge) #53 2
table(pd_merge$cell_type)
# Epi  IC 
# 11  42 

##################################################################################################
# irrelevant due to re-processing encode data and removing fibroblasts then
##################################################################################################

#remove the fibroblasts from beta_wbc_encode
# dim(beta_wbc_encode) #441401     60 (42 ic, 11 epi, 7 fib)
# colnames(beta_wbc_encode)
# beta_wbc_encode <- as.matrix(beta_wbc_encode)
# 
# beta_colnames <- colnames(beta_wbc_encode)
# 
# gsm <- grep("^GSM", beta_colnames, value = TRUE)
# gsm #length = 11 epi + 7 fib
# 
# #figure out which gsms are for fibroblasts
# #pd_encode only has epithelial cell types in it
# library(usefun)
# fib <- outersect(pd_encode$meth_id, gsm)
# length(fib) #7
# fib #"GSM999340" "GSM999342" "GSM999344" "GSM999345" "GSM999348" "GSM999350" "GSM999394"
# 
# #drop the fibroblast columns
# beta_wbc_epi = subset(beta_wbc_encode,
#                       select = -c(GSM999340,
#                                   GSM999342,
#                                   GSM999344,
#                                   GSM999345,
#                                   GSM999348,
#                                   GSM999350,
#                                   GSM999394))
# dim(beta_wbc_epi) #441401     53

beta_wbc_epi <- beta_wbc_encode
dim(beta_wbc_epi) #443574     53

#save the matrix that includes reinius and epithelial cells
setwd("~/Research 2020")
saveRDS(beta_wbc_epi, "07-09-20 beta_reinius_epi.rds")

################################### t test ######################################################

# Make dummy vars for cell types;
pd_merge$IC <-ifelse(pd_merge$cell_type =="IC", 1, 0)
pd_merge$Epi <-ifelse(pd_merge$cell_type =="Epi", 1, 0)
table(pd_merge$IC)

# Subset beta_final into cell types;
# check that order meth_ID in pd file = columns in beta matrix;
identical(colnames(beta_wbc_epi), pd_merge$meth_id) #TRUE

# Run row ttests;
library(genefilter)
#do a row ttest of one cell types against all other cell types within one probe;
ttest_ic <- rowttests(as.matrix(beta_wbc_epi), as.factor(pd_merge$cell_type), tstatOnly = FALSE)


# How many probes are significant at 1e-8 for rowttest of one cell vs other two;
## doing lower p value because larger chance of getting <0.05
length(which(ttest_ic$p.value<1e-8)) #64329 sites are different between encode and reinius
length(which(ttest_ic$p.value<0.05)) #301127 - 3/4 of the total probes
