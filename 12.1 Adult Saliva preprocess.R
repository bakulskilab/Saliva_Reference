#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##############################  PREPROCESS THE DNAm DATA AND CREATE BETA MATRIX  ##############################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose:  This code file processes the raw data from GSE111631
#          
# Inputs:   RGset          - object containing DNAm data compiled using minfi 
# 
#           meth           - list containing DNAm data compiled using ewastools
#           
#           pd_adult_clean - dataframe containing information about the saliva samples
#
# Outputs:  "w - adult_saliva_global_environment" - includes beta matrix and pd file

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#############################################################################################################
###################################### Load Libraries And Set Directory #####################################
#############################################################################################################

library(tidyverse)
library(minfi)
library(ewastools)
library(magrittr)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

current_directory <- "C:/Users/T7920/Desktop/Lauren/Saliva - Adult dataset"
control_metrics_directory <- "C:/Users/T7920/Desktop/Lauren/Saliva - Adult dataset/Control Metrics"
setwd(current_directory)

#run 12.0 code first

#############################################################################################################
################################# Calculate QC metrics based on Illumina Guide ##############################
#############################################################################################################

setwd(control_metrics_directory)

#this has all of the control metrics for the array
ctrl.metrics <- control_metrics(meth)

#make a PDF of the control metric graphs
pdf("12-03-20 control_metrics.pdf")
for(metric in names(ctrl.metrics)) {
  stripchart(ctrl.metrics[[metric]],
             method = "jitter",
             pch = 4, #symbol for points
             xlab=metric)
  abline(v = attr(ctrl.metrics[[metric]],
                  "threshold"),
         col=2, #this makes the line red
         lty=3) } #this makes the line type dotted
dev.off()


#plot only one control metric at a time and set xlim
#https://hhhh5.github.io/ewastools/articles/exemplary_ewas.html
#staining red metric
stripchart(ctrl.metrics$"Staining Red",
           method="jitter",
           pch = 4,
           xlab = "Staining Red",
           xlim = c(0, 150))
abline(v = 5, col = 2, lty = 3) #lty = line type

#staining green metric
stripchart(ctrl.metrics$"Staining Green",
           method="jitter",
           pch = 4,
           xlab = "Staining Green",
           xlim = c(0, 2500))
abline(v = 5, col = 2, lty = 3)

#extension green metric
stripchart(ctrl.metrics$"Extension Green",
           method = "jitter",
           pch = 4,
           xlab = "Extension Green",
           xlim = c(0, 100))
abline(v = 5, col = 2, lty = 3)

#############################################################################################################
############################### Flag ctrl.metrics For Any Metric Threshold Fail #############################
#############################################################################################################

#flag failures as TRUE if they fail on any of the 17 control metrics
fail.flag <- sample_failure(ctrl.metrics)
class(fail.flag)
#logical
table(fail.flag)
#13 failed, 155 passed

#############################################################################################################
############################### Create a Wide Dataset To Understand the Metrics #############################
#############################################################################################################

class(ctrl.metrics)
#list

#unlist
ctrl_metrics_unlist <- unlist(ctrl.metrics)

#change into a dataframe
ctrl_metrics_df <- as.data.frame(ctrl_metrics_unlist)

#change the rownames into a column
ctrl_metrics_tidy <- rownames_to_column(ctrl_metrics_df,
                                        var = "metric")

#add a sample id
ctrl_metrics_tidy$sample_id <- 1:168

#remove the sample ids from the metrics
ctrl_metrics_tidy_clean <- ctrl_metrics_tidy
ctrl_metrics_tidy_clean$metric <- gsub('[[:digit:]]+',
                                      '',
                                      ctrl_metrics_tidy$metric)

#Target Removal also has numbers in it (1, 2) so I need to put those back in
metric_list <- c("Restoration",
                 "Staining Green",
                 "Staining Red",
                 "Extension Green",
                 "Extension Red",
                 "Hybridization High/Medium",
                 "Hybridization Medium/Low",
                 "Target Removal 1",
                 "Target Removal 2",
                 "Bisulfite Conversion I Green",
                 "Bisulfite Conversion I Red",
                 "Bisulfite Conversion II",
                 "Specificity I Green",
                 "Specificity I Red",
                 "Specificity II",
                 "Non-polymorphic Green",
                 "Non-polymorphic Red"  )

#create a repeating list of the metrics
metrics <- rep(metric_list, each = 168)

#add that list as a column and remove the old column
ctrl_target_fix <- ctrl_metrics_tidy_clean %>%
  mutate(metrics = metrics) %>%
  select(-metric)

#create a wide dataset
metrics_wide <- pivot_wider(data = ctrl_target_fix,
                            id_cols = "sample_id",
                            names_from = "metrics",
                            values_from = "ctrl_metrics_unlist"
                            )

#############################################################################################################
################################## Which Samples Failed And On Which Metrics ################################
#############################################################################################################

#checking through pdf generated for figures that have samples below the threshold:
# Staining Red, #staining efficiency (background to noise?) - 5
# Extension Red, #efficiency of adding ACTG nucleotides - 5
# Bisulfite Conversion I Green, - 1
# Bisulfite Conversion I Red, - 1
# Bisulfite Conversion II, - 1
# Non−polymorphic Red #overall assay performance - 5

staining_fail <- which(x = metrics_wide$`Staining Red` <= 5)
extension_fail <- which(x = metrics_wide$`Extension Red` <= 5)
bisulfite_1_green <- which(x = metrics_wide$`Bisulfite Conversion I Green` <= 1)
bisulfite_1_red <- which(x = metrics_wide$`Bisulfite Conversion I Red` <= 1)
bisulfite_2 <- which(x = metrics_wide$`Bisulfite Conversion II` <= 1)
nonpolymorphic_red <- which(x = metrics_wide$`Non-polymorphic Red` <= 5)

failed_samples <- c(staining_fail,
                    extension_fail,
                    bisulfite_1_green,
                    bisulfite_1_red,
                    bisulfite_2,
                    nonpolymorphic_red)
table(failed_samples)
# 3   6  71  84  88  98 109 114 135 137 138 145 164 
# 1   1   2   1   2   3   2   3   3   2   1   3   1

#Summary: samples 98, 114, 135, and 145 failed on 3 metrics

#############################################################################################################
#################################### Reorder the pd file to match the meth ##################################
#############################################################################################################

#make a meth dataset that has sample_id with only the GSM###### as the ID
colnames(meth)
#NULL
rownames(pd_adult_clean)
#GSMs
head(meth$meta$sample_id)
#"GSM3035989_201530450003_R01C01"

#remove the row and column info from the sample id in the meth dataset to get the gsm
meth_gsm <- gsub('_.*', '', meth$meta$sample_id)
head(meth_gsm)

#add meth_gsm into meth dataset
meth$meta$gsm <- meth_gsm

#reorder the datasets so the orders match up by GSM
pd.flag <- pd_adult_clean[match(meth$meta$gsm, pd_adult_clean$geo_accession), ]
identical(pd.flag$geo_accession, meth$meta$gsm)
#TRUE

pd.flag$control_probe_flag <- fail.flag

#############################################################################################################
################################## Alternative Way to ID Which Samples Fail #################################
#############################################################################################################

#Which specific metrics fail
failed <- sapply(ctrl.metrics, function(metric) {
  metric < attr(metric, "threshold") 
})

#replace characters like spaces with underscores in failed
colnames(failed) %>% gsub(' |/|-', '_', .)
head(failed) #now R won't break from spaces in column labels

#combine pd.flag and failed by columns
head(pd.flag)
head(failed)
pd.flag <- cbind(pd.flag, failed)
head(pd.flag)

#turned failed dataset into a data frame
class(failed) #df
failed <- data.frame(failed)
class(failed) #data.frame

#add fail.flag as a new column called Any into failed dataset
failed$Any <- fail.flag

#count the column sums from failed
counts <- colSums(failed, na.rm = TRUE)
counts

#make a csv file of the control flag counts - metrics that failed
write.csv(counts, file = "12-03-20 adult control.flag.counts.csv")

#Any means 13 samples failed
# Restoration               Staining.Green                 Staining.Red              Extension.Green 
# 0                            0                            1                            0 
# Extension.Red    Hybridization.High.Medium     Hybridization.Medium.Low             Target.Removal.1 
# 1                            0                            0                            0 
# Target.Removal.2 Bisulfite.Conversion.I.Green   Bisulfite.Conversion.I.Red      Bisulfite.Conversion.II 
# 0                            5                            1                            8 
# Specificity.I.Green            Specificity.I.Red               Specificity.II        Non.polymorphic.Green 
# 0                                 0                            0                            0 
# Non.polymorphic.Red                          Any 
# 9                                             13

# Which samples were flagged
#from pd.flag, combine geo_accession with failed
# samp.flag <- pd.flag[pd.flag$control_probe_flag,c('Sample.ID',colnames(failed)[-18])]
colnames(failed) #18 doesn't exist so I'm going to ignore that piece

#colnames of failed are the control metrics
samp.flag <- pd.flag[pd.flag$control_probe_flag,]
#samples GSM999339 and GSM999368 were flagged - hepatocyte and prostate adenocarcinoma

samp.flag <- t(samp.flag)
samp.flag

write.csv(samp.flag, file = "12-03-20 adult control_probe_flagged_samples.csv")

#############################################################################################################
############################################ Correct For Dye Bias ###########################################
#############################################################################################################

meth <- correct_dye_bias(meth)

#check pd file lines up with meth object again
identical(rownames(pd.flag), meth$meta$gsm) #TRUE

#############################################################################################################
################################################# Sex Checks ################################################
#############################################################################################################

#X and Y normalized intensities
pd.sex <- pd.flag

#make 2 new columns for X and Y
pd.sex[ ,c('X', 'Y')] <- check_sex(meth)
#this gives X and Y median intensities
head(pd.sex[ ,c('X', 'Y')])

#fix the pd column for sex
pd_adult_clean <- pd_adult_clean %>%
  dplyr::rename(sex = characteristics_ch1.4)

pd_adult_clean$sex <- gsub("gender: ",
                           "",
                           pd_adult_clean$sex)

#sexes of samples based on pd.flag data
sex <- pd_adult_clean$sex
  
pd.sex$sex <- sex
colnames(pd.sex)

#predicted sex - lowercase
pd.sex$predicted_sex <- predict_sex(X = pd.sex$X, Y = pd.sex$Y,
                                    male = which(pd.sex$sex == "M"),
                                    female = which(pd.sex$sex == "F"))

#based on the pd file information
table(pd.sex$sex)
# F   M 
# 24 144

sex_predictions <- table(pd.sex$sex, pd.sex$predicted_sex,
                         useNA = "always")
#        m   f  <NA> #predicted
# F      1  23    0
# M    142   0    2
# <NA>   0   0    0
write.csv(sex_predictions, "12-03-20 sex predictions vs reported.csv")

#summary: one person who was reported female predicts as male, two people reported male don't predict


#############################################################################################################
############################################ Plot Predicted Sexes ###########################################
#############################################################################################################

#red means predicted doesn't match with annotated including unknowns(NA)
#if F, plot a circle etc - this is setting up the shapes
pd.sex_plot <- pd.sex %>%
  mutate(sex =
           case_when(sex == "M" ~ "m",
                     sex == "F" ~ "f")
         )

#red means predicted doesn't match with annotated including unknowns(NA)
#if F, plot a circle etc - this is setting up the shapes
tch <- ifelse(is.na(pd.sex_plot$predicted_sex), 0, #if NA, plot a square
              ifelse(pd.sex_plot$predicted_sex == "f", 1, 4)) #if F, plot circle; if M, plot X

pdf("12-03-20 sex_plot_ewastools.pdf")
#plot X vs Y
plot(Y ~ X, data = pd.sex_plot,
     pch = tch, asp = 1,
     xlab = "Normalized X chromosome intensities",
     ylab = "Normalized Y chromosome intensities",
     ylim = c(0, 1.0),
     cex.lab = 1.3)
#diff = NA or mismatch
diff <- pd.sex_plot[is.na(pd.sex_plot$predicted_sex) | pd.sex_plot$sex != pd.sex_plot$predicted_sex,]
#plots red mismatches or NA on top of previous points
tch.diff <- ifelse(is.na(diff$predicted_sex), 0, 
                   ifelse(diff$predicted_sex == "f", 1, 4))
points(Y ~ X, data = diff,
       pch = tch.diff,
       col = 2) #this makes it red ontop of plot
legend("topright", title = "Predicted Sex", pch = c(1,4,0),legend = c("Female", "Male", "Unassigned"))
legend("topleft", c("Sex match", "Sex mismatch"), fill = c("black", "red"))
#text labels the points with their GSM ID
#text(Y ~ X, labels = pd.sex_plot$geo_accession, data = pd.sex_plot)
dev.off()

#############################################################################################################
################################################ Detection p ################################################
#############################################################################################################

#"This function identifies failed positions defined as both the methylated and unmethylated channel
# reporting background signal levels."
detP <- ewastools::detectionP(meth)

#make a dataset called ewastools.detP from detP on meth dataset
ewastools.detP <- detP$detP

#make a dataset called samps from the gsm names in meth
#this is the names of the samples
samps <- meth$meta$gsm
samps

#make a dataset called probes from the manifest$probe_id in detP/meth
probes <- detP$manifest$probe_id

colnames(ewastools.detP) <- samps
rownames(ewastools.detP) <- probes
print(ewastools.detP[1:6, 1:6])


#this is supposed to check if the mismatches change the graph
#0.05 or 0.01 is a good cutoff for detP based on this graph
pdf("12-03-20 adult detP_eval_ewastools.pdf")
#reported
eval_detP_cutoffs(detP,
                  males = which(pd.sex$sex == "M"),
                  females = which(pd.sex$sex == "F"))
#predicted
eval_detP_cutoffs(detP,
                  males = which(pd.sex$predicted_sex == "m"),
                  females = which(pd.sex$predicted_sex == "f"))
dev.off()

#############################################################################################################
############################################ Detp Failed Samples ############################################
#############################################################################################################
#detection p: measure of the quality of the sites? smaller value is better quality

# Failed samples based on detP
failedP <- ewastools.detP > 0.05

#Fraction of failed positions per sample
per.samp <- colMeans(failedP, na.rm=T) 
summary(per.samp)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.007473 0.034141 0.063358 0.087919 0.105369 0.553346 

#Fraction of failed samples per position
per.probe <- rowMeans(failedP, na.rm=T) 
summary(per.probe)
# Min. 1st Qu.     Median    Mean    3rd Qu.  Max    NA's 
# 0.00000 0.00000 0.01786 0.08790 0.05357 1.00000       3 

#How many samples had more than 10% of sites fail?
sum(per.samp > 0.1)
#46
# length(per.samp) #168
colnames(per.samp)
#NULL
#since there are no sample names in per.samp, assign names from "samps"
names(per.samp) <- samps
per.samp

#who are you
fail_ten_perc <- which(per.samp > 0.1)
fail_10_perc_samples <- per.samp[fail_ten_perc]
#46/168 GSMs listed (child saliva: 17/60)
write.csv(fail_10_perc_samples, "12-03-20 samples that fail at over ten perc sites.csv")

#How many positions failed >5% samples
sum(per.probe > 0.05, na.rm=T)
# 226,853 (ENCODE paper had about 6,000 sites fail, child saliva: 179,136)

pdf("12-03-20 adult saliva detP-hist.pdf")
#Fraction of Failed Positions Per Sample
hist(per.samp,breaks=20,
     xlab='Fraction of Failed Positions Per Sample',
     main='All Samples')
#Samples with <5% Failed Probes - 
hist(per.samp[per.samp<0.05],
     breaks=40,
     xlab='Fraction of Failed Positions Per Sample',
     main='Samples with <5% Failed Probes')
#Fraction of Failed Samples Per Position
hist(per.probe,breaks=20,
     xlab='Fraction of Failed Samples Per Position',
     main='All Probes')
#Fraction of Failed Samples Per Position < 30% - zooms into left side
hist(per.probe[per.probe<0.3],
     breaks=60,
     xlab='Fraction of Failed Samples Per Position',
     main='Probes with <30% Sampled Failed')
dev.off()

which(per.samp > 0.5)
# GSM3036118 GSM3036150 
# 114        135

pd.detp <- pd.sex

#this is adding a new var of % probes fail per sample onto pd.detp
#this is adding the previous info from the section above into pd.detp
#rownames makes sure per.samp will be ordered the same as pd.detp
pd.detp$ew_probe_fail_pct <- per.samp[rownames(pd.detp)]
rownames(pd.detp)

#############################################################################################################
############################################# Genotype Calling ##############################################
#############################################################################################################

#don't normalize the beta values from meth
#simplifies meth dataset - speeds up the code
beta <- dont_normalize(meth)

#snp probes labeled as rs
snps <- meth$manifest[probe_type == 'rs', index]

#dataset of only snp probes
snps <- beta[snps, ]

#this says which bases correspond
genotypes <- call_genotypes(snps, learn=FALSE)

#############################################################################################################
############################################# Genotype Checking #############################################
#############################################################################################################

#double check ordering
identical(colnames(genotypes$snps), rownames(pd.detp))
#FALSE - they're both GSMs

#remove the chip position from the column names
colnames(genotypes$snps) <- substring(colnames(genotypes$snps), 1, 10)
identical(colnames(genotypes$snps), pd.detp$geo_accession)
#TRUE

# gsm_order <- colnames(genotypes$snps)
# 
# #reorder genotypes so it matches the order of the pd
# pd_reorder <- pd.detp %>%
#   arrange(gsm_order)
# 
# identical(colnames(genotypes$snps), (pd_reorder$geo_accession))
# #TRUE

#check agreement of samples from same person based on snps
check.snp <- check_snp_agreement(genotypes,
                                 donor_ids = pd.detp$geo_accession,
                                 sample_ids = pd.detp$geo_accession)

#CHECK.SNP GIVES NULL (github says NULL means no conflicts found)


pdf("12-03-120 adult SNP Agreement.pdf")
ewastools:::agreement_(genotypes, pd.detp$geo_accession, pd.detp$geo_accession)
dev.off()
# saveRDS(genotypes, file = "adult_genotypes.rds")

#############################################################################################################
############################################# Checks With Minfi #############################################
#############################################################################################################

#############################################################################################################
######################################### Overall Intensity: M vs. U ########################################
#############################################################################################################

#pulling methylated and unmeth intensities from red and green intensities
rawMSet <- preprocessRaw(RGset)

dim(rawMSet)
#865859    168
#rownames are probes, colnames are samples

#M signal per probe, per sample
Meth <- getMeth(rawMSet)
Meth[1:5,1:5]

#U signal per probe, per sample
Unmeth <- getUnmeth(rawMSet)
Unmeth[1:5,1:5]

MQC <- log2(colMedians(Meth))
UQC <- log2(colMedians(Unmeth))

intensities <- data.frame(meth_id = colnames(Meth), MQC=MQC, UQC=UQC)
head(intensities)
#need to shorten meth_id to match gsm
head(pd.detp)
#has geo_accession

#cut off the ending of meth_id to make gsm in intensities
intensities_gsm <- substring(intensities$meth_id, 1, 10)
head(intensities_gsm)
intensities$meth_id <- intensities_gsm

#merge intensities into pd.detp
pd.beta <- merge(pd.detp, intensities, by.x = "geo_accession", by.y = "meth_id")
head(pd.beta)

#############################################################################################################
############################################## Plot of M vs U ###############################################
#############################################################################################################

library(RColorBrewer)
myColors <- brewer.pal(3, "Set1")

pdf("12-03-20 M vs U plot - all samples.pdf")
palette(myColors)
plot(pd.beta$UQC,
     pd.beta$MQC,
     main = "M vs. U QC",
     pch=16,
     xlab="Log2 Median Unmethylated Intensity",
     ylab="Log2 Median Methylated Intensity",
     cex.lab=1.2,
     cex.main=2)
text(pd.beta$UQC, pd.beta$MQC, labels = pd.beta$geo_accession, pos=4, cex=0.6)
dev.off()

pdf("12-03-20 M vs U plot - predicted_sex.pdf")
palette(myColors)
plot(pd.beta$UQC,
     pd.beta$MQC,
     col = pd.beta$predicted_sex,
     main = "M vs. U QC by Predicted Sex",
     pch=16,
     xlab="Log2 Median Unmethylated Intensity",
     ylab="Log2 Median Methylated Intensity",
     cex.lab=1.2,
     cex.main=2)
legend("topleft", levels(pd.beta$predicted_sex), fill = myColors)
dev.off()

#############################################################################################################
############################################# Raw Density Plot ##############################################
#############################################################################################################

#Raw density plot
beta.raw <- getBeta(rawMSet)
# saveRDS(beta.raw, "raw_beta_matrix.rds")

pdf("12-03-20 Density plot - Raw.pdf")
densityPlot(beta.raw,
            sampGroups = pd.beta$characteristics_ch1.4,
            main="Raw Beta by Reported Gender",
            xlim = c(0.0, 1.2),
            ylim = c(0.0, 5))
dev.off()

#############################################################################################################
############################################### Density Plots ###############################################
#############################################################################################################

#samples where more than 50% of sites fail
# GSM3036118 GSM3036150

beta.raw <- as.data.frame(beta.raw)

beta_bad_samps <- beta.raw %>%
  select(GSM3036118_201172570014_R02C01,
         GSM3036150_201172570024_R07C01)

beta_bad_samps <- as.matrix(beta_bad_samps)
beta.raw <- as.matrix(beta.raw)

pdf("12-03-20 Density plot - Raw bad samples.pdf")
densityPlot(beta_bad_samps,
            sampGroups = pd.beta$characteristics_ch1.4,
            main="Raw Beta Bad Samples",
            xlim = c(0.0, 1.2),
            ylim = c(0.0, 5))
dev.off()

pdf("12-03-20 Density plot - Raw bad samples.pdf")
densityPlot(beta.raw,
            sampGroups = pd.beta$geo_accession == "GSM3036150",
            main="Raw Beta Bad Sample1",
            xlim = c(0.0, 1.2),
            ylim = c(0.0, 5))
densityPlot(beta.raw,
            sampGroups = pd.beta$geo_accession == "GSM3036118",
            main="Raw Beta Bad Sample2",
            xlim = c(0.0, 1.2),
            ylim = c(0.0, 5))
dev.off()

#############################################################################################################
############################################ Run Noob Correction ############################################
#############################################################################################################

#noob correction
noob <- preprocessNoob(RGset, offset=15, dyeCorr=TRUE, verbose = TRUE)

#turn it into a beta matrix
#this has all the probes and all samples
beta <- getBeta(noob)
dim(beta) #865859    168

#############################################################################################################
############################################## Drop Bad Samples #############################################
#############################################################################################################

# Drop low detection p samples: this is the >3%

#keep samples that have a fail rate of less than 0.03
low_pdrop <- pd.beta[pd.beta$ew_probe_fail_pct <0.03, ]
#check the dimensions
dim(low_pdrop)
#keeping 34 samples

#fix the beta column names
colnames(beta) <- substring(colnames(beta), 1, 10)

#drop the samples from the beta based on the low_pdrop samples
beta_drop_samps <- beta[, low_pdrop$geo_accession]

#############################################################################################################
############################################## Drop Bad Probes ##############################################
#############################################################################################################

#############################################################################################################
################################################## Set Up ###################################################
#############################################################################################################

setwd(current_directory)

cross_probes <- read.csv("cross_probes_info_EPIC_Pidsley_2016.csv")

data(Locations) #load locations data from illumina;
head(Locations) #X and Y chromosomes are called "chrX" or "chrY"

#fix the beta column names - they're fine
# colnames(beta_drop_samps) <- substring(colnames(beta_drop_samps), 1, 10)

#############################################################################################################
############################################# Drop Failed Probes ############################################
#############################################################################################################

failedP <- detP$detP > 0.05
failedP[1:5, 1:5] #false = not failed

#add the column names (samples) onto detP
head(detP$meta$gsm)
head(detP$manifest$probe_id)

#assign gsm to column names of failedP and probe_id to rows
colnames(failedP) <- detP$meta$gsm
rownames(failedP) <- detP$manifest$probe_id
failedP[1:5, 1:5]

#recalculate per.probe (probe fails in 1+ samples)
per.probe <- rowMeans(failedP, na.rm = TRUE)
table(per.probe >0.05, exclude = NULL)
#    FALSE   TRUE   <NA> 
#   639062 226853      3

#this has only the failed probes- if a probe fails in >5% it gets dropped
detp_fail <- per.probe[per.probe >0.05]
#still 3 missings

#checking the number is the same;
length(detp_fail) #226856 (includes 3 missing)
head(detp_fail)

#get the probe names
detp_fail <- names(detp_fail)

#starting with this many probes
dim(beta_drop_samps) #865859    34

#keep anything that is not the failed row names;
beta_nodetp <- beta_drop_samps[!rownames(beta_drop_samps) %in% detp_fail, ]
dim(beta_nodetp) #639023     34

#############################################################################################################
######################################### Drop Cross Reactive Probes ########################################
#############################################################################################################

head(cross_probes)
#Probe is an actual column here as opposed to row names;

#Number of cross reactive probes
dim(cross_probes) #29233

#this drops all cross reactive probes
beta_nocross <- beta_nodetp[!rownames(beta_nodetp) %in% cross_probes$"ï..TargetID", ]
dim(beta_nocross) #639023    34

#############################################################################################################
######################################### Drop X/Y Chromosome Probes ########################################
#############################################################################################################

#make a dataset of only X and Y;
sex_probes <- Locations[Locations$chr == "chrX" | Locations$chr == "chrY", ] 
dim(sex_probes) #19627

#check first 5 rows and columns of beta_nocross dataset
beta_nocross[1:5, 1:5] 

#keeping rownames that is not in x and y probes;
beta_noXY <- beta_nocross[!rownames(beta_nocross) %in% rownames(sex_probes), ]

#check dimensions of dataset without x and y probes;
dim(beta_noXY) #628879    34

#############################################################################################################
####################################### Order The pd And beta Datasets ######################################
#############################################################################################################

#we already did the subset in the same order so this is essentially just renaming beta_noXY
beta_adult <- beta_noXY[ , low_pdrop$geo_accession]
dim(beta_adult) #628879    34
head(beta_adult)


#save the beta matrix - all probe info, rows = probes, columns = samples, data = % methylation
saveRDS(beta_adult, file = "12-03-20 beta_adult.rds")

#save the pd file (demographic data)
saveRDS(low_pdrop, file = "12-03-20 pd_adult.rds")

#############################################################################################################
################################### Remake Density Plot With Good Samples ###################################
#############################################################################################################

setwd("C:/Users/HP/Google Drive/Colacino Lab/Saliva/Adult Saliva/Control Metrics")

pdf("12-09-20 Density plot - Clean.pdf")
densityPlot(beta_adult,
            sampGroups = low_pdrop$characteristics_ch1.4,
            main="Beta by Reported Gender",
            xlim = c(0.0, 1.2),
            ylim = c(0.0, 5))
dev.off()

#############################################################################################################
##################################### Put Together The Global Environment ###################################
#############################################################################################################

#directories
current_directory <- "C:/Users/T7920/Desktop/Lauren/Saliva - Adult dataset"
control_metrics_directory <- "C:/Users/T7920/Desktop/Lauren/Saliva - Adult dataset/Control Metrics"

#pd file
pd_adult <- readRDS("12-03-20 pd_adult.rds")

#beta matrix
# beta.raw <- readRDS("raw_beta_matrix.rds")
beta_adult <- readRDS("12-03-20 beta_adult.rds")