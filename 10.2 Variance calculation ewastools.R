# VARIANCE ESTIMATE FOR EWASTOOLS BASED ON QC BETA MATRIX

setwd("C:/Users/HP/Documents/Research 2019")

##################################################################################################
################################# INSTALL PACKAGES AND LOAD THEM #################################
##################################################################################################


# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("minfi")

library(minfi)
#library(BeadSorted.Saliva.EPIC)
library(tidyverse)
load("centEpiFibIC.m.rda")
library(broom)


##################################################################################################
#################### LOAD AND SUBSET BETA MATRIX TO ONLY INCLUDE WHOLE SAMPLES ###################
##################################################################################################


#Load the beta matrix QC'd by Lauren and John
beta_final_all_samps <- readRDS("06-11-19 beta_final_all_samps.rds")
# str(beta_final_all_samps) #rows are cgs, columns are meth_ids


#Subset the matrix down to only the Whole samples using the meth_ids
beta_samps <- beta_final_all_samps[ , c("203282450208_R02C01",
                                        "203282450208_R05C01",
                                        "203282450209_R08C01",
                                        "203286230053_R03C01",
                                        "203286230053_R08C01",
                                        "203286230103_R01C01",
                                        "203286230103_R02C01",
                                        "203286230120_R02C01",
                                        "203286230120_R05C01",
                                        "203286230120_R06C01",
                                        "203287590025_R03C01",
                                        "203287590048_R05C01",
                                        "203287590228_R05C01",
                                        "203287590243_R01C01",
                                        "203287590243_R05C01",
                                        "203293440009_R03C01",
                                        "203293440009_R04C01",
                                        "203293440009_R07C01")]
dim(beta_samps) #795694     18


##################################################################################################
############################ USE EWASTOOLS TO ESTIMATE THE CELL PROPORTIONS ######################
##################################################################################################


tools = estimateLC(beta_samps, ref = "saliva")


##################################################################################################
############################### RUN THE MODEL TO GENERATE THE R2S ################################
##################################################################################################


start_time <- print(Sys.time())


counter <<- 1
dish.out <- apply(beta_samps[1:5000, ], 1, function(probs){
  print(counter)
  df <- data.frame(methyl = probs, tools)
  # print(str(df))
  model <- lm(methyl ~ ., data = df)
  #print(class(model))
  counter <<- counter + 1
  return(model)
})
dim_beta_samps <- dim(beta_samps)
print(Sys.time())

#Compile the results of the regression
varestimates <- list(ewastools_saliva = dish.out)

beta_subset <- beta_samps


##################################################################################################
################################## PULL THE R2S OUT OF THE STATS #################################
##################################################################################################


run_regression_model_on_cell_types <- function(probs
                                               , stats_of_interest)
{
  probs_updated <- probs %>%
    select(-prob_names) %>% #drops the probe names column so the dimensions match
    t(.) #transforms the dataset 90 degrees
  # print(probs_updated)
  
  df_tools <- data.frame(methyl = probs_updated, tools) #makes this a data frame
  # print(df_tools)
  model <- lm(methyl ~ ., data = df_tools)
  
  #Separate out the coefficients from the other statistics
  if(stats_of_interest == "coefficients")
  {
    df_interest <- tidy(model)
  } else if(stats_of_interest == "other_regression_stats") {
    df_interest <- glance(model)
  }
  print(df_interest)
  return(df_interest)
}
print(Sys.time())


#Summaries is a list of stats for each probe including R2
summaries <- lapply(varestimates, function(x){
  print(counter)
  # print(names(x))
  lapply(x, function(y){
    # print(names(y))
    summary(y)
  })
  #counter <<- counter + 1
})
print(Sys.time())


#pull the r squareds into one dataset
rsquareds <- lapply(summaries, function(method){
  lapply(method, function(probe){
    probe$r.squared
    # print(probe$r.squared)
  }) %>% unlist
}) %>% bind_cols
#rsquareds
print(Sys.time())


##################################################################################################
############################### USE THE R2 VALUES TO CALCULATE STATS #############################
##################################################################################################


rsquareds <- as.data.frame(rsquareds)

#Save the R2 dataset
write.csv(rsquareds, "07-03-20 rsquared values.csv")

#Calculate some stats
mean <- mean(rsquareds$ewastools_saliva)
median <- median(rsquareds$ewastools_saliva)

#Make a histogram of the R2 values
pdf("07-03-20 rsquared histogram.pdf")
hist(rsquareds$ewastools_saliva)
dev.off()


end_time <- print(Sys.time())