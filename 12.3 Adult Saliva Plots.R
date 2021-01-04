#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
######################################  PLOT ADULT SAMPLES AS BY SAMPLE  ######################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose:  This code file estimates the proportion of cells in the adult saliva
#          
# Inputs:   adult_long - Long dataframe of cell proportion estimates from adult saliva
#
# Outputs:  Plots of cell proportions from adult samples

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#############################################################################################################
###################################### Load Libraries And Set Directory #####################################
#############################################################################################################

#NOTE: load w - adult_saliva_global_environment.RData to run this code

library(tidyverse)
library(ggplot2)

# setwd(current_directory)

#############################################################################################################
###################################### Spaghetti Plot Colored By Sample #####################################
#############################################################################################################

setwd("C:/Users/HP/Google Drive/Colacino Lab/Saliva/Adult Saliva/Spaghetti Plots")

png("12-04-20 saliva - adult saliva spaghetti plot.png", width = 500, height = 500)
ggplot(data = adult_long_saliva, aes(x = celltype, y = percent, group = GSM))+
  # geom_point(aes(color = GSM), size = 3)+
  # geom_jitter(aes(color = GSM), size = 3,
  #             position=position_jitter(0.2))+
  geom_line(aes(color = GSM,
                group = GSM), size = 1,
            position=position_jitter(0.2))+
  # geom_point(aes(color = GSM), size = 3,
  #             position=position_jitter(0.2))+
  xlab("")+
  ylab("Percent Cell Type (%)")+
  scale_y_continuous(breaks = seq(0,100,10),
                     limits = c(-1,100))+ #set the y axis from 0-100 labeled every 10 units
  theme(text = element_text(size=30))+ #makes text bigger
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))+ #change axis text color
  theme(panel.background = element_rect(fill = "white", #this is the background
                                        colour = "black", #this is the border
                                        size = 0.1, linetype = "solid"))+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(legend.position = "none")
dev.off()

png("12-04-20 encode - adult saliva spaghetti plot.png", width = 500, height = 500)
ggplot(data = adult_long_encode, aes(x = celltype, y = percent, group = GSM))+
  geom_point(aes(color = GSM), size = 3)+
  geom_line(aes(color = GSM), size = 1)+
  xlab("")+
  ylab("Percent Cell Type (%)")+
  scale_y_continuous(breaks = seq(0,100,10),
                     limits = c(-1,100))+ #set the y axis from 0-100 labeled every 10 units
  theme(text = element_text(size=30))+ #makes text bigger
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))+ #change axis text color
  theme(panel.background = element_rect(fill = "white", #this is the background
                                        colour = "black", #this is the border
                                        size = 0.1, linetype = "solid"))+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(legend.position = "none")
dev.off()

#############################################################################################################
####################################### Spaghetti Plot Colored By Sex #######################################
#############################################################################################################

setwd("C:/Users/HP/Google Drive/Colacino Lab/Saliva/Adult Saliva/Spaghetti Plots")

png("12-04-20 adult saliva spaghetti plot_by_sex.png", width = 500, height = 700)
ggplot(data = adult_long, aes(x = celltype, y = percent, group = celltype))+
  geom_point(aes(color = sex), size = 3)+
  geom_line(aes(color = sex), size = 1)+
  xlab("")+
  ylab("Percent Cell Type (%)")+
  scale_y_continuous(breaks = seq(0,120,10),
                     limits = c(-1,120))+ #set the y axis from 0-100 labeled every 10 units
  theme(text = element_text(size=30))+ #makes text bigger
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))+ #change axis text color
  theme(panel.background = element_rect(fill = "white", #this is the background
                                        colour = "black", #this is the border
                                        size = 0.1, linetype = "solid"))+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(legend.position = "none")
dev.off()

#############################################################################################################
######################################### Violin Plot By Epith/Immune #######################################
#############################################################################################################

setwd("C:/Users/HP/Google Drive/Colacino Lab/Saliva/Adult Saliva/Violin Plots")

png("12-04-20 saliva - adult saliva violin plot.png", width = 500, height = 700)
ggplot(data = adult_long_saliva, aes(x = celltype, y = percent, fill = celltype))+
  geom_violin()+
  geom_point(size = 2)+
  xlab("")+
  ylab("Percent Cell Type (%)")+
  # labs(title = "")+
  scale_y_continuous(breaks = seq(0,120,10),
                     limits = c(-1,120))+ #set the y axis from 0-120 labeled every 10 units
  theme(text = element_text(size=30))+ #makes text bigger
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))+ #change axis text color
  theme(panel.background = element_rect(fill = "white", #this is the background
                                        colour = "black", #this is the border
                                        size = 0.1, linetype = "solid"))+
  theme(legend.position="bottom")+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(legend.position = "none")
dev.off()

png("12-04-20 encode - adult saliva violin plot.png", width = 500, height = 700)
ggplot(data = adult_long_encode, aes(x = celltype, y = percent, fill = celltype))+
  geom_violin()+
  geom_point(size = 2)+
  xlab("")+
  ylab("Percent Cell Type (%)")+
  # labs(title = "")+
  scale_y_continuous(breaks = seq(0,120,10),
                     limits = c(-1,120))+ #set the y axis from 0-120 labeled every 10 units
  theme(text = element_text(size=30))+ #makes text bigger
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))+ #change axis text color
  theme(panel.background = element_rect(fill = "white", #this is the background
                                        colour = "black", #this is the border
                                        size = 0.1, linetype = "solid"))+
  theme(legend.position="bottom")+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(legend.position = "none")
dev.off()