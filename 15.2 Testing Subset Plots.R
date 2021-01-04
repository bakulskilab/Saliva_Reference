#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#######################################  PLOT SORTED SAMPLES BY SAMPLE  #######################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose:  This code file estimates the proportion of cells in the sorted testing saliva samples
#          
# Inputs:   long_saliva - Long dataframe of cell proportion estimates from sorted saliva - saliva ref
#           long_encode - Long dataframe of cell proportion estimates from sorted saliva - encode ref
#
# Outputs:  Plots of cell proportions from sorted samples

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#############################################################################################################
###################################### Load Libraries And Set Directory #####################################
#############################################################################################################

#NOTE: load w - training_testing_global_environment to run this code

library(tidyverse)
library(ggplot2)

setwd(current_directory)

#############################################################################################################
######################################### Violin Plot By Epith/Immune #######################################
#############################################################################################################

long_saliva$est_celltype <- factor(long_saliva$est_celltype,
                                   levels = c("Immune", "Epithelial"))
long_encode$est_celltype <- factor(long_encode$est_celltype,
                                   levels = c("Immune", "Epithelial"))

png("12-11-20 saliva - testing subset saliva violin plot.png", width = 700, height = 700)
ggplot(data = long_saliva, aes(x = est_celltype,
                               y = percent,
                               fill = est_celltype))+
  geom_violin()+
  geom_point(size = 2)+
  xlab("")+
  ylab("Percent Cell Type (%)")+
  # labs(title = "")+
  scale_y_continuous(breaks = seq(0,100,10),
                     limits = c(-1,100))+ #set the y axis from 0-120 labeled every 10 units
  # scale_x_reverse()+
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
  theme(legend.position = "none")+
  facet_wrap(vars(celltype))
dev.off()

png("12-11-20 encode - testing subset saliva violin plot.png", width = 700, height = 700)
ggplot(data = long_encode, aes(x = est_celltype,
                               y = percent,
                               fill = est_celltype))+
  geom_violin()+
  geom_point(size = 2)+
  xlab("")+
  ylab("Percent Cell Type (%)")+
  # labs(title = "")+
  scale_y_continuous(breaks = seq(0,100,10),
                     limits = c(-1,100))+ #set the y axis from 0-120 labeled every 10 units
  # scale_x_reverse()+
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
  theme(legend.position = "none")+
  facet_wrap(vars(celltype))
dev.off()