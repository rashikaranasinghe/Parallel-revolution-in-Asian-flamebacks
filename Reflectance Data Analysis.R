#####################################################################
## Reflectance_data_analysis_Asian_Flamebacks_script.R
# started by Rashika W. Ranasinghe 17 Sept 2024

# This file contains R scripts for analyzing reflectance spectrometric data, including 
# 1. generating reflectance spectra, 
# 2. extracting colorimetric variables,
# 3. visualizing colorimetric variables, and 
# 4. performing statistical analysis on these variables, 
# as used in the manuscript titled: '*********'."


############################################################################
################## Plot reflectance spectra ############################
############################################################################
## set the working directory
setwd("/OneDrive-UBC/Irwin Lab/My Research/carotenoid analysis_I/HPLC/Flameback-reflectance data/Reflectance analysis project")

## Load the libraries 
library(pavo)
library(stringr)

# I took three replicate reflectance readings from each sample. All the data from mantle feathers were organized into separate folders labeled by species. Be sure to maintain separate folders for crown feathers and mantle feathers.

################################################################################
######## Make the reflectance spectrum for mantle feathers ################

########  C.stricklandi #########
x <- getspec("C.stricklandi_reflectance-mantle/", ext = "txt", decimal = ".")
x_mean <- aggspec(x,FUN = mean) # To get the mean of replicates
quartz(height = 3.5, width=4)
aggplot(x, col = "red3",FUN.error = function(x) sd(x) / sqrt(length(x))) # just to check the plots
title(expression(italic("C. stricklandi") ~ "- mantle"))

########  C.haematribon #########
x <- getspec("C.haematribon_reflectance-mantle/", ext = "txt", decimal = ".")
x_mean <- aggspec(x,FUN = mean)
quartz(height = 3.5, width=4)
aggplot(x, col = "red3",FUN.error = function(x) sd(x) / sqrt(length(x))) # just to check the plots
title(expression(italic("C.haematribon")~ "- mantle"))

########  C.xanthocephalus #########
x <- getspec("C.xanthocephalus_reflectance-mantle/", ext = "txt", decimal = ".")
x_mean <- aggspec(x,FUN = mean)
quartz(height = 3.5, width=4)
aggplot(x, col = "red3",FUN.error = function(x) sd(x) / sqrt(length(x))) # just to check the plots
title(expression(italic("C.xanthocephalus")~ "- mantle"))

########  C.l.rufopunctatus #########
x <- getspec("C.l.rufopunctatus_reflectance_mantle/", ext = "txt", decimal = ".")
x_mean <- aggspec(x,FUN = mean)
quartz(height = 3.5, width=4)
aggplot(x, col = "red3",FUN.error = function(x) sd(x) / sqrt(length(x))) # just to check the plots
title(expression(italic("C.l.rufopunctatus")~ "- mantle"))

########  D.psarodes #########
x <- getspec("D.psarodes_Allo_R_txt/", ext = "txt", decimal = ".")

x_mean <- aggspec(x,FUN = mean)
quartz(height = 3.5, width=4)
aggplot(x, col = "red3",FUN.error = function(x) sd(x) / sqrt(length(x))) # just to check the plots
title(expression(italic("D.psarodes")~ "- mantle"))
######## END Making the reflectance spectra for mantle  feathers ################
###################################################################################


###############################################################################
######## Make the reflectance spectrum for crown feathers ################
        
########  C.stricklandi #########
x <- getspec("C.stricklandi_crown_txt/", ext = "txt", decimal = ".")
x_mean <- aggspec(x,FUN = mean) # to get the mean of replicates
quartz(height = 3.5, width=4)
aggplot(x, col = "red3",FUN.error = function(x) sd(x) / sqrt(length(x))) 
title(expression(italic("C. stricklandi") ~ "- crown"))

########  C.haematribon #########
x <- getspec("C.haematribon_crown/", ext = "txt", decimal = ".")
x_mean <- aggspec(x,FUN = mean)
quartz(height = 3.5, width=4)
aggplot(x, col = "red3",FUN.error = function(x) sd(x) / sqrt(length(x))) # just to check the plots
title(expression(italic("C.haematribon")~ "- crown"))

########  C.xanthocephalus #########
x <- getspec("C.xanthocephalus_crown/", ext = "txt", decimal = ".")
x_mean <- aggspec(x,FUN = mean)
quartz(height = 3.5, width=4)
aggplot(x, col = "red3",FUN.error = function(x) sd(x) / sqrt(length(x))) # just to check the plots
title(expression(italic("C.xanthocephalus")~ "- crown"))

########  C.l.rufopunctatus #########
x <- getspec("C.l.rufopunctatus_crown/", ext = "txt", decimal = ".")
x_mean <- aggspec(x,FUN = mean)
quartz(height = 3.5, width=4)
aggplot(x, col = "red3",FUN.error = function(x) sd(x) / sqrt(length(x))) # just to check the plots
title(expression(italic("C.l.rufopunctatus")~ "- crown"))

########  C.l.rufopunctatus #########
x <- getspec("C.l.rufopunctatus_reflectance_crown/", ext = "txt", decimal = ".")
x_mean <- aggspec(x,FUN = mean)
quartz(height = 3.5, width=4)
aggplot(x, col = "red3",FUN.error = function(x) sd(x) / sqrt(length(x))) # just to check the plots
title(expression(italic("C.l.rufopunctatus")~ "- crown"))

######## D.psarodes #########
x <- getspec("D.psarodes_crown/", ext = "txt", decimal = ".")
x_mean <- aggspec(x,FUN = mean)
quartz(height = 3.5, width=4)
aggplot(x, col = "red3",FUN.error = function(x) sd(x) / sqrt(length(x))) # just to check the plots
title(expression(italic("D.psarodes")~ "- crown"))

######## END making the reflectance spectrum for crown feathers ################
###############################################################################

############################################################################
############### END plotting reflectance spectra ############################
############################################################################



########################################################################################
############### Calculate colorimatric variables using pavo ############################
########################################################################################
# For this analysis, ensure that all reflectance data from all species are placed in a single folder. Keep separate folders for crown and mantle feathers. Each sample has three replicates, and the mean of the replicates was taken before extracting the colorimetric data.

########################### For Mantle feathers #####################
# read the data file 
mantle_Dino_Chrys <- getspec("Reflectance_red_plumage_birds-mantle/", ext = "txt", decimal = ".")

# replace all negative values with zeros
mantle_Dino_Chrys_nonNeg <- procspec(mantle_Dino_Chrys, fixneg = "zero") 
plot(mantle_Dino_Chrys_nonNeg)

# take the average of every 3 spectrum because I have taken 3 replicates from each sample.
dim(mantle_Dino_Chrys)
mantle_Dino_Chrys_mean <- aggspec(mantle_Dino_Chrys_nonNeg, by = 3, FUN = mean)
dim(mantle_Dino_Chrys_mean)

# Extract the colorimatric data and save it as a .csv file
write.csv(summary(mantle_Dino_mean), file = "Colorimatric.Data-Dinopium_mantle_from_pavo.csv")



######################### for crown feathers ######################
# red teh data file
mantle_Dino_Chrys <- getspec("Reflectance_red_plumage_birds-crown/", ext = "txt", decimal = ".")
        
# replace all negative values with zeros
mantle_Dino_Chrys_nonNeg <- procspec(mantle_Dino_Chrys, fixneg = "zero") plot(mantle_Dino_Chrys_nonNeg)

# take the average of every 3 spectrum because I have taken 3 replicates from each sample.
dim(mantle_Dino_Chrys)
mantle_Dino_Chrys_mean <- aggspec(mantle_Dino_Chrys_nonNeg, by = 3, FUN = mean)
dim(mantle_Dino_Chrys_mean)

# Extract the colorimatric data and save it as a .csv file
write.csv(summary(mantle_Dino_mean), file = "Colorimatric.Data-Dinopium_crown_from_pavo.csv")
########################################################################################
############### Calculate colorimatric variables using pavo ############################
########################################################################################





########################################################################################################
############### Visualize colorimatric parameters (H3, catotenoid chroma and brightness) ###############
########################################################################################################
# set the working directory
setwd("/Users/rashikaranasinghe/Library/CloudStorage/OneDrive-UBC/Irwin Lab/My Research/carotenoid analysis_I/HPLC/Flameback-reflectance data/Reflectance analysis project")

# Lodad the libraries
library(tidyverse)
library(tidyr)

# Read the data file 
colorimatrics <- read.csv("Colorimatric.Data-Dino_Chrys_crown_and_mantle_from_pavo.csv", stringsAsFactors = T)
dim(colorimatrics)
names(colorimatrics)

# Manipulate the data
Data_long <- colorimatrics %>%
  pivot_longer(cols = c(B2,S9, Reflectance.average.Lambda.R50),
               names_to = "Variable",
               values_to = "Values")
levels(Data_long$Scientific.names)
Data_long$Species <- factor(Data_long$Scientific.names, levels = c("C.haematribon", "C.stricklandi", "C.l.rufopunctatus", "C.xanthocephalus","D.psarodes"))
Data_long$Variable <- factor(Data_long$Variable, levels = c("Reflectance.average.Lambda.R50", "B2","S9"))

# make the plot
quartz(height = 3.5, width = 9)
ggplot(Data_long, aes(y = Values, x = Scientific.names, fill=Scientific.names)) +
  geom_jitter(position = position_jitter(0.3), size = 2, shape = 21, color = "black", alpha = 0.6) +
  stat_summary(fun.data = "mean_cl_normal",  fun.args = list(mult = 2.5), 
               geom = "pointrange",  size = 0.3, shape = 23, fill = "black")+
  facet_wrap(~feathers+Variable, ncol = 3, scales = "free",
             labeller = labeller(Variable=c("Reflectance.average.Lambda.R50"= "Hue", 
                                            "B2"="Mean brightness",
                                            "S9"="Carotenoid Chroma"))) +
  scale_fill_manual(values = c("cyan3", "green4", "blue", "red", "purple"),
                    name = NULL,
  ) +
  ylab("Intensity") +
  xlab("Species")+
  scale_x_discrete(labels = NULL) +
  theme_linedraw() +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "white"),
    strip.background = element_rect(color = "grey50", fill = "grey50"),  
    strip.text = element_text(face = "bold",  margin = margin(1, 10, 1, 10)),
    panel.spacing = unit(0.2, "lines"),
    legend.position = "top",
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold")) +
  ggtitle("Colorimatric variables (for 3 ovelapped feathers)") 
########################################################################################################
############## END visualizing colorimatric parameters (H3, catotenoid chroma and brightness) ############
########################################################################################################


###############################################################################################
########################## Statistical analysis of colorimatric variables #####################
###############################################################################################
## Conduct a Fisher-Pitman permutation test to assess whether there are significant differences in each variable between species. If the permutation test indicates significant differences among groups, perform post-hoc analysis with pairwise tests to identify which group pairs show significant differences.

# Load the lirbaries
library(coin)
library(rcompanion)
        
## Fisher-Pitman permutation test
oneway_test(S9~ Scientific.names, data=filter(colorimatrics, feathers == "mantle"), distribution= approximate(nresample = 1000000))
#####################
# data:  B2 - mantle by
# Scientific.names (C.haematribon, C.l.rufopunctatus, C.stricklandi, C.xanthocephalus, D.psarodes)
# chi-squared = 14.968, p-value = 0.000856.
#####################

## post-hoc analysis
PT <- pairwisePermutationTest(B2 ~ Scientific.names, 
                              data=filter(colorimatrics, feathers == "crown"),teststat = "quadratic",method   = "fdr")
PT
###############################################################################################
####################### END Statistical analysis of colorimatric variables ####################
###############################################################################################
