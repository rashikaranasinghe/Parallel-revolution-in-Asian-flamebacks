#####################################################################
## Reflectance_data_analysis_Asian_Flamebacks_script.R
# started by Rashika W. Ranasinghe 17 Sept 2024

# This file includes R scripts for analysis reflectance spectrometric data () used in the manuscript titled: "*********".


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
