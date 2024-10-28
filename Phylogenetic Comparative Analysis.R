##################################################################################
############### Phylogenetic generalized least squares (PGLS) ####################
#  to remove the effect of the evolutionary relationships of species when fitting a regression between two variables
# Depending on the model of evolution of the characters, the covariance matrix can be scaled using different approaches. For instance, one might assume that the character evolves under the Brownian motion model, or under an Ornstein-Uhlenbeck model where the co-variance between two species decreases exponentially according to a parameter alpha

library(nlme)

# Using BM model
## Read the molecular phylogenetic tree
tree <- read.nexus("Dino.Chrys.n60_40412SNPs_l10000.tree.nexusfile")
tree <- drop.tip(tree,"P.nahrattemsis") # Remove species in the tree that are not in the data matrix
## read the phenotype data 
p <- read.csv("Dino.Chrys.n60_40412SNPs_l10000.tree.pheno_Data.csv",row.names = 1)
pheno <- p[-3,]

# Get the correlation structure of the tree
bm.corr <- corBrownian(phy=tree)

# PGLS: Hue_mantle ~ C4-keto carotenoid
bm.pgls <- gls(Hue_mantle ~ Red_pig_conten_percentage , data = pheno, correlation = bm.corr)
summary(bm.pgls)

# PGLS: Hue_mantle ~ avg_4_keto_groups
bm.pgls <- gls(Hue_mantle ~ avg_4_keto_groups , data = pheno, correlation = bm.corr)
summary(bm.pgls)

# PGLS: Hue_mantle ~ No of epsilon rings
bm.pgls <- gls(Hue_mantle ~ avg_epsilon_rings , data = pheno, correlation = bm.corr)
summary(bm.pgls)

# PGLS: Hue_mantle ~ C3-oxygenated groups
bm.pgls <- gls(Hue_mantle ~ avg_3_hydroxylation , data = pheno, correlation = bm.corr)
summary(bm.pgls)

############### END Phylogenetic generalized least squares (PGLS) ####################
##################################################################################



###################################################################################
################################## Phylogenetic PCA ##############################
## Read the molecular phylogenetic tree
tree <- read.nexus("Dino.Chrys.n60_40412SNPs_l10000.tree.nexusfile")
tree <- drop.tip(tree,"P.nahrattemsis") # Remove species in the tree that are not in the data matrix
## read the phenotype data 
pheno <- read.csv("Dino_Chrys_Carotenod_pigments.csv",row.names = 1)

# PCA without accounting for phylogetic relationships
pca <- prcomp(pheno[,c(1:12)], scale. = T)
summary(pca)
pc_values <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2])
pc_values
# make the plot
colPalette <- c("black", "black", "red", "red", "black", "red", "red",  "black", "red", "black")
spp.names <- c("C.eryth", "C.fest", "C.haem", "C.l.ruf", "C.l.luc", "C.xanth", "C.stric", "D.b.jaff", "D.psar","D.b.punct")
quartz(height = 4.3, width = 4.3)
plot(pca$x,type="n", xlim=c(-3, 4.3), main="Standard PCA", 
     xlab="PC1 (59.7%)", ylab="PC2 (19.98%)")
text(pca$x,labels=spp.names,cex=0.8,col=colPalette)


# phylogenetic PCA
tree.pca <- phyl.pca(tree, pheno[,c(1:12)], method = "BM", mode = "corr")
summary(tree.pca)
tree_pc_values <- data.frame(PC1 = tree.pca$S[, 1], PC2 = tree.pca$S[, 2])
tree_pc_values
## Make the plot 
quartz(height = 4.3, width = 4.3)
plot(tree.pca$S,type="n", xlim=c(-0.6, 0.95),main="Phylogenetic PCA",
     xlab="PC1 (95.3%)", ylab="PC2 (3.8%)")
text(tree.pca$S,labels=spp.names,cex=0.8,col=colPalette)

################################## Phylogenetic PCA ##############################
###################################################################################


#################################################################################
##################### Ancestral states reconstruction ############################
#################################################################################

#################################################################################
######## Ancestral states reconstruction of Hue ################

# set the working direcotry
setwd("/Users/rashikaranasinghe/Library/CloudStorage/OneDrive-UBC/Dinopium_GBS_on_laptop/012NA_150_samples/Ancestor_reconstruction")

# load the packages
library(phytools)
library(ape)

## Read the molecular phylogenetic tree
tree <- read.nexus("Dino.Chrys.n60_40412SNPs_l10000.tree.nexusfile")
tree <- drop.tip(tree,"P.nahrattemsis") # Remove species in the tree that are not in the data matrix

## read the phenotype data 
pheno <- read.csv("Dino.Chrys.n60_40412SNPs_l10000.tree.pheno_Data.csv",row.names = 1)
# remove the pygmy woodpecker form teh data set because I don't have phenotype data for them.
pheno_new <- pheno[-c(3),]
pheno = pheno_new 

## extract character of interest
hue<-setNames(pheno$Hue_mantle,rownames(pheno))

## estimate ancestral state under BM model
fit.BM<-anc.ML(tree, hue)

## Breakdown continuous trait in categories
hue.cat <- cut(hue, breaks = 4, labels=FALSE) # for tip labels
hue.node.cat <- cut(fit.BM$ace, breaks = 4, labels=FALSE) # for nodes 

# Create a color palette from yellow to red
color_bins <- c("#FFFF00", "#FFDA00", "#FFB600", "#FF0000")

# Plot the data
quartz()
plot(tree, type="p", use.edge.length = T, label.offset = 0.8, cex = 0.9)
tiplabels(round(hue, 0.1), bg=color_bins[hue.cat], frame="circle", cex=0.4)
nodelabels(round(fit.BM$ace, 0.1), bg=color_bins[hue.node.cat], frame="circle", cex=0.4)

######## END Ancestral states reconstruction of Hue ################
#################################################################################


#########################################################################################
######## Tile plot showing percentages of each red pigment in for each individual ######
library(dplyr)
library(tidyr)

# read the data 
pheno <- read.csv("Dino.Chrys.n60_40412SNPs_l10000.tree.pheno_Data.csv", stringsAsFactors = T)

# Transform the data to a long format
data_long <- pheno %>%
  pivot_longer(cols = c(Astaxanthin_content_percentage:papilioerythrinone_percentage),
               names_to = "Pigment", 
               values_to = "Percentage") 
data_long$Tip_name <- factor(data_long$Tip_name, levels = c("C.erythrocephalus","C.festivus", "P.nahrattemsis", "C.haematribon", "C.l.rufopunctatus", "C.l.lucidus", "C.xanthocephalus","C.stricklandi", "D.b.jaffnense", "D.psarodes", "D.b.puncticolle"))

data_long$Pigment <- factor(data_long$Pigment, levels = c("Astaxanthin_content_percentage", "Adonirubin_percentage","Canthaxanthin_percentage","Alpa.doradexanthin_percentage","papilioerythrinone_percentage"))

# Define unique colors for each pigment
pigment_colors <- c("Canthaxanthin_percentage"="red3",
                    "Adonirubin_percentage"="red3",
                    "Astaxanthin_content_percentage"="red3",  
                    "Alpa.doradexanthin_percentage"="red3", 
                    "papilioerythrinone_percentage"="red3")

# Create the tile plot
quartz()
ggplot(data_long, aes(x = Pigment, y = Tip_name, fill = Pigment)) +
  geom_tile(color = "black", aes(alpha = Percentage)) +
  scale_fill_manual(values = pigment_colors) +
  scale_alpha(range = c(0, 3)) + # Vary transparency based on concentration
  geom_text(aes(label = round(Percentage, 2)), size = 3.5, color="grey2", face="bold") + 
  scale_x_discrete(labels = c("Canthaxanthin_percentage"="Canthaxanthin",
                              "Adonirubin_percentage"="Adonirubin",
                              "Astaxanthin_content_percentage"="Astaxanthin",  
                              "Alpa.doradexanthin_percentage"="Alpha-doradexanthin",
                              "papilioerythrinone_percentage"="Papilioerythrinone"))+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1,size=9),
        axis.text.y = element_blank(),
        axis.title = element_blank(), # to remove both x and y axis titles
        legend.position = "none",
        panel.grid.major = element_line(color = "white", linetype = "dashed")
  )

##### END Tile plot showing percentages of each red pigment in for each individual ######
#########################################################################################

#################################################################################
##################### END Ancestral states reconstruction ############################
#################################################################################

