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


###################################################################################################
####################### Phylogenetic independent contrasts (PIC) ################################
# Testing for statistical association between traitson the phylogenetic tree

## Read the molecular phylogenetic tree
tree <- read.nexus("Dino.Chrys.n60_40412SNPs_l10000.tree.nexusfile")
tree <- drop.tip(tree,"P.nahrattemsis") # Remove species in the tree that are not in the data matrix
tree$tip.label <- c("C.eryth", "C.fest", "C.haem", "C.l.ruf", "C.l.luc", "C.xanth", "C.stric", "D.b.jaff", "D.psar","D.b.punct")
## read the phenotype data
p <- read.csv("Dino.Chrys.n60_40412SNPs_l10000.tree.pheno_Data.csv",row.names = 1)
pheno <- p[-3,]
rownames(pheno) <- c("C.eryth", "C.fest", "C.haem", "C.l.ruf", "C.l.luc", "C.xanth", "C.stric", "D.b.jaff", "D.psar","D.b.punct")

# Order of the data in the trait file should be in the same order as the tip.label of the tree. 

#### Calculate contrasts
# Calculate the contrasts for each trait that have been scaled using the expected variances
C4_ketocarotenoids.contrast <- pic(pheno$Red_pig_conten_percentage, tree, scaled = T)
C4ketolation.contrast <- pic(pheno$avg_4_keto_groups, tree, scaled = T)
epsilon_rings.contrast <- pic(pheno$avg_epsilon_rings, tree, scaled = T)
C3oxygenation.contrast <- pic(pheno$avg_3_hydroxylation, tree, scaled = T)
hue.contrast <- pic(pheno$Hue_mantle, tree, scaled = T)

# Merge the contrasts into a dataframe
contrast_df <- data.frame(C4_ketocarotenoids.contrast, C4ketolation.contrast, epsilon_rings.contrast, C3oxygenation.contrast, hue.contrast)

################### Correlation between mantle hue and C4 ketocarotenoid percentage ####################
## Spearman's rank correlation without accounting for the phylogenetic assciations
cor.test(pheno$Hue_mantle, pheno$Red_pig_conten_percentage, method = "spearman")
#############
# data:  pheno$Hue_mantle and pheno$Red_pig_conten_percentage
# S = 28, p-value = 0.005557, rho =0.830303 

## Spearman's rank correlation accounting for the phylogenetic assciations
cor.test(contrast_df$hue.contrast, contrast_df$C4_ketocarotenoids.contrast, method = "spearman")
#############
# data:  contrast_df$hue.contrast and contrast_df$C4_ketocarotenoids.contrast
# S = 6, p-value = 0.0003527, rho = 0.95 

###### Visualize the data in a 2D plot 
## extract character of interest
C4_ketocarotenoids<-setNames(pheno$Red_pig_conten_percentage,rownames(pheno))
C4_ketocarotenoids.cat <- cut(C4_ketocarotenoids, breaks = 3, labels=FALSE) 

# Create a color palette from yellow to red
ColorPalette <- c("yellow","#FD8D3C","#E31A1C")

quartz(height = 4.5, width = 5)
phylomorphospace(tree, pheno[,c(11, 6)], xlim=c(520,611), ylim=c(-17,112),
                 xlab = "Mantle hue (nm)", ylab = "C4-keto-carotenoids (%)") 
points(pheno[,c(11, 6)], pch=21, bg=ColorPalette[C4_ketocarotenoids.cat],col="black",cex=1.2,adj=1)



################### Correlation between mantle hue and C4 ketoc groups  ####################
## Spearman's rank correlation without accounting for the phylogenetic associations
cor.test(pheno$Hue_mantle, pheno$avg_4_keto_groups, method = "spearman")
#############
# data:  pheno$Hue_mantle and pheno$avg_4_keto_groups
# S = 26, p-value = 0.004459, rho = 0.8424242

## Spearman's rank correlation accounting for the phylogenetic assciations
cor.test(contrast_df$hue.contrast, contrast_df$C4ketolation.contrast, method = "spearman")
#############
# data:  contrast_df$hue.contrast and contrast_df$C4ketolation.contrast
# S = 6, p-value = 0.0003527,  rho =0.95  

###### Visualize the data in a 2D plot 
## extract character of interest
C4_ketolation<-setNames(pheno$avg_4_keto_groups,rownames(pheno))
C4_ketolation.cat <- cut(C4_ketolation, breaks = 3, labels=FALSE) 

# Create a color palette from yellow to red
ColorPalette <- c("yellow","#FD8D3C","#E31A1C")

quartz(height = 4.5, width = 5)
phylomorphospace(tree, pheno[,c(11, 8)], xlim=c(520,611), ylim=c(-0.3,2.1), #  ylim=c(-0.3,2.1),
                 xlab = "Mantle hue (nm)", ylab = "No. of C4(4')-keto groups") 
points(pheno[,c(11, 8)], pch=21, bg=ColorPalette[C4_ketolation.cat],col="black",cex=1.2,adj=1)



################### Correlation between mantle hue and C3 oxygenated groups  ####################
## Spearman's rank correlation without accounting for the phylogenetic associations
cor.test(pheno$Hue_mantle, pheno$avg_3_hydroxylation, method = "spearman")
#############
# data:  pheno$Hue_mantle and pheno$avg_3_hydroxylation
# S = 266, p-value = 0.06647,  rho = -0.6121212

## Spearman's rank correlation accounting for the phylogenetic assciations
cor.test(contrast_df$hue.contrast, contrast_df$C3oxygenation.contrast, method = "spearman")
#############
# data:  contrast_df$hue.contrast and contrast_df$C3oxygenation.contrast
# S = 202, p-value = 0.05032,  rho = -0.6833333  

###### Visualize the data in a 2D plot 
# Create a color palette from yellow to red
ColorPalette <- c("yellow", "yellow","red", "red", "#FD8D3C","red", "red", "yellow","#E31A1C")
quartz(height = 4.5, width = 5)
phylomorphospace(tree, pheno[,c(11, 9)], xlim=c(520,613), ylim=c(0,2), #ylim=c(1.2,2.03),
                 xlab = "Mantle hue (nm)", ylab = "No. of C3(3')-oxygenated groups") 
points(pheno[,c(11, 9)], pch=21, bg=ColorPalette,col="black",cex=1.2,adj=1)


################### Correlation between mantle hue and epsilon rings  ####################
## Spearman's rank correlation without accounting for the phylogenetic associations
cor.test(pheno$Hue_mantle, pheno$avg_epsilon_rings, method = "spearman")
#############
# data:  pheno$Hue_mantle and pheno$avg_epsilon_rings
# S = 116, p-value = 0.407,  rho =0.2969697 

## Spearman's rank correlation accounting for the phylogenetic assciations
cor.test(contrast_df$hue.contrast, contrast_df$epsilon_rings.contrast, method = "spearman")
#############
# data:  contrast_df$hue.contrast and contrast_df$epsilon_rings.contrast
# S = 114, p-value = 0.9116, rho =0.05 

###### Visualize the data in a 2D plot 
# Create a color palette from yellow to red
ColorPalette <- c("yellow", "yellow","red", "red", "#FD8D3C","red", "red", "yellow","#E31A1C")
quartz(height = 4.5, width = 5)
phylomorphospace(tree, pheno[,c(11, 10)], xlim=c(520,613), ylim=c(0,2), #ylim=c(0.13,0.51),
                 xlab = "Mantle hue (nm)", ylab = "No. of epsilon rings") 
points(pheno[,c(11, 10)], pch=21, bg=ColorPalette,col="black",cex=1.2,adj=1)


####################### Phylogenetic independent contrasts (PIC) ################################
###################################################################################################



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
setwd("Dinopium_GBS_on_laptop/012NA_150_samples/Ancestor_reconstruction")

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

###########################################################################################
############## Ancestral state recostruction of other variables  #########################

# set the working direcotry
setwd("/Dinopium_GBS_on_laptop/012NA_150_samples/Ancestor reconstruction")

# load the packages
library(phytools)
library(ape)

## Read the molecular phylogenetic tree
tree.cladogram <- read.nexus("Dino.Chrys.n60_40412SNPs_l10000.tree_Edited_nexusfile.nex")
tree.cladogram <- drop.tip(tree,"P.nahrattemsis") # Remove species in the tree that are not in the data matrix
## read the phenotype data 
pheno <- read.csv("Dino.Chrys.n60_40412SNPs_l10000.tree.pheno_Data.csv",row.names = 1)
# remove the pygmy woodpecker form teh data set because I don't have phenotype data for them.
pheno_new <- pheno[-c(3),]
pheno = pheno_new 


############### C4-keto carotenoid content 
## extract character of interest
C4_ketocarotenoids <-setNames(pheno$Red_pig_conten_percentage,rownames(pheno))

## find the best fit model to estimate ancestral state for C4 ketolation
fit.BM<-anc.ML(tree.cladogram, C4_ketocarotenoids, model = "BM")
fit.EB<-anc.ML(tree.cladogram, C4_ketocarotenoids, model = "EB")

AIC(fit.BM, fit.EB)

## create "contMap" object
Dino.contMap<-contMap(tree.cladogram,C4_ketocarotenoids,plot=FALSE,res=200)
## change color scheme
Dino.contMap<-setMap(Dino.contMap,
                     c("white","grey98", "#FFCCCC", "red"))

quartz(width = 4.3, height = 4)
plot(Dino.contMap,fsize=c(0.7,0.8), offset=1.5,leg.txt="4-keto groups")
tiplabels(round(avg_4_keto_groups, 0.1), bg="white", frame="circle", cex=0.7, offset=1.5)
#nodelabels(round(fit.BM$ace, 0.1), bg="white", frame="circle")


############### C4 ketolation 
## extract character of interest
avg_4_keto_groups<-setNames(pheno$avg_4_keto_groups,rownames(pheno))

## find the best fit model to estimate ancestral state for C4 ketolation
fit.BM<-anc.ML(tree.cladogram, avg_4_keto_groups, model = "BM")
fit.EB<-anc.ML(tree.cladogram, avg_4_keto_groups, model = "EB")

AIC(fit.BM, fit.EB)

## create "contMap" object
Dino.contMap<-contMap(tree.cladogram,avg_4_keto_groups,plot=FALSE,res=200)
## change color scheme
Dino.contMap<-setMap(Dino.contMap,
                     c("white","grey95", "grey90", "brown", "brown4"))

quartz(width = 4.3, height = 4)
plot(Dino.contMap,fsize=c(0.7,0.8), offset=1.5,leg.txt="4-keto groups")
tiplabels(round(avg_4_keto_groups, 0.1), bg="white", frame="circle", offset=1.4)
#nodelabels(round(fit.BM$ace, 0.1), bg="white", frame="circle")

############### Number of hydroxy group composition  along the tree ##################
## extract character of interest
avg_3_hydroxylation<-setNames(pheno$avg_3_hydroxylation,rownames(pheno))

## find the best fit model to estimate ancestral state for C4 ketolation
fit.BM<-anc.ML(tree.cladogram, avg_3_hydroxylation, model = "BM")
fit.EB<-anc.ML(tree.cladogram, avg_3_hydroxylation, model = "EB")

AIC(fit.BM, fit.EB)

## create "contMap" object
Dino.contMap<-contMap(tree.cladogram,avg_3_hydroxylation,plot=FALSE,res=200)
## change color scheme
Dino.contMap<-setMap(Dino.contMap,
                     c("white","grey90", "grey", "grey30", "grey20"))

quartz(width = 4.3, height = 4)
plot(Dino.contMap,fsize=c(0.7,0.8), offset=1.4,leg.txt="3-hydoxyl groups")
tiplabels(round(avg_3_hydroxylation, 1), bg="white", frame="circle", cex=0.6,offset=1.5)
#nodelabels(round(fit.BM$ace, 1), bg="white", frame="circle", cex=0.6)
#title("avg_3_hydroxylation")


############### Number of epsilon rings composition  along the tree ##################
## extract character of interest
avg_epsilon_rings<-setNames(pheno$avg_epsilon_rings,rownames(pheno))

## find the best fit model to estimate ancestral state for C4 ketolation
fit.BM<-anc.ML(tree.cladogram, avg_3_hydroxylation, model = "BM")
fit.EB<-anc.ML(tree.cladogram, avg_3_hydroxylation, model = "EB")

AIC(fit.BM, fit.EB)

## create "contMap" object
Dino.contMap<-contMap(tree.cladogram,avg_epsilon_rings,plot=FALSE,res=200)
## change color scheme
Dino.contMap<-setMap(Dino.contMap,
                     c("white","grey95", "grey90", "yellow"))

quartz(width = 4.3, height = 4)
plot(Dino.contMap,fsize=c(0.7,0.8), offset=1.35,
     leg.txt="epsilon rings")
tiplabels(round(avg_epsilon_rings, 1), bg="white", frame="circle", cex=0.7, offset=1.2)
#nodelabels(round(fit.BM$ace, 0.1), bg="white", frame="circle")
############################ END Plot functional groups  #################################
###########################################################################################

#################################################################################
##################### END Ancestral states reconstruction ############################
#################################################################################

