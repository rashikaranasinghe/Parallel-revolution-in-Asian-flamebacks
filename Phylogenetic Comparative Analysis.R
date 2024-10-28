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
