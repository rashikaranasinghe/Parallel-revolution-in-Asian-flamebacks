##################################################################################
########## Compare red pigments among red flamebacks for corwn and  mantle feathers #########
##################################################################################

############################# Make the pigmnet composition plot for mantle feathers #########################
## set the working directory
setwd("My Research/carotenoid analysis_I/HPLC/carotenoid analysis")

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)


###### Makd plots for mantle data
# Read teh data file
parallel.red <- read.csv("HPLC_Dino_Chrys_Red_birds_mantle_Sep2024.csv")
levels(parallel.red$Scientific.name)

# wide to long format
Pigment_content <- parallel.red %>%
  # filter(Scientific.name=="C. haematribon" | Scientific.name=="C. l. rufopunctatus" |Scientific.name=="C. stricklandi " |Scientific.name=="C. xanthocephalus" |Scientific.name=="D. psarodes  " & feathers == "mantle") %>%
  pivot_longer(cols = c(beta.carotene_concentration:lutein.ident__concentration,canthaxanthin_concentration:papilioerythrinone_concentration),
               names_to = "Variable",
               values_to = "Values")
levels(as.factor(Pigment_content$Variable))
Pigment_content$Variable <- factor(Pigment_content$Variable, levels = c("beta.carotene_concentration" , "beta.cryptoxanthin_concentration", "zeaxanthin._concentration", "lutein.ident__concentration","canthaxanthin_concentration", "adonirubin_concentration", "astaxanthin_concentration", "alpha.doradexanthin_concentration" , "papilioerythrinone_concentration"))
levels(as.factor(Pigment_content$Scientific.name))
Pigment_content$Scientific.name <- factor(Pigment_content$Scientific.name, levels = c("C. haematribon", "C. stricklandi ", "C. l. rufopunctatus",  "C. xanthocephalus", "D. psarodes  "))

## Make the plot
quartz(height=10, width = 5)
ggplot(Pigment_content, aes(x = Variable, y = Values, fill = Variable)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.6) +
  scale_fill_manual(values = c("beta.carotene_concentration"="yellow",
                               "beta.cryptoxanthin_concentration" = "yellow3",
                               "zeaxanthin._concentration"="gold",
                               "lutein.ident__concentration" = "gold3", 
                               "canthaxanthin_concentration" = "red3", 
                               "adonirubin_concentration" = "red4", 
                               "astaxanthin_concentration" = "red1", 
                               "alpha.doradexanthin_concentration" = "indianred3", 
                               "papilioerythrinone_concentration" = "firebrick"),
                    labels = c("beta.carotene_concentration"= "beta-carotene",
                               "beta.cryptoxanthin_concentration" = "Beta-cryptoxanthin", 
                               "zeaxanthin._concentration"="Zeaxanthin",
                               "lutein.ident__concentration" = "Lutein", 
                               "canthaxanthin_concentration" = "Canthaxanthin", 
                               "adonirubin_concentration" = "Adonirubin", 
                               "astaxanthin_concentration" = "Astaxanthin", 
                               "alpha.doradexanthin_concentration" = "Alpha-doradexanthin", 
                               "papilioerythrinone_concentration" = "Papilioerythrinone"))+
  facet_wrap(~Scientific.name, ncol = 1
             ,scales = "fixed")+
  ylab("Concentration (mg/g)") +
  scale_x_discrete(labels = c("beta-carotene", "beta-cryptoxanthin", "Zeaxanthin", "Lutein", "Canthaxanthin", "Adonirubin", "Astaxanthin", "Alpha-doradexanthin", "Papilioerythrinone")) +
 ylim(c(0,60))+
  theme_linedraw() +
  scale_x_discrete(NULL)+  
  theme(
    axis.text.x = element_blank(),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white"),
    strip.background = element_rect(color = "white", fill = "white", size = 0.1),
    strip.text = element_text(face = "bold.italic", color = "black"),
    legend.position = "right",
    legend.title = element_blank(),
    plot.title = element_text(size = 7, face = "bold")
  ) +
  guides(fill = guide_legend(nrow = 9)) + 
  ggtitle("Red pigment concentrations in mantle among red flamebacks")


###### Makd plots for crown data
# read the file 
parallel.red <- read.csv("HPLC_Dino_Chrys_Red_birds_crown_Sep2024.csv")
levels(parallel.red$Scientific.name)

# wide to land format 
Pigment_content <- parallel.red %>%
  # filter(Scientific.name=="C. haematribon" | Scientific.name=="C. l. rufopunctatus" |Scientific.name=="C. stricklandi " |Scientific.name=="C. xanthocephalus" |Scientific.name=="D. psarodes  " & feathers == "crown" & color == "red") %>%
  pivot_longer(cols = c(beta.carotene_concentration:lutein.ident__concentration,canthaxanthin_concentration:papilioerythrinone_concentration),
               names_to = "Variable",
               values_to = "Values")
levels(as.factor(Pigment_content$Variable))
Pigment_content$Variable <- factor(Pigment_content$Variable, levels = c("beta.carotene_concentration" , "beta.cryptoxanthin_concentration", "zeaxanthin._concentration", "lutein.ident__concentration","canthaxanthin_concentration", "adonirubin_concentration", "astaxanthin_concentration", "alpha.doradexanthin_concentration" , "papilioerythrinone_concentration"))
levels(as.factor(Pigment_content$Scientific.name))
Pigment_content$Scientific.name <- factor(Pigment_content$Scientific.name, levels = c("C. haematribon", "C. stricklandi ", "C. l. rufopunctatus",  "C. xanthocephalus", "D. psarodes  "))
levels(as.factor(Pigment_content$color))

# make the plot
quartz(height=10, width = 5)
ggplot(Pigment_content, aes(x = Variable, y = Values, fill = Variable)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.6) +
  scale_fill_manual(values = c("beta.carotene_concentration"="yellow",
                               "beta.cryptoxanthin_concentration" = "yellow3",
                               "zeaxanthin._concentration"="gold",
                               "lutein.ident__concentration" = "gold3", 
                               "canthaxanthin_concentration" = "red3", 
                               "adonirubin_concentration" = "red4", 
                               "astaxanthin_concentration" = "red1", 
                               "alpha.doradexanthin_concentration" = "indianred3", 
                               "papilioerythrinone_concentration" = "firebrick"),
                    labels = c("beta.carotene_concentration"= "beta-carotene",
                               "beta.cryptoxanthin_concentration" = "Beta-cryptoxanthin", 
                               "zeaxanthin._concentration"="Zeaxanthin",
                               "lutein.ident__concentration" = "Lutein", 
                               "canthaxanthin_concentration" = "Canthaxanthin", 
                               "adonirubin_concentration" = "Adonirubin", 
                               "astaxanthin_concentration" = "Astaxanthin", 
                               "alpha.doradexanthin_concentration" = "Alpha-doradexanthin", 
                               "papilioerythrinone_concentration" = "Papilioerythrinone"))+
  facet_wrap(~Scientific.name, ncol = 1
             , scales = "fixed")+
  ylab("Concentration (mg/g)") +
  #ylim(c(0,60))+
  scale_x_discrete(labels = c("beta-carotene", "beta-cryptoxanthin", "Zeaxanthin", "Lutein", "Canthaxanthin", "Adonirubin", "Astaxanthin", "Alpha-doradexanthin", "Papilioerythrinone")) +
  theme_linedraw() +
  scale_x_discrete(NULL)+  
  theme(
    axis.text.x = element_blank(),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white"),
    strip.background = element_rect(color = "white", fill = "white", size = 0.1),
    strip.text = element_text(face = "bold.italic", color = "black"),
    legend.position = "right",
    legend.title = element_blank(),
    plot.title = element_text(size = 7, face = "bold")
  ) +
  guides(fill = guide_legend(nrow = 9)) + 
  ggtitle("Red pigment concentrations in crown among red flamebacks")


################################################################################################
########## Compare functional group compositions among red flamebacks ##########################

## I have to make the plots for crown and mantle seperately because i have to exclude all the female birds in C. xanthocephalus when making the plot for crown feathers. 

#################### Make the plot for mantle feathers ################
## set the working directory
setwd("/My Research/carotenoid analysis_I/HPLC/carotenoid analysis")

library(dplyr)
library(tidyr)
library(ggplot2)


# Read teh data file
#parallel.red <- read.csv("HPLC_Dino_Chrys_mantle_crown_Sep2024.csv", stringsAsFactors = T)
parallel.red <- read.csv("HPLC_Dino_Chrys_Red_birds_mantle_Sep2024.csv", stringsAsFactors = T)
levels(parallel.red$Scientific.name)
names(parallel.red)


Pigment_content <- parallel.red %>%
  # filter(Scientific.name=="C. haematribon" | Scientific.name=="C. l. rufopunctatus" |Scientific.name=="C. stricklandi " |Scientific.name=="C. xanthocephalus" |Scientific.name=="D. psarodes  " & feathers == "mantle") %>%
  pivot_longer(cols = c(avg.4.ketonization,avg.hydroxylation_JH,epsilon_rings),
               names_to = "Variable",
               values_to = "Values")
levels(as.factor(Pigment_content$Variable))
Pigment_content$Variable <- factor(Pigment_content$Variable, levels = c("epsilon_rings","avg.hydroxylation_JH" ,"avg.4.ketonization"))
levels(as.factor(Pigment_content$Scientific.name))
Pigment_content$Scientific.name <- factor(Pigment_content$Scientific.name, levels = c("C. haematribon", "C. stricklandi ", "C. l. rufopunctatus",  "C. xanthocephalus", "D. psarodes  "))

quartz(height=3.7, width = 8)
ggplot(Pigment_content, aes(x = Scientific.name, y = Values, fill = Variable)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.6, alpha = 0.75, size = 0.2) +
  scale_fill_manual(values = c("epsilon_rings"="yellow",
                               "avg.hydroxylation" ="grey",
                               "avg.4.ketonization"= "brown"),
                    labels = c("epsilon_rings"="Epsilon rings",
                               "avg.hydroxylation"="C3(3')-hydroxyl groups" ,
                               "avg.4.ketonization"="C4(4')-keto groups"))+
  ylab("Number of functional groups") +
  theme_linedraw() +
  scale_x_discrete(NULL)+  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face="italic"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey80"),
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(size = 7, face = "bold")
  ) +
  guides(fill = guide_legend(nrow = 1)) + 
  ggtitle("Functional group composition for mantle feathers")
#################### END Make the plot for mantle feathers ################

#################### Make the plot for crown feathers ################
# read the file 
parallel.red <- read.csv("HPLC_Dino_Chrys_Red_birds_crown_Sep2024.csv", stringsAsFactors = T)

Pigment_content <- parallel.red %>%
  # filter(Scientific.name=="C. haematribon" | Scientific.name=="C. l. rufopunctatus" |Scientific.name=="C. stricklandi " |Scientific.name=="C. xanthocephalus" |Scientific.name=="D. psarodes  " & feathers == "crown" & color == "red")%>%
  pivot_longer(cols = c(avg.4.ketonization,avg.hydroxylation,epsilon_rings),
               names_to = "Variable",
               values_to = "Values")
levels(as.factor(Pigment_content$Variable))
Pigment_content$Variable <- factor(Pigment_content$Variable, levels = c("epsilon_rings","avg.hydroxylation" ,"avg.4.ketonization"))
levels(as.factor(Pigment_content$Scientific.name))
Pigment_content$Scientific.name <- factor(Pigment_content$Scientific.name, levels = c("C. haematribon", "C. stricklandi ", "C. l. rufopunctatus",  "C. xanthocephalus", "D. psarodes  "))

quartz(height=3.7, width = 8)
ggplot(Pigment_content, aes(x = Scientific.name, y = Values, fill = Variable)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.6, alpha = 0.75, size=0.2) +
  scale_fill_manual(values = c("epsilon_rings"="yellow",
                               "avg.hydroxylation" ="grey",
                               "avg.4.ketonization"= "brown"),
                    labels = c("epsilon_rings"="Epsilon rings",
                               "avg.hydroxylation"="C3(3')-hydroxyl groups" ,
                               "avg.4.ketonization"="C4(4')-keto groups"))+
  ylab("Number of functional groups") +
  theme_linedraw() +
  scale_x_discrete(NULL)+  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face="italic"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey80"),
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(size = 7, face = "bold")
  ) +
  guides(fill = guide_legend(nrow = 1)) + 
  ggtitle("Functional group composition for crown feathers")
#################### END Make the plot for crown feathers ################


################# Sumarry statistics of functional groups #####################
## set the working directory
setwd("My Research/carotenoid analysis_I/HPLC/carotenoid analysis")

library(dplyr)
library(tidyr)
library(ggplot2)

# Read teh data file
parallel.red <- read.csv("HPLC_Dino_Chrys_mantle_crown_Sep2024.csv", stringsAsFactors = T)

######## Mantle #############
# Calculate mean, sd, LL of 95%CI and UL of 95%CI for all the variables for crown and mantle 
summary_stats <- parallel.red %>%
  filter(Scientific.name=="C. haematribon" | Scientific.name=="C. l. rufopunctatus" |Scientific.name=="C. stricklandi " |Scientific.name=="C. xanthocephalus" |Scientific.name=="D. psarodes  " & feathers == "mantle")%>%
  group_by(Scientific.name) %>%
  summarise(across(avg.4.ketonization:Total_red_pigment_concentration, 
                   list(
                     mean = ~mean(.x, na.rm = TRUE), 
                     sd = ~sd(.x, na.rm = TRUE),
                     CI_lower = ~mean(.x, na.rm = TRUE) - (1.96 * sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))),
                     CI_upper = ~mean(.x, na.rm = TRUE) + (1.96 * sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))))
                   ),
                   .names = "{col}_{fn}"))
# save the data
#write.csv(summary_stats, file = "../../../../Dinopium_GBS_on_laptop/012NA_150_samples/parallel_evolution/summary_stats_mantle.csv")


######## END Mantle #############

######## Crown #############
# Calculate mean, sd, LL of 95%CI and UL of 95%CI for all the variables for crown and mantle 
summary_stats <- parallel.red %>%
  filter(Scientific.name=="C. haematribon" | Scientific.name=="C. l. rufopunctatus" |Scientific.name=="C. stricklandi " |Scientific.name=="C. xanthocephalus" |Scientific.name=="D. psarodes  " & feathers == "crown" & color == "red")%>%
  group_by(Scientific.name) %>%
  summarise(across(avg.4.ketonization:Total_red_pigment_concentration, 
                   list(
                     mean = ~mean(.x, na.rm = TRUE), 
                     sd = ~sd(.x, na.rm = TRUE),
                     CI_lower = ~mean(.x, na.rm = TRUE) - (1.96 * sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))),
                     CI_upper = ~mean(.x, na.rm = TRUE) + (1.96 * sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))))
                   ),
                   .names = "{col}_{fn}"))
# save the data
#write.csv(summary_stats, file = "../../../../Dinopium_GBS_on_laptop/012NA_150_samples/parallel_evolution/summary_stats_crown.csv")

######## END Crown #############

################# END Sumarry statistics of functional groups #####################


######################################################################################
####################### All functional group vairables on the same plot ###############
# Load the ggplot2 library
library(ggplot2)

# Read teh data file
parallel.red <- read.csv("HPLC_Dino_Chrys_mantle_crown_Sep2024.csv", stringsAsFactors = T)

# Create the scatter plot
quartz(height = 3.8, width = 5.8)
ggplot(parallel.red, aes(x = avg.4.ketonization, y = epsilon_rings, fill = avg.hydroxylation_JH, shape = Scientific.name)) +
  geom_point(alpha = 0.7,size= 2.5) +
  scale_fill_gradient2(low = "white", mid = "#ADD8E6", high = "#00008B", midpoint = mean(parallel.red$avg.hydroxylation_JH, na.rm = TRUE)) +
  scale_shape_manual(values = c(25, 22, 23, 24, 21), labels = c(expression(italic("C. haematribon")), 
                                                              expression(italic("C. l. rufopunctatus")), 
                                                              expression(italic("C. xanthocephalus")), 
                                                              expression(italic("C. stricklandi")), 
                                                              expression(italic("D. psarodes")))) +
  labs(title = "Scatter Plot of functional groups by species",
       x = "C4(4')-keto groups",
       y = "Epsilon rings",
       fill = "C3(3')-oxygenated groups",
       size = "Species",
       shape = "Species") +
  theme_test()+
  theme(panel.background = element_blank(),  # Remove panel background
        plot.background = element_blank(),   # Remove plot background
        legend.background = element_blank(), # Remove legend background
        panel.grid.major = element_blank(),   # Optional: remove major grid lines
        panel.grid.minor = element_blank()) 

####################### All functional group vairables on the same plot ##############
######################################################################################

