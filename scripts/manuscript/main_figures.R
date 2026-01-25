setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)




################ ################ ################ ################ ################ 
################ FIGURE 1 ################ ################ ################ ################ 

## Combined data frames (feature importance metrics, posterior means and gene binnings) 
merged_df_Mouse <- read_delim("../../intermediate_data/Merged_df_feature_posterior_Mouse", 
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE)

merged_df_Dmel <- read_delim("../../intermediate_data/Merged_df_feature_posterior_Mouse", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)

merged_df_Yeast <- read_delim("../../intermediate_data/Merged_df_feature_posterior_Mouse", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)

merged_df_Mouse$shet_bin <- factor(merged_df_Mouse$shet_bin,
                             levels= c("High constraint", "Moderate-high constraint", "Moderate-low constraint", "Low constraint"))
merged_df_Dmel$shet_bin <- factor(merged_df_Dmel$shet_bin,
                             levels= c("High constraint", "Moderate-high constraint", "Moderate-low constraint", "Low constraint"))
merged_df_Yeast$shet_bin <- factor(merged_df_Yeast$shet_bin,
                             levels= c("High constraint", "Moderate-high constraint", "Moderate-low constraint", "Low constraint"))

### Ex. to create boxplots

ggplot(merged_df_Mouse, aes(x = shet_bin, y = mean_phastCons)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6, notch = TRUE) +
  theme_minimal() +
  labs(x = "Shet bin (post_mean)", y = "Mean phastCons") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+guides(colour="none")+xlab("")


################ ################ ################ ################ ################ 
################ FIGURE 2 ################ ################ ################ ################ 

Bestfits_all_paramsMLE <- read_delim("../../intermediate_data/Bestfits_DFE_all_paramsMLE", 
                                     delim = "\t", escape_double = FALSE, 
                                     trim_ws = TRUE)

Bestfits_all_paramsMLE$Species<- factor(Bestfits_all_paramsMLE$Species, levels=c("Mouse","Dmel","Yeast"),labels = c("Mouse","Fruit fly","Yeast"))


Bestfits_all_paramsMLE_ss<-filter(Bestfits_all_paramsMLE, ! Parameter %in% c("theta","eps"))

Bestfits_all_paramsMLE_ss$inference <- factor(Bestfits_all_paramsMLE_ss$inference,
                                              levels= c("full","High constraint", "Moderate-high constraint", "Moderate-low constraint", "Low constraint"))

### 2A
pdf(file = "Fig2A.pdf",   width = 8, height = 4) 

ggplot(Bestfits_all_paramsMLE_ss[Bestfits_all_paramsMLE_ss$Parameter=="S_d",],aes(x=inference,y=log10(abs(MLE)),fill=inference))+geom_bar(position = "dodge", stat = "identity")+
  theme(panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"), legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size=7))+ylab("MLE")+xlab("")+
  scale_fill_viridis(discrete = TRUE)+  facet_wrap(Species~ ., nrow = 1) +
  geom_errorbar(aes(ymax=log10(abs(Upper_Bound)), ymin=log10(abs(Lower_Bound))))+ ylab(expression(log10(S[d])))

dev.off()

### 2B
pdf(file = "Fig2B.pdf",   width = 8, height = 4) 

ggplot(Bestfits_all_paramsMLE_ss[Bestfits_all_paramsMLE_ss$Parameter=="b",],aes(x=inference,y=abs(MLE),fill=inference))+geom_bar(position = "dodge", stat = "identity")+
  theme(panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"), legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size=7))+ylab("MLE")+xlab("")+
  scale_fill_viridis(discrete = TRUE)+  facet_wrap(~ Species,  nrow = 1) +
  geom_errorbar(aes(ymax=abs(Upper_Bound), ymin=abs(Lower_Bound)))+ ylab(expression(b))

dev.off()


### 2C

pdf(file = "Fig2C.pdf",   width = 8, height = 4) 

ggplot(Bestfits_all_paramsMLE_ss[Bestfits_all_paramsMLE_ss$Parameter=="S_b",],aes(x=inference,y=log10(abs(MLE)),fill=inference))+geom_bar(position = "dodge", stat = "identity")+
  theme(panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"), legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size=7))+ylab("MLE")+xlab("")+
  scale_fill_viridis(discrete = TRUE)+  facet_wrap(~ Species, nrow = 1)+ 
  geom_errorbar(aes(ymax=log10(abs(Upper_Bound)), ymin=log10(abs(Lower_Bound))))+ylab(expression(log10(S[b])))

dev.off()

### 2D

pdf(file = "Fig2D.pdf",   width = 8, height = 4) 

ggplot(Bestfits_all_paramsMLE_ss[Bestfits_all_paramsMLE_ss$Parameter=="p_b",],aes(x=inference,y=abs(MLE),fill=inference))+geom_bar(position = "dodge", stat = "identity")+
  theme(panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"), legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size=7))+ylab("MLE")+xlab("")+
  scale_fill_viridis(discrete = TRUE)+  facet_wrap(~ Species, nrow = 1) +
  geom_errorbar(aes(ymax=abs(Upper_Bound), ymin=abs(Lower_Bound)))+ ylab("proportion of beneficial mutations")

dev.off()

################ ################ ################ ################ ################ 
################ FIGURE 3 ################ ################ ################ ################ 

Bestfits_all_paramsMLE_alphathr <- read_delim("../../intermediate_data/Bestfits_DFE_alphathr_all_paramsMLE", 
                                     delim = "\t", escape_double = FALSE, 
                                     trim_ws = TRUE)

Bestfits_all_paramsMLE_alphathr$Species<- factor(Bestfits_all_paramsMLE_alphathr$Species, levels=c("Mouse","Dmel","Yeast"),labels = c("Mouse","Fruit fly","Yeast"))


Bestfits_all_paramsMLE_alphathr_ss<-filter(Bestfits_all_paramsMLE_alphathr, ! Parameter %in% c("theta","eps"))

Bestfits_all_paramsMLE_alphathr_ss$inference <- factor(Bestfits_all_paramsMLE_alphathr_ss$inference,
                                              levels= c("full","High constraint", "Moderate-high constraint", "Moderate-low constraint", "Low constraint"))

### Upper panel


pdf(file = "Fig3_upper.pdf",   width = 10, height = 4) 

ggplot(Bestfits_all_paramsMLE_ss[Bestfits_all_paramsMLE_ss$Parameter=="alpha",],aes(x=inference,y=abs(MLE),fill=inference))+geom_bar(position = "dodge", stat = "identity")+
  theme(panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"), legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size=7))+ylab("MLE")+xlab("")+
  scale_fill_viridis(discrete = TRUE)+  facet_wrap(~ Species, nrow = 1) +
  geom_errorbar(aes(ymax=abs(Upper_Bound), ymin=abs(Lower_Bound)))+ ylab(expression(alpha))

dev.off()

### Lower panel

pdf(file = "Fig3_lower.pdf",   width = 10, height = 4) 

ggplot(Bestfits_all_paramsMLE_alphathr_ss[Bestfits_all_paramsMLE_alphathr_ss$Parameter=="alpha",],aes(x=inference,y=abs(MLE),fill=inference))+geom_bar(position = "dodge", stat = "identity")+
  theme(panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"), legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size=7))+ylab("MLE")+xlab("")+
  scale_fill_viridis(discrete = TRUE)+  facet_wrap(~ Species, nrow = 1) +
  geom_errorbar(aes(ymax=abs(Upper_Bound), ymin=abs(Lower_Bound)))+ ylab(expression(alpha))

dev.off()

################ ################ ################ ################ ################ 
################ FIGURE 4 ################ ################ ################ ################ 

all_bingenes <- read_delim("../../intermediate_data/Propgenes_binned_all", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)

all_bingenes$Species<- factor(all_bingenes$Species, levels=c("Mouse","Dmel","Yeast"),labels = c("Mouse","Fruit fly","Yeast"))
all_bingenes$Bin <- factor(all_bingenes$Bin,
                           levels = c("High constraint", "Moderate-high constraint", "Moderate-low constraint", "Low constraint"))


### 4A 


pdf(file = "Fig4A.pdf",   width = 8, height = 4) 

ggplot(all_bingenes, aes(x = Bin, y = Percentage, fill=Species, colour = Species)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+xlab("")

dev.off()

## 4B

all_bingenes$Bin2<-rep(c(rep("Highly deleterious",2), rep("Weakly deleterious",2)), 3)


pdf(file = "Fig4B.pdf",   width = 8, height = 4) 

ggplot(all_bingenes, aes(x = Bin2, y = Percentage, fill=Species, colour = Species)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+xlab("")

dev.off()