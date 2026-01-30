setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(heatmap3)
library(ggnewscale)

################ ################ ################ ################ ################ 
################ FIGURE S1 ################ ################ ################ ################ 


full_all_feature_Mouse <- read_delim("../../GeneBayes_out/Mouse_fullfeature_best_model.feature_importance.tsv", 
                               delim = "\t", escape_double = FALSE, 
                               trim_ws = TRUE)



full_all_feature_Mouse$Category6<- factor(full_all_feature_Mouse$Category6, levels=c("gene_structure","aa_composition",
                                                           "GO_terms",
                                                           "Gene_expression",
                                                           "ppi",
                                                           "Conservation"))

df_long <- full_all_feature_Mouse %>%
  dplyr::select(feature, Category6,
                param0_importance, param1_importance) %>%
  tidyr::pivot_longer(
    cols = c(param0_importance, param1_importance),
    names_to = "parameter",
    values_to = "importance"
  ) %>%
  dplyr::mutate(
    parameter = recode(parameter,
                       param0_importance = "Prior mean",
                       param1_importance = "Prior variance")
  )


gene_expression <- df_long %>%
  dplyr::filter(Category6 == "Gene_expression")

gene_expression <- gene_expression %>%
  dplyr::mutate(importance = ifelse(importance == 0, NA, importance))


ggplot(gene_expression,
       aes(x = 1, y = feature, fill = importance)) +
  geom_tile() +
  scale_fill_viridis_c(
    trans = "sqrt",
    na.value = "grey90",
    name = "Relative contribution"
  ) +
  facet_wrap(~ parameter, nrow = 1) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

GO_terms <- df_long %>%
  dplyr::filter(Category6 == "GO_terms")

GO_terms <- GO_terms %>%
  dplyr::mutate(importance = ifelse(importance == 0, NA, importance))


ggplot(GO_terms,
       aes(x = 1, y = feature, fill = importance)) +
  geom_tile() +
  scale_fill_viridis_c(
    trans = "sqrt",
    na.value = "grey90",
    name = "Relative contribution"
  ) +
  facet_wrap(~ parameter, nrow = 1) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

# Rest of the categories
sel_cats <- c("gene_structure", "aa_composition","Conservation", "ppi")

# Define colors for categories
cat_cols <- c(
  Conservation   = "#1b9e77",
  gene_structure = "#d95f02",
  aa_composition = "#7570b3",
  ppi            = "#e7298a"
)

# Prepare data
df_plot <- df_long %>%
  filter(Category6 %in% sel_cats) %>%
  mutate(
    importance = ifelse(importance == 0, NA, importance),
    x_main = 1,    # main heatmap column
    x_cat  = 0     # category annotation column
  ) %>%
  arrange(Category6, feature) %>%
  mutate(feature = factor(feature, levels = unique(feature)))

# Compute category boundaries for labels / separators
cat_bounds <- df_plot %>%
  distinct(feature, Category6) %>%
  mutate(y = as.numeric(feature)) %>%
  group_by(Category6) %>%
  summarise(
    y_min = min(y) - 0.5,
    y_max = max(y) + 0.5,
    y_mid = mean(c(min(y), max(y)))
  ) %>%
  mutate(x_cat = 0)   # add x_cat so labels can reference it



# Plot
ggplot() +
  
  # Category side stripe
  geom_tile(
    data = df_plot,
    aes(x = x_cat, y = feature, fill = Category6),
    width = 0.08
  ) +
  scale_fill_manual(values = cat_cols, guide = "none") +
  
  ggnewscale::new_scale_fill() +
  
  # Main heatmap
  geom_tile(
    data = df_plot,
    aes(x = x_main, y = feature, fill = importance),
    width = 0.9
  ) +
  scale_fill_viridis_c(
    trans = "sqrt",
    na.value = "grey90",
    name = "Relative contribution"
  ) +
  
  # Horizontal lines separating categories
  geom_hline(
    data = cat_bounds,
    aes(yintercept = y_min),
    color = "grey40",
    linewidth = 0.4
  ) +
  
  # Category labels at left of stripe
  geom_text(
    data = cat_bounds,
    aes(x = x_cat - 0.1, y = y_mid, label = Category6),
    hjust = 1,
    size = 3
  ) +
  
  # Facet by prior mean / variance
  facet_wrap(~ parameter, nrow = 1) +
  
  scale_x_continuous(expand = c(0, 0)) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank()
  )


################ ################ ################ ################ ################ 
################ FIGURE S2 ################ ################ ################ ################ 


full_all_feature_Dmel <- read_delim("../../GeneBayes_out/Dmel_fullfeature_best_model.feature_importance.tsv", 
                                     delim = "\t", escape_double = FALSE, 
                                     trim_ws = TRUE)


full_all_feature_Dmel$Category6<- factor(full_all_feature_Dmel$Category6, levels=c("gene_structure","aa_composition",
                                                           "pathway",
                                                           "Gene_expression",
                                                           "interactions",
                                                           "Conservation"))


### Figure S1 with full_all_feature_Dmel

################ ################ ################ ################ ################ 
################ FIGURE S3 ################ ################ ################ ################ 

full_all_feature_Yeast <- read_delim("../../GeneBayes_out/Yeast_fullfeature_best_model.feature_importance.tsv", 
                                    delim = "\t", escape_double = FALSE, 
                                    trim_ws = TRUE)


full_all_feature_Yeast$Category6<- factor(full_all_feature_Yeast$Category6, levels=c("gene_structure","aa_composition",
                                                           "GO_terms",
                                                           "Gene_expression",
                                                           "ppi",
                                                           "Conservation"))

### Figure S1 with full_all_feature_Yeast

################ ################ ################ ################ ################ 
################ FIGURE S4 ################ ################ ################ ################ 

plot_spearman_heatmap <- function(df_subset, rho_limit) {
  
  df_plot <- df_subset %>%
    arrange(category, feature) %>%
    mutate(
      feature = factor(feature, levels = unique(feature)),
      x_main = 1,
      x_cat  = 0,
      # Mask non-significant correlations
      spearman_rho_plot = ifelse(FDR <= 0.05, spearman_rho, NA)
    )
  
  cat_bounds <- df_plot %>%
    distinct(feature, category) %>%
    mutate(y = as.numeric(feature)) %>%
    group_by(category) %>%
    summarise(
      y_min = min(y) - 0.5,
      y_max = max(y) + 0.5,
      y_mid = mean(c(min(y), max(y)))
    ) %>%
    mutate(x_cat = 0)
  
  ggplot() +
    # Category stripe
    geom_tile(
      data = df_plot,
      aes(x = x_cat, y = feature, fill = category),
      width = 0.08
    ) +
    scale_fill_manual(values = cat_cols, guide = "none") +
    
    ggnewscale::new_scale_fill() +
    
    # Main heatmap
    geom_tile(
      data = df_plot,
      aes(x = x_main, y = feature, fill = spearman_rho_plot),
      width = 0.9
    ) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0,
      limits = c(-rho_limit, rho_limit),
      na.value = "grey80",   # non-significant correlations appear gray
      name = "Spearman rho"
    ) +
    
    # Horizontal separators
    geom_hline(
      data = cat_bounds,
      aes(yintercept = y_min),
      color = "grey40",
      linewidth = 0.4
    ) +
    
    # Category labels
    geom_text(
      data = cat_bounds,
      aes(x = x_cat - 0.1, y = y_mid, label = category),
      hjust = 1,
      size = 3
    ) +
    
    scale_x_continuous(expand = c(0, 0)) +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )
}



postmean_feature_cors_Mouse <- read_delim("../../intermediate_data/postmeanBO_feature_cors_Mouse", 
                                          delim = "\t", escape_double = FALSE, 
                                          trim_ws = TRUE)

postmean_feature_cors_Dmel <- read_delim("../../intermediate_data/postmeanBO_feature_cors_Dmel", 
                                         delim = "\t", escape_double = FALSE, 
                                         trim_ws = TRUE)
postmean_feature_cors_Yeast <- read_delim("../../intermediate_data/postmeanBO_feature_cors_Yeast", 
                                          delim = "\t", escape_double = FALSE, 
                                          trim_ws = TRUE)

postmean_feature_cors_all<- rbind(postmean_feature_cors_Mouse, postmean_feature_cors_Dmel,postmean_feature_cors_Yeast)

# Global min and max Spearman rho across all features
rho_min <- min(postmean_feature_cors_all$spearman_rho, na.rm = TRUE)
rho_max <- max(postmean_feature_cors_all$spearman_rho, na.rm = TRUE)

# Symmetric around zero (common for correlations)
rho_limit <- max(abs(rho_min), abs(rho_max))

### S4A

### Mouse ### 

cat_cols <- c(
  Conservation     = "#1b9e77",
  gene_structure   = "#d95f02",
  aa_composition   = "#7570b3",
  ppi              = "#e7298a",
  GO_terms         = "#66a61e",
  Gene_expression  = "#e6ab02"
)

group1 <- c("GO_terms")
group2 <- c("Gene_expression")
group3 <- c("Conservation", "gene_structure", "aa_composition", "ppi")
group4 <- c("GO_terms","Gene_expression","Conservation", "gene_structure", "aa_composition", "ppi")

# Gene_expression

plot_spearman_heatmap(postmean_feature_cors_Mouse %>% filter(category %in% group2),rho_limit)

# GO_terms

plot_spearman_heatmap(postmean_feature_cors_Mouse %>% filter(category == "GO_terms"), rho_limit)

# Rest together

plot_spearman_heatmap(postmean_feature_cors_Mouse %>% filter(category %in% group3),rho_limit)

### S4B

### Same with postmean_feature_cors_Dmel

### S4C

### Same with postmean_feature_cors_Yeast



################ ################ ################ ################ ################ 
################ FIGURE S5 ################ ################ ################ ################ 

Pivalues_polarized_GBbinned_all<- read_delim("../../intermediate_data/Pivalues_GBbinned_BO_all", 
                                             delim = "\t", escape_double = FALSE, 
                                             trim_ws = TRUE)

Pivalues_polarized_GBbinned_all$Species<- factor(Pivalues_polarized_GBbinned_all$Species, levels=c("Mouse","Dmel","Yeast"),labels = c("Mouse","Fruit fly","Yeast"))

Pivalues_polarized_GBbinned_all$Shetbin <- factor(Pivalues_polarized_GBbinned_all$Shetbin,
                                                  levels =  c("High constraint", "Moderate-high constraint", "Moderate-low constraint", "Low constraint"))

## 5A 

ggplot(Pivalues_polarized_GBbinned_all[Pivalues_polarized_GBbinned_all$Degeneracy=="4fold",], aes(x = Shetbin, y = Pi)) +
  geom_bar(position = "dodge", stat = "identity")+
  theme(panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"),
        strip.text.x = element_text(size = 10),axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = expression(pi[4]))+geom_errorbar(aes(ymax=BS_UB, ymin=BS_LB))+facet_wrap(. ~ Species) 

## 5B 

ggplot(Pivalues_polarized_GBbinned_all[Pivalues_polarized_GBbinned_all$Degeneracy=="0fold",], aes(x = Shetbin, y = Pi)) +
  geom_bar(position = "dodge", stat = "identity")+
  theme(panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"),
        strip.text.x = element_text(size = 10),axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = expression(pi[0]))+geom_errorbar(aes(ymax=BS_UB, ymin=BS_LB))+facet_wrap(. ~ Species) 

################ ################ ################ ################ ################ 
################ FIGURE S6 ################ ################ ################ ################ 

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
################ FIGURE S7 ################ ################ ################ ################ 

### Ex. to create barplots

go_summary <- merged_df_Mouse %>%
  group_by(shet_bin) %>%
  summarise(
    count_total = n(),
    count_with_go = sum(`GO:0006357`, na.rm = TRUE),
    proportion_with_go = count_with_go / count_total
  )

ggplot(go_summary, aes(x = shet_bin, y = proportion_with_go)) +
  geom_col(width = 0.6) +
  theme_minimal() +
  labs(x = "Shet bin", y = "Proportion of genes with GO:0006357") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")+xlab("")

################ ################ ################ ################ ################ 
################ FIGURE S8 ################ ################ ################ ################ 

Bestfits_all_discretizedDFE <- read_delim("../../intermediate_data/Bestfits_DFE_all_discretizedDFE", 
                                          delim = "\t", escape_double = FALSE, 
                                          trim_ws = TRUE)

Bestfits_all_discretizedDFE$Species<- factor(Bestfits_all_discretizedDFE$Species, levels=c("Mouse","Dmel","Yeast"),labels = c("Mouse","Fruit fly","Yeast"))


Bestfits_all_discretizedDFE$inference <- factor(Bestfits_all_discretizedDFE$inference,
                                                levels= c("full","High constraint", "Moderate-high constraint", "Moderate-low constraint", "Low constraint"))

ggplot(Bestfits_all_discretizedDFE,aes(x=Species,y=Value,fill=Bin))+geom_bar(stat="identity",position="dodge")+
  theme(panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"), legend.title = element_blank(),
        strip.text.y = element_text(size = 10),strip.text.x = element_text(size = 10))+ylab("fraction")+xlab("S")+
  scale_fill_brewer(palette = "RdGy")+facet_grid(.~inference)
