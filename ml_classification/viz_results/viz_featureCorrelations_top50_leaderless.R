
####### To evaluate correlations among top important features
####### for log phase, hypoxia, half-life fold change in hypoxia
####### feature selected by log phase in combined selected model
####### leaderless transcript

library(tidyverse)
library(RColorBrewer)
library(scales)
library(GGally)

###### log phase
##### get feature values
logPhase_feature_value <- read.csv('../../feature/FeatureTables/featureTable_combinedSelected_halfLifeClass/HalfLifeCls_logPhase_combinedSelected_leaderless_byLeaderlessLogPhase_jointlyByLogPhaseHypoxia.csv')

##### get feautre Gini rankings
logPhase_feature_gini <- read.table('../ml_features_halfLifeClass/HalfLifeCls_logPhase_combinedSelected_leaderless_byLeaderlessLogPhase_jointlyByLogPhaseHypoxia_feature.txt',
                                    col.names = c("feature", "Gini_mean", "Gini_sd"))
logPhase_feature_gini_sort <- logPhase_feature_gini %>% 
  arrange(desc(Gini_mean)) %>% 
  head(50) %>% 
  pull(feature)

##### get feautre value of top50 important features
logPhase_feature_value_top50 <- select(logPhase_feature_value, one_of(logPhase_feature_gini_sort))

##### viz of top50 important features correlation
ggcorr(logPhase_feature_value_top50, hjust = 1, size = 3, 
       # color = "white",
       # legend.size = 40,
       layout.exp = 7,
       legend.position = "None",
       color = "grey50",
       # label = TRUE, label_size = 5, label_color = "white",
       method = c("all.obs", "spearman"))


###### hypoxia
##### get feature values
hypoxia_feature_value <- read.csv('../../feature/FeatureTables/featureTable_combinedSelected_halfLifeClass/HalfLifeCls_hypoxia_combinedSelected_leaderless_byLeaderlessLogPhase_jointlyByLogPhaseHypoxia.csv')

##### get feautre Gini rankings
hypoxia_feature_gini <- read.table('../ml_features_halfLifeClass/HalfLifeCls_hypoxia_combinedSelected_leaderless_byLeaderlessLogPhase_jointlyByLogPhaseHypoxia_feature.txt',
                                   col.names = c("feature", "Gini_mean", "Gini_sd"))
hypoxia_feature_gini_sort <- hypoxia_feature_gini %>% 
  arrange(desc(Gini_mean)) %>% 
  head(50) %>% 
  pull(feature)

##### get feautre value of top50 important features
hypoxia_feature_value_top50 <- select(hypoxia_feature_value, one_of(hypoxia_feature_gini_sort))

##### viz of top50 important features correlation
ggcorr(hypoxia_feature_value_top50, hjust = 1, size = 3, 
       # color = "white",
       # legend.size = 40,
       layout.exp = 7,
       legend.position = "None",
       color = "grey50",
       # label = TRUE, label_size = 5, label_color = "white",
       method = c("all.obs", "spearman"))


###### fold change in hypoxia
##### get feature values
hypoxia2LogPhase_feature_value <- read.csv('../../feature/FeatureTables/featureTable_combinedSelected_halfLifeFcClass/HalfLifeFcCls_hypoxiaToLogPhase_combinedSelected_leaderless_byLeaderlessLogPhase_jointlyByLogPhaseFoldChange.csv')

##### get feautre Gini rankings
hypoxia2LogPhase_feature_gini <- read.table('../ml_features_halfLifeFcClass/HalfLifeFcCls_hypoxiaToLogPhase_combinedSelected_leaderless_byLeaderlessLogPhase_jointlyByLogPhaseFoldChange_feature.txt',
                                            col.names = c("feature", "Gini_mean", "Gini_sd"))
hypoxia2LogPhase_feature_gini_sort <- hypoxia2LogPhase_feature_gini %>% 
  arrange(desc(Gini_mean)) %>% 
  head(50) %>% 
  pull(feature)

##### get feautre value of top50 important features
hypoxia2LogPhase_feature_value_top50 <- select(hypoxia2LogPhase_feature_value, one_of(hypoxia2LogPhase_feature_gini_sort))

##### viz of top50 important features correlation
ggcorr(hypoxia2LogPhase_feature_value_top50, hjust = 1, size = 3, 
       # color = "white",
       # legend.size = 40,
       layout.exp = 7,
       legend.position = "None",
       color = "grey50",
       # label = TRUE, label_size = 5, label_color = "white",
       method = c("all.obs", "spearman"))


