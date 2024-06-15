
####### To compare the relative ranking of top20 important features in each of the three conditions 
####### log phase, hypoxia and fold change in hypoxia
####### feature selected by log phase in combined selected model
####### leadered transcripts

library(tidyverse)
library(RColorBrewer)
library(scales)

####### get color scheme for ranking percentage
rankingPercent_palette <- c("#a63603", "#e6550d", "#fd8d3c", "#fdbe85", "#feedde")

####### get feature type annotation
###### 5'UTR
path = '../../feature/FeatureSet_byType_5pUTR_UTRrelated_leadered/'
files = list.files(path, pattern = '.txt')
setwd('../../feature/FeatureSet_byType_5pUTR_UTRrelated_leadered/')
table_list = lapply(files, read.table, header = TRUE)
features_UTR = do.call(cbind, table_list)
setwd('../../ml_classification/viz_results/')

###### CDS nucleotide
path = '../../feature/FeatureSet_byType_CDSnucleotide/'
files = list.files(path, pattern = '.txt')
setwd('../../feature/FeatureSet_byType_CDSnucleotide/')
table_list = lapply(files, read.table, header = TRUE)
features_CDSnucleotide = do.call(cbind, table_list)
setwd('../../ml_classification/viz_results/')

###### CDS MFE
path = '../../feature/FeatureSet_byType_CDSmfe/'
files = list.files(path, pattern = '.txt')
setwd('../../feature/FeatureSet_byType_CDSmfe/')
table_list = lapply(files, read.table, header = TRUE)
features_CDSmfe = do.call(cbind, table_list)
setwd('../../ml_classification/viz_results/')

###### Codon
path = '../../feature/FeatureSet_byType_Codon/'
files = list.files(path, pattern = '.txt')
setwd('../../feature/FeatureSet_byType_Codon/')
table_list = lapply(files, read.table, header = TRUE)
features_Codon = do.call(cbind, table_list)
setwd('../../ml_classification/viz_results/')

###### Translation
path = '../../feature/FeatureSet_byType_Translation_leadered/'
files = list.files(path, pattern = '.txt')
setwd('../../feature/FeatureSet_byType_Translation_leadered/')
table_list = lapply(files, read.table, header = TRUE)
features_Translation = do.call(cbind, table_list)
setwd('../../ml_classification/viz_results/')

###### Others
##### log phase
path = '../../feature/FeatureSet_byType_Others_logPhase/'
files = list.files(path, pattern = '.txt')
setwd('../../feature/FeatureSet_byType_Others_logPhase/')
table_list = lapply(files, read.table, header = TRUE)
features_Others_logPhase = do.call(cbind, table_list)
setwd('../../ml_classification/viz_results/')

##### hypoxia
path = '../../feature/FeatureSet_byType_Others_hypoxia/'
files = list.files(path, pattern = '.txt')
setwd('../../feature/FeatureSet_byType_Others_hypoxia/')
table_list = lapply(files, read.table, header = TRUE)
features_Others_hypoxia = do.call(cbind, table_list)
setwd('../../ml_classification/viz_results/')

UTR <- data.frame(colnames(features_UTR), "UTR")
colnames(UTR) <- c("feature", "type")
CDSnucleotide <- data.frame(colnames(features_CDSnucleotide), "CDSnucleotide")
colnames(CDSnucleotide) <- c("feature", "type")
CDSmfe <- data.frame(colnames(features_CDSmfe), "CDSmfe")
colnames(CDSmfe) <- c("feature", "type")
Codon <- data.frame(colnames(features_Codon), "Codon")
colnames(Codon) <- c("feature", "type")
Translation <- data.frame(colnames(features_Translation), "Translation")
colnames(Translation) <- c("feature", "type")
other_p1 <- data.frame(colnames(features_Others_logPhase), "other")
colnames(other_p1) <- c("feature", "type")
other_p2 <- data.frame(colnames(features_Others_hypoxia), "other")
colnames(other_p2) <- c("feature", "type")

###### combine all feature type
feature_completeList <- bind_rows(UTR, CDSnucleotide, CDSmfe, Codon, Translation, other_p1, other_p2) %>% distinct()

####### get important features
###### log phase
logPhase_leadered_featureSort <- read.table('../ml_features_halfLifeClass/HalfLifeCls_logPhase_combinedSelected_leadered_byLeaderedLogPhase_jointlyByLogPhaseHypoxia_feature.txt',
                                            col.names = c("feature", "Gini_mean", "Gini_sd"))
logPhase_leadered_featureSort_type <- logPhase_leadered_featureSort %>% 
  inner_join(feature_completeList, by = "feature") %>% 
  arrange(desc(Gini_mean)) %>% 
  mutate(ranking_p1 = 1:n()) %>% 
  mutate(ranking = (ranking_p1 / nrow(logPhase_leadered_featureSort)) * 100) %>% 
  select(-Gini_mean, -Gini_sd, -ranking_p1)

###### hypoxia
hypoxia_leadered_featureSort <- read.table('../ml_features_halfLifeClass/HalfLifeCls_hypoxia_combinedSelected_leadered_byLeaderedLogPhase_jointlyByLogPhaseHypoxia_feature.txt',
                                           col.names = c("feature", "Gini_mean", "Gini_sd"))
hypoxia_leadered_featureSort_type <- hypoxia_leadered_featureSort %>% 
  inner_join(feature_completeList, by = "feature") %>% 
  arrange(desc(Gini_mean)) %>% 
  mutate(ranking_p1 = 1:n()) %>% 
  mutate(ranking = (ranking_p1 / nrow(hypoxia_leadered_featureSort)) * 100) %>% 
  select(-Gini_mean, -Gini_sd, -ranking_p1)

###### fold change in hypoxia
hypoxia2LogPhase_leadered_featureSort <- read.table('../ml_features_halfLifeFcClass/HalfLifeFcCls_hypoxiaToLogPhase_combinedSelected_leadered_byLeaderedLogPhase_jointlyByLogPhaseFoldChange_feature.txt',
                                                    col.names = c("feature", "Gini_mean", "Gini_sd"))

hypoxia2LogPhase_leadered_featureSort_type <- hypoxia2LogPhase_leadered_featureSort %>% 
  inner_join(feature_completeList, by = "feature") %>% 
  arrange(desc(Gini_mean)) %>% 
  mutate(ranking_p1 = 1:n()) %>% 
  mutate(ranking = (ranking_p1 / nrow(hypoxia2LogPhase_leadered_featureSort)) * 100) %>% 
  select(-Gini_mean, -Gini_sd, -ranking_p1)

####### combine important features
###### combined with top 20 features (horizontal figures)
feature_toCompareTop20 <- bind_rows(logPhase_leadered_featureSort_type[1:20, ], hypoxia_leadered_featureSort_type[1:20, ], hypoxia2LogPhase_leadered_featureSort_type[1:20, ])
feature_toCompareTop20_unique <- feature_toCompareTop20 %>% 
  select(-ranking) %>% 
  distinct() %>% 
  filter(feature != "initialAbundance_hypoxia" & feature != "initialAbundance_logPhase")

feature_toCompareTop20_unique_initialAbundance <- logPhase_leadered_featureSort_type %>% 
  filter(feature == "initialAbundance_logPhase") %>% 
  rename(feature_logPhase = feature, logPhase_ranking = ranking) %>% 
  select(-type) %>% 
  mutate(feature = "initialAbundance_hypoxia") %>% 
  inner_join(hypoxia_leadered_featureSort_type, by = "feature") %>% 
  rename(hypoxia_ranking = ranking) %>% 
  mutate(feature = "initialAbundance_logPhase") %>% 
  inner_join(hypoxia2LogPhase_leadered_featureSort_type, by = "feature") %>% 
  rename(hypoxia2LogPhase_ranking = ranking) %>% 
  mutate(feature = "initialAbundance") %>% 
  select(-feature_logPhase) %>% 
  mutate(type = coalesce(type.x, type.y)) %>% 
  select(feature, type, logPhase_ranking, hypoxia_ranking, hypoxia2LogPhase_ranking)

feature_toCompareTop20_unique_withRanking_p1 <- feature_toCompareTop20_unique %>% 
  left_join(logPhase_leadered_featureSort_type, by = "feature") %>% 
  rename(logPhase_ranking = ranking) %>% 
  mutate(type_p1 = coalesce(type.x, type.y)) %>% 
  left_join(hypoxia_leadered_featureSort_type, by = "feature") %>%
  rename(hypoxia_ranking = ranking) %>% 
  mutate(type_p2 = coalesce(type_p1, type)) %>% 
  left_join(hypoxia2LogPhase_leadered_featureSort_type, by = "feature") %>% 
  rename(hypoxia2LogPhase_ranking = ranking) %>% 
  select(feature, type_p2, logPhase_ranking, hypoxia_ranking, hypoxia2LogPhase_ranking) %>% 
  rename(type = type_p2)

##### IF needed, adding fake features to match the number of feature for leaderless
#### for viz purpose only, to make same size square with same number of features
### depending on the number of unique features in each run, it might not be necessary OR the number of fake features might be changed
feature_toCompareTop20_unique_withRanking_p2 <- data.frame(feature = c("initialAbundance_fake1", "initialAbundance_fake2"),
                                                           type = "other", logPhase_ranking = 2.888087,
                                                           hypoxia_ranking = 1.083032, hypoxia2LogPhase_ranking = 2.166065)

feature_toCompareTop20_unique_withRanking <- feature_toCompareTop20_unique_withRanking_p1 %>%
  bind_rows(feature_toCompareTop20_unique_initialAbundance) %>%
  bind_rows(feature_toCompareTop20_unique_withRanking_p2) %>%
  mutate(type = factor(type, levels = c("UTR", "CDSnucleotide", "CDSmfe", "Codon", "Translation", "other"))) %>%
  arrange(type, feature) %>%
  mutate(feature = factor(feature, levels = feature)) %>%
  gather("condition", "ranking", -type, -feature) %>%
  mutate(condition = factor(condition, levels = c("hypoxia2LogPhase_ranking", "hypoxia_ranking", "logPhase_ranking")))

####### viz of relative feature ranking comparison
###### combined with top 20 features
ggplot(feature_toCompareTop20_unique_withRanking, aes(feature, condition)) + 
  geom_tile(aes(fill = ranking), colour = "white", size = 0.7) + 
  scale_fill_gradientn(colours = rankingPercent_palette) + 
  theme_minimal() + 
  theme(
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
    axis.title = element_blank(), 
    panel.grid = element_blank()) + 
  coord_equal()
ggsave("./viz_featureRanks_toCompareCondition/featureRankingTop20combined_toCompareCondition_combinedSelectedByLogPhase_leadered_jointlyByLogPhaseHypoxiaFoldChange.png", width = 70, height = 20, units = "cm", dpi = 600)
