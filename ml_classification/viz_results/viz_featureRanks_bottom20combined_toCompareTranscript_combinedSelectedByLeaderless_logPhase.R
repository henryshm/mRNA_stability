
####### To compare the relative ranking of bottom 20 features for the two transcript types
####### leadered and leaderless transcripts
####### for log phase, hypoxia and fold change in hypoxia
####### feature selected by leaderless in combined selected model
####### log phase

library(tidyverse)
library(RColorBrewer)
library(scales)

####### get color scheme for ranking percentage
rankingPercent_palette <- c("#a63603", "#e6550d", "#fd8d3c", "#fdbe85", "#feedde")

####### get feature type annotation
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
path = '../../feature/FeatureSet_byType_Translation_leaderless/'
files = list.files(path, pattern = '.txt')
setwd('../../feature/FeatureSet_byType_Translation_leaderless/')
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
feature_completeList <- bind_rows(CDSnucleotide, CDSmfe, Codon, Translation, other_p1, other_p2) %>% distinct()

####### get important features
###### leadered
logPhase_leadered_featureSort <- read.table('../ml_features_halfLifeClass/HalfLifeCls_logPhase_combinedSelected_leadered_byLeaderlessLogPhase_jointlyByLeaderedLeaderless_feature.txt',
                                            col.names = c("feature", "Gini_mean", "Gini_sd"))
logPhase_leadered_featureSort_type <- logPhase_leadered_featureSort %>% 
  inner_join(feature_completeList, by = "feature") %>% 
  arrange(desc(Gini_mean)) %>% 
  mutate(ranking_p1 = 1:n()) %>% 
  mutate(ranking = (ranking_p1 / nrow(logPhase_leadered_featureSort)) * 100) %>% 
  select(-Gini_mean, -Gini_sd, -ranking_p1)

###### leaderless
logPhase_leaderless_featureSort <- read.table('../ml_features_halfLifeClass/HalfLifeCls_logPhase_combinedSelected_leaderless_byLeaderlessLogPhase_jointlyByLeaderedLeaderless_feature.txt',
                                              col.names = c("feature", "Gini_mean", "Gini_sd"))
logPhase_leaderless_featureSort_type <- logPhase_leaderless_featureSort %>% 
  inner_join(feature_completeList, by = "feature") %>% 
  arrange(desc(Gini_mean)) %>% 
  mutate(ranking_p1 = 1:n()) %>% 
  mutate(ranking = (ranking_p1 / nrow(logPhase_leaderless_featureSort)) * 100) %>% 
  select(-Gini_mean, -Gini_sd, -ranking_p1)

####### combine important features
###### combined with bottom 20 features (horizontal figures)
feature_toCompareBottom20 <- bind_rows(logPhase_leadered_featureSort_type[(nrow(logPhase_leadered_featureSort_type) - 19) : nrow(logPhase_leadered_featureSort_type), ], 
                                    logPhase_leaderless_featureSort_type[(nrow(logPhase_leaderless_featureSort_type) -19) : nrow(logPhase_leaderless_featureSort_type), ])

feature_toCompareBottom20_unique <- feature_toCompareBottom20 %>% 
  select(-ranking) %>% 
  distinct() 

feature_toCompareBottom20_unique_withRanking_p1 <- feature_toCompareBottom20_unique %>% 
  left_join(logPhase_leadered_featureSort_type, by = "feature") %>% 
  rename(logPhase_leadered_ranking = ranking) %>% 
  mutate(type_p1 = coalesce(type.x, type.y)) %>% 
  left_join(logPhase_leaderless_featureSort_type, by = "feature") %>%
  rename(logPhase_leaderless_ranking = ranking)  %>% 
  select(feature, type, logPhase_leadered_ranking, logPhase_leaderless_ranking) %>%
  mutate(type = factor(type, levels = c("other", "Translation", "Codon", "CDSmfe", "CDSnucleotide"))) %>%
  arrange(desc(type), desc(feature))

##### IF needed, adding fake features to match the number of feature for hypoxia and log phase to hypoxia fold change
#### for viz purpose only, to make same size square with same number of features
### depending on the number of unique features in each run, it might not be necessary OR the number of fake features might be changed
## also added top 1 feature in top 20 to adjust color
feature_toCompareBottom20_unique_withRanking_p2 <- data.frame(feature = c("GGC_CDSnstop_fake1", "GGC_CDSnstop_fake2", "GGC_CDSnstop_fake3",
                                                                       "GGC_CDSnstop_fake4", "GGC_CDSnstop_fake5", "GGC_CDSnstop_fake6",
                                                                       "GGC_CDSnstop_fake7"),
                                                           type = "other", logPhase_leadered_ranking = 0.4098361,
                                                           logPhase_leaderless_ranking = 0.4098361)

feature_toCompareBottom20_unique_withRanking_p3 <- rbind(feature_toCompareBottom20_unique_withRanking_p1, feature_toCompareBottom20_unique_withRanking_p2)


feature_toCompareBottom20_unique_withRanking <- feature_toCompareBottom20_unique_withRanking_p3 %>%
  mutate(type = factor(type, levels = c("CDSnucleotide", "CDSmfe", "Codon", "Translation", "other"))) %>%
  arrange(type, feature) %>%
  mutate(feature = factor(feature, levels = feature)) %>%
  gather("condition", "ranking", -type, -feature) %>%
  mutate(condition = factor(condition, levels = c("logPhase_leaderless_ranking", "logPhase_leadered_ranking")))

####### viz of relative feature ranking comparison
###### combined with top 20 features
##### with y axis labels
ggplot(feature_toCompareBottom20_unique_withRanking, aes(feature, condition)) + 
  geom_tile(aes(fill = ranking), colour = "white", size = 0.7) + 
  scale_fill_gradientn(colours = rankingPercent_palette) + 
  theme_minimal() + 
  theme(
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
    axis.title = element_blank(), 
    panel.grid = element_blank()) +
  coord_equal()
ggsave("./viz_featureRanks_toCompareTranscript/featureRankingBottom20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_withYLabels.png", width = 70, height = 20, units = "cm", dpi = 600)

##### without y axis labels
ggplot(feature_toCompareBottom20_unique_withRanking, aes(feature, condition)) + 
  geom_tile(aes(fill = ranking), colour = "white", size = 0.7) + 
  scale_fill_gradientn(colours = rankingPercent_palette) + 
  theme_minimal() + 
  theme(
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    axis.title = element_blank(), 
    panel.grid = element_blank()) +
  coord_equal()
ggsave("./viz_featureRanks_toCompareTranscript/featureRankingBottom20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless.png", width = 70, height = 20, units = "cm", dpi = 600)
