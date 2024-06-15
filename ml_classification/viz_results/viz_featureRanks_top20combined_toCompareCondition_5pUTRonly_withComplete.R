
####### To compare the relative ranking of top20 important features in each of the three conditions 
####### log phase, hypoxia and fold change in hypoxia
####### using only 5'UTR features for leadered transcripts
####### with complete 5'UTR features

library(tidyverse)
library(RColorBrewer)
library(scales)

####### get color scheme for ranking percentage
rankingPercent_palette <- c("#a63603", "#e6550d", "#fd8d3c", "#fdbe85", "#feedde")

####### get feature type annotation
###### 5'UTR complete
path = '../../feature/FeatureSet_byType_5pUTR_leadered/'
files = list.files(path, pattern = '.txt')
setwd('../../feature/FeatureSet_byType_5pUTR_leadered/')
table_list = lapply(files, read.table, header = TRUE)
features_UTRcomplete = do.call(cbind, table_list)
setwd('../../ml_classification/viz_results/')

UTRcomplete <- data.frame(colnames(features_UTRcomplete), "UTRcomplete")
colnames(UTRcomplete) <- c("feature", "type")

####### get important features
###### log phase
logPhase_UTRcomplete_featureSort <- read.table('../ml_features_halfLifeClass/HalfLifeCls_logPhase_nonSelected_leadered_byType_5pUTR_complete_feature.txt',
                                            col.names = c("feature", "Gini_mean", "Gini_sd"))
logPhase_UTRcomplete_featureSort_type <- logPhase_UTRcomplete_featureSort %>% 
  inner_join(UTRcomplete, by = "feature") %>% 
  arrange(desc(Gini_mean)) %>% 
  mutate(ranking_p1 = 1:n()) %>% 
  mutate(ranking = (ranking_p1 / nrow(logPhase_UTRcomplete_featureSort)) * 100) %>% 
  select(-Gini_mean, -Gini_sd, -ranking_p1)

###### hypoxia
hypoxia_UTRcomplete_featureSort <- read.table('../ml_features_halfLifeClass/HalfLifeCls_hypoxia_nonSelected_leadered_byType_5pUTR_complete_feature.txt',
                                               col.names = c("feature", "Gini_mean", "Gini_sd"))
hypoxia_UTRcomplete_featureSort_type <- hypoxia_UTRcomplete_featureSort %>% 
  inner_join(UTRcomplete, by = "feature") %>% 
  arrange(desc(Gini_mean)) %>% 
  mutate(ranking_p1 = 1:n()) %>% 
  mutate(ranking = (ranking_p1 / nrow(hypoxia_UTRcomplete_featureSort)) * 100) %>% 
  select(-Gini_mean, -Gini_sd, -ranking_p1)

###### fold change in hypoxia
hypoxia2LogPhase_UTRcomplete_featureSort <- read.table('../ml_features_halfLifeFcClass/HalfLifeFcCls_hypoxiaToLogPhase_nonSelected_leadered_byType_5pUTR_complete_feature.txt',
                                               col.names = c("feature", "Gini_mean", "Gini_sd"))
hypoxia2LogPhase_UTRcomplete_featureSort_type <- hypoxia2LogPhase_UTRcomplete_featureSort %>% 
  inner_join(UTRcomplete, by = "feature") %>% 
  arrange(desc(Gini_mean)) %>% 
  mutate(ranking_p1 = 1:n()) %>% 
  mutate(ranking = (ranking_p1 / nrow(hypoxia2LogPhase_UTRcomplete_featureSort)) * 100) %>% 
  select(-Gini_mean, -Gini_sd, -ranking_p1)

####### combine important features
###### combined with top 20 features (horizontal figures)
feature_toCompareTop20 <- bind_rows(logPhase_UTRcomplete_featureSort_type[1:20, ], hypoxia_UTRcomplete_featureSort_type[1:20, ], hypoxia2LogPhase_UTRcomplete_featureSort_type[1:20, ])
feature_toCompareTop20_unique <- feature_toCompareTop20 %>% 
  select(-ranking) %>% 
  distinct()

feature_toCompareTop20_unique_withRanking_p1 <- feature_toCompareTop20_unique %>% 
  left_join(logPhase_UTRcomplete_featureSort_type, by = "feature") %>% 
  rename(logPhase_ranking = ranking) %>% 
  mutate(type_p1 = coalesce(type.x, type.y)) %>% 
  left_join(hypoxia_UTRcomplete_featureSort_type, by = "feature") %>%
  rename(hypoxia_ranking = ranking) %>% 
  mutate(type_p2 = coalesce(type_p1, type)) %>% 
  left_join(hypoxia2LogPhase_UTRcomplete_featureSort_type, by = "feature") %>% 
  rename(hypoxia2LogPhase_ranking = ranking) %>% 
  select(feature, type_p2, logPhase_ranking, hypoxia_ranking, hypoxia2LogPhase_ranking) %>% 
  rename(type = type_p2)

##### In UTR related only features code, if necessary, added fake features in models to match the number of feature here for 5'UTR complete features
feature_toCompareTop20_unique_withRanking <- feature_toCompareTop20_unique_withRanking_p1 %>%
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
ggsave("./viz_featureRanks_toCompareCondition/featureRankingTop20combined_toCompareCondition_nonSelected_leadered_byType_5pUTRcomplete.png", width = 70, height = 20, units = "cm", dpi = 600)
