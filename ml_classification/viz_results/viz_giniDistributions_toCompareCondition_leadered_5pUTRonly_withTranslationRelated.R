
####### To compare the Gini score distributions within 5'UTR
####### secondary structure features vs. (SD) sequence related features
####### Gini scores from model using only translation related features in log phase

library(tidyverse)
library(ggbeeswarm)
library(RColorBrewer)
library(scales)

####### get important features
logPhase_TranslationRelated_featureSort <- read.table('../ml_features_halfLifeClass/HalfLifeCls_logPhase_nonSelected_leadered_byType_5pUTR_TranslationRelated_feature.txt',
                                              col.names = c("feature", "Gini_mean", "Gini_sd"))

####### split by secondary structure vs. SD sequence related
logPhase_TranslationRelated_featureSort_secondaryStructure <- logPhase_TranslationRelated_featureSort %>% 
  filter(grepl('UnpairedProb', feature)) %>% 
  mutate(type = 'SecondaryStructure')

logPhase_TranslationRelated_featureSort_SDsequence <- logPhase_TranslationRelated_featureSort %>% 
  anti_join(logPhase_TranslationRelated_featureSort_secondaryStructure, join_by(feature)) %>% 
  mutate(type = 'SDsequence')

wilcox.test(logPhase_TranslationRelated_featureSort_secondaryStructure$Gini_mean, 
       logPhase_TranslationRelated_featureSort_SDsequence$Gini_mean)

TranslationRelated_gini <- logPhase_TranslationRelated_featureSort_secondaryStructure %>% 
  bind_rows(logPhase_TranslationRelated_featureSort_SDsequence) %>% 
  mutate(type = factor(type))

####### viz of Gini score distribution
set.seed(7)
gini_median <- TranslationRelated_gini %>% 
  group_by(type) %>% 
  summarise(median_val = median(Gini_mean))

TranslationRelated_gini %>%
  ggplot(aes(type, Gini_mean)) +
  # geom_bar(stat = "summary", fill = NA, color = "#96999A",
  #          alpha = 0.8, fun = median, width = 0.5, size = 0.8) +
  # geom_hline(data = gini_median, aes(yintercept = median_val, col = type)) +
  geom_crossbar(data = gini_median, aes(x = type, y = median_val, ymin = median_val, ymax = median_val),
                width = 0.5, colour = "#96999A", size = 0.5) +
  geom_quasirandom(shape = 21, size = 2.5, color = "white",
                   fill = "#DC9d39", width = 0.2) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
    axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_blank(),
  )
ggsave("./viz_giniDistributions_toCompareFeatures_5pUTRonly/giniDistributions_avgs_5pUTRonly_TranslationRelated_logPhase.png", width = 7, height = 10, units = "cm", dpi = 600)
