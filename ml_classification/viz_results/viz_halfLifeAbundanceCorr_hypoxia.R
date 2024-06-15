
####### viz correlation between half-life and abundance for leadered and leaderless transcript
####### both half-life and abundance are in log2 scale
####### hypxoa

library(tidyverse)

####### get half-life values
hl_CDS_complete_class <- read.table('../../degradation_classes/half_life/CDS_halfLife_completeClass.txt', header = TRUE)

hl_CDS_hypoxia <- hl_CDS_complete_class %>% 
  select(gene, HalfLife_vals_hypoxia) %>% 
  drop_na() %>% 
  rename(HalfLife_vals = HalfLife_vals_hypoxia)

####### leadered transcripts
###### get feature table
featureTable_leadered <- read.csv('../../feature/FeatureTables/featureTable_combinedSelected_halfLifeClass/HalfLifeCls_hypoxia_combinedSelected_leadered_byLeaderedLogPhase_jointlyByLogPhaseHypoxia.csv')

###### viz correlation of half-life values and initial abundance
hl_CDS_hypoxia_initialAbundance_leadered <- hl_CDS_hypoxia %>% 
  inner_join(featureTable_leadered, by = "gene") %>% 
  select(HalfLife_vals, initialAbundance_hypoxia) %>% 
  mutate(initialAbundance_hypoxia_log2 = log2(initialAbundance_hypoxia)) %>% 
  mutate(HalfLife_vals_log2 = log2(HalfLife_vals))
cor.test(hl_CDS_hypoxia_initialAbundance_leadered$initialAbundance_hypoxia_log2, hl_CDS_hypoxia_initialAbundance_leadered$HalfLife_vals_log2, method = "spearman")

hl_CDS_hypoxia_initialAbundance_leadered %>% 
  ggplot(aes(x = HalfLife_vals_log2, y = initialAbundance_hypoxia_log2)) + 
  geom_point(size = 4, fill = "#83B7DF",
             color = "white",
             shape = 21) +
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE, colour = "#EB636A") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.text = element_text(size = 25, face = "bold"),
        axis.title =  element_text(size = 15, face = "bold"),
        legend.position = "None") +
  xlim(2.5, 12) +
  ylim(5, 22)
ggsave("./viz_halfLifeAbundanceCorr/featureCorrelationWithHalfLife_toCompareCondition_HalfLifeCls_hypoxia_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseHypoxia_leadered_initialAbundance_hypoxia.png", width = 15, height = 10, units = "cm", dpi = 600)

model <- lm(hl_CDS_hypoxia_initialAbundance_leadered$initialAbundance_hypoxia_log2 ~ hl_CDS_hypoxia_initialAbundance_leadered$HalfLife_vals_log2)
summary(model)
summary(model)$r.squared

####### leaderless transcripts
###### get feature table
featureTable_leaderless <- read.csv('../../feature/FeatureTables/featureTable_combinedSelected_halfLifeClass/HalfLifeCls_hypoxia_combinedSelected_leaderless_byLeaderlessLogPhase_jointlyByLogPhaseHypoxia.csv')

###### viz correlation of half-life values and initial abundance
hl_CDS_hypoxia_initialAbundance_leaderless <- hl_CDS_hypoxia %>% 
  inner_join(featureTable_leaderless, by = "gene") %>% 
  select(HalfLife_vals, initialAbundance_hypoxia) %>% 
  mutate(initialAbundance_hypoxia_log2 = log2(initialAbundance_hypoxia)) %>% 
  mutate(HalfLife_vals_log2 = log2(HalfLife_vals))
cor.test(hl_CDS_hypoxia_initialAbundance_leaderless$initialAbundance_hypoxia_log2, hl_CDS_hypoxia_initialAbundance_leaderless$HalfLife_vals_log2, method = "spearman")

hl_CDS_hypoxia_initialAbundance_leaderless %>% 
  ggplot(aes(x = HalfLife_vals_log2, y = initialAbundance_hypoxia_log2)) + 
  geom_point(size = 4, fill = "#83B7DF",
             color = "white",
             shape = 21) +
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE, colour = "#EB636A") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.text = element_text(size = 25, face = "bold"),
        axis.title =  element_text(size = 15, face = "bold"),
        legend.position = "None") +
  xlim(2.5, 12) +
  ylim(5, 22)
ggsave("./viz_halfLifeAbundanceCorr/featureCorrelationWithHalfLife_toCompareCondition_HalfLifeCls_hypoxia_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseHypoxia_leaderless_initialAbundance_hypoxia.png", width = 15, height = 10, units = "cm", dpi = 600)

model <- lm(hl_CDS_hypoxia_initialAbundance_leaderless$initialAbundance_hypoxia_log2 ~ hl_CDS_hypoxia_initialAbundance_leaderless$HalfLife_vals_log2)
summary(model)
summary(model)$r.squared
