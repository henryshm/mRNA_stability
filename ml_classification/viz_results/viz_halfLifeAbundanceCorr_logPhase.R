
####### viz correlation between half-life and abundance for leadered and leaderless transcript
####### both half-life and abundance are in log2 scale
####### log phase

library(tidyverse)

####### get half-life values
hl_CDS_complete_class <- read.table('../../degradation_classes/half_life/CDS_halfLife_completeClass.txt', header = TRUE)

hl_CDS_logPhase <- hl_CDS_complete_class %>% 
  select(gene, HalfLife_vals_logPhase) %>% 
  drop_na() %>% 
  rename(HalfLife_vals = HalfLife_vals_logPhase)

####### leadered transcripts
###### get feature table
featureTable_leadered <- read.csv('../../feature/FeatureTables/featureTable_combinedSelected_halfLifeClass/HalfLifeCls_logPhase_combinedSelected_leadered_byLeaderedLogPhase_jointlyByLogPhaseHypoxia.csv')

###### viz correlation of half-life values and initial abundance
hl_CDS_logPhase_initialAbundance_leadered <- hl_CDS_logPhase %>% 
  inner_join(featureTable_leadered, by = "gene") %>% 
  select(HalfLife_vals, initialAbundance_logPhase) %>% 
  mutate(initialAbundance_logPhase_log2 = log2(initialAbundance_logPhase)) %>% 
  mutate(HalfLife_vals_log2 = log2(HalfLife_vals))
cor.test(hl_CDS_logPhase_initialAbundance_leadered$initialAbundance_logPhase_log2, hl_CDS_logPhase_initialAbundance_leadered$HalfLife_vals_log2, method = "spearman")

hl_CDS_logPhase_initialAbundance_leadered %>% 
  ggplot(aes(x = HalfLife_vals_log2, y = initialAbundance_logPhase_log2)) + 
  geom_point(size = 4, fill = "#036D9C",
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
  xlim(-1.5, 3) +
  ylim(5, 22)
ggsave("./viz_halfLifeAbundanceCorr/featureCorrelationWithHalfLife_toCompareCondition_HalfLifeCls_logPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseHypoxia_leadered_initialAbundance_logPhase.png", width = 15, height = 10, units = "cm", dpi = 600)

model <- lm(hl_CDS_logPhase_initialAbundance_leadered$initialAbundance_logPhase_log2 ~ hl_CDS_logPhase_initialAbundance_leadered$HalfLife_vals_log2)
summary(model)
summary(model)$r.squared

####### leaderless transcripts
###### get feature table
featureTable_leaderless <- read.csv('../../feature/FeatureTables/featureTable_combinedSelected_halfLifeClass/HalfLifeCls_logPhase_combinedSelected_leaderless_byLeaderlessLogPhase_jointlyByLogPhaseHypoxia.csv')

###### viz correlation of half-life values and initial abundance
hl_CDS_logPhase_initialAbundance_leaderless <- hl_CDS_logPhase %>% 
  inner_join(featureTable_leaderless, by = "gene") %>% 
  select(HalfLife_vals, initialAbundance_logPhase) %>% 
  mutate(initialAbundance_logPhase_log2 = log2(initialAbundance_logPhase)) %>% 
  mutate(HalfLife_vals_log2 = log2(HalfLife_vals))
cor.test(hl_CDS_logPhase_initialAbundance_leaderless$initialAbundance_logPhase_log2, hl_CDS_logPhase_initialAbundance_leaderless$HalfLife_vals_log2, method = "spearman")

hl_CDS_logPhase_initialAbundance_leaderless %>% 
  ggplot(aes(x = HalfLife_vals_log2, y = initialAbundance_logPhase_log2)) + 
  geom_point(size = 4, fill = "#036D9C",
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
  xlim(-1.5, 3) +
  ylim(5, 22)
ggsave("./viz_halfLifeAbundanceCorr/featureCorrelationWithHalfLife_toCompareCondition_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseHypoxia_leaderless_initialAbundance_logPhase.png", width = 15, height = 10, units = "cm", dpi = 600)

model <- lm(hl_CDS_logPhase_initialAbundance_leaderless$initialAbundance_logPhase_log2 ~ hl_CDS_logPhase_initialAbundance_leaderless$HalfLife_vals_log2)
summary(model)
summary(model)$r.squared
