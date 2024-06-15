
####### To compare models between log phase and hypoxia
####### models including individual property and combined selected models
####### feature selected by log phase in combined selected model
####### statistical tests using average fscores
####### leadered transcript

library(tidyverse)
library(RColorBrewer)
library(scales)

####### get color scheme for condition
condition_type <- c("#036D9C", "#83B7DF")

####### get averaged fscores
###### 5' UTR
fscore_avg_HalfLife_logPhase_5pUTR_complete <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_logPhase_nonSelected_leadered_byType_5pUTR_complete_fscore_avg.txt')
fscore_avg_HalfLife_hypoxia_5pUTR_complete <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_hypoxia_nonSelected_leadered_byType_5pUTR_complete_fscore_avg.txt')

fscore_avg_10repeats_logPhase_5pUTR_complete <- fscore_avg_HalfLife_logPhase_5pUTR_complete[, c(6:15)]
delta_logPhase_5pUTR_complete_p1 <- fscore_avg_10repeats_logPhase_5pUTR_complete[1, ] - fscore_avg_10repeats_logPhase_5pUTR_complete[2, ]
delta_logPhase_5pUTR_complete <- delta_logPhase_5pUTR_complete_p1 %>% 
  gather(repeats, fscore_diff) %>% 
  mutate(type = "logPhase")

fscore_avg_10repeats_hypoxia_5pUTR_complete <- fscore_avg_HalfLife_hypoxia_5pUTR_complete[, c(6:15)]
delta_hypoxia_5pUTR_complete_p1 <- fscore_avg_10repeats_hypoxia_5pUTR_complete[1, ] - fscore_avg_10repeats_hypoxia_5pUTR_complete[2, ]
delta_hypoxia_5pUTR_complete <- delta_hypoxia_5pUTR_complete_p1 %>% 
  gather(repeats, fscore_diff) %>% 
  mutate(type = "hypoxia")

wilcox.test(delta_logPhase_5pUTR_complete$fscore_diff, delta_hypoxia_5pUTR_complete$fscore_diff, paired = FALSE)

delta_5pUTR_complete <- bind_rows(delta_logPhase_5pUTR_complete, delta_hypoxia_5pUTR_complete)

delta_5pUTR_complete %>% 
  mutate(type = factor(type, level = c("logPhase", "hypoxia"))) %>% 
  group_by(type) %>% 
  summarise(
    sd = sd(fscore_diff), 
    delta = mean(fscore_diff)
  ) %>% 
  ggplot(aes(type, delta)) + 
  scale_colour_manual(values = condition_type) +
  geom_pointrange(aes(ymin = delta - sd, ymax = delta + sd, color = type),
                  linewidth = 2.3, size = 2.7) + 
  ylim(-0.05, 0.2) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "#E8EBEC", linewidth = 1.5),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.text = element_text(size = 10, color = "#96999A"),
    axis.text.x = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.key.size = unit(1, 'cm')
  )
ggsave("./viz_deltaFscores_toCompareCondition/deltaFscores_avgs_toCompareConditionType_5pUTR_complete_leadered.png", width = 8, height = 20, units = "cm", dpi = 600)

###### CDS nucleotide
fscore_avg_HalfLife_logPhase_CDSnucleotide <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_logPhase_nonSelected_leadered_byType_CDSnucleotide_fscore_avg.txt')
fscore_avg_HalfLife_hypoxia_CDSnucleotide <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_hypoxia_nonSelected_leadered_byType_CDSnucleotide_fscore_avg.txt')

fscore_avg_10repeats_logPhase_CDSnucleotide <- fscore_avg_HalfLife_logPhase_CDSnucleotide[, c(6:15)]
delta_logPhase_CDSnucleotide_p1 <- fscore_avg_10repeats_logPhase_CDSnucleotide[1, ] - fscore_avg_10repeats_logPhase_CDSnucleotide[2, ]
delta_logPhase_CDSnucleotide <- delta_logPhase_CDSnucleotide_p1 %>% 
  gather(repeats, fscore_diff) %>% 
  mutate(type = "logPhase")

fscore_avg_10repeats_hypoxia_CDSnucleotide <- fscore_avg_HalfLife_hypoxia_CDSnucleotide[, c(6:15)]
delta_hypoxia_CDSnucleotide_p1 <- fscore_avg_10repeats_hypoxia_CDSnucleotide[1, ] - fscore_avg_10repeats_hypoxia_CDSnucleotide[2, ]
delta_hypoxia_CDSnucleotide <- delta_hypoxia_CDSnucleotide_p1 %>% 
  gather(repeats, fscore_diff) %>% 
  mutate(type = "hypoxia")

wilcox.test(delta_logPhase_CDSnucleotide$fscore_diff, delta_hypoxia_CDSnucleotide$fscore_diff, paired = FALSE)

delta_CDSnucleotide <- bind_rows(delta_logPhase_CDSnucleotide, delta_hypoxia_CDSnucleotide)

delta_CDSnucleotide %>% 
  mutate(type = factor(type, level = c("logPhase", "hypoxia"))) %>% 
  group_by(type) %>% 
  summarise(
    sd = sd(fscore_diff), 
    delta = mean(fscore_diff)
  ) %>% 
  ggplot(aes(type, delta)) + 
  scale_colour_manual(values = condition_type) +
  geom_pointrange(aes(ymin = delta - sd, ymax = delta + sd, color = type),
                  linewidth = 2.3, size = 2.7) + 
  ylim(-0.05, 0.2) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "#E8EBEC", linewidth = 1.5),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.text = element_text(size = 10, color = "#96999A"),
    axis.text.x = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.key.size = unit(1, 'cm')
  )
ggsave("./viz_deltaFscores_toCompareCondition/deltaFscores_avgs_toCompareConditionType_CDSnucleotide_leadered.png", width = 8, height = 20, units = "cm", dpi = 600)

###### CDS MFE
fscore_avg_HalfLife_logPhase_CDSmfe <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_logPhase_nonSelected_leadered_byType_CDSmfe_fscore_avg.txt')
fscore_avg_HalfLife_hypoxia_CDSmfe <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_hypoxia_nonSelected_leadered_byType_CDSmfe_fscore_avg.txt')

fscore_avg_10repeats_logPhase_CDSmfe <- fscore_avg_HalfLife_logPhase_CDSmfe[, c(6:15)]
delta_logPhase_CDSmfe_p1 <- fscore_avg_10repeats_logPhase_CDSmfe[1, ] - fscore_avg_10repeats_logPhase_CDSmfe[2, ]
delta_logPhase_CDSmfe <- delta_logPhase_CDSmfe_p1 %>% 
  gather(repeats, fscore_diff) %>% 
  mutate(type = "logPhase")

fscore_avg_10repeats_hypoxia_CDSmfe <- fscore_avg_HalfLife_hypoxia_CDSmfe[, c(6:15)]
delta_hypoxia_CDSmfe_p1 <- fscore_avg_10repeats_hypoxia_CDSmfe[1, ] - fscore_avg_10repeats_hypoxia_CDSmfe[2, ]
delta_hypoxia_CDSmfe <- delta_hypoxia_CDSmfe_p1 %>% 
  gather(repeats, fscore_diff) %>% 
  mutate(type = "hypoxia")

wilcox.test(delta_logPhase_CDSmfe$fscore_diff, delta_hypoxia_CDSmfe$fscore_diff, paired = FALSE)

delta_CDSmfe <- bind_rows(delta_logPhase_CDSmfe, delta_hypoxia_CDSmfe)

delta_CDSmfe %>% 
  mutate(type = factor(type, level = c("logPhase", "hypoxia"))) %>% 
  group_by(type) %>% 
  summarise(
    sd = sd(fscore_diff), 
    delta = mean(fscore_diff)
  ) %>% 
  ggplot(aes(type, delta)) + 
  scale_colour_manual(values = condition_type) +
  geom_pointrange(aes(ymin = delta - sd, ymax = delta + sd, color = type),
                  linewidth = 2.3, size = 2.7) + 
  ylim(-0.05, 0.2) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "#E8EBEC", linewidth = 1.5),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.text = element_text(size = 10, color = "#96999A"),
    axis.text.x = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.key.size = unit(1, 'cm')
  )
ggsave("./viz_deltaFscores_toCompareCondition/deltaFscores_avgs_toCompareConditionType_CDSmfe_leadered.png", width = 8, height = 20, units = "cm", dpi = 600)

###### Codon
fscore_avg_HalfLife_logPhase_Codon <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_logPhase_nonSelected_leadered_byType_Codon_fscore_avg.txt')
fscore_avg_HalfLife_hypoxia_Codon <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_hypoxia_nonSelected_leadered_byType_Codon_fscore_avg.txt')

fscore_avg_10repeats_logPhase_Codon <- fscore_avg_HalfLife_logPhase_Codon[, c(6:15)]
delta_logPhase_Codon_p1 <- fscore_avg_10repeats_logPhase_Codon[1, ] - fscore_avg_10repeats_logPhase_Codon[2, ]
delta_logPhase_Codon <- delta_logPhase_Codon_p1 %>% 
  gather(repeats, fscore_diff) %>% 
  mutate(type = "logPhase")

fscore_avg_10repeats_hypoxia_Codon <- fscore_avg_HalfLife_hypoxia_Codon[, c(6:15)]
delta_hypoxia_Codon_p1 <- fscore_avg_10repeats_hypoxia_Codon[1, ] - fscore_avg_10repeats_hypoxia_Codon[2, ]
delta_hypoxia_Codon <- delta_hypoxia_Codon_p1 %>% 
  gather(repeats, fscore_diff) %>% 
  mutate(type = "hypoxia")

wilcox.test(delta_logPhase_Codon$fscore_diff, delta_hypoxia_Codon$fscore_diff, paired = FALSE)

delta_Codon <- bind_rows(delta_logPhase_Codon, delta_hypoxia_Codon)

delta_Codon %>% 
  mutate(type = factor(type, level = c("logPhase", "hypoxia"))) %>% 
  group_by(type) %>% 
  summarise(
    sd = sd(fscore_diff), 
    delta = mean(fscore_diff)
  ) %>% 
  ggplot(aes(type, delta)) + 
  scale_colour_manual(values = condition_type) +
  geom_pointrange(aes(ymin = delta - sd, ymax = delta + sd, color = type),
                  linewidth = 2.3, size = 2.7) + 
  ylim(-0.05, 0.2) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "#E8EBEC", linewidth = 1.5),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.text = element_text(size = 10, color = "#96999A"),
    axis.text.x = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.key.size = unit(1, 'cm')
  )
ggsave("./viz_deltaFscores_toCompareCondition/deltaFscores_avgs_toCompareConditionType_Codon_leadered.png", width = 8, height = 20, units = "cm", dpi = 600)

###### Translation
fscore_avg_HalfLife_logPhase_Translation <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_logPhase_nonSelected_leadered_byType_TranslationWTribo_fscore_avg.txt')
fscore_avg_HalfLife_hypoxia_Translation <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_hypoxia_nonSelected_leadered_byType_TranslationWTribo_fscore_avg.txt')

fscore_avg_10repeats_logPhase_Translation <- fscore_avg_HalfLife_logPhase_Translation[, c(6:15)]
delta_logPhase_Translation_p1 <- fscore_avg_10repeats_logPhase_Translation[1, ] - fscore_avg_10repeats_logPhase_Translation[2, ]
delta_logPhase_Translation <- delta_logPhase_Translation_p1 %>% 
  gather(repeats, fscore_diff) %>% 
  mutate(type = "logPhase")

fscore_avg_10repeats_hypoxia_Translation <- fscore_avg_HalfLife_hypoxia_Translation[, c(6:15)]
delta_hypoxia_Translation_p1 <- fscore_avg_10repeats_hypoxia_Translation[1, ] - fscore_avg_10repeats_hypoxia_Translation[2, ]
delta_hypoxia_Translation <- delta_hypoxia_Translation_p1 %>% 
  gather(repeats, fscore_diff) %>% 
  mutate(type = "hypoxia")

wilcox.test(delta_logPhase_Translation$fscore_diff, delta_hypoxia_Translation$fscore_diff, paired = FALSE)

delta_Translation <- bind_rows(delta_logPhase_Translation, delta_hypoxia_Translation)

delta_Translation %>% 
  mutate(type = factor(type, level = c("logPhase", "hypoxia"))) %>% 
  group_by(type) %>% 
  summarise(
    sd = sd(fscore_diff), 
    delta = mean(fscore_diff)
  ) %>% 
  ggplot(aes(type, delta)) + 
  scale_colour_manual(values = condition_type) +
  geom_pointrange(aes(ymin = delta - sd, ymax = delta + sd, color = type),
                  linewidth = 2.3, size = 2.7) + 
  ylim(-0.05, 0.2) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "#E8EBEC", linewidth = 1.5),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.text = element_text(size = 10, color = "#96999A"),
    axis.text.x = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.key.size = unit(1, 'cm')
  )
ggsave("./viz_deltaFscores_toCompareCondition/deltaFscores_avgs_toCompareConditionType_TranslationWTribo_leadered.png", width = 8, height = 20, units = "cm", dpi = 600)

###### Others
fscore_avg_HalfLife_logPhase_Others <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_logPhase_nonSelected_leadered_byType_Others_fscore_avg.txt')
fscore_avg_HalfLife_hypoxia_Others <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_hypoxia_nonSelected_leadered_byType_Others_fscore_avg.txt')

fscore_avg_10repeats_logPhase_Others <- fscore_avg_HalfLife_logPhase_Others[, c(6:15)]
delta_logPhase_Others_p1 <- fscore_avg_10repeats_logPhase_Others[1, ] - fscore_avg_10repeats_logPhase_Others[2, ]
delta_logPhase_Others <- delta_logPhase_Others_p1 %>% 
  gather(repeats, fscore_diff) %>% 
  mutate(type = "logPhase")

fscore_avg_10repeats_hypoxia_Others <- fscore_avg_HalfLife_hypoxia_Others[, c(6:15)]
delta_hypoxia_Others_p1 <- fscore_avg_10repeats_hypoxia_Others[1, ] - fscore_avg_10repeats_hypoxia_Others[2, ]
delta_hypoxia_Others <- delta_hypoxia_Others_p1 %>% 
  gather(repeats, fscore_diff) %>% 
  mutate(type = "hypoxia")

wilcox.test(delta_logPhase_Others$fscore_diff, delta_hypoxia_Others$fscore_diff, paired = FALSE)

delta_Others <- bind_rows(delta_logPhase_Others, delta_hypoxia_Others)

delta_Others %>% 
  mutate(type = factor(type, level = c("logPhase", "hypoxia"))) %>% 
  group_by(type) %>% 
  summarise(
    sd = sd(fscore_diff), 
    delta = mean(fscore_diff)
  ) %>% 
  ggplot(aes(type, delta)) + 
  scale_colour_manual(values = condition_type) +
  geom_pointrange(aes(ymin = delta - sd, ymax = delta + sd, color = type),
                  linewidth = 2.3, size = 2.7) + 
  ylim(-0.05, 0.2) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "#E8EBEC", linewidth = 1.5),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.text = element_text(size = 10, color = "#96999A"),
    axis.text.x = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.key.size = unit(1, 'cm')
  )
ggsave("./viz_deltaFscores_toCompareCondition/deltaFscores_avgs_toCompareConditionType_Others_leadered.png", width = 8, height = 20, units = "cm", dpi = 600)

###### combined selected
fscore_avg_HalfLife_logPhase_CombinedSelected <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_logPhase_combinedSelected_leadered_byLeaderedLogPhase_jointlyByLogPhaseHypoxia_fscore_avg.txt')
fscore_avg_HalfLife_hypoxia_CombinedSelected <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_hypoxia_combinedSelected_leadered_byLeaderedLogPhase_jointlyByLogPhaseHypoxia_fscore_avg.txt')

fscore_avg_10repeats_logPhase_CombinedSelected <- fscore_avg_HalfLife_logPhase_CombinedSelected[, c(6:15)]
delta_logPhase_CombinedSelected_p1 <- fscore_avg_10repeats_logPhase_CombinedSelected[1, ] - fscore_avg_10repeats_logPhase_CombinedSelected[2, ]
delta_logPhase_CombinedSelected <- delta_logPhase_CombinedSelected_p1 %>% 
  gather(repeats, fscore_diff) %>% 
  mutate(type = "logPhase")

fscore_avg_10repeats_hypoxia_CombinedSelected <- fscore_avg_HalfLife_hypoxia_CombinedSelected[, c(6:15)]
delta_hypoxia_CombinedSelected_p1 <- fscore_avg_10repeats_hypoxia_CombinedSelected[1, ] - fscore_avg_10repeats_hypoxia_CombinedSelected[2, ]
delta_hypoxia_CombinedSelected <- delta_hypoxia_CombinedSelected_p1 %>% 
  gather(repeats, fscore_diff) %>% 
  mutate(type = "hypoxia")

wilcox.test(delta_logPhase_CombinedSelected$fscore_diff, delta_hypoxia_CombinedSelected$fscore_diff, paired = FALSE)

delta_CombinedSelected <- bind_rows(delta_logPhase_CombinedSelected, delta_hypoxia_CombinedSelected)

delta_CombinedSelected %>% 
  mutate(type = factor(type, level = c("logPhase", "hypoxia"))) %>% 
  group_by(type) %>% 
  summarise(
    sd = sd(fscore_diff), 
    delta = mean(fscore_diff)
  ) %>% 
  ggplot(aes(type, delta)) + 
  scale_colour_manual(values = condition_type) +
  geom_pointrange(aes(ymin = delta - sd, ymax = delta + sd, color = type),
                  linewidth = 2.3, size = 2.7) + 
  ylim(-0.05, 0.2) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "#E8EBEC", linewidth = 1.5),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.text = element_text(size = 10, color = "#96999A"),
    axis.text.x = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.key.size = unit(1, 'cm')
  )
ggsave("./viz_deltaFscores_toCompareCondition/deltaFscores_avgs_toCompareConditionType_Combined_CombinedSelectedByLogPhase_leadered.png", width = 8, height = 20, units = "cm", dpi = 600)
