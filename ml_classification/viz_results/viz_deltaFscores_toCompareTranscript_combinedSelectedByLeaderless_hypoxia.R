
####### To compare models between leadered and leaderless transcript
####### models including individual property and combined selected models
####### feature selected by leaderless in combined selected model
####### statistical tests using average fscores
####### hypoxia

library(tidyverse)
library(RColorBrewer)
library(scales)

####### get color scheme for transcript type
transcript_type <- c("#EB636A", "#F1C36B")

####### get averaged fscores
###### CDS nucleotide
fscore_avg_HalfLife_leadered_CDSnucleotide <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_hypoxia_nonSelected_leadered_byType_CDSnucleotide_fscore_avg.txt')
fscore_avg_HalfLife_leaderless_CDSnucleotide <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_hypoxia_nonSelected_leaderless_byType_CDSnucleotide_fscore_avg.txt')

fscore_avg_10repeats_leadered_CDSnucleotide <- fscore_avg_HalfLife_leadered_CDSnucleotide[, c(6:15)]
delta_leadered_CDSnucleotide_p1 <- fscore_avg_10repeats_leadered_CDSnucleotide[1, ] - fscore_avg_10repeats_leadered_CDSnucleotide[2, ]
delta_leadered_CDSnucleotide <- delta_leadered_CDSnucleotide_p1 %>% 
  gather(repeats, fscore_diff) %>% 
  mutate(type = "leadered")

fscore_avg_10repeats_leaderless_CDSnucleotide <- fscore_avg_HalfLife_leaderless_CDSnucleotide[, c(6:15)]
delta_leaderless_CDSnucleotide_p1 <- fscore_avg_10repeats_leaderless_CDSnucleotide[1, ] - fscore_avg_10repeats_leaderless_CDSnucleotide[2, ]
delta_leaderless_CDSnucleotide <- delta_leaderless_CDSnucleotide_p1 %>% 
  gather(repeats, fscore_diff) %>% 
  mutate(type = "leaderless")

wilcox.test(delta_leadered_CDSnucleotide$fscore_diff, delta_leaderless_CDSnucleotide$fscore_diff, paired = FALSE)

delta_CDSnucleotide <- bind_rows(delta_leadered_CDSnucleotide, delta_leaderless_CDSnucleotide)

delta_CDSnucleotide %>% 
  mutate(type = factor(type, level = c("leadered", "leaderless"))) %>% 
  group_by(type) %>% 
  summarise(
    sd = sd(fscore_diff), 
    delta = mean(fscore_diff)
  ) %>% 
  ggplot(aes(type, delta)) + 
  scale_colour_manual(values = transcript_type) +
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
ggsave("./viz_deltaFscores_toCompareTranscript/deltaFscores_avgs_toCompareTranscriptType_CDSnucleotide_hypoxia.png", width = 8, height = 20, units = "cm", dpi = 600)

###### CDS MFE
fscore_avg_HalfLife_leadered_CDSmfe <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_hypoxia_nonSelected_leadered_byType_CDSmfe_fscore_avg.txt')
fscore_avg_HalfLife_leaderless_CDSmfe <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_hypoxia_nonSelected_leaderless_byType_CDSmfe_fscore_avg.txt')

fscore_avg_10repeats_leadered_CDSmfe <- fscore_avg_HalfLife_leadered_CDSmfe[, c(6:15)]
delta_leadered_CDSmfe_p1 <- fscore_avg_10repeats_leadered_CDSmfe[1, ] - fscore_avg_10repeats_leadered_CDSmfe[2, ]
delta_leadered_CDSmfe <- delta_leadered_CDSmfe_p1 %>% 
  gather(repeats, fscore_diff) %>% 
  mutate(type = "leadered")

fscore_avg_10repeats_leaderless_CDSmfe <- fscore_avg_HalfLife_leaderless_CDSmfe[, c(6:15)]
delta_leaderless_CDSmfe_p1 <- fscore_avg_10repeats_leaderless_CDSmfe[1, ] - fscore_avg_10repeats_leaderless_CDSmfe[2, ]
delta_leaderless_CDSmfe <- delta_leaderless_CDSmfe_p1 %>% 
  gather(repeats, fscore_diff) %>% 
  mutate(type = "leaderless")

wilcox.test(delta_leadered_CDSmfe$fscore_diff, delta_leaderless_CDSmfe$fscore_diff, paired = FALSE)

delta_CDSmfe <- bind_rows(delta_leadered_CDSmfe, delta_leaderless_CDSmfe)

delta_CDSmfe %>% 
  mutate(type = factor(type, level = c("leadered", "leaderless"))) %>% 
  group_by(type) %>% 
  summarise(
    sd = sd(fscore_diff), 
    delta = mean(fscore_diff)
  ) %>% 
  ggplot(aes(type, delta)) + 
  scale_colour_manual(values = transcript_type) +
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
ggsave("./viz_deltaFscores_toCompareTranscript/deltaFscores_avgs_toCompareTranscriptType_CDSmfe_hypoxia.png", width = 8, height = 20, units = "cm", dpi = 600)

###### Codon
fscore_avg_HalfLife_leadered_Codon <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_hypoxia_nonSelected_leadered_byType_Codon_fscore_avg.txt')
fscore_avg_HalfLife_leaderless_Codon <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_hypoxia_nonSelected_leaderless_byType_Codon_fscore_avg.txt')

fscore_avg_10repeats_leadered_Codon <- fscore_avg_HalfLife_leadered_Codon[, c(6:15)]
delta_leadered_Codon_p1 <- fscore_avg_10repeats_leadered_Codon[1, ] - fscore_avg_10repeats_leadered_Codon[2, ]
delta_leadered_Codon <- delta_leadered_Codon_p1 %>% 
  gather(repeats, fscore_diff) %>% 
  mutate(type = "leadered")

fscore_avg_10repeats_leaderless_Codon <- fscore_avg_HalfLife_leaderless_Codon[, c(6:15)]
delta_leaderless_Codon_p1 <- fscore_avg_10repeats_leaderless_Codon[1, ] - fscore_avg_10repeats_leaderless_Codon[2, ]
delta_leaderless_Codon <- delta_leaderless_Codon_p1 %>% 
  gather(repeats, fscore_diff) %>% 
  mutate(type = "leaderless")

wilcox.test(delta_leadered_Codon$fscore_diff, delta_leaderless_Codon$fscore_diff, paired = FALSE)

delta_Codon <- bind_rows(delta_leadered_Codon, delta_leaderless_Codon)

delta_Codon %>% 
  mutate(type = factor(type, level = c("leadered", "leaderless"))) %>% 
  group_by(type) %>% 
  summarise(
    sd = sd(fscore_diff), 
    delta = mean(fscore_diff)
  ) %>% 
  ggplot(aes(type, delta)) + 
  scale_colour_manual(values = transcript_type) +
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
ggsave("./viz_deltaFscores_toCompareTranscript/deltaFscores_avgs_toCompareTranscriptType_Codon_hypoxia.png", width = 8, height = 20, units = "cm", dpi = 600)

###### Translation
fscore_avg_HalfLife_leadered_Translation <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_hypoxia_nonSelected_leadered_byType_Translation_fscore_avg.txt')
fscore_avg_HalfLife_leaderless_Translation <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_hypoxia_nonSelected_leaderless_byType_Translation_fscore_avg.txt')

fscore_avg_10repeats_leadered_Translation <- fscore_avg_HalfLife_leadered_Translation[, c(6:15)]
delta_leadered_Translation_p1 <- fscore_avg_10repeats_leadered_Translation[1, ] - fscore_avg_10repeats_leadered_Translation[2, ]
delta_leadered_Translation <- delta_leadered_Translation_p1 %>% 
  gather(repeats, fscore_diff) %>% 
  mutate(type = "leadered")

fscore_avg_10repeats_leaderless_Translation <- fscore_avg_HalfLife_leaderless_Translation[, c(6:15)]
delta_leaderless_Translation_p1 <- fscore_avg_10repeats_leaderless_Translation[1, ] - fscore_avg_10repeats_leaderless_Translation[2, ]
delta_leaderless_Translation <- delta_leaderless_Translation_p1 %>% 
  gather(repeats, fscore_diff) %>% 
  mutate(type = "leaderless")

wilcox.test(delta_leadered_Translation$fscore_diff, delta_leaderless_Translation$fscore_diff, paired = FALSE)

delta_Translation <- bind_rows(delta_leadered_Translation, delta_leaderless_Translation)

delta_Translation %>% 
  mutate(type = factor(type, level = c("leadered", "leaderless"))) %>% 
  group_by(type) %>% 
  summarise(
    sd = sd(fscore_diff), 
    delta = mean(fscore_diff)
  ) %>% 
  ggplot(aes(type, delta)) + 
  scale_colour_manual(values = transcript_type) +
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
ggsave("./viz_deltaFscores_toCompareTranscript/deltaFscores_avgs_toCompareTranscriptType_Translation_hypoxia.png", width = 8, height = 20, units = "cm", dpi = 600)

###### Others
fscore_avg_HalfLife_leadered_Others <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_hypoxia_nonSelected_leadered_byType_Others_fscore_avg.txt')
fscore_avg_HalfLife_leaderless_Others <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_hypoxia_nonSelected_leaderless_byType_Others_fscore_avg.txt')

fscore_avg_10repeats_leadered_Others <- fscore_avg_HalfLife_leadered_Others[, c(6:15)]
delta_leadered_Others_p1 <- fscore_avg_10repeats_leadered_Others[1, ] - fscore_avg_10repeats_leadered_Others[2, ]
delta_leadered_Others <- delta_leadered_Others_p1 %>% 
  gather(repeats, fscore_diff) %>% 
  mutate(type = "leadered")

fscore_avg_10repeats_leaderless_Others <- fscore_avg_HalfLife_leaderless_Others[, c(6:15)]
delta_leaderless_Others_p1 <- fscore_avg_10repeats_leaderless_Others[1, ] - fscore_avg_10repeats_leaderless_Others[2, ]
delta_leaderless_Others <- delta_leaderless_Others_p1 %>% 
  gather(repeats, fscore_diff) %>% 
  mutate(type = "leaderless")

wilcox.test(delta_leadered_Others$fscore_diff, delta_leaderless_Others$fscore_diff, paired = FALSE)

delta_Others <- bind_rows(delta_leadered_Others, delta_leaderless_Others)

delta_Others %>% 
  mutate(type = factor(type, level = c("leadered", "leaderless"))) %>% 
  group_by(type) %>% 
  summarise(
    sd = sd(fscore_diff), 
    delta = mean(fscore_diff)
  ) %>% 
  ggplot(aes(type, delta)) + 
  scale_colour_manual(values = transcript_type) +
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
ggsave("./viz_deltaFscores_toCompareTranscript/deltaFscores_avgs_toCompareTranscriptType_Others_hypoxia.png", width = 8, height = 20, units = "cm", dpi = 600)

###### combined selected
fscore_avg_HalfLife_leadered_CombinedSelected <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_hypoxia_combinedSelected_leadered_byLeaderlessHypoxia_jointlyByLeaderedLeaderless_fscore_avg.txt')
fscore_avg_HalfLife_leaderless_CombinedSelected <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_hypoxia_combinedSelected_leaderless_byLeaderlessHypoxia_jointlyByLeaderedLeaderless_fscore_avg.txt')

fscore_avg_10repeats_leadered_CombinedSelected <- fscore_avg_HalfLife_leadered_CombinedSelected[, c(6:15)]
delta_leadered_CombinedSelected_p1 <- fscore_avg_10repeats_leadered_CombinedSelected[1, ] - fscore_avg_10repeats_leadered_CombinedSelected[2, ]
delta_leadered_CombinedSelected <- delta_leadered_CombinedSelected_p1 %>% 
  gather(repeats, fscore_diff) %>% 
  mutate(type = "leadered")

fscore_avg_10repeats_leaderless_CombinedSelected <- fscore_avg_HalfLife_leaderless_CombinedSelected[, c(6:15)]
delta_leaderless_CombinedSelected_p1 <- fscore_avg_10repeats_leaderless_CombinedSelected[1, ] - fscore_avg_10repeats_leaderless_CombinedSelected[2, ]
delta_leaderless_CombinedSelected <- delta_leaderless_CombinedSelected_p1 %>% 
  gather(repeats, fscore_diff) %>% 
  mutate(type = "leaderless")

wilcox.test(delta_leadered_CombinedSelected$fscore_diff, delta_leaderless_CombinedSelected$fscore_diff, paired = FALSE)

delta_CombinedSelected <- bind_rows(delta_leadered_CombinedSelected, delta_leaderless_CombinedSelected)

delta_CombinedSelected %>% 
  mutate(type = factor(type, level = c("leadered", "leaderless"))) %>% 
  group_by(type) %>% 
  summarise(
    sd = sd(fscore_diff), 
    delta = mean(fscore_diff)
  ) %>% 
  ggplot(aes(type, delta)) + 
  scale_colour_manual(values = transcript_type) +
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
ggsave("./viz_deltaFscores_toCompareTranscript/deltaFscores_avgs_toCompareTranscriptType_Combined_CombinedSelectedByLeaderless_hypoxia.png", width = 8, height = 20, units = "cm", dpi = 600)
