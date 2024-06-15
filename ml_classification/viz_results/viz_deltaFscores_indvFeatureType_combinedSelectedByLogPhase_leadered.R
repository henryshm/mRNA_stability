
####### To compare training with individual feature type with combined features
####### comparing log phase, hypoxia, half-life fold change in hypoxia
####### feature selected by log phase in combined feature model
####### for leadered genes

library(tidyverse)
library(RColorBrewer)
library(scales)

# display.brewer.all(type = 'qual')
# display.brewer.pal(n = 12, name = 'Paired')
# brewer.pal(n = 12, name = 'Paired')

####### get fscore_avgs
###### 5'UTR
fscore_avg_HalfLife_logPhase_5pUTR_complete <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_logPhase_nonSelected_leadered_byType_5pUTR_complete_fscore_avg.txt')
fscore_avg_HalfLife_hypoxia_5pUTR_complete <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_hypoxia_nonSelected_leadered_byType_5pUTR_complete_fscore_avg.txt')
fscore_avg_HalfLife_hypoxia2LogPhase_5pUTR_complete <- read.table('../ml_metrics_halfLifeFcClass/HalfLifeFcCls_hypoxiaToLogPhase_nonSelected_leadered_byType_5pUTR_complete_fscore_avg.txt')

fscore_avg_10repeats_logPhase_5pUTR_complete <- fscore_avg_HalfLife_logPhase_5pUTR_complete[, c(6:15)]
delta_logPhase_5pUTR_complete <- fscore_avg_10repeats_logPhase_5pUTR_complete[1, ] - fscore_avg_10repeats_logPhase_5pUTR_complete[2, ]

fscore_avg_10repeats_hypoxia_5pUTR_complete <- fscore_avg_HalfLife_hypoxia_5pUTR_complete[, c(6:15)]
delta_hypoxia_5pUTR_complete <- fscore_avg_10repeats_hypoxia_5pUTR_complete[1, ] - fscore_avg_10repeats_hypoxia_5pUTR_complete[2, ]

fscore_avg_10repeats_hypoxia2LogPhase_5pUTR_complete <- fscore_avg_HalfLife_hypoxia2LogPhase_5pUTR_complete[, c(6:15)]
delta_hypoxia2LogPhase_5pUTR_complete <- fscore_avg_10repeats_hypoxia2LogPhase_5pUTR_complete[1, ] - fscore_avg_10repeats_hypoxia2LogPhase_5pUTR_complete[2, ]

###### CDS nucleotide
fscore_avg_HalfLife_logPhase_CDSnucleotide <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_logPhase_nonSelected_leadered_byType_CDSnucleotide_fscore_avg.txt')
fscore_avg_HalfLife_hypoxia_CDSnucleotide <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_hypoxia_nonSelected_leadered_byType_CDSnucleotide_fscore_avg.txt')
fscore_avg_HalfLife_hypoxia2LogPhase_CDSnucleotide <- read.table('../ml_metrics_halfLifeFcClass/HalfLifeFcCls_hypoxiaToLogPhase_nonSelected_leadered_byType_CDSnucleotide_fscore_avg.txt')

fscore_avg_10repeats_logPhase_CDSnucleotide <- fscore_avg_HalfLife_logPhase_CDSnucleotide[, c(6:15)]
delta_logPhase_CDSnucleotide <- fscore_avg_10repeats_logPhase_CDSnucleotide[1, ] - fscore_avg_10repeats_logPhase_CDSnucleotide[2, ]

fscore_avg_10repeats_hypoxia_CDSnucleotide <- fscore_avg_HalfLife_hypoxia_CDSnucleotide[, c(6:15)]
delta_hypoxia_CDSnucleotide <- fscore_avg_10repeats_hypoxia_CDSnucleotide[1, ] - fscore_avg_10repeats_hypoxia_CDSnucleotide[2, ]

fscore_avg_10repeats_hypoxia2LogPhase_CDSnucleotide <- fscore_avg_HalfLife_hypoxia2LogPhase_CDSnucleotide[, c(6:15)]
delta_hypoxia2LogPhase_CDSnucleotide <- fscore_avg_10repeats_hypoxia2LogPhase_CDSnucleotide[1, ] - fscore_avg_10repeats_hypoxia2LogPhase_CDSnucleotide[2, ]

###### CDS MFE
fscore_avg_HalfLife_logPhase_CDSmfe <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_logPhase_nonSelected_leadered_byType_CDSmfe_fscore_avg.txt')
fscore_avg_HalfLife_hypoxia_CDSmfe <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_hypoxia_nonSelected_leadered_byType_CDSmfe_fscore_avg.txt')
fscore_avg_HalfLife_hypoxia2LogPhase_CDSmfe <- read.table('../ml_metrics_halfLifeFcClass/HalfLifeFcCls_hypoxiaToLogPhase_nonSelected_leadered_byType_CDSmfe_fscore_avg.txt')

fscore_avg_10repeats_logPhase_CDSmfe <- fscore_avg_HalfLife_logPhase_CDSmfe[, c(6:15)]
delta_logPhase_CDSmfe <- fscore_avg_10repeats_logPhase_CDSmfe[1, ] - fscore_avg_10repeats_logPhase_CDSmfe[2, ]

fscore_avg_10repeats_hypoxia_CDSmfe <- fscore_avg_HalfLife_hypoxia_CDSmfe[, c(6:15)]
delta_hypoxia_CDSmfe <- fscore_avg_10repeats_hypoxia_CDSmfe[1, ] - fscore_avg_10repeats_hypoxia_CDSmfe[2, ]

fscore_avg_10repeats_hypoxia2LogPhase_CDSmfe <- fscore_avg_HalfLife_hypoxia2LogPhase_CDSmfe[, c(6:15)]
delta_hypoxia2LogPhase_CDSmfe <- fscore_avg_10repeats_hypoxia2LogPhase_CDSmfe[1, ] - fscore_avg_10repeats_hypoxia2LogPhase_CDSmfe[2, ]

###### Codon
fscore_avg_HalfLife_logPhase_Codon <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_logPhase_nonSelected_leadered_byType_Codon_fscore_avg.txt')
fscore_avg_HalfLife_hypoxia_Codon <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_hypoxia_nonSelected_leadered_byType_Codon_fscore_avg.txt')
fscore_avg_HalfLife_hypoxia2LogPhase_Codon <- read.table('../ml_metrics_halfLifeFcClass/HalfLifeFcCls_hypoxiaToLogPhase_nonSelected_leadered_byType_Codon_fscore_avg.txt')

fscore_avg_10repeats_logPhase_Codon <- fscore_avg_HalfLife_logPhase_Codon[, c(6:15)]
delta_logPhase_Codon <- fscore_avg_10repeats_logPhase_Codon[1, ] - fscore_avg_10repeats_logPhase_Codon[2, ]

fscore_avg_10repeats_hypoxia_Codon <- fscore_avg_HalfLife_hypoxia_Codon[, c(6:15)]
delta_hypoxia_Codon <- fscore_avg_10repeats_hypoxia_Codon[1, ] - fscore_avg_10repeats_hypoxia_Codon[2, ]

fscore_avg_10repeats_hypoxia2LogPhase_Codon <- fscore_avg_HalfLife_hypoxia2LogPhase_Codon[, c(6:15)]
delta_hypoxia2LogPhase_Codon <- fscore_avg_10repeats_hypoxia2LogPhase_Codon[1, ] - fscore_avg_10repeats_hypoxia2LogPhase_Codon[2, ]

###### Translation
fscore_avg_HalfLife_logPhase_Translation <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_logPhase_nonSelected_leadered_byType_Translation_fscore_avg.txt')
fscore_avg_HalfLife_hypoxia_Translation <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_hypoxia_nonSelected_leadered_byType_Translation_fscore_avg.txt')
fscore_avg_HalfLife_hypoxia2LogPhase_Translation <- read.table('../ml_metrics_halfLifeFcClass/HalfLifeFcCls_hypoxiaToLogPhase_nonSelected_leadered_byType_Translation_fscore_avg.txt')

fscore_avg_10repeats_logPhase_Translation <- fscore_avg_HalfLife_logPhase_Translation[, c(6:15)]
delta_logPhase_Translation <- fscore_avg_10repeats_logPhase_Translation[1, ] - fscore_avg_10repeats_logPhase_Translation[2, ]

fscore_avg_10repeats_hypoxia_Translation <- fscore_avg_HalfLife_hypoxia_Translation[, c(6:15)]
delta_hypoxia_Translation <- fscore_avg_10repeats_hypoxia_Translation[1, ] - fscore_avg_10repeats_hypoxia_Translation[2, ]

fscore_avg_10repeats_hypoxia2LogPhase_Translation <- fscore_avg_HalfLife_hypoxia2LogPhase_Translation[, c(6:15)]
delta_hypoxia2LogPhase_Translation <- fscore_avg_10repeats_hypoxia2LogPhase_Translation[1, ] - fscore_avg_10repeats_hypoxia2LogPhase_Translation[2, ]

###### Others
fscore_avg_HalfLife_logPhase_Others <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_logPhase_nonSelected_leadered_byType_Others_fscore_avg.txt')
fscore_avg_HalfLife_hypoxia_Others <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_hypoxia_nonSelected_leadered_byType_Others_fscore_avg.txt')
fscore_avg_HalfLife_hypoxia2LogPhase_Others <- read.table('../ml_metrics_halfLifeFcClass/HalfLifeFcCls_hypoxiaToLogPhase_nonSelected_leadered_byType_Others_fscore_avg.txt')

fscore_avg_10repeats_logPhase_Others <- fscore_avg_HalfLife_logPhase_Others[, c(6:15)]
delta_logPhase_Others <- fscore_avg_10repeats_logPhase_Others[1, ] - fscore_avg_10repeats_logPhase_Others[2, ]

fscore_avg_10repeats_hypoxia_Others <- fscore_avg_HalfLife_hypoxia_Others[, c(6:15)]
delta_hypoxia_Others <- fscore_avg_10repeats_hypoxia_Others[1, ] - fscore_avg_10repeats_hypoxia_Others[2, ]

fscore_avg_10repeats_hypoxia2LogPhase_Others <- fscore_avg_HalfLife_hypoxia2LogPhase_Others[, c(6:15)]
delta_hypoxia2LogPhase_Others <- fscore_avg_10repeats_hypoxia2LogPhase_Others[1, ] - fscore_avg_10repeats_hypoxia2LogPhase_Others[2, ]

##### combined selected features
fscore_avg_HalfLife_logPhase_CombinedSelected_leadered <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_logPhase_combinedSelected_leadered_byLeaderedLogPhase_jointlyByLogPhaseHypoxia_fscore_avg.txt')
fscore_avg_HalfLife_hypoxia_CombinedSelected_leadered <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_hypoxia_combinedSelected_leadered_byLeaderedLogPhase_jointlyByLogPhaseHypoxia_fscore_avg.txt')
fscore_avg_HalfLife_hypoxia2LogPhase_CombinedSelected_leadered <- read.table('../ml_metrics_halfLifeFcClass/HalfLifeFcCls_hypoxiaToLogPhase_combinedSelected_leadered_byLeaderedLogPhase_jointlyByLogPhaseFoldChange_fscore_avg.txt')

fscore_avg_10repeats_logPhase_CombinedSelected_leadered <- fscore_avg_HalfLife_logPhase_CombinedSelected_leadered[, c(6:15)]
delta_logPhase_CombinedSelected_leadered <- fscore_avg_10repeats_logPhase_CombinedSelected_leadered[1, ] - fscore_avg_10repeats_logPhase_CombinedSelected_leadered[2, ]

fscore_avg_10repeats_hypoxia_CombinedSelected_leadered <- fscore_avg_HalfLife_hypoxia_CombinedSelected_leadered[, c(6:15)]
delta_hypoxia_CombinedSelected_leadered <- fscore_avg_10repeats_hypoxia_CombinedSelected_leadered[1, ] - fscore_avg_10repeats_hypoxia_CombinedSelected_leadered[2, ]

fscore_avg_10repeats_hypoxia2LogPhase_CombinedSelected_leadered <- fscore_avg_HalfLife_hypoxia2LogPhase_CombinedSelected_leadered[, c(6:15)]
delta_hypoxia2LogPhase_CombinedSelected_leadered <- fscore_avg_10repeats_hypoxia2LogPhase_CombinedSelected_leadered[1, ] - fscore_avg_10repeats_hypoxia2LogPhase_CombinedSelected_leadered[2, ]


####### delta fscore_avg of each individual feature model
###### log phase
fscore_avg_delta_logPhase <- rbind(delta_logPhase_5pUTR_complete, 
                               delta_logPhase_CDSnucleotide, 
                               delta_logPhase_CDSmfe, 
                               delta_logPhase_Codon, 
                               delta_logPhase_Translation, 
                               delta_logPhase_Others,
                               delta_logPhase_CombinedSelected_leadered)

###### adding feature type labels
fscore_avg_delta_logPhase$FeatureType <- c('5pUTR', 'CDSnucleotide', 'CDSmfe', 'Codon', 'Translation', 'Others', 'CombinedSelected')
fscore_avg_delta_logPhase$FeatureType <- factor(fscore_avg_delta_logPhase$FeatureType,
                                            level = c('CombinedSelected','Others', 'Translation',
                                                      'Codon', 'CDSmfe', 'CDSnucleotide', '5pUTR'))

###### viz of delta fscore_avg
fscore_avg_delta_logPhase %>% 
  gather(repeat_10, delta, -FeatureType) %>%
  group_by(FeatureType) %>% 
  summarise(
    sd = sd(delta), 
    delta = mean(delta)
  ) %>% 
  ggplot(aes(FeatureType, delta)) + 
  geom_pointrange(aes(ymin = delta - sd, ymax = delta + sd), 
                  color = "#036D9C",
                  size = 0.6) +
  ylim(-0.05, 0.15) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
    axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 3, face = "bold"),
    axis.title = element_blank()
    ) +
  coord_flip()
ggsave("./viz_deltaFscores_toCompareAllmodels/deltaFscores_avgs_toCompareAllmodels_CombinedSelectedByLogPhase_logPhase_leadered.png", width = 7, height = 8, units = "cm", dpi = 600)

###### hypoxia
fscore_avg_delta_hypoxia <- rbind(delta_hypoxia_5pUTR_complete, 
                              delta_hypoxia_CDSnucleotide, 
                              delta_hypoxia_CDSmfe,
                              delta_hypoxia_Codon,
                              delta_hypoxia_Translation,
                              delta_hypoxia_Others,
                              delta_hypoxia_CombinedSelected_leadered)

###### adding feature type labels
fscore_avg_delta_hypoxia$FeatureType <- c('5pUTR', 'CDSnucleotide', 'CDSmfe', 'Codon', 'Translation', 'Others', 'CombinedSelected')
fscore_avg_delta_hypoxia$FeatureType <- factor(fscore_avg_delta_hypoxia$FeatureType,
                                           level = c('CombinedSelected','Others', 'Translation',
                                                     'Codon', 'CDSmfe', 'CDSnucleotide', '5pUTR'))

###### viz of delta fscore_avg
fscore_avg_delta_hypoxia %>% 
  gather(repeat_10, delta, -FeatureType) %>%
  group_by(FeatureType) %>% 
  summarise(
    sd = sd(delta), 
    delta = mean(delta)
  ) %>% 
  ggplot(aes(FeatureType, delta)) + 
  geom_pointrange(aes(ymin = delta - sd, ymax = delta + sd), 
                  color = "#83B7DF",
                  size = 0.6) +
  ylim(-0.05, 0.15) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
    axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 3, face = "bold"),
    axis.title = element_blank()
    ) +
  coord_flip()
ggsave("./viz_deltaFscores_toCompareAllmodels/deltaFscores_avgs_toCompareAllmodels_CombinedSelectedByLogPhase_hypoxia_leadered.png", width = 7, height = 8, units = "cm", dpi = 600)

###### fold change in hypoxia
fscore_avg_delta_hypoxia2LogPhase <- rbind(delta_hypoxia2LogPhase_5pUTR_complete,
                                       delta_hypoxia2LogPhase_CDSnucleotide,
                                       delta_hypoxia2LogPhase_CDSmfe,
                                       delta_hypoxia2LogPhase_Codon,
                                       delta_hypoxia2LogPhase_Translation,
                                       delta_hypoxia2LogPhase_Others,
                                       delta_hypoxia2LogPhase_CombinedSelected_leadered)

###### adding feature type labels
fscore_avg_delta_hypoxia2LogPhase$FeatureType <- c('5pUTR', 'CDSnucleotide', 'CDSmfe', 'Codon', 'Translation', 'Others', 'CombinedSelected')
fscore_avg_delta_hypoxia2LogPhase$FeatureType <- factor(fscore_avg_delta_hypoxia2LogPhase$FeatureType,
                                                    level = c('CombinedSelected','Others', 'Translation',
                                                              'Codon', 'CDSmfe', 'CDSnucleotide', '5pUTR'))

###### viz of delta fscore_avg
fscore_avg_delta_hypoxia2LogPhase %>% 
  gather(repeat_10, delta, -FeatureType) %>%
  group_by(FeatureType) %>% 
  summarise(
    sd = sd(delta), 
    delta = mean(delta)
  ) %>% 
  ggplot(aes(FeatureType, delta)) + 
  geom_pointrange(aes(ymin = delta - sd, ymax = delta + sd), 
                  color = "#0D95D0",
                  size = 0.6) +
  ylim(-0.05, 0.15) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
    axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 3, face = "bold"),
    axis.title = element_blank()
    ) +
  coord_flip()
ggsave("./viz_deltaFscores_toCompareAllmodels/deltaFscores_avgs_toCompareAllmodels_CombinedSelectedByLogPhase_hypoxiaToLogPhase_leadered.png", width = 7, height = 8, units = "cm", dpi = 600)
