
####### To compare fscore_avgs of models using only 5'UTR features
####### 5'UTR features are also split by translation related and UTR (not translation) related
####### (final version remove the viz for completed feature models) for log phase, hypoxia, half-life fold change in hypoxia and their fscore_avgs with combined features

library(tidyverse)
library(RColorBrewer)
library(scales)

# display.brewer.all(type = 'qual')
# display.brewer.pal(n = 12, name = 'Set3')
# brewer.pal(n = 12, name = 'Set3')

####### get fscore_avgs
###### complete 5'UTR features
fscore_avg_HalfLife_logPhase_5pUTR_complete <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_logPhase_nonSelected_leadered_byType_5pUTR_complete_fscore_avg.txt')
fscore_avg_HalfLife_hypoxia_5pUTR_complete <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_hypoxia_nonSelected_leadered_byType_5pUTR_complete_fscore_avg.txt')
fscore_avg_HalfLifeFc_hypoxia2LogPhase_5pUTR_complete <- read.table('../ml_metrics_halfLifeFcClass/HalfLifeFcCls_hypoxiaToLogPhase_nonSelected_leadered_byType_5pUTR_complete_fscore_avg.txt')

###### translation related 5'UTR features
fscore_avg_HalfLife_logPhase_5pUTR_TranslationRelated <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_logPhase_nonSelected_leadered_byType_5pUTR_TranslationRelated_fscore_avg.txt')
fscore_avg_HalfLife_hypoxia_5pUTR_TranslationRelated <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_hypoxia_nonSelected_leadered_byType_5pUTR_TranslationRelated_fscore_avg.txt')
fscore_avg_HalfLifeFc_hypoxia2LogPhase_5pUTR_TranslationRelated <- read.table('../ml_metrics_halfLifeFcClass/HalfLifeFcCls_hypoxiaToLogPhase_nonSelected_leadered_byType_5pUTR_TranslationRelated_fscore_avg.txt')

###### UTR (not translation) related 5'UTR features
fscore_avg_HalfLife_logPhase_5pUTR_UTRrelated <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_logPhase_nonSelected_leadered_byType_5pUTR_UTRrelated_fscore_avg.txt')
fscore_avg_HalfLife_hypoxia_5pUTR_UTRrelated <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_hypoxia_nonSelected_leadered_byType_5pUTR_UTRrelated_fscore_avg.txt')
fscore_avg_HalfLifeFc_hypoxia2LogPhase_5pUTR_UTRrelated <- read.table('../ml_metrics_halfLifeFcClass/HalfLifeFcCls_hypoxiaToLogPhase_nonSelected_leadered_byType_5pUTR_UTRrelated_fscore_avg.txt')

###### using combined features
##### final version remove this part of viz
# fscore_avg_HalfLife_logPhase_CombinedSelected_leadered <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_logPhase_combinedSelected_leadered_byLeaderedLogPhase_jointlyByLogPhaseHypoxia_fscore_avg.txt')
# fscore_avg_HalfLife_hypoxia_CombinedSelected_leadered <- read.table('../ml_metrics_halfLifeClass/HalfLifeCls_hypoxia_combinedSelected_leadered_byLeaderedLogPhase_jointlyByLogPhaseHypoxia_fscore_avg.txt')
# fscore_avg_HalfLifeFc_hypoxia2LogPhase_CombinedSelected_leadered <- read.table('../ml_metrics_halfLifeFcClass/HalfLifeFcCls_hypoxiaToLogPhase_combinedSelected_leadered_byLeaderedLogPhase_jointlyByLogPhaseFoldChange_fscore_avg.txt')

####### get fscore_avgs relative to random prediction
###### complete 5'UTR features
fscore_avg_10repeats_logPhase_5pUTR_complete <- fscore_avg_HalfLife_logPhase_5pUTR_complete[, c(6:15)]
delta_logPhase_5pUTR_complete <- fscore_avg_10repeats_logPhase_5pUTR_complete[1, ] - fscore_avg_10repeats_logPhase_5pUTR_complete[2, ]

fscore_avg_10repeats_hypoxia_5pUTR_complete <- fscore_avg_HalfLife_hypoxia_5pUTR_complete[, c(6:15)]
delta_hypoxia_5pUTR_complete <- fscore_avg_10repeats_hypoxia_5pUTR_complete[1, ] - fscore_avg_10repeats_hypoxia_5pUTR_complete[2, ]

fscore_avg_10repeats_hypoxia2LogPhase_5pUTR_complete <- fscore_avg_HalfLifeFc_hypoxia2LogPhase_5pUTR_complete[, c(6:15)]
delta_hypoxia2LogPhase_5pUTR_complete <- fscore_avg_10repeats_hypoxia2LogPhase_5pUTR_complete[1, ] - fscore_avg_10repeats_hypoxia2LogPhase_5pUTR_complete[2, ]

###### translation related 5'UTR features
fscore_avg_10repeats_logPhase_5pUTR_TranslationRelated <- fscore_avg_HalfLife_logPhase_5pUTR_TranslationRelated[, c(6:15)]
delta_logPhase_5pUTR_TranslationRelated <- fscore_avg_10repeats_logPhase_5pUTR_TranslationRelated[1, ] - fscore_avg_10repeats_logPhase_5pUTR_TranslationRelated[2, ]

delta_logPhase_5pUTR_TranslationRelated_forTest <- delta_logPhase_5pUTR_TranslationRelated %>% 
  gather(repeats, fscore_diff) %>% 
  mutate(type = "logPhase")

fscore_avg_10repeats_hypoxia_5pUTR_TranslationRelated <- fscore_avg_HalfLife_hypoxia_5pUTR_TranslationRelated[, c(6:15)]
delta_hypoxia_5pUTR_TranslationRelated <- fscore_avg_10repeats_hypoxia_5pUTR_TranslationRelated[1, ] - fscore_avg_10repeats_hypoxia_5pUTR_TranslationRelated[2, ]

delta_hypoxia_5pUTR_TranslationRelated_forTest <- delta_hypoxia_5pUTR_TranslationRelated %>% 
  gather(repeats, fscore_diff) %>% 
  mutate(type = "hypoxia")

wilcox.test(delta_logPhase_5pUTR_TranslationRelated_forTest$fscore_diff, delta_hypoxia_5pUTR_TranslationRelated_forTest$fscore_diff, paired = FALSE)

fscore_avg_10repeats_hypoxia2LogPhase_5pUTR_TranslationRelated <- fscore_avg_HalfLifeFc_hypoxia2LogPhase_5pUTR_TranslationRelated[, c(6:15)]
delta_hypoxia2LogPhase_5pUTR_TranslationRelated <- fscore_avg_10repeats_hypoxia2LogPhase_5pUTR_TranslationRelated[1, ] - fscore_avg_10repeats_hypoxia2LogPhase_5pUTR_TranslationRelated[2, ]

###### UTR (not translation) related 5'UTR features
fscore_avg_10repeats_logPhase_5pUTR_UTRrelated <- fscore_avg_HalfLife_logPhase_5pUTR_UTRrelated[, c(6:15)]
delta_logPhase_5pUTR_UTRrelated <- fscore_avg_10repeats_logPhase_5pUTR_UTRrelated[1, ] - fscore_avg_10repeats_logPhase_5pUTR_UTRrelated[2, ]

fscore_avg_10repeats_hypoxia_5pUTR_UTRrelated <- fscore_avg_HalfLife_hypoxia_5pUTR_UTRrelated[, c(6:15)]
delta_hypoxia_5pUTR_UTRrelated <- fscore_avg_10repeats_hypoxia_5pUTR_UTRrelated[1, ] - fscore_avg_10repeats_hypoxia_5pUTR_UTRrelated[2, ]

fscore_avg_10repeats_hypoxia2LogPhase_5pUTR_UTRrelated <- fscore_avg_HalfLifeFc_hypoxia2LogPhase_5pUTR_UTRrelated[, c(6:15)]
delta_hypoxia2LogPhase_5pUTR_UTRrelated <- fscore_avg_10repeats_hypoxia2LogPhase_5pUTR_UTRrelated[1, ] - fscore_avg_10repeats_hypoxia2LogPhase_5pUTR_UTRrelated[2, ]

###### using combined features
##### previous version included models with combined features, final version removed it
# fscore_avg_10repeats_logPhase_CombinedSelected_leadered <- fscore_avg_HalfLife_logPhase_CombinedSelected_leadered[, c(6:15)]
# delta_logPhase_CombinedSelected_leadered <- fscore_avg_10repeats_logPhase_CombinedSelected_leadered[1, ] - fscore_avg_10repeats_logPhase_CombinedSelected_leadered[2, ]
# 
# fscore_avg_10repeats_hypoxia_CombinedSelected_leadered <- fscore_avg_HalfLife_hypoxia_CombinedSelected_leadered[, c(6:15)]
# delta_hypoxia_CombinedSelected_leadered <- fscore_avg_10repeats_hypoxia_CombinedSelected_leadered[1, ] - fscore_avg_10repeats_hypoxia_CombinedSelected_leadered[2, ]
# 
# fscore_avg_10repeats_hypoxia2LogPhase_CombinedSelected_leadered <- fscore_avg_HalfLifeFc_hypoxia2LogPhase_CombinedSelected_leadered[, c(6:15)]
# delta_hypoxia2LogPhase_CombinedSelected_leadered <- fscore_avg_10repeats_hypoxia2LogPhase_CombinedSelected_leadered[1, ] - fscore_avg_10repeats_hypoxia2LogPhase_CombinedSelected_leadered[2, ]

###### combine delta fscore_avg
##### previous version included models with combined features, final version removed it
# fscore_avg_delta <- rbind(delta_logPhase_5pUTR_complete, delta_hypoxia_5pUTR_complete, delta_hypoxia2LogPhase_5pUTR_complete,
#                       delta_logPhase_5pUTR_TranslationRelated, delta_hypoxia_5pUTR_TranslationRelated, delta_hypoxia2LogPhase_5pUTR_TranslationRelated,
#                       delta_logPhase_5pUTR_UTRrelated, delta_hypoxia_5pUTR_UTRrelated, delta_hypoxia2LogPhase_5pUTR_UTRrelated,
#                       delta_logPhase_CombinedSelected_leadered, delta_hypoxia_CombinedSelected_leadered, 
#                       delta_hypoxia2LogPhase_CombinedSelected_leadered)
# fscore_avg_delta$model <- c('logPhase_5pUTR_complete', 'hypoxia_5pUTR_complete', 'hypoxia2LogPhase_5pUTR_complete', 
#                         'logPhase_5pUTR_TranslationRelated', 'hypoxia_5pUTR_TranslationRelated', 'hypoxia2LogPhase_5pUTR_TranslationRelated',
#                         'logPhase_5pUTR_UTRrelated', 'hypoxia_5pUTR_UTRrelated', 'hypoxia2LogPhase_5pUTR_UTRrelated',
#                         'logPhase_CombinedSelected', 'hypoxia_CombinedSelected', 'hypoxia2LogPhase_CombinedSelected')
# fscore_avg_delta$model <- factor(fscore_avg_delta$model, levels = c('hypoxia2LogPhase_5pUTR_UTRrelated', 'hypoxia_5pUTR_UTRrelated', 'logPhase_5pUTR_UTRrelated', 
#                                                             'hypoxia2LogPhase_5pUTR_TranslationRelated', 'hypoxia_5pUTR_TranslationRelated', 'logPhase_5pUTR_TranslationRelated',
#                                                             'hypoxia2LogPhase_5pUTR_complete', 'hypoxia_5pUTR_complete', 'logPhase_5pUTR_complete',
#                                                             'hypoxia2LogPhase_CombinedSelected', 'hypoxia_CombinedSelected', 'logPhase_CombinedSelected'))

fscore_avg_delta <- rbind(delta_logPhase_5pUTR_complete, delta_hypoxia_5pUTR_complete, delta_hypoxia2LogPhase_5pUTR_complete,
                          delta_logPhase_5pUTR_TranslationRelated, delta_hypoxia_5pUTR_TranslationRelated, delta_hypoxia2LogPhase_5pUTR_TranslationRelated,
                          delta_logPhase_5pUTR_UTRrelated, delta_hypoxia_5pUTR_UTRrelated, delta_hypoxia2LogPhase_5pUTR_UTRrelated)
fscore_avg_delta$model <- c('logPhase_5pUTR_complete', 'hypoxia_5pUTR_complete', 'hypoxia2LogPhase_5pUTR_complete', 
                            'logPhase_5pUTR_TranslationRelated', 'hypoxia_5pUTR_TranslationRelated', 'hypoxia2LogPhase_5pUTR_TranslationRelated',
                            'logPhase_5pUTR_UTRrelated', 'hypoxia_5pUTR_UTRrelated', 'hypoxia2LogPhase_5pUTR_UTRrelated')
fscore_avg_delta$model <- factor(fscore_avg_delta$model, levels = c('hypoxia2LogPhase_5pUTR_UTRrelated', 'hypoxia_5pUTR_UTRrelated', 'logPhase_5pUTR_UTRrelated', 
                                                                    'hypoxia2LogPhase_5pUTR_TranslationRelated', 'hypoxia_5pUTR_TranslationRelated', 'logPhase_5pUTR_TranslationRelated',
                                                                    'hypoxia2LogPhase_5pUTR_complete', 'hypoxia_5pUTR_complete', 'logPhase_5pUTR_complete'))

###### viz of delta fscore_avg
fscore_avg_delta %>% 
  gather(repeat_10, delta, -model) %>%
  group_by(model) %>%
  summarise(
    sd = sd(delta), 
    delta = mean(delta)
  ) %>% 
  ggplot(aes(model, delta)) + 
  geom_pointrange(aes(ymin = delta - sd, ymax = delta + sd, color = model), 
                  size = 0.8) +
  scale_color_manual(values = c("#0D95D0", "#83B7DF", "#036D9C", "#0D95D0", "#83B7DF", "#036D9C",
                                "#0D95D0", "#83B7DF", "#036D9C", "#0D95D0", "#83B7DF", "#036D9C")) +
  ylim(-0.05, 0.15) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
    axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 5, face = "bold"),
    axis.title = element_blank(),
    legend.key.size = unit(1, 'cm')
  ) +
  coord_flip()
ggsave("./viz_deltaFscores_toCompare5pUTR/deltaFscores_avgs_toCompare5pUTR.png", width = 20, height = 10, units = "cm", dpi = 600)
