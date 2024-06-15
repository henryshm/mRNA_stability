
####### viz of feature distributions
####### To compare features of leadered and leaderless transcripts 
####### in log phase, hypoxia and fold change in hypoxia
####### features are the top20 important features for the two transcript types
####### feature selected by leaderless in combined selected model
####### log phase
####### leadered transcript

library(tidyverse)
library(RColorBrewer)
library(PupillometryR)
library(scales)

####### get color scheme for feature distributions
col_hl_cls <- c("#7DC462", "#0D95D0", "#774FA0", "#E72F52")

####### get feature table
featureTable <- read.csv('../../feature/FeatureTables/featureTable_combinedSelected_halfLifeClass/HalfLifeCls_logPhase_combinedSelected_leadered_byLeaderlessLogPhase_jointlyByLeaderedLeaderless.csv')

####### viz of feature distributions
###### adja_CG_CDS
feature_i = "adja_CG_CDS"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  group_by(HalfLife_cls_logPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_cls_logPhase)) +
  scale_colour_manual(values = col_hl_cls) +
  scale_fill_manual(values = col_hl_cls) +
  stat_summary(mapping = aes(HalfLife_cls_logPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_cls_logPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_cls_logPhase, feature_toViz, fill = HalfLife_cls_logPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
  ylim(0.09, 0.2) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareTranscript_logPhase/featureDistributionTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_adja_CG_CDS.png", width = 35, height = 18, units = "cm", dpi = 600)

###### adja_GU_CDS
feature_i = "adja_GU_CDS"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  group_by(HalfLife_cls_logPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_cls_logPhase)) +
  scale_colour_manual(values = col_hl_cls) +
  scale_fill_manual(values = col_hl_cls) +
  stat_summary(mapping = aes(HalfLife_cls_logPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_cls_logPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_cls_logPhase, feature_toViz, fill = HalfLife_cls_logPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
  ylim(0.02, 0.1) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareTranscript_logPhase/featureDistributionTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_adja_GU_CDS.png", width = 35, height = 18, units = "cm", dpi = 600)

###### adja_UA_CDS
feature_i = "adja_UA_CDS"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  group_by(HalfLife_cls_logPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_cls_logPhase)) +
  scale_colour_manual(values = col_hl_cls) +
  scale_fill_manual(values = col_hl_cls) +
  stat_summary(mapping = aes(HalfLife_cls_logPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_cls_logPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_cls_logPhase, feature_toViz, fill = HalfLife_cls_logPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
  ylim(-0.01, 0.025) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareTranscript_logPhase/featureDistributionTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_adja_UA_CDS.png", width = 35, height = 18, units = "cm", dpi = 600)

###### CDS_MFE_100_50nt_3p
feature_i = "CDS_MFE_100_50nt_3p"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  group_by(HalfLife_cls_logPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_cls_logPhase)) +
  scale_colour_manual(values = col_hl_cls) +
  scale_fill_manual(values = col_hl_cls) +
  stat_summary(mapping = aes(HalfLife_cls_logPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_cls_logPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_cls_logPhase, feature_toViz, fill = HalfLife_cls_logPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
  ylim(-55, -10) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareTranscript_logPhase/featureDistributionTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_CDS_MFE_100_50nt_3p.png", width = 35, height = 18, units = "cm", dpi = 600)

###### CDS_MFE_100_50nt_5p
feature_i = "CDS_MFE_100_50nt_5p"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  group_by(HalfLife_cls_logPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_cls_logPhase)) +
  scale_colour_manual(values = col_hl_cls) +
  scale_fill_manual(values = col_hl_cls) +
  stat_summary(mapping = aes(HalfLife_cls_logPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_cls_logPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_cls_logPhase, feature_toViz, fill = HalfLife_cls_logPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
  ylim(-55, -20) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareTranscript_logPhase/featureDistributionTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_CDS_MFE_100_50nt_5p.png", width = 35, height = 18, units = "cm", dpi = 600)

###### CDS_MFE_20_10nt_mid
feature_i = "CDS_MFE_20_10nt_mid"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  group_by(HalfLife_cls_logPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_cls_logPhase)) +
  scale_colour_manual(values = col_hl_cls) +
  scale_fill_manual(values = col_hl_cls) +
  stat_summary(mapping = aes(HalfLife_cls_logPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_cls_logPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_cls_logPhase, feature_toViz, fill = HalfLife_cls_logPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
  ylim(-6, -1) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareTranscript_logPhase/featureDistributionTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_CDS_MFE_20_10nt_mid.png", width = 35, height = 18, units = "cm", dpi = 600)

###### CDS_MFE_50_25nt_3p
feature_i = "CDS_MFE_50_25nt_3p"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  group_by(HalfLife_cls_logPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_cls_logPhase)) +
  scale_colour_manual(values = col_hl_cls) +
  scale_fill_manual(values = col_hl_cls) +
  stat_summary(mapping = aes(HalfLife_cls_logPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_cls_logPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_cls_logPhase, feature_toViz, fill = HalfLife_cls_logPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
  ylim(-25, -5) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareTranscript_logPhase/featureDistributionTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_CDS_MFE_50_25nt_3p.png", width = 35, height = 18, units = "cm", dpi = 600)

###### CDS_MFE_50_25nt_5p
feature_i = "CDS_MFE_50_25nt_5p"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  group_by(HalfLife_cls_logPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_cls_logPhase)) +
  scale_colour_manual(values = col_hl_cls) +
  scale_fill_manual(values = col_hl_cls) +
  stat_summary(mapping = aes(HalfLife_cls_logPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_cls_logPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_cls_logPhase, feature_toViz, fill = HalfLife_cls_logPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
  ylim(-25, -10) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareTranscript_logPhase/featureDistributionTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_CDS_MFE_50_25nt_5p.png", width = 35, height = 18, units = "cm", dpi = 600)

###### AAA_CDSnstop
feature_i = "AAA_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  group_by(HalfLife_cls_logPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_cls_logPhase)) +
  scale_colour_manual(values = col_hl_cls) +
  scale_fill_manual(values = col_hl_cls) +
  stat_summary(mapping = aes(HalfLife_cls_logPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_cls_logPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_cls_logPhase, feature_toViz, fill = HalfLife_cls_logPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
  ylim(-0.01, 0.02) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareTranscript_logPhase/featureDistributionTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_AAA_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### CGC_CDSnstop
feature_i = "CGC_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  group_by(HalfLife_cls_logPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_cls_logPhase)) +
  scale_colour_manual(values = col_hl_cls) +
  scale_fill_manual(values = col_hl_cls) +
  stat_summary(mapping = aes(HalfLife_cls_logPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_cls_logPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_cls_logPhase, feature_toViz, fill = HalfLife_cls_logPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
  ylim(-0.01, 0.08) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareTranscript_logPhase/featureDistributionTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_CGC_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### CGG_CDSnstop
feature_i = "CGG_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  group_by(HalfLife_cls_logPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_cls_logPhase)) +
  scale_colour_manual(values = col_hl_cls) +
  scale_fill_manual(values = col_hl_cls) +
  stat_summary(mapping = aes(HalfLife_cls_logPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_cls_logPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_cls_logPhase, feature_toViz, fill = HalfLife_cls_logPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
  ylim(-0.01, 0.06) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareTranscript_logPhase/featureDistributionTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_CGG_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### CGU_CDSnstop
feature_i = "CGU_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  group_by(HalfLife_cls_logPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_cls_logPhase)) +
  scale_colour_manual(values = col_hl_cls) +
  scale_fill_manual(values = col_hl_cls) +
  stat_summary(mapping = aes(HalfLife_cls_logPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_cls_logPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_cls_logPhase, feature_toViz, fill = HalfLife_cls_logPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
  ylim(-0.01, 0.035) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareTranscript_logPhase/featureDistributionTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_CGU_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### CUG_CDSnstop
feature_i = "CUG_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  group_by(HalfLife_cls_logPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_cls_logPhase)) +
  scale_colour_manual(values = col_hl_cls) +
  scale_fill_manual(values = col_hl_cls) +
  stat_summary(mapping = aes(HalfLife_cls_logPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_cls_logPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_cls_logPhase, feature_toViz, fill = HalfLife_cls_logPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
  ylim(-0.01, 0.12) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareTranscript_logPhase/featureDistributionTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_CUG_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### GCC_CDSnstop
feature_i = "GCC_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  group_by(HalfLife_cls_logPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_cls_logPhase)) +
  scale_colour_manual(values = col_hl_cls) +
  scale_fill_manual(values = col_hl_cls) +
  stat_summary(mapping = aes(HalfLife_cls_logPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_cls_logPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_cls_logPhase, feature_toViz, fill = HalfLife_cls_logPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
  ylim(-0.01, 0.15) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareTranscript_logPhase/featureDistributionTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_GCC_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### GGC_CDSnstop
feature_i = "GGC_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  group_by(HalfLife_cls_logPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_cls_logPhase)) +
  scale_colour_manual(values = col_hl_cls) +
  scale_fill_manual(values = col_hl_cls) +
  stat_summary(mapping = aes(HalfLife_cls_logPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_cls_logPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_cls_logPhase, feature_toViz, fill = HalfLife_cls_logPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
  ylim(-0.01, 0.1) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareTranscript_logPhase/featureDistributionTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_GGC_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### GGG_CDSnstop
feature_i = "GGG_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  group_by(HalfLife_cls_logPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_cls_logPhase)) +
  scale_colour_manual(values = col_hl_cls) +
  scale_fill_manual(values = col_hl_cls) +
  stat_summary(mapping = aes(HalfLife_cls_logPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_cls_logPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_cls_logPhase, feature_toViz, fill = HalfLife_cls_logPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
  ylim(-0.01, 0.05) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareTranscript_logPhase/featureDistributionTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_GGG_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### GGU_CDSnstop
feature_i = "GGU_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  group_by(HalfLife_cls_logPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_cls_logPhase)) +
  scale_colour_manual(values = col_hl_cls) +
  scale_fill_manual(values = col_hl_cls) +
  stat_summary(mapping = aes(HalfLife_cls_logPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_cls_logPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_cls_logPhase, feature_toViz, fill = HalfLife_cls_logPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
  ylim(-0.01, 0.05) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareTranscript_logPhase/featureDistributionTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_GGU_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### GUG_CDSnstop
feature_i = "GUG_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  group_by(HalfLife_cls_logPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_cls_logPhase)) +
  scale_colour_manual(values = col_hl_cls) +
  scale_fill_manual(values = col_hl_cls) +
  stat_summary(mapping = aes(HalfLife_cls_logPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_cls_logPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_cls_logPhase, feature_toViz, fill = HalfLife_cls_logPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
  ylim(-0.01, 0.1) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareTranscript_logPhase/featureDistributionTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_GUG_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### UAU_CDSnstop
feature_i = "UAU_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  group_by(HalfLife_cls_logPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_cls_logPhase)) +
  scale_colour_manual(values = col_hl_cls) +
  scale_fill_manual(values = col_hl_cls) +
  stat_summary(mapping = aes(HalfLife_cls_logPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_cls_logPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_cls_logPhase, feature_toViz, fill = HalfLife_cls_logPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
  ylim(-0.005, 0.015) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareTranscript_logPhase/featureDistributionTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_UAU_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### UCC_CDSnstop
feature_i = "UCC_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  group_by(HalfLife_cls_logPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_cls_logPhase)) +
  scale_colour_manual(values = col_hl_cls) +
  scale_fill_manual(values = col_hl_cls) +
  stat_summary(mapping = aes(HalfLife_cls_logPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_cls_logPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_cls_logPhase, feature_toViz, fill = HalfLife_cls_logPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
  ylim(-0.01, 0.04) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareTranscript_logPhase/featureDistributionTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_UCC_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### UCG_CDSnstop
feature_i = "UCG_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  group_by(HalfLife_cls_logPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_cls_logPhase)) +
  scale_colour_manual(values = col_hl_cls) +
  scale_fill_manual(values = col_hl_cls) +
  stat_summary(mapping = aes(HalfLife_cls_logPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_cls_logPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_cls_logPhase, feature_toViz, fill = HalfLife_cls_logPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
  ylim(-0.01, 0.06) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareTranscript_logPhase/featureDistributionTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_UCG_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### UGG_CDSnstop
feature_i = "UGG_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  group_by(HalfLife_cls_logPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_cls_logPhase)) +
  scale_colour_manual(values = col_hl_cls) +
  scale_fill_manual(values = col_hl_cls) +
  stat_summary(mapping = aes(HalfLife_cls_logPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_cls_logPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_cls_logPhase, feature_toViz, fill = HalfLife_cls_logPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
  ylim(-0.01, 0.05) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareTranscript_logPhase/featureDistributionTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_UGG_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### UUG_CDSnstop
feature_i = "UUG_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  group_by(HalfLife_cls_logPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_cls_logPhase)) +
  scale_colour_manual(values = col_hl_cls) +
  scale_fill_manual(values = col_hl_cls) +
  stat_summary(mapping = aes(HalfLife_cls_logPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_cls_logPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_cls_logPhase, feature_toViz, fill = HalfLife_cls_logPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
  ylim(-0.01, 0.03) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareTranscript_logPhase/featureDistributionTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_UUG_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### GC_CDS_5p18nt
feature_i = "GC_CDS_5p18nt"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  group_by(HalfLife_cls_logPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_cls_logPhase)) +
  scale_colour_manual(values = col_hl_cls) +
  scale_fill_manual(values = col_hl_cls) +
  stat_summary(mapping = aes(HalfLife_cls_logPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_cls_logPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_cls_logPhase, feature_toViz, fill = HalfLife_cls_logPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
  ylim(0.25, 0.8) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareTranscript_logPhase/featureDistributionTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_GC_CDS_5p18nt.png", width = 35, height = 18, units = "cm", dpi = 600)

###### MFE_unfold
feature_i = "MFE_unfold"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  group_by(HalfLife_cls_logPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_cls_logPhase)) +
  scale_colour_manual(values = col_hl_cls) +
  scale_fill_manual(values = col_hl_cls) +
  stat_summary(mapping = aes(HalfLife_cls_logPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_cls_logPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_cls_logPhase, feature_toViz, fill = HalfLife_cls_logPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
  ylim(-5, 20) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareTranscript_logPhase/featureDistributionTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_MFE_unfold.png", width = 35, height = 18, units = "cm", dpi = 600)

###### Nucl_A_CDS_5p18nt
feature_i = "Nucl_A_CDS_5p18nt"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  group_by(HalfLife_cls_logPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_cls_logPhase)) +
  scale_colour_manual(values = col_hl_cls) +
  scale_fill_manual(values = col_hl_cls) +
  stat_summary(mapping = aes(HalfLife_cls_logPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_cls_logPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_cls_logPhase, feature_toViz, fill = HalfLife_cls_logPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
  ylim(-0.01, 0.5) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareTranscript_logPhase/featureDistributionTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_Nucl_A_CDS_5p18nt.png", width = 35, height = 18, units = "cm", dpi = 600)

###### ribo_5p
feature_i = "ribo_5p"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  group_by(HalfLife_cls_logPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_cls_logPhase)) +
  scale_colour_manual(values = col_hl_cls) +
  scale_fill_manual(values = col_hl_cls) +
  stat_summary(mapping = aes(HalfLife_cls_logPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_cls_logPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_cls_logPhase, feature_toViz, fill = HalfLife_cls_logPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
  ylim(-1, 6) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareTranscript_logPhase/featureDistributionTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_ribo_5p.png", width = 35, height = 18, units = "cm", dpi = 600)

###### ribo_5p_excl
feature_i = "ribo_5p_excl"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  group_by(HalfLife_cls_logPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_cls_logPhase)) +
  scale_colour_manual(values = col_hl_cls) +
  scale_fill_manual(values = col_hl_cls) +
  stat_summary(mapping = aes(HalfLife_cls_logPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_cls_logPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_cls_logPhase, feature_toViz, fill = HalfLife_cls_logPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
  ylim(0, 1.5) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareTranscript_logPhase/featureDistributionTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_ribo_5p_excl.png", width = 35, height = 18, units = "cm", dpi = 600)

###### CDS_length
feature_i = "CDS_length"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  group_by(HalfLife_cls_logPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_cls_logPhase)) +
  scale_colour_manual(values = col_hl_cls) +
  scale_fill_manual(values = col_hl_cls) +
  stat_summary(mapping = aes(HalfLife_cls_logPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_cls_logPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_cls_logPhase, feature_toViz, fill = HalfLife_cls_logPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
  ylim(-500, 2500) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareTranscript_logPhase/featureDistributionTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_CDS_length.png", width = 35, height = 18, units = "cm", dpi = 600)

###### initialAbundance_logPhase
feature_i = "initialAbundance_logPhase"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  rename(raw_value = all_of(feature_i)) %>%
  mutate(initialAbundance_logPhase = log2(raw_value)) %>%
  select(c("HalfLife_cls_logPhase", "initialAbundance_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  group_by(HalfLife_cls_logPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  rename(raw_value = all_of(feature_i)) %>%
  mutate(initialAbundance_logPhase = log2(raw_value)) %>%
  select(c("HalfLife_cls_logPhase", "initialAbundance_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_cls_logPhase)) +
  scale_colour_manual(values = col_hl_cls) +
  scale_fill_manual(values = col_hl_cls) +
  stat_summary(mapping = aes(HalfLife_cls_logPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_cls_logPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_cls_logPhase, feature_toViz, fill = HalfLife_cls_logPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Fast', 'Med-fast', 'Med-slow', 'Slow')) +
  ylim(8, 20) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareTranscript_logPhase/featureDistributionTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_initialAbundance_logPhase.png", width = 35, height = 18, units = "cm", dpi = 600)
