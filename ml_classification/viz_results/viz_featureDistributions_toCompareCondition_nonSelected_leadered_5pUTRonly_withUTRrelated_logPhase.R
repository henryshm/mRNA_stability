
####### viz of feature distributions
####### To compare features of log phase, hypoxia and fold change in hypoxia
####### for leadered transcripts non-selected feature models with only 5'UTR features
####### features are the top20 important features for the three conditions
####### with only UTR related features
####### log phase

library(tidyverse)
library(RColorBrewer)
library(PupillometryR)
library(scales)

####### get color scheme for feature distributions
col_hl_cls <- c("#7DC462", "#0D95D0", "#774FA0", "#E72F52")

####### get feature table
featureTable <- read.csv('../../feature/FeatureTables/featureTable_separatedByTypeNonSelected_halfLifeClass/HalfLifeCls_logPhase_nonSelected_leadered_5pUTR_UTRrelated.csv')

####### viz of feature distributions
###### Nucl_G_5pUTR
feature_i = "Nucl_G_5pUTR"

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
  ylim(0.1, 0.5) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered_5pUTRonly_withUTRrelated/featureDistributionTop20combined_toCompareCondition_HalfLifeCls_logPhase_nonSelected_leadered_5pUTRonly_withUTRrelated_Nucl_G_5pUTR.png", width = 35, height = 18, units = "cm", dpi = 600)

###### adja_UG_5pUTR
feature_i = "adja_UG_5pUTR"

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
  ylim(-0.02, 0.15) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered_5pUTRonly_withUTRrelated/featureDistributionTop20combined_toCompareCondition_HalfLifeCls_logPhase_nonSelected_leadered_5pUTRonly_withUTRrelated_adja_UG_5pUTR.png", width = 35, height = 18, units = "cm", dpi = 600)

###### adja_GC_5pUTR
feature_i = "adja_GC_5pUTR"

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
  ylim(-0.02, 0.23) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered_5pUTRonly_withUTRrelated/featureDistributionTop20combined_toCompareCondition_HalfLifeCls_logPhase_nonSelected_leadered_5pUTRonly_withUTRrelated_adja_GC_5pUTR.png", width = 35, height = 18, units = "cm", dpi = 600)

###### transcript_5p20nt_Avg_First5nt_UnpairedProb
feature_i = "transcript_5p20nt_Avg_First5nt_UnpairedProb"

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
  ylim(-0.2, 1.3) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered_5pUTRonly_withUTRrelated/featureDistributionTop20combined_toCompareCondition_HalfLifeCls_logPhase_nonSelected_leadered_5pUTRonly_withUTRrelated_transcript_5p20nt_Avg_First5nt_UnpairedProb.png", width = 35, height = 18, units = "cm", dpi = 600)

###### adja_CA_5pUTR
feature_i = "adja_CA_5pUTR"

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
  ylim(-0.02, 0.17) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered_5pUTRonly_withUTRrelated/featureDistributionTop20combined_toCompareCondition_HalfLifeCls_logPhase_nonSelected_leadered_5pUTRonly_withUTRrelated_adja_CA_5pUTR.png", width = 35, height = 18, units = "cm", dpi = 600)

###### adja_UC_5pUTR
feature_i = "adja_UC_5pUTR"

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
  ylim(-0.02, 0.17) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered_5pUTRonly_withUTRrelated/featureDistributionTop20combined_toCompareCondition_HalfLifeCls_logPhase_nonSelected_leadered_5pUTRonly_withUTRrelated_adja_UC_5pUTR.png", width = 35, height = 18, units = "cm", dpi = 600)

###### fprUTR_MFE_20_10nt_mid
feature_i = "fprUTR_MFE_20_10nt_mid"

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
  ylim(-10, 3) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered_5pUTRonly_withUTRrelated/featureDistributionTop20combined_toCompareCondition_HalfLifeCls_logPhase_nonSelected_leadered_5pUTRonly_withUTRrelated_fprUTR_MFE_20_10nt_mid.png", width = 35, height = 18, units = "cm", dpi = 600)

###### adja_CC_5pUTR
feature_i = "adja_CC_5pUTR"

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
  ylim(-0.03, 0.23) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered_5pUTRonly_withUTRrelated/featureDistributionTop20combined_toCompareCondition_HalfLifeCls_logPhase_nonSelected_leadered_5pUTRonly_withUTRrelated_adja_CC_5pUTR.png", width = 35, height = 18, units = "cm", dpi = 600)

###### adja_GU_5pUTR
feature_i = "adja_GU_5pUTR"

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
  ylim(-0.02, 0.15) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered_5pUTRonly_withUTRrelated/featureDistributionTop20combined_toCompareCondition_HalfLifeCls_logPhase_nonSelected_leadered_5pUTRonly_withUTRrelated_adja_GU_5pUTR.png", width = 35, height = 18, units = "cm", dpi = 600)

###### adja_CU_5pUTR
feature_i = "adja_CU_5pUTR"

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
  ylim(-0.03, 0.13) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered_5pUTRonly_withUTRrelated/featureDistributionTop20combined_toCompareCondition_HalfLifeCls_logPhase_nonSelected_leadered_5pUTRonly_withUTRrelated_adja_CU_5pUTR.png", width = 35, height = 18, units = "cm", dpi = 600)

###### transcript_5p20nt_all_First3nt_UnpairedProb
feature_i = "transcript_5p20nt_all_First3nt_UnpairedProb"

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
  ylim(-0.4, 1.4) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered_5pUTRonly_withUTRrelated/featureDistributionTop20combined_toCompareCondition_HalfLifeCls_logPhase_nonSelected_leadered_5pUTRonly_withUTRrelated_transcript_5p20nt_all_First3nt_UnpairedProb.png", width = 35, height = 18, units = "cm", dpi = 600)

###### Nucl_A_5pUTR
feature_i = "Nucl_A_5pUTR"

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
  ylim(-0.05, 0.45) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered_5pUTRonly_withUTRrelated/featureDistributionTop20combined_toCompareCondition_HalfLifeCls_logPhase_nonSelected_leadered_5pUTRonly_withUTRrelated_Nucl_A_5pUTR.png", width = 35, height = 18, units = "cm", dpi = 600)

###### GC_5pUTR
feature_i = "GC_5pUTR"

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
  ylim(0.38, 0.84) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered_5pUTRonly_withUTRrelated/featureDistributionTop20combined_toCompareCondition_HalfLifeCls_logPhase_nonSelected_leadered_5pUTRonly_withUTRrelated_GC_5pUTR.png", width = 35, height = 18, units = "cm", dpi = 600)

###### adja_GA_5pUTR
feature_i = "adja_GA_5pUTR"

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
  ylim(-0.05, 0.22) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered_5pUTRonly_withUTRrelated/featureDistributionTop20combined_toCompareCondition_HalfLifeCls_logPhase_nonSelected_leadered_5pUTRonly_withUTRrelated_adja_GA_5pUTR.png", width = 35, height = 18, units = "cm", dpi = 600)

###### adja_GG_5pUTR
feature_i = "adja_GG_5pUTR"

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
  ylim(-0.05, 0.22) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered_5pUTRonly_withUTRrelated/featureDistributionTop20combined_toCompareCondition_HalfLifeCls_logPhase_nonSelected_leadered_5pUTRonly_withUTRrelated_adja_GG_5pUTR.png", width = 35, height = 18, units = "cm", dpi = 600)

###### Nucl_C_5pUTR
feature_i = "Nucl_C_5pUTR"

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
  ylim(0.1, 0.5) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered_5pUTRonly_withUTRrelated/featureDistributionTop20combined_toCompareCondition_HalfLifeCls_logPhase_nonSelected_leadered_5pUTRonly_withUTRrelated_Nucl_C_5pUTR.png", width = 35, height = 18, units = "cm", dpi = 600)

###### transcript_5p20nt_Avg_First3nt_UnpairedProb
feature_i = "transcript_5p20nt_Avg_First3nt_UnpairedProb"

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
  ylim(-0.4, 1.4) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered_5pUTRonly_withUTRrelated/featureDistributionTop20combined_toCompareCondition_HalfLifeCls_logPhase_nonSelected_leadered_5pUTRonly_withUTRrelated_transcript_5p20nt_Avg_First3nt_UnpairedProb.png", width = 35, height = 18, units = "cm", dpi = 600)

###### transcript_5p20nt_all_First5nt_UnpairedProb
feature_i = "transcript_5p20nt_all_First5nt_UnpairedProb"

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
  ylim(-0.4, 1.4) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered_5pUTRonly_withUTRrelated/featureDistributionTop20combined_toCompareCondition_HalfLifeCls_logPhase_nonSelected_leadered_5pUTRonly_withUTRrelated_transcript_5p20nt_all_First5nt_UnpairedProb.png", width = 35, height = 18, units = "cm", dpi = 600)

###### adja_AG_5pUTR
feature_i = "adja_AG_5pUTR"

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
  ylim(-0.02, 0.17) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered_5pUTRonly_withUTRrelated/featureDistributionTop20combined_toCompareCondition_HalfLifeCls_logPhase_nonSelected_leadered_5pUTRonly_withUTRrelated_adja_AG_5pUTR.png", width = 35, height = 18, units = "cm", dpi = 600)

###### adja_AA_5pUTR
feature_i = "adja_AA_5pUTR"

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
  ylim(-0.03, 0.15) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered_5pUTRonly_withUTRrelated/featureDistributionTop20combined_toCompareCondition_HalfLifeCls_logPhase_nonSelected_leadered_5pUTRonly_withUTRrelated_adja_AA_5pUTR.png", width = 35, height = 18, units = "cm", dpi = 600)

###### fpr_UTR_length
feature_i = "fpr_UTR_length"

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
  ylim(-50, 375) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered_5pUTRonly_withUTRrelated/featureDistributionTop20combined_toCompareCondition_HalfLifeCls_logPhase_nonSelected_leadered_5pUTRonly_withUTRrelated_fpr_UTR_length.png", width = 35, height = 18, units = "cm", dpi = 600)

###### adja_UU_5pUTR
feature_i = "adja_UU_5pUTR"

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
  ylim(-0.03, 0.11) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered_5pUTRonly_withUTRrelated/featureDistributionTop20combined_toCompareCondition_HalfLifeCls_logPhase_nonSelected_leadered_5pUTRonly_withUTRrelated_adja_UU_5pUTR.png", width = 35, height = 18, units = "cm", dpi = 600)

###### Nucl_T_5pUTR
feature_i = "Nucl_T_5pUTR"

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
  ylim(0, 0.3) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered_5pUTRonly_withUTRrelated/featureDistributionTop20combined_toCompareCondition_HalfLifeCls_logPhase_nonSelected_leadered_5pUTRonly_withUTRrelated_Nucl_T_5pUTR.png", width = 35, height = 18, units = "cm", dpi = 600)

###### adja_AU_5pUTR
feature_i = "adja_AU_5pUTR"

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
  ylim(-0.03, 0.1) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered_5pUTRonly_withUTRrelated/featureDistributionTop20combined_toCompareCondition_HalfLifeCls_logPhase_nonSelected_leadered_5pUTRonly_withUTRrelated_adja_AU_5pUTR.png", width = 35, height = 18, units = "cm", dpi = 600)

###### fprUTR_MFE_20_10nt_5pUTR
feature_i = "fprUTR_MFE_20_10nt_5pUTR"

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
  ylim(-7.5, 1) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered_5pUTRonly_withUTRrelated/featureDistributionTop20combined_toCompareCondition_HalfLifeCls_logPhase_nonSelected_leadered_5pUTRonly_withUTRrelated_fprUTR_MFE_20_10nt_5pUTR.png", width = 35, height = 18, units = "cm", dpi = 600)

###### adja_CG_5pUTR
feature_i = "adja_CG_5pUTR"

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
  ylim(-0.01, 0.23) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered_5pUTRonly_withUTRrelated/featureDistributionTop20combined_toCompareCondition_HalfLifeCls_logPhase_nonSelected_leadered_5pUTRonly_withUTRrelated_adja_CG_5pUTR.png", width = 35, height = 18, units = "cm", dpi = 600)

###### fprUTR_MFE_20_10nt_5p
feature_i = "fprUTR_MFE_20_10nt_5p"

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
  ylim(-10, 3) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered_5pUTRonly_withUTRrelated/featureDistributionTop20combined_toCompareCondition_HalfLifeCls_logPhase_nonSelected_leadered_5pUTRonly_withUTRrelated_fprUTR_MFE_20_10nt_5p.png", width = 35, height = 18, units = "cm", dpi = 600)

###### adja_AC_5pUTR
feature_i = "adja_AC_5pUTR"

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
  ylim(-0.02, 0.17) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered_5pUTRonly_withUTRrelated/featureDistributionTop20combined_toCompareCondition_HalfLifeCls_logPhase_nonSelected_leadered_5pUTRonly_withUTRrelated_adja_AC_5pUTR.png", width = 35, height = 18, units = "cm", dpi = 600)
