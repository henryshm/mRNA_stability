
####### viz of feature distributions
####### To compare features of log phase, hypoxia and fold change in hypoxia
####### for leadered and leaderless transcripts 
####### features are the top20 important features for the three conditions
####### feature selected by log phase in combined selected model
####### leaderless transcript
####### fold change in hypoxia

library(tidyverse)
library(RColorBrewer)
library(PupillometryR)
library(scales)

####### get color scheme for feature distributions
col_FC_hypoxiaToLogPhase <- c("#E72F52", "#774FA0", "#0D95D0", "#7DC462")

####### get feature table
featureTable <- read.csv('../../feature/FeatureTables/featureTable_combinedSelected_halfLifeFcClass/HalfLifeFcCls_hypoxiaToLogPhase_combinedSelected_leaderless_byLeaderlessLogPhase_jointlyByLogPhaseFoldChange.csv')

####### viz of feature distributions
###### adja_AU_CDS
feature_i = "adja_AU_CDS"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(0, 0.06) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_adja_AU_CDS.png", width = 35, height = 18, units = "cm", dpi = 600)

###### adja_CA_CDS
feature_i = "adja_CA_CDS"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(0.01, 0.1) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_adja_CA_CDS.png", width = 35, height = 18, units = "cm", dpi = 600)

###### adja_CG_CDS
feature_i = "adja_CG_CDS"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(0.09, 0.2) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_adja_CG_CDS.png", width = 35, height = 18, units = "cm", dpi = 600)

###### adja_CU_CDS
feature_i = "adja_CU_CDS"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(0.015, 0.081) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_adja_CU_CDS.png", width = 35, height = 18, units = "cm", dpi = 600)

###### adja_GU_CDS
feature_i = "adja_GU_CDS"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(0.02, 0.1) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_adja_GU_CDS.png", width = 35, height = 18, units = "cm", dpi = 600)

###### adja_UA_CDS
feature_i = "adja_UA_CDS"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-0.01, 0.025) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_adja_UA_CDS.png", width = 35, height = 18, units = "cm", dpi = 600)

###### CDS_MFE_100_50nt_3p
feature_i = "CDS_MFE_100_50nt_3p"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-55, -10) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_CDS_MFE_100_50nt_3p.png", width = 35, height = 18, units = "cm", dpi = 600)

###### CDS_MFE_100_50nt_5p
feature_i = "CDS_MFE_100_50nt_5p"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-55, -20) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_CDS_MFE_100_50nt_5p.png", width = 35, height = 18, units = "cm", dpi = 600)

###### CDS_MFE_20_10nt_mid
feature_i = "CDS_MFE_20_10nt_mid"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-6, -1) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_CDS_MFE_20_10nt_mid.png", width = 35, height = 18, units = "cm", dpi = 600)

###### CDS_MFE_50_25nt_3p
feature_i = "CDS_MFE_50_25nt_3p"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-25, -5) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_CDS_MFE_50_25nt_3p.png", width = 35, height = 18, units = "cm", dpi = 600)

###### CDS_MFE_50_25nt_5p
feature_i = "CDS_MFE_50_25nt_5p"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-25, -10) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_CDS_MFE_50_25nt_5p.png", width = 35, height = 18, units = "cm", dpi = 600)

###### AAA_CDSnstop
feature_i = "AAA_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-0.01, 0.02) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_AAA_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### ACA_CDSnstop
feature_i = "ACA_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-0.005, 0.025) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_ACA_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### AGC_CDSnstop
feature_i = "AGC_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-0.007, 0.047) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_AGC_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### AUC_CDSnstop
feature_i = "AUC_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-0.005, 0.11) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_AUC_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### CAA_CDSnstop
feature_i = "CAA_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-0.005, 0.022) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_CAA_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### CAG_CDSnstop
feature_i = "CAG_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-0.01, 0.07) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_CAG_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### CGC_CDSnstop
feature_i = "CGC_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-0.01, 0.08) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_CGC_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### CGG_CDSnstop
feature_i = "CGG_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-0.01, 0.06) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_CGG_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### CGU_CDSnstop
feature_i = "CGU_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-0.01, 0.04) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_CGU_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### CUC_CDSnstop
feature_i = "CUC_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-0.015, 0.08) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_CUC_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### CUG_CDSnstop
feature_i = "CUG_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-0.01, 0.12) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_CUG_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### CUU_CDSnstop
feature_i = "CUU_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-0.0035, 0.02) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_CUU_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### CodonPairBias
feature_i = "CodonPairBias"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-0.1, 0.2) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_CodonPairBias.png", width = 35, height = 18, units = "cm", dpi = 600)

###### GAA_CDSnstop
feature_i = "GAA_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-0.015, 0.06) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_GAA_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### GCC_CDSnstop
feature_i = "GCC_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-0.01, 0.15) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_GCC_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### GCG_CDSnstop
feature_i = "GCG_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-0.015, 0.15) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_GCG_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### GGC_CDSnstop
feature_i = "GGC_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-0.01, 0.1) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_GGC_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### GGG_CDSnstop
feature_i = "GGG_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-0.01, 0.05) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_GGG_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### GGU_CDSnstop
feature_i = "GGU_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-0.01, 0.05) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_GGU_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### GUC_CDSnstop
feature_i = "GUC_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-0.01, 0.1) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_GUC_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### GUG_CDSnstop
feature_i = "GUG_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-0.01, 0.1) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_GUG_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### UAU_CDSnstop
feature_i = "UAU_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-0.005, 0.015) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_UAU_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### UCC_CDSnstop
feature_i = "UCC_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-0.01, 0.04) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_UCC_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### UCG_CDSnstop
feature_i = "UCG_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-0.01, 0.06) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_UCG_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### UGC_CDSnstop
feature_i = "UGC_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-0.005, 0.035) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_UGC_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### UGG_CDSnstop
feature_i = "UGG_CDSnstop"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-0.01, 0.05) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_UGG_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### GC_CDS_5p18nt
feature_i = "GC_CDS_5p18nt"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(0.25, 0.8) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_GC_CDS_5p18nt.png", width = 35, height = 18, units = "cm", dpi = 600)

###### Nucl_A_CDS_5p18nt
feature_i = "Nucl_A_CDS_5p18nt"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-0.01, 0.5) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_Nucl_A_CDS_5p18nt.png", width = 35, height = 18, units = "cm", dpi = 600)

###### ribo_5p
feature_i = "ribo_5p"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-1, 6) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_ribo_5p.png", width = 35, height = 18, units = "cm", dpi = 600)

###### ribo_5p_excl
feature_i = "ribo_5p_excl"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(0, 1.5) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_ribo_5p_excl.png", width = 35, height = 18, units = "cm", dpi = 600)

###### CDS_length
feature_i = "CDS_length"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-500, 2500) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_CDS_length.png", width = 35, height = 18, units = "cm", dpi = 600)

###### initialAbundance
feature_i = "initialAbundance_logPhase"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  rename(raw_value = all_of(feature_i)) %>%
  mutate(initialAbundance_logPhase = log2(raw_value)) %>%
  select(c("HalfLife_FCcls_hypoxiaToLogPhase", "initialAbundance_logPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  rename(raw_value = all_of(feature_i)) %>%
  mutate(initialAbundance_logPhase = log2(raw_value)) %>%
  select(c("HalfLife_FCcls_hypoxiaToLogPhase", "initialAbundance_logPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(6, 22) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_initialAbundance_logPhase.png", width = 35, height = 18, units = "cm", dpi = 600)

###### thrprUTR_MFE_20_10nt_3p
feature_i = "thrprUTR_MFE_20_10nt_3p"

featureTable_den <- featureTable %>% 
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  group_by(HalfLife_FCcls_hypoxiaToLogPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  select(c(all_of(feature_i), "HalfLife_FCcls_hypoxiaToLogPhase")) %>%
  mutate_at("HalfLife_FCcls_hypoxiaToLogPhase", as.factor) %>%
  rename(feature_toViz = all_of(feature_i)) %>%
  ggplot(aes(colour = HalfLife_FCcls_hypoxiaToLogPhase)) +
  scale_colour_manual(values = col_FC_hypoxiaToLogPhase) +
  scale_fill_manual(values = col_FC_hypoxiaToLogPhase) +
  stat_summary(mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, group = 1),
               fun = "median", geom = "line", color = "#424242", linewidth = 6) +
  stat_summary(
    mapping = aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz),
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median, size = 5,
    linewidth = 5) +
  geom_flat_violin(data = featureTable_den, aes(HalfLife_FCcls_hypoxiaToLogPhase, feature_toViz, fill = HalfLife_FCcls_hypoxiaToLogPhase),
                   alpha = 0.5, trim = FALSE, adjust = 2, linetype = "blank") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#E8EBEC", linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 50, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.key.size = unit(2, 'cm'), legend.title = element_blank()) +
  scale_x_discrete(limits = c('Small', 'Med-small', 'Med-large', 'Large')) +
  ylim(-13, 5) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leaderless/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseFoldChange_leaderless_thrprUTR_MFE_20_10nt_3p.png", width = 35, height = 18, units = "cm", dpi = 600)
