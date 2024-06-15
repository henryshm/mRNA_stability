
####### viz of feature distributions
####### To compare features of log phase, hypoxia and fold change in hypoxia
####### for leadered and leaderless transcripts 
####### features are the top20 important features for the three conditions
####### feature selected by log phase in combined selected model
####### leadered transcript
####### fold change in hypoxia

library(tidyverse)
library(RColorBrewer)
library(PupillometryR)
library(scales)

####### get color scheme for feature distributions
col_FC_hypoxiaToLogPhase <- c("#E72F52", "#774FA0", "#0D95D0", "#7DC462")

####### get feature table
featureTable <- read.csv('../../feature/FeatureTables/featureTable_combinedSelected_halfLifeFcClass/HalfLifeFcCls_hypoxiaToLogPhase_combinedSelected_leadered_byLeaderedLogPhase_jointlyByLogPhaseFoldChange.csv')

####### viz of feature distributions
###### Nucl_C_5pUTR
feature_i = "Nucl_C_5pUTR"

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
  ylim(0.1, 0.5) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_Nucl_C_5pUTR.png", width = 35, height = 18, units = "cm", dpi = 600)

###### Nucl_G_5pUTR
feature_i = "Nucl_G_5pUTR"

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
  ylim(0.1, 0.5) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_Nucl_G_5pUTR.png", width = 35, height = 18, units = "cm", dpi = 600)

###### Nucl_T_5pUTR
feature_i = "Nucl_T_5pUTR"

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
  ylim(0, 0.3) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_Nucl_T_5pUTR.png", width = 35, height = 18, units = "cm", dpi = 600)

###### adja_AC_5pUTR
feature_i = "adja_AC_5pUTR"

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
  ylim(-0.02, 0.17) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_adja_AC_5pUTR.png", width = 35, height = 18, units = "cm", dpi = 600)

###### adja_CA_5pUTR
feature_i = "adja_CA_5pUTR"

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
  ylim(-0.02, 0.17) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_adja_CA_5pUTR.png", width = 35, height = 18, units = "cm", dpi = 600)

###### adja_GC_5pUTR
feature_i = "adja_GC_5pUTR"

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
  ylim(-0.02, 0.23) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_adja_GC_5pUTR.png", width = 35, height = 18, units = "cm", dpi = 600)

###### adja_GU_5pUTR
feature_i = "adja_GU_5pUTR"

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
  ylim(-0.02, 0.15) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_adja_GU_5pUTR.png", width = 35, height = 18, units = "cm", dpi = 600)

###### adja_UG_5pUTR
feature_i = "adja_UG_5pUTR"

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
  ylim(-0.02, 0.15) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_adja_UG_5pUTR.png", width = 35, height = 18, units = "cm", dpi = 600)

###### fpr_UTR_length
feature_i = "fpr_UTR_length"

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
  ylim(-50, 375) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_fpr_UTR_length.png", width = 35, height = 18, units = "cm", dpi = 600)

###### adja_GG_CDS
feature_i = "adja_GG_CDS"

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
  ylim(0.03, 0.15) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_adja_GG_CDS.png", width = 35, height = 18, units = "cm", dpi = 600)

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
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_adja_GU_CDS.png", width = 35, height = 18, units = "cm", dpi = 600)

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
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_adja_UA_CDS.png", width = 35, height = 18, units = "cm", dpi = 600)

###### CDS_MFE_20_10nt_5p
feature_i = "CDS_MFE_20_10nt_5p"

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
  ylim(-6, 0) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_CDS_MFE_20_10nt_5p.png", width = 35, height = 18, units = "cm", dpi = 600)

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
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_AAA_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### AAG_CDSnstop
feature_i = "AAG_CDSnstop"

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
  ylim(-0.025, 0.1) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_AAG_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### ACG_CDSnstop
feature_i = "ACG_CDSnstop"

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
  ylim(-0.015, 0.055) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_ACG_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

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
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_CGC_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

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
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_CGG_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

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
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_CGU_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

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
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_CUG_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

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
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_CodonPairBias.png", width = 35, height = 18, units = "cm", dpi = 600)

###### GAU_CDSnstop
feature_i = "GAU_CDSnstop"

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
  ylim(-0.015, 0.05) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_GAU_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### GCA_CDSnstop
feature_i = "GCA_CDSnstop"

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
  ylim(-0.015, 0.05) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_GCA_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

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
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_GCC_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

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
  ylim(-0.02, 0.13) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_GCG_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### GGA_CDSnstop
feature_i = "GGA_CDSnstop"

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
  ylim(-0.015, 0.04) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_GGA_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

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
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_GGC_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

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
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_GGG_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

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
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_GGU_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

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
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_UAU_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

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
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_UGG_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### UUG_CDSnstop
feature_i = "UUG_CDSnstop"

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
  ylim(-0.01, 0.03) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_UUG_CDSnstop.png", width = 35, height = 18, units = "cm", dpi = 600)

###### MFE_unfold
feature_i = "MFE_unfold"

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
  ylim(-5, 25) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_MFE_unfold.png", width = 35, height = 18, units = "cm", dpi = 600)

###### fprUTRlast30nt_plus20ntCDS_2ntAvg_WholeSeq_UnpairedProb
feature_i = "fprUTRlast30nt_plus20ntCDS_2ntAvg_WholeSeq_UnpairedProb"

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
  ylim(0, 0.8) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_fprUTRlast30nt_plus20ntCDS_2ntAvg_WholeSeq_UnpairedProb.png", width = 35, height = 18, units = "cm", dpi = 600)

###### fprUTRlast30nt_plus20ntCDS_all_StartCodon_UnpairedProb
feature_i = "fprUTRlast30nt_plus20ntCDS_all_StartCodon_UnpairedProb"

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
  ylim(-0.5, 1.2) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_fprUTRlast30nt_plus20ntCDS_all_StartCodon_UnpairedProb.png", width = 35, height = 18, units = "cm", dpi = 600)

###### fprUTRlast30nt_plusStartCodon_2ntAvg_SDseq_UnpairedProb
feature_i = "fprUTRlast30nt_plusStartCodon_2ntAvg_SDseq_UnpairedProb"

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
  ylim(-0.2, 1.2) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_fprUTRlast30nt_plusStartCodon_2ntAvg_SDseq_UnpairedProb.png", width = 35, height = 18, units = "cm", dpi = 600)

###### fprUTRlast30nt_plusStartCodon_2ntAvg_WholeSeq_UnpairedProb
feature_i = "fprUTRlast30nt_plusStartCodon_2ntAvg_WholeSeq_UnpairedProb"

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
  ylim(0, 1) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_fprUTRlast30nt_plusStartCodon_2ntAvg_WholeSeq_UnpairedProb.png", width = 35, height = 18, units = "cm", dpi = 600)

###### fprUTRlast30nt_plusStartCodon_all_StartCodon_UnpairedProb
feature_i = "fprUTRlast30nt_plusStartCodon_all_StartCodon_UnpairedProb"

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
  ylim(-0.3, 1.5) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_fprUTRlast30nt_plusStartCodon_all_StartCodon_UnpairedProb.png", width = 35, height = 18, units = "cm", dpi = 600)

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
  ylim(-0.5, 2.5) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_ribo_5p.png", width = 35, height = 18, units = "cm", dpi = 600)

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
  ylim(-0.7, 2) +
  coord_flip()
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_ribo_5p_excl.png", width = 35, height = 18, units = "cm", dpi = 600)

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
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_CDS_length.png", width = 35, height = 18, units = "cm", dpi = 600)

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
ggsave("./viz_featureDistributions_toCompareCondition_leadered/featureDistributionTop20combined_toCompareCondition_HalfLifeFcCls_hypoxiaToLogPhase_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseFoldChange_leadered_initialAbundance_logPhase.png", width = 35, height = 18, units = "cm", dpi = 600)
