
####### check correlation between CDS MFE with regional ribosome profiling in leaderless models
####### observed stronger/linear in logPhase, but not log phase slow class
####### to see if the slow half-life class in log phase is protected by ribosome binding
####### and therefore related with translation
####### using features in the top20 important features selected by log phase in combined selected model
####### ribosome profile using one third of CDS (5'end, mid and 3'end)
####### leaderless transcript
####### log phase

library(tidyverse)
library(RColorBrewer)
library(PupillometryR)
library(scales)

####### get color scheme for half-life classes
col_hl_cls <- c("#7DC462", "#0D95D0", "#774FA0", "#E72F52")

####### get half-life of entire gene
hl_CDS_complete_class <- read.table('../../degradation_classes/half_life/CDS_halfLife_completeClass.txt', header = TRUE)

hl_CDS_logPhase <- hl_CDS_complete_class %>% 
  select(gene, HalfLife_vals_logPhase) %>% 
  drop_na()

####### get feature table
featureTable <- read.csv('../../feature/FeatureTables/featureTable_combinedSelected_halfLifeClass/HalfLifeCls_logPhase_combinedSelected_leaderless_byLeaderlessLogPhase_jointlyByLogPhaseHypoxia.csv')

####### get regional ribosome profiles
ribo_oneThird <- read.table('../../feature/ribosome_profiling/ribo_normalizedBy_totalmRNA_oneThird.txt', header = FALSE, 
                            col.names = c("gene", "ribo_CDSoneThird_5p", "ribo_CDSoneThird_mid", "ribo_CDSoneThird_3p"))

####### distributions of regional ribosome profiles
###### one third 5'end
feature_i = "ribo_CDSoneThird_5p"

featureTable_den <- featureTable %>% 
  inner_join(ribo_oneThird, by = "gene") %>% 
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  group_by(HalfLife_cls_logPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  inner_join(ribo_oneThird, by = "gene") %>%
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
  ylim(0, 2) +
  coord_flip()
ggsave("./viz_CDSmfeCorr_withRegionalRibosome/featureDistribution_HalfLifeCls_logPhase_leaderless_ribo_CDSoneThird_5p.png", width = 35, height = 18, units = "cm", dpi = 600)

###### one third mid
feature_i = "ribo_CDSoneThird_mid"

featureTable_den <- featureTable %>% 
  inner_join(ribo_oneThird, by = "gene") %>% 
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  group_by(HalfLife_cls_logPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  inner_join(ribo_oneThird, by = "gene") %>%
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
  ylim(0, 2) +
  coord_flip()
ggsave("./viz_CDSmfeCorr_withRegionalRibosome/featureDistribution_HalfLifeCls_logPhase_leaderless_ribo_CDSoneThird_mid.png", width = 35, height = 18, units = "cm", dpi = 600)

###### one third 3'end
feature_i = "ribo_CDSoneThird_3p"

featureTable_den <- featureTable %>% 
  inner_join(ribo_oneThird, by = "gene") %>% 
  select(c(all_of(feature_i), "HalfLife_cls_logPhase")) %>%
  mutate_at("HalfLife_cls_logPhase", as.factor) %>%
  group_by(HalfLife_cls_logPhase) %>%
  rename(feature_toViz = all_of(feature_i))

featureTable %>%
  inner_join(ribo_oneThird, by = "gene") %>%
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
  ylim(0, 2) +
  coord_flip()
ggsave("./viz_CDSmfeCorr_withRegionalRibosome/featureDistribution_HalfLifeCls_logPhase_leaderless_ribo_CDSoneThird_3p.png", width = 35, height = 18, units = "cm", dpi = 600)

####### correlations between half-lives and regional ribosome profile within each MFE window split by MFE values
###### get CDS MFE and split into subgroups by MFE values
##### CDS_MFE_50_25nt_5p
feature_CDSmfe_50_5p <- featureTable %>% 
  select(gene, CDS_MFE_50_25nt_5p)

feature_CDSmfe_50_5p_subset <- feature_CDSmfe_50_5p %>% 
  mutate(logPhase_CDSmfe_50_5p_subset = case_when(CDS_MFE_50_25nt_5p <= -17.38 ~ "logPhase_CDSmfe_50_5p_G1",
                                                  CDS_MFE_50_25nt_5p > -17.38 & CDS_MFE_50_25nt_5p <= -16.09 ~ "logPhase_CDSmfe_50_5p_G2",
                                                  CDS_MFE_50_25nt_5p > -16.09 & CDS_MFE_50_25nt_5p <= -14.74 ~ "logPhase_CDSmfe_50_5p_G3",
                                                  CDS_MFE_50_25nt_5p > -14.74 ~ "logPhase_CDSmfe_50_5p_G4"))

feature_CDSmfe_50_5p_subset_ribo <- feature_CDSmfe_50_5p_subset %>% 
  inner_join(ribo_oneThird, by = "gene") %>% 
  inner_join(hl_CDS_logPhase, by = "gene")

feature_CDSmfe_50_5p_subset_ribo %>% 
  group_by(logPhase_CDSmfe_50_5p_subset) %>%
  summarize(cor = cor.test(HalfLife_vals_logPhase, ribo_CDSoneThird_5p, method = "spearman")$estimate,
            count = n())

feature_CDSmfe_50_5p_subset_ribo %>% 
  ggplot(aes(x = ribo_CDSoneThird_5p, y = HalfLife_vals_logPhase)) + 
  geom_point(size = 4, color = "#036D9C", 
             alpha = 0.3,
             stroke = 0) +
  xlim(0, 2) +
  facet_wrap(~ logPhase_CDSmfe_50_5p_subset) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.title =  element_text(size = 25, face = "bold"),
        legend.position = "None",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15, face = "bold")
  )

###### CDS_MFE_100_50nt_5p
feature_CDSmfe_100_5p <- featureTable %>% 
  select(gene, CDS_MFE_100_50nt_5p)

summary(feature_CDSmfe_100_5p$CDS_MFE_100_50nt_5p)

feature_CDSmfe_100_5p_subset <- feature_CDSmfe_100_5p %>% 
  mutate(logPhase_CDSmfe_100_5p_subset = case_when(CDS_MFE_100_50nt_5p <= -41.72 ~ "logPhase_CDSmfe_100_5p_G1",
                                                   CDS_MFE_100_50nt_5p > -41.72 & CDS_MFE_100_50nt_5p <= -39.13 ~ "logPhase_CDSmfe_100_5p_G2",
                                                   CDS_MFE_100_50nt_5p > -39.13 & CDS_MFE_100_50nt_5p <= -36.32 ~ "logPhase_CDSmfe_100_5p_G3",
                                                   CDS_MFE_100_50nt_5p > -36.32 ~ "logPhase_CDSmfe_100_5p_G4"))

feature_CDSmfe_100_5p_subset_ribo <- feature_CDSmfe_100_5p_subset %>% 
  inner_join(ribo_oneThird, by = "gene") %>% 
  inner_join(hl_CDS_logPhase, by = "gene")

feature_CDSmfe_100_5p_subset_ribo %>% 
  group_by(logPhase_CDSmfe_100_5p_subset) %>%
  summarize(cor = cor.test(HalfLife_vals_logPhase, ribo_CDSoneThird_5p, method = "spearman")$estimate,
            count = n())

feature_CDSmfe_100_5p_subset_ribo %>% 
  ggplot(aes(x = ribo_CDSoneThird_5p, y = HalfLife_vals_logPhase)) + 
  geom_point(size = 4, color = "#036D9C", 
             alpha = 0.3,
             stroke = 0) +
  xlim(0, 2) +
  facet_wrap(~ logPhase_CDSmfe_100_5p_subset) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.title =  element_text(size = 25, face = "bold"),
        legend.position = "None",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15, face = "bold")
  )

###### CDS_MFE_20_10nt_mid
feature_CDSmfe_20_mid <- featureTable %>% 
  select(gene, CDS_MFE_20_10nt_mid)

summary(feature_CDSmfe_20_mid$CDS_MFE_20_10nt_mid)

feature_CDSmfe_20_mid_subset <- feature_CDSmfe_20_mid %>% 
  mutate(logPhase_CDSmfe_20_mid_subset = case_when(CDS_MFE_20_10nt_mid <= -3.871 ~ "logPhase_CDSmfe_20_mid_G1",
                                                   CDS_MFE_20_10nt_mid > -3.871 & CDS_MFE_20_10nt_mid <= -3.428 ~ "logPhase_CDSmfe_20_mid_G2",
                                                   CDS_MFE_20_10nt_mid > -3.428 & CDS_MFE_20_10nt_mid <= -2.972 ~ "logPhase_CDSmfe_20_mid_G3",
                                                   CDS_MFE_20_10nt_mid > -2.972 ~ "logPhase_CDSmfe_20_mid_G4"))

feature_CDSmfe_20_mid_subset_ribo <- feature_CDSmfe_20_mid_subset %>% 
  inner_join(ribo_oneThird, by = "gene") %>% 
  inner_join(hl_CDS_logPhase, by = "gene")

feature_CDSmfe_20_mid_subset_ribo %>% 
  group_by(logPhase_CDSmfe_20_mid_subset) %>%
  summarize(cor = cor.test(HalfLife_vals_logPhase, ribo_CDSoneThird_mid, method = "spearman")$estimate,
            count = n())

feature_CDSmfe_20_mid_subset_ribo %>% 
  ggplot(aes(x = ribo_CDSoneThird_mid, y = HalfLife_vals_logPhase)) + 
  geom_point(size = 4, color = "#036D9C", 
             alpha = 0.3,
             stroke = 0) +
  xlim(0, 2) +
  facet_wrap(~ logPhase_CDSmfe_20_mid_subset) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.title =  element_text(size = 25, face = "bold"),
        legend.position = "None",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15, face = "bold")
  )

###### CDS_MFE_50_25nt_3p
feature_CDSmfe_50_3p <- featureTable %>% 
  select(gene, CDS_MFE_50_25nt_3p)

summary(feature_CDSmfe_50_3p$CDS_MFE_50_25nt_3p)

feature_CDSmfe_50_3p_subset <- feature_CDSmfe_50_3p %>% 
  mutate(logPhase_CDSmfe_50_3p_subset = case_when(CDS_MFE_50_25nt_3p <= -17.05 ~ "logPhase_CDSmfe_50_3p_G1",
                                                  CDS_MFE_50_25nt_3p > -17.05 & CDS_MFE_50_25nt_3p <= -15.66 ~ "logPhase_CDSmfe_50_3p_G2",
                                                  CDS_MFE_50_25nt_3p > -15.66 & CDS_MFE_50_25nt_3p <= -14.06 ~ "logPhase_CDSmfe_50_3p_G3",
                                                  CDS_MFE_50_25nt_3p > -14.06 ~ "logPhase_CDSmfe_50_3p_G4"))

feature_CDSmfe_50_3p_subset_ribo <- feature_CDSmfe_50_3p_subset %>% 
  inner_join(ribo_oneThird, by = "gene") %>% 
  inner_join(hl_CDS_logPhase, by = "gene")

feature_CDSmfe_50_3p_subset_ribo %>% 
  group_by(logPhase_CDSmfe_50_3p_subset) %>%
  summarize(cor = cor.test(HalfLife_vals_logPhase, ribo_CDSoneThird_3p, method = "spearman")$estimate,
            count = n())

feature_CDSmfe_50_3p_subset_ribo %>% 
  ggplot(aes(x = ribo_CDSoneThird_3p, y = HalfLife_vals_logPhase)) + 
  geom_point(size = 4, color = "#036D9C", 
             alpha = 0.3,
             stroke = 0) +
  xlim(0, 2) +
  facet_wrap(~ logPhase_CDSmfe_50_3p_subset) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.title =  element_text(size = 25, face = "bold"),
        legend.position = "None",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15, face = "bold")
  )

###### CDS_MFE_100_50nt_3p
feature_CDSmfe_100_3p <- featureTable %>% 
  select(gene, CDS_MFE_100_50nt_3p)

summary(feature_CDSmfe_100_3p$CDS_MFE_100_50nt_3p)

feature_CDSmfe_100_3p_subset <- feature_CDSmfe_100_3p %>% 
  mutate(logPhase_CDSmfe_100_3p_subset = case_when(CDS_MFE_100_50nt_3p <= -39.72 ~ "logPhase_CDSmfe_100_3p_G1",
                                                   CDS_MFE_100_50nt_3p > -39.72 & CDS_MFE_100_50nt_3p <= -36.16 ~ "logPhase_CDSmfe_100_3p_G2",
                                                   CDS_MFE_100_50nt_3p > -36.16 & CDS_MFE_100_50nt_3p <= -32.71 ~ "logPhase_CDSmfe_100_3p_G3",
                                                   CDS_MFE_100_50nt_3p > -32.71 ~ "logPhase_CDSmfe_100_3p_G4"))

feature_CDSmfe_100_3p_subset_ribo <- feature_CDSmfe_100_3p_subset %>% 
  inner_join(ribo_oneThird, by = "gene") %>% 
  inner_join(hl_CDS_logPhase, by = "gene")

feature_CDSmfe_100_3p_subset_ribo %>% 
  group_by(logPhase_CDSmfe_100_3p_subset) %>%
  summarize(cor = cor.test(HalfLife_vals_logPhase, ribo_CDSoneThird_3p, method = "spearman")$estimate,
            count = n())

feature_CDSmfe_100_3p_subset_ribo %>% 
  ggplot(aes(x = ribo_CDSoneThird_3p, y = HalfLife_vals_logPhase)) + 
  geom_point(size = 4, color = "#036D9C", 
             alpha = 0.3,
             stroke = 0) +
  xlim(0, 2) +
  facet_wrap(~ logPhase_CDSmfe_100_3p_subset) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.title =  element_text(size = 25, face = "bold"),
        legend.position = "None",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15, face = "bold")
  )

####### correlations between regional ribosome profile and CDS MFE 
###### CDS_MFE_50_25nt_5p
feature_CDSmfe_50_5p <- featureTable %>% 
  select(gene, CDS_MFE_50_25nt_5p)

feature_CDSmfe_50_5p_ribo <- feature_CDSmfe_50_5p %>% 
  inner_join(ribo_oneThird, by = "gene")

feature_CDSmfe_50_5p_ribo %>% 
  ggplot(aes(x = CDS_MFE_50_25nt_5p, y = ribo_CDSoneThird_5p)) + 
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
        axis.text = element_text(size = 20, face = "bold"),
        axis.title =  element_text(size = 25, face = "bold"),
        legend.position = "None")
cor.test(feature_CDSmfe_50_5p_ribo$CDS_MFE_50_25nt_5p, feature_CDSmfe_50_5p_ribo$ribo_CDSoneThird_5p, method = "spearman")

###### CDS_MFE_100_50nt_5p
feature_CDSmfe_100_5p <- featureTable %>% 
  select(gene, CDS_MFE_100_50nt_5p)

feature_CDSmfe_100_5p_ribo <- feature_CDSmfe_100_5p %>% 
  inner_join(ribo_oneThird, by = "gene")

feature_CDSmfe_100_5p_ribo %>% 
  ggplot(aes(x = CDS_MFE_100_50nt_5p, y = ribo_CDSoneThird_5p)) + 
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
        axis.text = element_text(size = 20, face = "bold"),
        axis.title =  element_text(size = 25, face = "bold"),
        legend.position = "None")
cor.test(feature_CDSmfe_100_5p_ribo$CDS_MFE_100_50nt_5p, feature_CDSmfe_100_5p_ribo$ribo_CDSoneThird_5p, method = "spearman")

###### CDS_MFE_20_10nt_mid
feature_CDSmfe_20_mid <- featureTable %>% 
  select(gene, CDS_MFE_20_10nt_mid)

feature_CDSmfe_20_mid_ribo <- feature_CDSmfe_20_mid %>% 
  inner_join(ribo_oneThird, by = "gene")

feature_CDSmfe_20_mid_ribo %>% 
  ggplot(aes(x = CDS_MFE_20_10nt_mid, y = ribo_CDSoneThird_mid)) + 
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
        axis.text = element_text(size = 20, face = "bold"),
        axis.title =  element_text(size = 25, face = "bold"),
        legend.position = "None")
cor.test(feature_CDSmfe_20_mid_ribo$CDS_MFE_20_10nt_mid, feature_CDSmfe_20_mid_ribo$ribo_CDSoneThird_mid, method = "spearman")

###### CDS_MFE_50_25nt_3p
feature_CDSmfe_50_3p <- featureTable %>% 
  select(gene, CDS_MFE_50_25nt_3p)

feature_CDSmfe_50_3p_ribo <- feature_CDSmfe_50_3p %>% 
  inner_join(ribo_oneThird, by = "gene")

feature_CDSmfe_50_3p_ribo %>% 
  ggplot(aes(x = CDS_MFE_50_25nt_3p, y = ribo_CDSoneThird_3p)) + 
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
        axis.text = element_text(size = 20, face = "bold"),
        axis.title =  element_text(size = 25, face = "bold"),
        legend.position = "None")
cor.test(feature_CDSmfe_50_3p_ribo$CDS_MFE_50_25nt_3p, feature_CDSmfe_50_3p_ribo$ribo_CDSoneThird_3p, method = "spearman")

###### CDS_MFE_100_50nt_3p
feature_CDSmfe_100_3p <- featureTable %>% 
  select(gene, CDS_MFE_100_50nt_3p)

feature_CDSmfe_100_3p_ribo <- feature_CDSmfe_100_3p %>% 
  inner_join(ribo_oneThird, by = "gene")

feature_CDSmfe_100_3p_ribo %>% 
  ggplot(aes(x = CDS_MFE_100_50nt_3p, y = ribo_CDSoneThird_3p)) + 
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
        axis.text = element_text(size = 20, face = "bold"),
        axis.title =  element_text(size = 25, face = "bold"),
        legend.position = "None")
cor.test(feature_CDSmfe_100_3p_ribo$CDS_MFE_100_50nt_3p, feature_CDSmfe_100_3p_ribo$ribo_CDSoneThird_3p, method = "spearman")

####### correlations between regional ribosome profile and CDS MFE seperated by half-life class
###### CDS_MFE_50_25nt_5p
feature_CDSmfe_50_5p <- featureTable %>% 
  select(gene, CDS_MFE_50_25nt_5p, HalfLife_cls_logPhase)

feature_CDSmfe_50_5p_ribo <- feature_CDSmfe_50_5p %>% 
  inner_join(ribo_oneThird, by = "gene")

feature_CDSmfe_50_5p_ribo %>% 
  ggplot(aes(x = CDS_MFE_50_25nt_5p, y = ribo_CDSoneThird_5p)) + 
  geom_point(size = 4, fill = "#036D9C",
             color = "white",
             shape = 21) +
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE, colour = "#EB636A") +
  xlim(-24, -8) +
  ylim(0, 3) +
  facet_wrap(~ HalfLife_cls_logPhase, scales = "free") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.text = element_text(size = 15, face = "bold"),
        axis.title =  element_text(size = 20, face = "bold"),
        legend.position = "None",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15, face = "bold"))
ggsave("./viz_CDSmfeCorr_withRegionalRibosome/ribo_CDSoneThird_5p_CDS_MFE_50_25nt_5p_byHalfLifeCls.png", width = 17, height = 12, units = "cm", dpi = 600)

feature_CDSmfe_50_5p_ribo %>% 
  group_by(HalfLife_cls_logPhase) %>%
  summarize(cor = cor.test(CDS_MFE_50_25nt_5p, ribo_CDSoneThird_5p, method = "spearman")$estimate,
            count = n(),
            pval = cor.test(CDS_MFE_50_25nt_5p, ribo_CDSoneThird_5p, method = "spearman")$p.value)

feature_CDSmfe_50_5p_ribo %>% 
  group_by(HalfLife_cls_logPhase) %>%
  group_modify(~ broom::tidy(lm(ribo_CDSoneThird_5p ~ CDS_MFE_50_25nt_5p, data = .)))

feature_CDSmfe_50_5p_ribo_fast <- feature_CDSmfe_50_5p_ribo %>% 
  filter(HalfLife_cls_logPhase == 'Fast')
model_feature_CDSmfe_50_5p_ribo_fast <- lm(ribo_CDSoneThird_5p ~ CDS_MFE_50_25nt_5p, data = feature_CDSmfe_50_5p_ribo_fast)
summary(model_feature_CDSmfe_50_5p_ribo_fast)
summary(model_feature_CDSmfe_50_5p_ribo_fast)$r.squared

feature_CDSmfe_50_5p_ribo_medFast <- feature_CDSmfe_50_5p_ribo %>% 
  filter(HalfLife_cls_logPhase == 'Med-fast')
model_feature_CDSmfe_50_5p_ribo_medFast <- lm(ribo_CDSoneThird_5p ~ CDS_MFE_50_25nt_5p, data = feature_CDSmfe_50_5p_ribo_medFast)
summary(model_feature_CDSmfe_50_5p_ribo_medFast)
summary(model_feature_CDSmfe_50_5p_ribo_medFast)$r.squared

feature_CDSmfe_50_5p_ribo_medSlow <- feature_CDSmfe_50_5p_ribo %>% 
  filter(HalfLife_cls_logPhase == 'Med-slow')
model_feature_CDSmfe_50_5p_ribo_medSlow <- lm(ribo_CDSoneThird_5p ~ CDS_MFE_50_25nt_5p, data = feature_CDSmfe_50_5p_ribo_medSlow)
summary(model_feature_CDSmfe_50_5p_ribo_medSlow)
summary(model_feature_CDSmfe_50_5p_ribo_medSlow)$r.squared

feature_CDSmfe_50_5p_ribo_slow <- feature_CDSmfe_50_5p_ribo %>% 
  filter(HalfLife_cls_logPhase == 'Slow')
model_feature_CDSmfe_50_5p_ribo_slow <- lm(ribo_CDSoneThird_5p ~ CDS_MFE_50_25nt_5p, data = feature_CDSmfe_50_5p_ribo_slow)
summary(model_feature_CDSmfe_50_5p_ribo_slow)
summary(model_feature_CDSmfe_50_5p_ribo_slow)$r.squared

############## test session ############## 
test <- feature_CDSmfe_50_5p_ribo %>% 
  filter(HalfLife_cls_logPhase == "Fast")
model <- lm(test$ribo_CDSoneThird_5p ~ test$CDS_MFE_50_25nt_5p)
summary(model)
############## end of test session ############## 

###### CDS_MFE_100_50nt_5p
feature_CDSmfe_100_5p <- featureTable %>% 
  select(gene, CDS_MFE_100_50nt_5p, HalfLife_cls_logPhase)

feature_CDSmfe_100_5p_ribo <- feature_CDSmfe_100_5p %>% 
  inner_join(ribo_oneThird, by = "gene")

feature_CDSmfe_100_5p_ribo %>% 
  ggplot(aes(x = CDS_MFE_100_50nt_5p, y = ribo_CDSoneThird_5p)) + 
  geom_point(size = 4, fill = "#036D9C",
             color = "white",
             shape = 21) +
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE, colour = "#EB636A") +
  xlim(-55, -25) +
  ylim(0, 3) +
  facet_wrap(~ HalfLife_cls_logPhase, scales = "free") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.text = element_text(size = 15, face = "bold"),
        axis.title =  element_text(size = 20, face = "bold"),
        legend.position = "None",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15, face = "bold"))
ggsave("./viz_CDSmfeCorr_withRegionalRibosome/ribo_CDSoneThird_5p_CDS_MFE_100_50nt_5p_byHalfLifeCls.png", width = 17, height = 12, units = "cm", dpi = 600)

feature_CDSmfe_100_5p_ribo %>% 
  group_by(HalfLife_cls_logPhase) %>%
  summarize(cor = cor.test(CDS_MFE_100_50nt_5p, ribo_CDSoneThird_5p, method = "spearman")$estimate,
            count = n(),
            pval = cor.test(CDS_MFE_100_50nt_5p, ribo_CDSoneThird_5p, method = "spearman")$p.value)

feature_CDSmfe_100_5p_ribo %>% 
  group_by(HalfLife_cls_logPhase) %>%
  group_modify(~ broom::tidy(lm(ribo_CDSoneThird_5p ~ CDS_MFE_100_50nt_5p, data = .)))

###### CDS_MFE_20_10nt_mid
feature_CDSmfe_20_mid <- featureTable %>% 
  select(gene, CDS_MFE_20_10nt_mid, HalfLife_cls_logPhase)

feature_CDSmfe_20_mid_ribo <- feature_CDSmfe_20_mid %>% 
  inner_join(ribo_oneThird, by = "gene")

feature_CDSmfe_20_mid_ribo %>% 
  ggplot(aes(x = CDS_MFE_20_10nt_mid, y = ribo_CDSoneThird_mid)) + 
  geom_point(size = 4, fill = "#036D9C",
             color = "white",
             shape = 21) +
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE, colour = "#EB636A") +
  xlim(-7, -1) +
  ylim(0, 3) +
  facet_wrap(~ HalfLife_cls_logPhase, scales = "free") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.text = element_text(size = 15, face = "bold"),
        axis.title =  element_text(size = 20, face = "bold"),
        legend.position = "None",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15, face = "bold"))
ggsave("./viz_CDSmfeCorr_withRegionalRibosome/ribo_CDSoneThird_mid_CDS_MFE_20_10nt_mid_byHalfLifeCls.png", width = 17, height = 12, units = "cm", dpi = 600)

feature_CDSmfe_20_mid_ribo %>% 
  group_by(HalfLife_cls_logPhase) %>%
  summarize(cor = cor.test(CDS_MFE_20_10nt_mid, ribo_CDSoneThird_mid, method = "spearman")$estimate,
            count = n(),
            pval = cor.test(CDS_MFE_20_10nt_mid, ribo_CDSoneThird_mid, method = "spearman")$p.value)

feature_CDSmfe_20_mid_ribo %>% 
  group_by(HalfLife_cls_logPhase) %>%
  group_modify(~ broom::tidy(lm(ribo_CDSoneThird_mid ~ CDS_MFE_20_10nt_mid, data = .)))

###### CDS_MFE_50_25nt_3p
feature_CDSmfe_50_3p <- featureTable %>% 
  select(gene, CDS_MFE_50_25nt_3p, HalfLife_cls_logPhase)

feature_CDSmfe_50_3p_ribo <- feature_CDSmfe_50_3p %>% 
  inner_join(ribo_oneThird, by = "gene")

feature_CDSmfe_50_3p_ribo %>% 
  ggplot(aes(x = CDS_MFE_50_25nt_3p, y = ribo_CDSoneThird_3p)) + 
  geom_point(size = 4, fill = "#036D9C",
             color = "white",
             shape = 21) +
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE, colour = "#EB636A") +
  xlim(-23, -4) +
  ylim(0, 3) +
  facet_wrap(~ HalfLife_cls_logPhase, scales = "free") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.text = element_text(size = 15, face = "bold"),
        axis.title =  element_text(size = 20, face = "bold"),
        legend.position = "None",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15, face = "bold"))
ggsave("./viz_CDSmfeCorr_withRegionalRibosome/ribo_CDSoneThird_3p_CDS_MFE_50_25nt_3p_byHalfLifeCls.png", width = 17, height = 12, units = "cm", dpi = 600)

feature_CDSmfe_50_3p_ribo %>% 
  group_by(HalfLife_cls_logPhase) %>%
  summarize(cor = cor.test(CDS_MFE_50_25nt_3p, ribo_CDSoneThird_3p, method = "spearman")$estimate,
            count = n(),
            pval = cor.test(CDS_MFE_50_25nt_3p, ribo_CDSoneThird_3p, method = "spearman")$p.value)

feature_CDSmfe_50_3p_ribo %>% 
  group_by(HalfLife_cls_logPhase) %>%
  group_modify(~ broom::tidy(lm(ribo_CDSoneThird_3p ~ CDS_MFE_50_25nt_3p, data = .)))

###### CDS_MFE_100_50nt_3p
feature_CDSmfe_100_3p <- featureTable %>% 
  select(gene, CDS_MFE_100_50nt_3p, HalfLife_cls_logPhase)

feature_CDSmfe_100_3p_ribo <- feature_CDSmfe_100_3p %>% 
  inner_join(ribo_oneThird, by = "gene")

feature_CDSmfe_100_3p_ribo %>% 
  ggplot(aes(x = CDS_MFE_100_50nt_3p, y = ribo_CDSoneThird_3p)) + 
  geom_point(size = 4, fill = "#036D9C",
             color = "white",
             shape = 21) +
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE, colour = "#EB636A") +
  xlim(-51, -7) +
  ylim(0, 3) +
  facet_wrap(~ HalfLife_cls_logPhase, scales = "free") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.text = element_text(size = 15, face = "bold"),
        axis.title =  element_text(size = 20, face = "bold"),
        legend.position = "None",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15, face = "bold"))
ggsave("./viz_CDSmfeCorr_withRegionalRibosome/ribo_CDSoneThird_3p_CDS_MFE_100_50nt_3p_byHalfLifeCls.png", width = 17, height = 12, units = "cm", dpi = 600)

feature_CDSmfe_100_3p_ribo %>% 
  group_by(HalfLife_cls_logPhase) %>%
  summarize(cor = cor.test(CDS_MFE_100_50nt_3p, ribo_CDSoneThird_3p, method = "spearman")$estimate,
            count = n(),
            pval = cor.test(CDS_MFE_100_50nt_3p, ribo_CDSoneThird_3p, method = "spearman")$p.value)

feature_CDSmfe_100_3p_ribo %>% 
  group_by(HalfLife_cls_logPhase) %>%
  group_modify(~ broom::tidy(lm(ribo_CDSoneThird_3p ~ CDS_MFE_100_50nt_3p, data = .)))
