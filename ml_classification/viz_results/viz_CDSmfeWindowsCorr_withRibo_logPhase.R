
####### check correlation between CDS MFE with ribosome profiling in leaderless models
####### observed stronger/linear in logPhase, but not log phase slow class
####### to see if the slow half-life class in log phase is protected by ribosome binding
####### and therefore related with translation
####### using features in the top20 important features selected by log phase in combined selected model
####### ribosome profile using ribo_5p_excl
####### leaderless transcript
####### log phase

library(tidyverse)

####### get feature table
featureTable <- read.csv('../../feature/FeatureTables/featureTable_combinedSelected_halfLifeClass/HalfLifeCls_logPhase_combinedSelected_leaderless_byLeaderlessLogPhase_jointlyByLogPhaseHypoxia.csv')

####### get half-life of entire gene
hl_CDS_complete_class <- read.table('../../degradation_classes/half_life/CDS_halfLife_completeClass.txt', header = TRUE)

hl_CDS_logPhase <- hl_CDS_complete_class %>% 
  select(gene, HalfLife_vals_logPhase) %>% 
  drop_na()

####### get ribosome profile
###### ribo_5p_excl
feature_ribo <- featureTable %>% 
  select(gene, ribo_5p_excl)

####### correlations between window MFE and ribo_5p_excl
####### get CDS MFE and split into subgroups by MFE values
###### CDS_MFE_50_25nt_5p
feature_CDSmfe_50_5p <- featureTable %>% 
  select(gene, CDS_MFE_50_25nt_5p)

summary(feature_CDSmfe_50_5p$CDS_MFE_50_25nt_5p)

feature_CDSmfe_50_5p_subset <- feature_CDSmfe_50_5p %>% 
  mutate(logPhase_CDSmfe_50_5p_subset = case_when(CDS_MFE_50_25nt_5p <= -17.38 ~ "logPhase_CDSmfe_50_5p_G1",
                                                  CDS_MFE_50_25nt_5p > -17.38 & CDS_MFE_50_25nt_5p <= -16.09 ~ "logPhase_CDSmfe_50_5p_G2",
                                                  CDS_MFE_50_25nt_5p > -16.09 & CDS_MFE_50_25nt_5p <= -14.74 ~ "logPhase_CDSmfe_50_5p_G3",
                                                  CDS_MFE_50_25nt_5p > -14.74 ~ "logPhase_CDSmfe_50_5p_G4"))

feature_CDSmfe_50_5p_subset_ribo <- feature_CDSmfe_50_5p_subset %>% 
  inner_join(feature_ribo, by = "gene")

feature_CDSmfe_50_5p_subset_ribo %>% 
  ggplot(aes(x = CDS_MFE_50_25nt_5p, y = ribo_5p_excl)) + 
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
ggsave("./viz_CDSmfeCorr_withRibosome/ribo_5pExcl_CDS_MFE_50_25nt_5p.png", width = 15, height = 12, units = "cm", dpi = 600)
cor.test(feature_CDSmfe_50_5p_subset_ribo$CDS_MFE_50_25nt_5p, feature_CDSmfe_50_5p_subset_ribo$ribo_5p_excl, method = "spearman")

feature_CDSmfe_50_5p_subset_ribo %>% 
  group_by(logPhase_CDSmfe_50_5p_subset) %>%
  summarize(cor = cor.test(CDS_MFE_50_25nt_5p, ribo_5p_excl, method = "spearman")$estimate,
            count = n())

feature_CDSmfe_50_5p_subset_ribo %>% 
  ggplot(aes(x = ribo_5p_excl, y = CDS_MFE_50_25nt_5p)) + 
  geom_point(size = 4, color = "#036D9C", 
             alpha = 0.3,
             stroke = 0) +
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
  inner_join(feature_ribo, by = "gene")

feature_CDSmfe_100_5p_subset_ribo %>% 
  ggplot(aes(x = CDS_MFE_100_50nt_5p, y = ribo_5p_excl)) + 
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
ggsave("./viz_CDSmfeCorr_withRibosome/ribo_5pExcl_CDS_MFE_100_50nt_5p.png", width = 15, height = 12, units = "cm", dpi = 600)
cor.test(feature_CDSmfe_100_5p_subset_ribo$CDS_MFE_100_50nt_5p, feature_CDSmfe_100_5p_subset_ribo$ribo_5p_excl, method = "spearman")

feature_CDSmfe_100_5p_subset_ribo %>% 
  group_by(logPhase_CDSmfe_100_5p_subset) %>%
  summarize(cor = cor.test(CDS_MFE_100_50nt_5p, ribo_5p_excl, method = "spearman")$estimate,
            count = n())

feature_CDSmfe_100_5p_subset_ribo %>% 
  ggplot(aes(x = ribo_5p_excl, y = CDS_MFE_100_50nt_5p)) + 
  geom_point(size = 4, color = "#036D9C", 
             alpha = 0.3,
             stroke = 0) +
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
  inner_join(feature_ribo, by = "gene")

feature_CDSmfe_20_mid_subset_ribo %>% 
  ggplot(aes(x = CDS_MFE_20_10nt_mid, y = ribo_5p_excl)) + 
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
ggsave("./viz_CDSmfeCorr_withRibosome/ribo_5pExcl_CDS_MFE_20_10nt_mid.png", width = 15, height = 12, units = "cm", dpi = 600)
cor.test(feature_CDSmfe_20_mid_subset_ribo$CDS_MFE_20_10nt_mid, feature_CDSmfe_20_mid_subset_ribo$ribo_5p_excl, method = "spearman")

feature_CDSmfe_20_mid_subset_ribo %>% 
  group_by(logPhase_CDSmfe_20_mid_subset) %>%
  summarize(cor = cor.test(CDS_MFE_20_10nt_mid, ribo_5p_excl, method = "spearman")$estimate,
            count = n())

feature_CDSmfe_20_mid_subset_ribo %>% 
  ggplot(aes(x = ribo_5p_excl, y = CDS_MFE_20_10nt_mid)) + 
  geom_point(size = 4, color = "#036D9C", 
             alpha = 0.3,
             stroke = 0) +
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
  inner_join(feature_ribo, by = "gene")

feature_CDSmfe_50_3p_subset_ribo %>% 
  ggplot(aes(x = CDS_MFE_50_25nt_3p, y = ribo_5p_excl)) + 
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
ggsave("./viz_CDSmfeCorr_withRibosome/ribo_5pExcl_CDS_MFE_50_25nt_3p.png", width = 15, height = 12, units = "cm", dpi = 600)
cor.test(feature_CDSmfe_50_3p_subset_ribo$CDS_MFE_50_25nt_3p, feature_CDSmfe_50_3p_subset_ribo$ribo_5p_excl, method = "spearman")

feature_CDSmfe_50_3p_subset_ribo %>% 
  group_by(logPhase_CDSmfe_50_3p_subset) %>%
  summarize(cor = cor.test(CDS_MFE_50_25nt_3p, ribo_5p_excl, method = "spearman")$estimate,
            count = n())

feature_CDSmfe_50_3p_subset_ribo %>% 
  ggplot(aes(x = ribo_5p_excl, y = CDS_MFE_50_25nt_3p)) + 
  geom_point(size = 4, color = "#036D9C", 
             alpha = 0.3,
             stroke = 0) +
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
  inner_join(feature_ribo, by = "gene")

feature_CDSmfe_100_3p_subset_ribo %>% 
  ggplot(aes(x = CDS_MFE_100_50nt_3p, y = ribo_5p_excl)) + 
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
ggsave("./viz_CDSmfeCorr_withRibosome/ribo_5pExcl_CDS_MFE_100_50nt_3p.png", width = 15, height = 12, units = "cm", dpi = 600)
cor.test(feature_CDSmfe_100_3p_subset_ribo$CDS_MFE_100_50nt_3p, feature_CDSmfe_100_3p_subset_ribo$ribo_5p_excl, method = "spearman")

feature_CDSmfe_100_3p_subset_ribo %>% 
  group_by(logPhase_CDSmfe_100_3p_subset) %>%
  summarize(cor = cor.test(CDS_MFE_100_50nt_3p, ribo_5p_excl, method = "spearman")$estimate,
            count = n())

feature_CDSmfe_100_3p_subset_ribo %>% 
  ggplot(aes(x = ribo_5p_excl, y = CDS_MFE_100_50nt_3p)) + 
  geom_point(size = 4, color = "#036D9C", 
             alpha = 0.3,
             stroke = 0) +
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

####### correlations between half-lives and ribo_5p_excl within each MFE window defined as above
####### get CDS MFE and split into subgroups by MFE values
###### CDS_MFE_50_25nt_5p
feature_CDSmfe_50_5p <- featureTable %>% 
  select(gene, CDS_MFE_50_25nt_5p)

summary(feature_CDSmfe_50_5p$CDS_MFE_50_25nt_5p)

feature_CDSmfe_50_5p_subset <- feature_CDSmfe_50_5p %>% 
  mutate(logPhase_CDSmfe_50_5p_subset = case_when(CDS_MFE_50_25nt_5p <= -17.38 ~ "logPhase_CDSmfe_50_5p_G1",
                                                  CDS_MFE_50_25nt_5p > -17.38 & CDS_MFE_50_25nt_5p <= -16.09 ~ "logPhase_CDSmfe_50_5p_G2",
                                                  CDS_MFE_50_25nt_5p > -16.09 & CDS_MFE_50_25nt_5p <= -14.74 ~ "logPhase_CDSmfe_50_5p_G3",
                                                  CDS_MFE_50_25nt_5p > -14.74 ~ "logPhase_CDSmfe_50_5p_G4"))

feature_CDSmfe_50_5p_subset_ribo <- feature_CDSmfe_50_5p_subset %>% 
  inner_join(feature_ribo, by = "gene") %>% 
  inner_join(hl_CDS_logPhase, by = "gene")

feature_CDSmfe_50_5p_subset_ribo %>% 
  group_by(logPhase_CDSmfe_50_5p_subset) %>%
  summarize(cor = cor.test(HalfLife_vals_logPhase, ribo_5p_excl, method = "spearman")$estimate,
            count = n())

feature_CDSmfe_50_5p_subset_ribo %>% 
  ggplot(aes(x = ribo_5p_excl, y = HalfLife_vals_logPhase)) + 
  geom_point(size = 4, color = "#036D9C", 
             alpha = 0.3,
             stroke = 0) +
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
  inner_join(feature_ribo, by = "gene") %>% 
  inner_join(hl_CDS_logPhase, by = "gene")

feature_CDSmfe_100_5p_subset_ribo %>% 
  group_by(logPhase_CDSmfe_100_5p_subset) %>%
  summarize(cor = cor.test(HalfLife_vals_logPhase, ribo_5p_excl, method = "spearman")$estimate,
            count = n())

feature_CDSmfe_100_5p_subset_ribo %>% 
  ggplot(aes(x = ribo_5p_excl, y = HalfLife_vals_logPhase)) + 
  geom_point(size = 4, color = "#036D9C", 
             alpha = 0.3,
             stroke = 0) +
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
  inner_join(feature_ribo, by = "gene") %>% 
  inner_join(hl_CDS_logPhase, by = "gene")

feature_CDSmfe_20_mid_subset_ribo %>% 
  group_by(logPhase_CDSmfe_20_mid_subset) %>%
  summarize(cor = cor.test(HalfLife_vals_logPhase, ribo_5p_excl, method = "spearman")$estimate,
            count = n())

feature_CDSmfe_20_mid_subset_ribo %>% 
  ggplot(aes(x = ribo_5p_excl, y = HalfLife_vals_logPhase)) + 
  geom_point(size = 4, color = "#036D9C", 
             alpha = 0.3,
             stroke = 0) +
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
  inner_join(feature_ribo, by = "gene") %>% 
  inner_join(hl_CDS_logPhase, by = "gene")

feature_CDSmfe_50_3p_subset_ribo %>% 
  group_by(logPhase_CDSmfe_50_3p_subset) %>%
  summarize(cor = cor.test(HalfLife_vals_logPhase, ribo_5p_excl, method = "spearman")$estimate,
            count = n())

feature_CDSmfe_50_3p_subset_ribo %>% 
  ggplot(aes(x = ribo_5p_excl, y = HalfLife_vals_logPhase)) + 
  geom_point(size = 4, color = "#036D9C", 
             alpha = 0.3,
             stroke = 0) +
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
feature_CDSmfe_100_3p = "CDS_MFE_100_50nt_3p"

feature_CDSmfe_100_3p <- featureTable %>% 
  select(gene, CDS_MFE_100_50nt_3p)

summary(feature_CDSmfe_100_3p$CDS_MFE_100_50nt_3p)

feature_CDSmfe_100_3p_subset <- feature_CDSmfe_100_3p %>% 
  mutate(logPhase_CDSmfe_100_3p_subset = case_when(CDS_MFE_100_50nt_3p <= -39.72 ~ "logPhase_CDSmfe_100_3p_G1",
                                                   CDS_MFE_100_50nt_3p > -39.72 & CDS_MFE_100_50nt_3p <= -36.16 ~ "logPhase_CDSmfe_100_3p_G2",
                                                   CDS_MFE_100_50nt_3p > -36.16 & CDS_MFE_100_50nt_3p <= -32.71 ~ "logPhase_CDSmfe_100_3p_G3",
                                                   CDS_MFE_100_50nt_3p > -32.71 ~ "logPhase_CDSmfe_100_3p_G4"))

feature_CDSmfe_100_3p_subset_ribo <- feature_CDSmfe_100_3p_subset %>% 
  inner_join(feature_ribo, by = "gene") %>% 
  inner_join(hl_CDS_logPhase, by = "gene")

feature_CDSmfe_100_3p_subset_ribo %>% 
  group_by(logPhase_CDSmfe_100_3p_subset) %>%
  summarize(cor = cor.test(HalfLife_vals_logPhase, ribo_5p_excl, method = "spearman")$estimate,
            count = n())

feature_CDSmfe_100_3p_subset_ribo %>% 
  ggplot(aes(x = ribo_5p_excl, y = HalfLife_vals_logPhase)) + 
  geom_point(size = 4, color = "#036D9C", 
             alpha = 0.3,
             stroke = 0) +
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
