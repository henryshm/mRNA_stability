
####### viz SHAP values of features
####### To compare features of log phase, hypoxia and fold change in hypoxia
####### for leadered and leaderless transcripts 
####### features are the top20 important features for the three conditions
####### feature selected by log phase in combined selected model
####### leaderless transcript
####### hypoxia

library(tidyverse)
library(RColorBrewer)
library(scales)
library(viridis)

####### get color scheme for SHAP values
color_SHAP = c("#900DA4", "#B12A90", "#CC4678", "#E16462", "#F1844B", "#FCA636", "#FCCE25")

####### get SHAP values for each class
SHAP_vals_Fast <- read.csv('../SHAP_values/SHAPvals_HalfLifeCls_hypoxia_combinedSelected_leaderless_byLeaderlessLogPhase_jointlyByLogPhaseHypoxia_Fast.csv') %>% 
  mutate(class = 'Fast')
SHAP_vals_MedFast <- read.csv('../SHAP_values/SHAPvals_HalfLifeCls_hypoxia_combinedSelected_leaderless_byLeaderlessLogPhase_jointlyByLogPhaseHypoxia_MedFast.csv') %>% 
  mutate(class = 'MedFast')
SHAP_vals_MedSlow <- read.csv('../SHAP_values/SHAPvals_HalfLifeCls_hypoxia_combinedSelected_leaderless_byLeaderlessLogPhase_jointlyByLogPhaseHypoxia_MedSlow.csv') %>% 
  mutate(class = 'MedSlow')
SHAP_vals_Slow <- read.csv('../SHAP_values/SHAPvals_HalfLifeCls_hypoxia_combinedSelected_leaderless_byLeaderlessLogPhase_jointlyByLogPhaseHypoxia_Slow.csv') %>% 
  mutate(class = 'Slow')

####### viz of SHAP values
###### Nucl_A_CDS_5p18nt
SHAP_vals_Fast_Nucl_A_CDS_5p18nt <- SHAP_vals_Fast %>% 
  select(Nucl_A_CDS_5p18nt, Nucl_A_CDS_5p18nt_SHAP, class)
SHAP_vals_MedFast_Nucl_A_CDS_5p18nt <- SHAP_vals_MedFast %>% 
  select(Nucl_A_CDS_5p18nt, Nucl_A_CDS_5p18nt_SHAP, class)
SHAP_vals_MedSlow_Nucl_A_CDS_5p18nt <- SHAP_vals_MedSlow %>% 
  select(Nucl_A_CDS_5p18nt, Nucl_A_CDS_5p18nt_SHAP, class)
SHAP_vals_Slow_Nucl_A_CDS_5p18nt <- SHAP_vals_Slow %>% 
  select(Nucl_A_CDS_5p18nt, Nucl_A_CDS_5p18nt_SHAP, class)

SHAP_Nucl_A_CDS_5p18nt <- bind_rows(SHAP_vals_Slow_Nucl_A_CDS_5p18nt, SHAP_vals_MedSlow_Nucl_A_CDS_5p18nt,
                                    SHAP_vals_MedFast_Nucl_A_CDS_5p18nt, SHAP_vals_Fast_Nucl_A_CDS_5p18nt) %>% 
  mutate(class = as.factor(class))

SHAP_Nucl_A_CDS_5p18nt %>% 
  ggplot(aes(class, Nucl_A_CDS_5p18nt_SHAP)) +
  geom_point(aes(color = Nucl_A_CDS_5p18nt), 
             size = 7,
             alpha = 0.7,
             stroke = 0, 
             position = position_jitter(height = 0, width = 0.2, seed = 7)) +
  scale_colour_gradientn(colours = color_SHAP) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey90", size = 1),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank()
  ) +
  coord_flip()
ggsave("./viz_SHAPvalues_toCompareCondition_leaderless/SHAPvalueTop20combined_toCompareCondition_HalfLifeCls_hypoxia_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseHypoxia_leaderless_Nucl_A_CDS_5p18nt.png", width = 55, height = 22, units = "cm", dpi = 600)

###### GC_CDS_5p18nt
SHAP_vals_Fast_GC_CDS_5p18nt <- SHAP_vals_Fast %>% 
  select(GC_CDS_5p18nt, GC_CDS_5p18nt_SHAP, class)
SHAP_vals_MedFast_GC_CDS_5p18nt <- SHAP_vals_MedFast %>% 
  select(GC_CDS_5p18nt, GC_CDS_5p18nt_SHAP, class)
SHAP_vals_MedSlow_GC_CDS_5p18nt <- SHAP_vals_MedSlow %>% 
  select(GC_CDS_5p18nt, GC_CDS_5p18nt_SHAP, class)
SHAP_vals_Slow_GC_CDS_5p18nt <- SHAP_vals_Slow %>% 
  select(GC_CDS_5p18nt, GC_CDS_5p18nt_SHAP, class)

SHAP_GC_CDS_5p18nt <- bind_rows(SHAP_vals_Slow_GC_CDS_5p18nt, SHAP_vals_MedSlow_GC_CDS_5p18nt,
                                SHAP_vals_MedFast_GC_CDS_5p18nt, SHAP_vals_Fast_GC_CDS_5p18nt) %>% 
  mutate(class = as.factor(class))

SHAP_GC_CDS_5p18nt %>% 
  ggplot(aes(class, GC_CDS_5p18nt_SHAP)) +
  geom_point(aes(color = GC_CDS_5p18nt), 
             size = 7,
             alpha = 0.7,
             stroke = 0, 
             position = position_jitter(height = 0, width = 0.2, seed = 7)) +
  scale_colour_gradientn(colours = color_SHAP) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey90", size = 1),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank()
  ) +
  coord_flip()
ggsave("./viz_SHAPvalues_toCompareCondition_leaderless/SHAPvalueTop20combined_toCompareCondition_HalfLifeCls_hypoxia_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseHypoxia_leaderless_GC_CDS_5p18nt.png", width = 55, height = 22, units = "cm", dpi = 600)

###### CDS_MFE_50_25nt_5p
SHAP_vals_Fast_CDS_MFE_50_25nt_5p <- SHAP_vals_Fast %>% 
  select(CDS_MFE_50_25nt_5p, CDS_MFE_50_25nt_5p_SHAP, class)
SHAP_vals_MedFast_CDS_MFE_50_25nt_5p <- SHAP_vals_MedFast %>% 
  select(CDS_MFE_50_25nt_5p, CDS_MFE_50_25nt_5p_SHAP, class)
SHAP_vals_MedSlow_CDS_MFE_50_25nt_5p <- SHAP_vals_MedSlow %>% 
  select(CDS_MFE_50_25nt_5p, CDS_MFE_50_25nt_5p_SHAP, class)
SHAP_vals_Slow_CDS_MFE_50_25nt_5p <- SHAP_vals_Slow %>% 
  select(CDS_MFE_50_25nt_5p, CDS_MFE_50_25nt_5p_SHAP, class)

SHAP_CDS_MFE_50_25nt_5p <- bind_rows(SHAP_vals_Slow_CDS_MFE_50_25nt_5p, SHAP_vals_MedSlow_CDS_MFE_50_25nt_5p,
                                     SHAP_vals_MedFast_CDS_MFE_50_25nt_5p, SHAP_vals_Fast_CDS_MFE_50_25nt_5p) %>% 
  mutate(class = as.factor(class))

SHAP_CDS_MFE_50_25nt_5p %>% 
  ggplot(aes(class, CDS_MFE_50_25nt_5p_SHAP)) +
  geom_point(aes(color = CDS_MFE_50_25nt_5p), 
             size = 7,
             alpha = 0.7,
             stroke = 0, 
             position = position_jitter(height = 0, width = 0.2, seed = 7)) +
  scale_colour_gradientn(colours = color_SHAP) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey90", size = 1),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank()
  ) +
  coord_flip()
ggsave("./viz_SHAPvalues_toCompareCondition_leaderless/SHAPvalueTop20combined_toCompareCondition_HalfLifeCls_hypoxia_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseHypoxia_leaderless_CDS_MFE_50_25nt_5p.png", width = 55, height = 22, units = "cm", dpi = 600)

###### CDS_MFE_20_10nt_mid
SHAP_vals_Fast_CDS_MFE_20_10nt_mid <- SHAP_vals_Fast %>% 
  select(CDS_MFE_20_10nt_mid, CDS_MFE_20_10nt_mid_SHAP, class)
SHAP_vals_MedFast_CDS_MFE_20_10nt_mid <- SHAP_vals_MedFast %>% 
  select(CDS_MFE_20_10nt_mid, CDS_MFE_20_10nt_mid_SHAP, class)
SHAP_vals_MedSlow_CDS_MFE_20_10nt_mid <- SHAP_vals_MedSlow %>% 
  select(CDS_MFE_20_10nt_mid, CDS_MFE_20_10nt_mid_SHAP, class)
SHAP_vals_Slow_CDS_MFE_20_10nt_mid <- SHAP_vals_Slow %>% 
  select(CDS_MFE_20_10nt_mid, CDS_MFE_20_10nt_mid_SHAP, class)

SHAP_CDS_MFE_20_10nt_mid <- bind_rows(SHAP_vals_Slow_CDS_MFE_20_10nt_mid, SHAP_vals_MedSlow_CDS_MFE_20_10nt_mid,
                                      SHAP_vals_MedFast_CDS_MFE_20_10nt_mid, SHAP_vals_Fast_CDS_MFE_20_10nt_mid) %>% 
  mutate(class = as.factor(class))

SHAP_CDS_MFE_20_10nt_mid %>% 
  ggplot(aes(class, CDS_MFE_20_10nt_mid_SHAP)) +
  geom_point(aes(color = CDS_MFE_20_10nt_mid), 
             size = 7,
             alpha = 0.7,
             stroke = 0, 
             position = position_jitter(height = 0, width = 0.2, seed = 7)) +
  scale_colour_gradientn(colours = color_SHAP) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey90", size = 1),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank()
  ) +
  coord_flip()
ggsave("./viz_SHAPvalues_toCompareCondition_leaderless/SHAPvalueTop20combined_toCompareCondition_HalfLifeCls_hypoxia_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseHypoxia_leaderless_CDS_MFE_20_10nt_mid.png", width = 55, height = 22, units = "cm", dpi = 600)

###### CDS_MFE_50_25nt_3p
SHAP_vals_Fast_CDS_MFE_50_25nt_3p <- SHAP_vals_Fast %>% 
  select(CDS_MFE_50_25nt_3p, CDS_MFE_50_25nt_3p_SHAP, class)
SHAP_vals_MedFast_CDS_MFE_50_25nt_3p <- SHAP_vals_MedFast %>% 
  select(CDS_MFE_50_25nt_3p, CDS_MFE_50_25nt_3p_SHAP, class)
SHAP_vals_MedSlow_CDS_MFE_50_25nt_3p <- SHAP_vals_MedSlow %>% 
  select(CDS_MFE_50_25nt_3p, CDS_MFE_50_25nt_3p_SHAP, class)
SHAP_vals_Slow_CDS_MFE_50_25nt_3p <- SHAP_vals_Slow %>% 
  select(CDS_MFE_50_25nt_3p, CDS_MFE_50_25nt_3p_SHAP, class)

SHAP_CDS_MFE_50_25nt_3p <- bind_rows(SHAP_vals_Slow_CDS_MFE_50_25nt_3p, SHAP_vals_MedSlow_CDS_MFE_50_25nt_3p,
                                     SHAP_vals_MedFast_CDS_MFE_50_25nt_3p, SHAP_vals_Fast_CDS_MFE_50_25nt_3p) %>% 
  mutate(class = as.factor(class))

SHAP_CDS_MFE_50_25nt_3p %>% 
  ggplot(aes(class, CDS_MFE_50_25nt_3p_SHAP)) +
  geom_point(aes(color = CDS_MFE_50_25nt_3p), 
             size = 7,
             alpha = 0.7,
             stroke = 0, 
             position = position_jitter(height = 0, width = 0.2, seed = 7)) +
  scale_colour_gradientn(colours = color_SHAP) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey90", size = 1),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank()
  ) +
  coord_flip()
ggsave("./viz_SHAPvalues_toCompareCondition_leaderless/SHAPvalueTop20combined_toCompareCondition_HalfLifeCls_hypoxia_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseHypoxia_leaderless_CDS_MFE_50_25nt_3p.png", width = 55, height = 22, units = "cm", dpi = 600)

###### adja_CG_CDS
SHAP_vals_Fast_adja_CG_CDS <- SHAP_vals_Fast %>% 
  select(adja_CG_CDS, adja_CG_CDS_SHAP, class)
SHAP_vals_MedFast_adja_CG_CDS <- SHAP_vals_MedFast %>% 
  select(adja_CG_CDS, adja_CG_CDS_SHAP, class)
SHAP_vals_MedSlow_adja_CG_CDS <- SHAP_vals_MedSlow %>% 
  select(adja_CG_CDS, adja_CG_CDS_SHAP, class)
SHAP_vals_Slow_adja_CG_CDS <- SHAP_vals_Slow %>% 
  select(adja_CG_CDS, adja_CG_CDS_SHAP, class)

SHAP_adja_CG_CDS <- bind_rows(SHAP_vals_Slow_adja_CG_CDS, SHAP_vals_MedSlow_adja_CG_CDS,
                              SHAP_vals_MedFast_adja_CG_CDS, SHAP_vals_Fast_adja_CG_CDS) %>% 
  mutate(class = as.factor(class))

SHAP_adja_CG_CDS %>% 
  ggplot(aes(class, adja_CG_CDS_SHAP)) +
  geom_point(aes(color = adja_CG_CDS), 
             size = 7,
             alpha = 0.7,
             stroke = 0, 
             position = position_jitter(height = 0, width = 0.2, seed = 7)) +
  scale_colour_gradientn(colours = color_SHAP) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey90", size = 1),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank()
  ) +
  coord_flip()
ggsave("./viz_SHAPvalues_toCompareCondition_leaderless/SHAPvalueTop20combined_toCompareCondition_HalfLifeCls_hypoxia_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseHypoxia_leaderless_adja_CG_CDS.png", width = 55, height = 22, units = "cm", dpi = 600)

###### adja_UA_CDS
SHAP_vals_Fast_adja_UA_CDS <- SHAP_vals_Fast %>% 
  select(adja_UA_CDS, adja_UA_CDS_SHAP, class)
SHAP_vals_MedFast_adja_UA_CDS <- SHAP_vals_MedFast %>% 
  select(adja_UA_CDS, adja_UA_CDS_SHAP, class)
SHAP_vals_MedSlow_adja_UA_CDS <- SHAP_vals_MedSlow %>% 
  select(adja_UA_CDS, adja_UA_CDS_SHAP, class)
SHAP_vals_Slow_adja_UA_CDS <- SHAP_vals_Slow %>% 
  select(adja_UA_CDS, adja_UA_CDS_SHAP, class)

SHAP_adja_UA_CDS <- bind_rows(SHAP_vals_Slow_adja_UA_CDS, SHAP_vals_MedSlow_adja_UA_CDS,
                              SHAP_vals_MedFast_adja_UA_CDS, SHAP_vals_Fast_adja_UA_CDS) %>% 
  mutate(class = as.factor(class))

SHAP_adja_UA_CDS %>% 
  ggplot(aes(class, adja_UA_CDS_SHAP)) +
  geom_point(aes(color = adja_UA_CDS), 
             size = 7,
             alpha = 0.7,
             stroke = 0, 
             position = position_jitter(height = 0, width = 0.2, seed = 7)) +
  scale_colour_gradientn(colours = color_SHAP) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey90", size = 1),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank()
  ) +
  coord_flip()
ggsave("./viz_SHAPvalues_toCompareCondition_leaderless/SHAPvalueTop20combined_toCompareCondition_HalfLifeCls_hypoxia_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseHypoxia_leaderless_adja_UA_CDS.png", width = 55, height = 22, units = "cm", dpi = 600)

###### initialAbundance_hypoxia
SHAP_vals_Fast_initialAbundance_hypoxia <- SHAP_vals_Fast %>% 
  select(initialAbundance_hypoxia, initialAbundance_hypoxia_SHAP, class)
SHAP_vals_MedFast_initialAbundance_hypoxia <- SHAP_vals_MedFast %>% 
  select(initialAbundance_hypoxia, initialAbundance_hypoxia_SHAP, class)
SHAP_vals_MedSlow_initialAbundance_hypoxia <- SHAP_vals_MedSlow %>% 
  select(initialAbundance_hypoxia, initialAbundance_hypoxia_SHAP, class)
SHAP_vals_Slow_initialAbundance_hypoxia <- SHAP_vals_Slow %>% 
  select(initialAbundance_hypoxia, initialAbundance_hypoxia_SHAP, class)

SHAP_initialAbundance_hypoxia <- bind_rows(SHAP_vals_Slow_initialAbundance_hypoxia, SHAP_vals_MedSlow_initialAbundance_hypoxia,
                                            SHAP_vals_MedFast_initialAbundance_hypoxia, SHAP_vals_Fast_initialAbundance_hypoxia) %>% 
  mutate(class = as.factor(class))

SHAP_initialAbundance_hypoxia %>% 
  ggplot(aes(class, initialAbundance_hypoxia_SHAP)) +
  geom_point(aes(color = initialAbundance_hypoxia), 
             size = 7,
             alpha = 0.7,
             stroke = 0, 
             position = position_jitter(height = 0, width = 0.2, seed = 7)) +
  scale_colour_gradientn(colours = color_SHAP) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey90", size = 1),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank()
  ) +
  coord_flip()
ggsave("./viz_SHAPvalues_toCompareCondition_leaderless/SHAPvalueTop20combined_toCompareCondition_HalfLifeCls_hypoxia_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseHypoxia_leaderless_initialAbundance_hypoxia.png", width = 55, height = 22, units = "cm", dpi = 600)

###### CDS_length
SHAP_vals_Fast_CDS_length <- SHAP_vals_Fast %>% 
  select(CDS_length, CDS_length_SHAP, class)
SHAP_vals_MedFast_CDS_length <- SHAP_vals_MedFast %>% 
  select(CDS_length, CDS_length_SHAP, class)
SHAP_vals_MedSlow_CDS_length <- SHAP_vals_MedSlow %>% 
  select(CDS_length, CDS_length_SHAP, class)
SHAP_vals_Slow_CDS_length <- SHAP_vals_Slow %>% 
  select(CDS_length, CDS_length_SHAP, class)

SHAP_CDS_length <- bind_rows(SHAP_vals_Slow_CDS_length, SHAP_vals_MedSlow_CDS_length,
                             SHAP_vals_MedFast_CDS_length, SHAP_vals_Fast_CDS_length) %>% 
  mutate(class = as.factor(class))

SHAP_CDS_length %>% 
  ggplot(aes(class, CDS_length_SHAP)) +
  geom_point(aes(color = CDS_length), 
             size = 7,
             alpha = 0.7,
             stroke = 0, 
             position = position_jitter(height = 0, width = 0.2, seed = 7)) +
  scale_colour_gradientn(colours = color_SHAP) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey90", size = 1),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.ticks = element_blank(), axis.title = element_blank()
  ) +
  coord_flip()
ggsave("./viz_SHAPvalues_toCompareCondition_leaderless/SHAPvalueTop20combined_toCompareCondition_HalfLifeCls_hypoxia_combinedSelectedByLeaderlessLogPhase_jointlyByLogPhaseHypoxia_leaderless_CDS_length.png", width = 55, height = 22, units = "cm", dpi = 600)
