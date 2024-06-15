
####### viz SHAP values of features
####### To compare features of leadered and leaderless transcripts 
####### in log phase, hypoxia and fold change in hypoxia
####### features are the top20 important features for the two transcript types
####### feature selected by leaderless in combined selected model
####### log phase
####### leadered transcript

library(tidyverse)
library(RColorBrewer)
library(scales)
library(viridis)

####### get color scheme for SHAP values
color_SHAP = c("#900DA4", "#B12A90", "#CC4678", "#E16462", "#F1844B", "#FCA636", "#FCCE25")

####### get SHAP values for each class
SHAP_vals_Fast <- read.csv('../SHAP_values/SHAPvals_HalfLifeCls_logPhase_combinedSelected_leadered_byLeaderlessLogPhase_jointlyByLeaderedLeaderless_Fast.csv') %>% 
  mutate(class = 'Fast')
SHAP_vals_MedFast <- read.csv('../SHAP_values/SHAPvals_HalfLifeCls_logPhase_combinedSelected_leadered_byLeaderlessLogPhase_jointlyByLeaderedLeaderless_MedFast.csv') %>% 
  mutate(class = 'MedFast')
SHAP_vals_MedSlow <- read.csv('../SHAP_values/SHAPvals_HalfLifeCls_logPhase_combinedSelected_leadered_byLeaderlessLogPhase_jointlyByLeaderedLeaderless_MedSlow.csv') %>% 
  mutate(class = 'MedSlow')
SHAP_vals_Slow <- read.csv('../SHAP_values/SHAPvals_HalfLifeCls_logPhase_combinedSelected_leadered_byLeaderlessLogPhase_jointlyByLeaderedLeaderless_Slow.csv') %>% 
  mutate(class = 'Slow')

####### viz of SHAP values
###### CGC_CDSnstop
SHAP_vals_Fast_CGC_CDSnstop <- SHAP_vals_Fast %>% 
  select(CGC_CDSnstop, CGC_CDSnstop_SHAP, class)
SHAP_vals_MedFast_CGC_CDSnstop <- SHAP_vals_MedFast %>% 
  select(CGC_CDSnstop, CGC_CDSnstop_SHAP, class)
SHAP_vals_MedSlow_CGC_CDSnstop <- SHAP_vals_MedSlow %>% 
  select(CGC_CDSnstop, CGC_CDSnstop_SHAP, class)
SHAP_vals_Slow_CGC_CDSnstop <- SHAP_vals_Slow %>% 
  select(CGC_CDSnstop, CGC_CDSnstop_SHAP, class)

SHAP_CGC_CDSnstop <- bind_rows(SHAP_vals_Slow_CGC_CDSnstop, SHAP_vals_MedSlow_CGC_CDSnstop, 
                               SHAP_vals_MedFast_CGC_CDSnstop, SHAP_vals_Fast_CGC_CDSnstop) %>% 
  mutate(class = as.factor(class))

SHAP_CGC_CDSnstop %>% 
  ggplot(aes(class, CGC_CDSnstop_SHAP)) +
  geom_point(aes(color = CGC_CDSnstop), 
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
ggsave("./viz_SHAPvalues_toCompareTranscript_logPhase/SHAPvalueTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_CGC_CDSnstop.png", width = 55, height = 22, units = "cm", dpi = 600)

###### CGG_CDSnstop
SHAP_vals_Fast_CGG_CDSnstop <- SHAP_vals_Fast %>% 
  select(CGG_CDSnstop, CGG_CDSnstop_SHAP, class)
SHAP_vals_MedFast_CGG_CDSnstop <- SHAP_vals_MedFast %>% 
  select(CGG_CDSnstop, CGG_CDSnstop_SHAP, class)
SHAP_vals_MedSlow_CGG_CDSnstop <- SHAP_vals_MedSlow %>% 
  select(CGG_CDSnstop, CGG_CDSnstop_SHAP, class)
SHAP_vals_Slow_CGG_CDSnstop <- SHAP_vals_Slow %>% 
  select(CGG_CDSnstop, CGG_CDSnstop_SHAP, class)

SHAP_CGG_CDSnstop <- bind_rows(SHAP_vals_Slow_CGG_CDSnstop, SHAP_vals_MedSlow_CGG_CDSnstop, 
                               SHAP_vals_MedFast_CGG_CDSnstop, SHAP_vals_Fast_CGG_CDSnstop) %>% 
  mutate(class = as.factor(class))

SHAP_CGG_CDSnstop %>% 
  ggplot(aes(class, CGG_CDSnstop_SHAP)) +
  geom_point(aes(color = CGG_CDSnstop), 
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
ggsave("./viz_SHAPvalues_toCompareTranscript_logPhase/SHAPvalueTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_CGG_CDSnstop.png", width = 55, height = 22, units = "cm", dpi = 600)

###### UUG_CDSnstop
SHAP_vals_Fast_UUG_CDSnstop <- SHAP_vals_Fast %>% 
  select(UUG_CDSnstop, UUG_CDSnstop_SHAP, class)
SHAP_vals_MedFast_UUG_CDSnstop <- SHAP_vals_MedFast %>% 
  select(UUG_CDSnstop, UUG_CDSnstop_SHAP, class)
SHAP_vals_MedSlow_UUG_CDSnstop <- SHAP_vals_MedSlow %>% 
  select(UUG_CDSnstop, UUG_CDSnstop_SHAP, class)
SHAP_vals_Slow_UUG_CDSnstop <- SHAP_vals_Slow %>% 
  select(UUG_CDSnstop, UUG_CDSnstop_SHAP, class)

SHAP_UUG_CDSnstop <- bind_rows(SHAP_vals_Slow_UUG_CDSnstop, SHAP_vals_MedSlow_UUG_CDSnstop, 
                               SHAP_vals_MedFast_UUG_CDSnstop, SHAP_vals_Fast_UUG_CDSnstop) %>% 
  mutate(class = as.factor(class))

SHAP_UUG_CDSnstop %>% 
  ggplot(aes(class, UUG_CDSnstop_SHAP)) +
  geom_point(aes(color = UUG_CDSnstop), 
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
ggsave("./viz_SHAPvalues_toCompareTranscript_logPhase/SHAPvalueTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_UUG_CDSnstop.png", width = 55, height = 22, units = "cm", dpi = 600)

###### CGU_CDSnstop
SHAP_vals_Fast_CGU_CDSnstop <- SHAP_vals_Fast %>% 
  select(CGU_CDSnstop, CGU_CDSnstop_SHAP, class)
SHAP_vals_MedFast_CGU_CDSnstop <- SHAP_vals_MedFast %>% 
  select(CGU_CDSnstop, CGU_CDSnstop_SHAP, class)
SHAP_vals_MedSlow_CGU_CDSnstop <- SHAP_vals_MedSlow %>% 
  select(CGU_CDSnstop, CGU_CDSnstop_SHAP, class)
SHAP_vals_Slow_CGU_CDSnstop <- SHAP_vals_Slow %>% 
  select(CGU_CDSnstop, CGU_CDSnstop_SHAP, class)

SHAP_CGU_CDSnstop <- bind_rows(SHAP_vals_Slow_CGU_CDSnstop, SHAP_vals_MedSlow_CGU_CDSnstop, 
                               SHAP_vals_MedFast_CGU_CDSnstop, SHAP_vals_Fast_CGU_CDSnstop) %>% 
  mutate(class = as.factor(class))

SHAP_CGU_CDSnstop %>% 
  ggplot(aes(class, CGU_CDSnstop_SHAP)) +
  geom_point(aes(color = CGU_CDSnstop), 
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
ggsave("./viz_SHAPvalues_toCompareTranscript_logPhase/SHAPvalueTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_CGU_CDSnstop.png", width = 55, height = 22, units = "cm", dpi = 600)

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
ggsave("./viz_SHAPvalues_toCompareTranscript_logPhase/SHAPvalueTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_CDS_MFE_20_10nt_mid.png", width = 55, height = 22, units = "cm", dpi = 600)

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
ggsave("./viz_SHAPvalues_toCompareTranscript_logPhase/SHAPvalueTop20combined_toCompareTranscript_HalfLifeCls_logPhase_combinedSelectedByLeaderlessLogPhase_jointlyByLeaderedLeaderless_leadered_adja_CG_CDS.png", width = 55, height = 22, units = "cm", dpi = 600)
