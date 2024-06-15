
####### viz SHAP values of features
####### To compare features of log phase, hypoxia and fold change in hypoxia
####### for leadered and leaderless transcripts 
####### features are the top20 important features for the three conditions
####### feature selected by log phase in combined selected model
####### leadered transcript
####### hypoxia

library(tidyverse)
library(RColorBrewer)
library(scales)
library(viridis)

####### get color scheme for SHAP values
color_SHAP = c("#900DA4", "#B12A90", "#CC4678", "#E16462", "#F1844B", "#FCA636", "#FCCE25")

####### get SHAP values for each class
SHAP_vals_Fast <- read.csv('../SHAP_values/SHAPvals_HalfLifeCls_hypoxia_combinedSelected_leadered_byLeaderedLogPhase_jointlyByLogPhaseHypoxia_Fast.csv') %>% 
  mutate(class = 'Fast')
SHAP_vals_MedFast <- read.csv('../SHAP_values/SHAPvals_HalfLifeCls_hypoxia_combinedSelected_leadered_byLeaderedLogPhase_jointlyByLogPhaseHypoxia_MedFast.csv') %>% 
  mutate(class = 'MedFast')
SHAP_vals_MedSlow <- read.csv('../SHAP_values/SHAPvals_HalfLifeCls_hypoxia_combinedSelected_leadered_byLeaderedLogPhase_jointlyByLogPhaseHypoxia_MedSlow.csv') %>% 
  mutate(class = 'MedSlow')
SHAP_vals_Slow <- read.csv('../SHAP_values/SHAPvals_HalfLifeCls_hypoxia_combinedSelected_leadered_byLeaderedLogPhase_jointlyByLogPhaseHypoxia_Slow.csv') %>% 
  mutate(class = 'Slow')

####### viz of SHAP values
###### MFE_unfold
SHAP_vals_Fast_MFE_unfold <- SHAP_vals_Fast %>% 
  select(MFE_unfold, MFE_unfold_SHAP, class)
SHAP_vals_MedFast_MFE_unfold <- SHAP_vals_MedFast %>% 
  select(MFE_unfold, MFE_unfold_SHAP, class)
SHAP_vals_MedSlow_MFE_unfold <- SHAP_vals_MedSlow %>% 
  select(MFE_unfold, MFE_unfold_SHAP, class)
SHAP_vals_Slow_MFE_unfold <- SHAP_vals_Slow %>% 
  select(MFE_unfold, MFE_unfold_SHAP, class)

SHAP_MFE_unfold <- bind_rows(SHAP_vals_Slow_MFE_unfold, SHAP_vals_MedSlow_MFE_unfold,
                             SHAP_vals_MedFast_MFE_unfold, SHAP_vals_Fast_MFE_unfold) %>% 
  mutate(class = as.factor(class))

SHAP_MFE_unfold %>% 
  ggplot(aes(class, MFE_unfold_SHAP)) +
  geom_point(aes(color = MFE_unfold), 
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
ggsave("./viz_SHAPvalues_toCompareCondition_leadered/SHAPvalueTop20combined_toCompareCondition_HalfLifeCls_hypoxia_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseHypoxia_leadered_MFE_unfold.png", width = 55, height = 22, units = "cm", dpi = 600)

###### AAA_CDSnstop
SHAP_vals_Fast_AAA_CDSnstop <- SHAP_vals_Fast %>% 
  select(AAA_CDSnstop, AAA_CDSnstop_SHAP, class)
SHAP_vals_MedFast_AAA_CDSnstop <- SHAP_vals_MedFast %>% 
  select(AAA_CDSnstop, AAA_CDSnstop_SHAP, class)
SHAP_vals_MedSlow_AAA_CDSnstop <- SHAP_vals_MedSlow %>% 
  select(AAA_CDSnstop, AAA_CDSnstop_SHAP, class)
SHAP_vals_Slow_AAA_CDSnstop <- SHAP_vals_Slow %>% 
  select(AAA_CDSnstop, AAA_CDSnstop_SHAP, class)

SHAP_AAA_CDSnstop <- bind_rows(SHAP_vals_Slow_AAA_CDSnstop, SHAP_vals_MedSlow_AAA_CDSnstop,
                               SHAP_vals_MedFast_AAA_CDSnstop, SHAP_vals_Fast_AAA_CDSnstop) %>% 
  mutate(class = as.factor(class))

SHAP_AAA_CDSnstop %>% 
  ggplot(aes(class, AAA_CDSnstop_SHAP)) +
  geom_point(aes(color = AAA_CDSnstop), 
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
ggsave("./viz_SHAPvalues_toCompareCondition_leadered/SHAPvalueTop20combined_toCompareCondition_HalfLifeCls_hypoxia_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseHypoxia_leadered_AAA_CDSnstop.png", width = 55, height = 22, units = "cm", dpi = 600)

###### ACG_CDSnstop
SHAP_vals_Fast_ACG_CDSnstop <- SHAP_vals_Fast %>% 
  select(ACG_CDSnstop, ACG_CDSnstop_SHAP, class)
SHAP_vals_MedFast_ACG_CDSnstop <- SHAP_vals_MedFast %>% 
  select(ACG_CDSnstop, ACG_CDSnstop_SHAP, class)
SHAP_vals_MedSlow_ACG_CDSnstop <- SHAP_vals_MedSlow %>% 
  select(ACG_CDSnstop, ACG_CDSnstop_SHAP, class)
SHAP_vals_Slow_ACG_CDSnstop <- SHAP_vals_Slow %>% 
  select(ACG_CDSnstop, ACG_CDSnstop_SHAP, class)

SHAP_ACG_CDSnstop <- bind_rows(SHAP_vals_Slow_ACG_CDSnstop, SHAP_vals_MedSlow_ACG_CDSnstop,
                               SHAP_vals_MedFast_ACG_CDSnstop, SHAP_vals_Fast_ACG_CDSnstop) %>% 
  mutate(class = as.factor(class))

SHAP_ACG_CDSnstop %>% 
  ggplot(aes(class, ACG_CDSnstop_SHAP)) +
  geom_point(aes(color = ACG_CDSnstop), 
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
ggsave("./viz_SHAPvalues_toCompareCondition_leadered/SHAPvalueTop20combined_toCompareCondition_HalfLifeCls_hypoxia_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseHypoxia_leadered_ACG_CDSnstop.png", width = 55, height = 22, units = "cm", dpi = 600)

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
ggsave("./viz_SHAPvalues_toCompareCondition_leadered/SHAPvalueTop20combined_toCompareCondition_HalfLifeCls_hypoxia_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseHypoxia_leadered_Nucl_A_CDS_5p18nt.png", width = 55, height = 22, units = "cm", dpi = 600)

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
ggsave("./viz_SHAPvalues_toCompareCondition_leadered/SHAPvalueTop20combined_toCompareCondition_HalfLifeCls_hypoxia_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseHypoxia_leadered_GC_CDS_5p18nt.png", width = 55, height = 22, units = "cm", dpi = 600)

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
ggsave("./viz_SHAPvalues_toCompareCondition_leadered/SHAPvalueTop20combined_toCompareCondition_HalfLifeCls_hypoxia_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseHypoxia_leadered_initialAbundance_hypoxia.png", width = 55, height = 22, units = "cm", dpi = 600)

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
ggsave("./viz_SHAPvalues_toCompareCondition_leadered/SHAPvalueTop20combined_toCompareCondition_HalfLifeCls_hypoxia_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseHypoxia_leadered_CDS_length.png", width = 55, height = 22, units = "cm", dpi = 600)

###### fpr_UTR_length
SHAP_vals_Fast_fpr_UTR_length <- SHAP_vals_Fast %>% 
  select(fpr_UTR_length, fpr_UTR_length_SHAP, class)
SHAP_vals_MedFast_fpr_UTR_length <- SHAP_vals_MedFast %>% 
  select(fpr_UTR_length, fpr_UTR_length_SHAP, class)
SHAP_vals_MedSlow_fpr_UTR_length <- SHAP_vals_MedSlow %>% 
  select(fpr_UTR_length, fpr_UTR_length_SHAP, class)
SHAP_vals_Slow_fpr_UTR_length <- SHAP_vals_Slow %>% 
  select(fpr_UTR_length, fpr_UTR_length_SHAP, class)

SHAP_fpr_UTR_length <- bind_rows(SHAP_vals_Slow_fpr_UTR_length, SHAP_vals_MedSlow_fpr_UTR_length, 
                             SHAP_vals_MedFast_fpr_UTR_length, SHAP_vals_Fast_fpr_UTR_length) %>% 
  mutate(class = as.factor(class))

SHAP_fpr_UTR_length %>% 
  ggplot(aes(class, fpr_UTR_length_SHAP)) +
  geom_point(aes(color = fpr_UTR_length), 
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
ggsave("./viz_SHAPvalues_toCompareCondition_leadered/SHAPvalueTop20combined_toCompareCondition_HalfLifeCls_hypoxia_combinedSelectedByLeaderedLogPhase_jointlyByLogPhaseHypoxia_leadered_fpr_UTR_length.png", width = 55, height = 22, units = "cm", dpi = 600)
