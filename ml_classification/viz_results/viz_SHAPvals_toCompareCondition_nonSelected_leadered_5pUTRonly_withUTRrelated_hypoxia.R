
####### viz SHAP values of features
####### To compare features of log phase, hypoxia and fold change in hypoxia 
####### for leadered transcripts non-selected feature models with only 5'UTR features
####### features are the top20 important features for the three conditions
####### with only UTR related features
####### hypoxia

library(tidyverse)
library(RColorBrewer)
library(scales)
library(viridis)

####### get color scheme for SHAP values
color_SHAP = c("#900DA4", "#B12A90", "#CC4678", "#E16462", "#F1844B", "#FCA636", "#FCCE25")

####### get SHAP values for each class
SHAP_vals_Fast <- read.csv('../SHAP_values/SHAPvals_HalfLifeCls_hypoxia_nonSelected_leadered_byType_5pUTR_UTRrelated_Fast.csv') %>% 
  mutate(class = 'Fast')
SHAP_vals_MedFast <- read.csv('../SHAP_values/SHAPvals_HalfLifeCls_hypoxia_nonSelected_leadered_byType_5pUTR_UTRrelated_MedFast.csv') %>% 
  mutate(class = 'MedFast')
SHAP_vals_MedSlow <- read.csv('../SHAP_values/SHAPvals_HalfLifeCls_hypoxia_nonSelected_leadered_byType_5pUTR_UTRrelated_MedSlow.csv') %>% 
  mutate(class = 'MedSlow')
SHAP_vals_Slow <- read.csv('../SHAP_values/SHAPvals_HalfLifeCls_hypoxia_nonSelected_leadered_byType_5pUTR_UTRrelated_Slow.csv') %>% 
  mutate(class = 'Slow')

####### viz of SHAP values
###### adja_GC_5pUTR
SHAP_vals_Fast_adja_GC_5pUTR <- SHAP_vals_Fast %>% 
  select(adja_GC_5pUTR, adja_GC_5pUTR_SHAP, class)
SHAP_vals_MedFast_adja_GC_5pUTR <- SHAP_vals_MedFast %>% 
  select(adja_GC_5pUTR, adja_GC_5pUTR_SHAP, class)
SHAP_vals_MedSlow_adja_GC_5pUTR <- SHAP_vals_MedSlow %>% 
  select(adja_GC_5pUTR, adja_GC_5pUTR_SHAP, class)
SHAP_vals_Slow_adja_GC_5pUTR <- SHAP_vals_Slow %>% 
  select(adja_GC_5pUTR, adja_GC_5pUTR_SHAP, class)

SHAP_adja_GC_5pUTR <- bind_rows(SHAP_vals_Slow_adja_GC_5pUTR, SHAP_vals_MedSlow_adja_GC_5pUTR,
                                SHAP_vals_MedFast_adja_GC_5pUTR, SHAP_vals_Fast_adja_GC_5pUTR) %>% 
  mutate(class = as.factor(class))

SHAP_adja_GC_5pUTR %>% 
  ggplot(aes(class, adja_GC_5pUTR_SHAP)) +
  geom_point(aes(color = adja_GC_5pUTR), 
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
ggsave("./viz_SHAPvalues_toCompareCondition_leadered_5pUTRonly_withUTRrelated/SHAPvalueTop20combined_toCompareCondition_HalfLifeCls_hypoxia_nonSelected_leadered_5pUTRonly_withUTRrelated_adja_GC_5pUTR.png", width = 55, height = 22, units = "cm", dpi = 600)

###### fprUTR_MFE_20_10nt_mid
SHAP_vals_Fast_fprUTR_MFE_20_10nt_mid <- SHAP_vals_Fast %>% 
  select(fprUTR_MFE_20_10nt_mid, fprUTR_MFE_20_10nt_mid_SHAP, class)
SHAP_vals_MedFast_fprUTR_MFE_20_10nt_mid <- SHAP_vals_MedFast %>% 
  select(fprUTR_MFE_20_10nt_mid, fprUTR_MFE_20_10nt_mid_SHAP, class)
SHAP_vals_MedSlow_fprUTR_MFE_20_10nt_mid <- SHAP_vals_MedSlow %>% 
  select(fprUTR_MFE_20_10nt_mid, fprUTR_MFE_20_10nt_mid_SHAP, class)
SHAP_vals_Slow_fprUTR_MFE_20_10nt_mid <- SHAP_vals_Slow %>% 
  select(fprUTR_MFE_20_10nt_mid, fprUTR_MFE_20_10nt_mid_SHAP, class)

SHAP_fprUTR_MFE_20_10nt_mid <- bind_rows(SHAP_vals_Slow_fprUTR_MFE_20_10nt_mid, SHAP_vals_MedSlow_fprUTR_MFE_20_10nt_mid,
                                         SHAP_vals_MedFast_fprUTR_MFE_20_10nt_mid, SHAP_vals_Fast_fprUTR_MFE_20_10nt_mid) %>% 
  mutate(class = as.factor(class))

SHAP_fprUTR_MFE_20_10nt_mid %>% 
  ggplot(aes(class, fprUTR_MFE_20_10nt_mid_SHAP)) +
  geom_point(aes(color = fprUTR_MFE_20_10nt_mid), 
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
ggsave("./viz_SHAPvalues_toCompareCondition_leadered_5pUTRonly_withUTRrelated/SHAPvalueTop20combined_toCompareCondition_HalfLifeCls_hypoxia_nonSelected_leadered_5pUTRonly_withUTRrelated_fprUTR_MFE_20_10nt_mid.png", width = 55, height = 22, units = "cm", dpi = 600)

###### GC_5pUTR
SHAP_vals_Fast_GC_5pUTR <- SHAP_vals_Fast %>% 
  select(GC_5pUTR, GC_5pUTR_SHAP, class)
SHAP_vals_MedFast_GC_5pUTR <- SHAP_vals_MedFast %>% 
  select(GC_5pUTR, GC_5pUTR_SHAP, class)
SHAP_vals_MedSlow_GC_5pUTR <- SHAP_vals_MedSlow %>% 
  select(GC_5pUTR, GC_5pUTR_SHAP, class)
SHAP_vals_Slow_GC_5pUTR <- SHAP_vals_Slow %>% 
  select(GC_5pUTR, GC_5pUTR_SHAP, class)

SHAP_GC_5pUTR <- bind_rows(SHAP_vals_Slow_GC_5pUTR, SHAP_vals_MedSlow_GC_5pUTR, 
                           SHAP_vals_MedFast_GC_5pUTR, SHAP_vals_Fast_GC_5pUTR) %>% 
  mutate(class = as.factor(class))

SHAP_GC_5pUTR %>% 
  ggplot(aes(class, GC_5pUTR_SHAP)) +
  geom_point(aes(color = GC_5pUTR), 
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
ggsave("./viz_SHAPvalues_toCompareCondition_leadered_5pUTRonly_withUTRrelated/SHAPvalueTop20combined_toCompareCondition_HalfLifeCls_hypoxia_nonSelected_leadered_5pUTRonly_withUTRrelated_GC_5pUTR.png", width = 55, height = 22, units = "cm", dpi = 600)

###### Nucl_T_5pUTR
SHAP_vals_Fast_Nucl_T_5pUTR <- SHAP_vals_Fast %>% 
  select(Nucl_T_5pUTR, Nucl_T_5pUTR_SHAP, class)
SHAP_vals_MedFast_Nucl_T_5pUTR <- SHAP_vals_MedFast %>% 
  select(Nucl_T_5pUTR, Nucl_T_5pUTR_SHAP, class)
SHAP_vals_MedSlow_Nucl_T_5pUTR <- SHAP_vals_MedSlow %>% 
  select(Nucl_T_5pUTR, Nucl_T_5pUTR_SHAP, class)
SHAP_vals_Slow_Nucl_T_5pUTR <- SHAP_vals_Slow %>% 
  select(Nucl_T_5pUTR, Nucl_T_5pUTR_SHAP, class)

SHAP_Nucl_T_5pUTR <- bind_rows(SHAP_vals_Slow_Nucl_T_5pUTR, SHAP_vals_MedSlow_Nucl_T_5pUTR,
                               SHAP_vals_MedFast_Nucl_T_5pUTR, SHAP_vals_Fast_Nucl_T_5pUTR) %>% 
  mutate(class = as.factor(class))

SHAP_Nucl_T_5pUTR %>% 
  ggplot(aes(class, Nucl_T_5pUTR_SHAP)) +
  geom_point(aes(color = Nucl_T_5pUTR), 
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
ggsave("./viz_SHAPvalues_toCompareCondition_leadered_5pUTRonly_withUTRrelated/SHAPvalueTop20combined_toCompareCondition_HalfLifeCls_hypoxia_nonSelected_leadered_5pUTRonly_withUTRrelated_Nucl_T_5pUTR.png", width = 55, height = 22, units = "cm", dpi = 600)
