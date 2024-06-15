
####### To get the feature type for features that are excluded from selection
####### so that, to quantify the number of features being selected in each type
####### for both leadered and leaderless transcripts
####### for selection in log phase, hypoxia and fold change in hypoxia

library(tidyverse)

####### get feature type annotation
###### 5'UTR
path = '../../feature/FeatureSet_byType_5pUTR_UTRrelated_leadered/'
files = list.files(path, pattern = '.txt')
setwd('../../feature/FeatureSet_byType_5pUTR_UTRrelated_leadered/')
table_list = lapply(files, read.table, header = TRUE)
features_UTR = do.call(cbind, table_list)
setwd('../../ml_classification/viz_results/')

###### CDS nucleotide
path = '../../feature/FeatureSet_byType_CDSnucleotide/'
files = list.files(path, pattern = '.txt')
setwd('../../feature/FeatureSet_byType_CDSnucleotide/')
table_list = lapply(files, read.table, header = TRUE)
features_CDSnucleotide = do.call(cbind, table_list)
setwd('../../ml_classification/viz_results/')

###### CDS MFE
path = '../../feature/FeatureSet_byType_CDSmfe/'
files = list.files(path, pattern = '.txt')
setwd('../../feature/FeatureSet_byType_CDSmfe/')
table_list = lapply(files, read.table, header = TRUE)
features_CDSmfe = do.call(cbind, table_list)
setwd('../../ml_classification/viz_results/')

###### Codon
path = '../../feature/FeatureSet_byType_Codon/'
files = list.files(path, pattern = '.txt')
setwd('../../feature/FeatureSet_byType_Codon/')
table_list = lapply(files, read.table, header = TRUE)
features_Codon = do.call(cbind, table_list)
setwd('../../ml_classification/viz_results/')

###### Translation
path = '../../feature/FeatureSet_byType_Translation_leadered/'
files = list.files(path, pattern = '.txt')
setwd('../../feature/FeatureSet_byType_Translation_leadered/')
table_list = lapply(files, read.table, header = TRUE)
features_Translation = do.call(cbind, table_list)
setwd('../../ml_classification/viz_results/')

###### Others
##### log phase
path = '../../feature/FeatureSet_byType_Others_logPhase/'
files = list.files(path, pattern = '.txt')
setwd('../../feature/FeatureSet_byType_Others_logPhase/')
table_list = lapply(files, read.table, header = TRUE)
features_Others_logPhase = do.call(cbind, table_list)
setwd('../../ml_classification/viz_results/')

##### hypoxia
path = '../../feature/FeatureSet_byType_Others_hypoxia/'
files = list.files(path, pattern = '.txt')
setwd('../../feature/FeatureSet_byType_Others_hypoxia/')
table_list = lapply(files, read.table, header = TRUE)
features_Others_hypoxia = do.call(cbind, table_list)
setwd('../../ml_classification/viz_results/')

UTR <- data.frame(colnames(features_UTR), "UTR")
colnames(UTR) <- c("feature", "type")
CDSnucleotide <- data.frame(colnames(features_CDSnucleotide), "CDSnucleotide")
colnames(CDSnucleotide) <- c("feature", "type")
CDSmfe <- data.frame(colnames(features_CDSmfe), "CDSmfe")
colnames(CDSmfe) <- c("feature", "type")
Codon <- data.frame(colnames(features_Codon), "Codon")
colnames(Codon) <- c("feature", "type")
Translation <- data.frame(colnames(features_Translation), "Translation")
colnames(Translation) <- c("feature", "type")
other_p1 <- data.frame(colnames(features_Others_logPhase), "other")
colnames(other_p1) <- c("feature", "type")
other_p2 <- data.frame(colnames(features_Others_hypoxia), "other")
colnames(other_p2) <- c("feature", "type")

###### combine all feature type
feature_completeList <- bind_rows(UTR, CDSnucleotide, CDSmfe, Codon, Translation, other_p1, other_p2) %>% distinct()

####### get excluded features
###### leadered 
feature_halfLifeCls_logPhase_excludedByLeadered <- read.table('../../feature/feature_selection/feature_halfLifeCls_logPhase_excludedByLeadered.txt', col.names = "feature")
feature_halfLifeCls_logPhase_excludedByLeadered_withType <- feature_halfLifeCls_logPhase_excludedByLeadered %>% 
  inner_join(feature_completeList, by = "feature") %>% 
  mutate(type = factor(type))
summary(feature_halfLifeCls_logPhase_excludedByLeadered_withType$type)
feature_halfLifeCls_logPhase_excludedByLeadered_withType %>% count(type)

feature_halfLifeCls_hypoxia_excludedByLeadered <- read.table('../../feature/feature_selection/feature_halfLifeCls_hypoxia_excludedByLeadered.txt', col.names = "feature")
feature_halfLifeCls_hypoxia_excludedByLeadered_withType <- feature_halfLifeCls_hypoxia_excludedByLeadered %>% 
  inner_join(feature_completeList, by = "feature") %>%
  mutate(type = factor(type))
summary(feature_halfLifeCls_hypoxia_excludedByLeadered_withType$type)
feature_halfLifeCls_hypoxia_excludedByLeadered_withType %>% count(type)

feature_halfLifeFcCls_hypoxiaToLogPhase_excludedByLeadered <- read.table('../../feature/feature_selection/feature_halfLifeFcCls_hypoxiaToLogPhase_excludedByLeadered.txt', col.names = "feature")
feature_halfLifeFcCls_hypoxiaToLogPhase_excludedByLeadered_withType <- feature_halfLifeFcCls_hypoxiaToLogPhase_excludedByLeadered %>% 
  inner_join(feature_completeList, by = "feature") %>%
  mutate(type = factor(type))
summary(feature_halfLifeFcCls_hypoxiaToLogPhase_excludedByLeadered_withType$type)
feature_halfLifeFcCls_hypoxiaToLogPhase_excludedByLeadered_withType %>% count(type)

###### leaderless
feature_halfLifeCls_logPhase_excludedByLeaderless <- read.table('../../feature/feature_selection/feature_halfLifeCls_logPhase_excludedByLeaderless.txt', col.names = "feature")
feature_halfLifeCls_logPhase_excludedByLeaderless_withType <- feature_halfLifeCls_logPhase_excludedByLeaderless %>% 
  inner_join(feature_completeList, by = "feature") %>%
  mutate(type = factor(type))
summary(feature_halfLifeCls_logPhase_excludedByLeaderless_withType$type)
feature_halfLifeCls_logPhase_excludedByLeaderless_withType %>% count(type)

feature_halfLifeCls_hypoxia_excludedByLeaderless <- read.table('../../feature/feature_selection/feature_halfLifeCls_hypoxia_excludedByLeaderless.txt', col.names = "feature")
feature_halfLifeCls_hypoxia_excludedByLeaderless_withType <- feature_halfLifeCls_hypoxia_excludedByLeaderless %>% 
  inner_join(feature_completeList, by = "feature") %>%
  mutate(type = factor(type))
summary(feature_halfLifeCls_hypoxia_excludedByLeaderless_withType$type)
feature_halfLifeCls_hypoxia_excludedByLeaderless_withType %>% count(type)

feature_halfLifeFcCls_hypoxiaToLogPhase_excludedByLeaderless <- read.table('../../feature/feature_selection/feature_halfLifeFcCls_hypoxiaToLogPhase_excludedByLeaderless.txt', col.names = "feature")
feature_halfLifeFcCls_hypoxiaToLogPhase_excludedByLeaderless_withType <- feature_halfLifeFcCls_hypoxiaToLogPhase_excludedByLeaderless %>%
  inner_join(feature_completeList, by = "feature") %>%
  mutate(type = factor(type))
summary(feature_halfLifeFcCls_hypoxiaToLogPhase_excludedByLeaderless_withType$type)
feature_halfLifeFcCls_hypoxiaToLogPhase_excludedByLeaderless_withType %>% count(type)
