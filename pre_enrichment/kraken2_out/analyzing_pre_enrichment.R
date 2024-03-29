#####import libraries and data

library(tidyverse)
library(zoo)
library(ggpubr)
library(gridExtra)
library(stringr)


#Run script from directory where outputs are located
setwd("C:/Users/fanta/Dropbox/PC (2)/Documents/OBrienMetagenome/pre_enrichment/kraken2_out/")
merged_pre_reports_full <- read.csv("merged_pre_reports_full.csv")
######

merged_pre_reports_full <- merged_pre_reports_full %>%
  filter(Treatment %in% c("FA", "RA", "FB", "RB", "UA", "UB", "UD", "RD", "FD")) %>%
  mutate(fraction = substr(Treatment, 1, 1)) %>%
  mutate(fraction = factor(fraction, labels = 
                             c("Filtrate", "Retentate", "Unfiltered"))) %>%
  mutate(concentration_type = substr(Treatment, 2, 2)) %>%
  mutate(concentration_type = factor(concentration_type, labels = 
                                       c("A", "A&B", "None"))) %>%
  mutate(LIMS_ID = factor(LIMS_ID)) %>%
  mutate(site = factor(site))

###Add other column for families
merged_pre_reports_full <- merged_pre_reports_full %>%
  mutate(f_other = factor(F)) %>%
  mutate(f_other = fct_lump_min(f_other, min = 0.0001, w = RA, other_level = "Other families"))

ggplot(filter(merged_pre_reports_full, R == "root" & 
                SampleID != "001-36397-FA-RP_S1_L006_k2report.tsv"), 
       aes(x = Treatment, y = RA, fill = D)) + 
  geom_bar(position = 'stack', stat = 'identity') +
  ylab("Relative abundance") +
  xlab("Preparation method") +
  theme_bw() +
  labs(fill = "Domains") + #legend title
  theme(text = element_text(size = 15)) +
  facet_wrap(~LIMS_ID, scales = "free_y") +
  ggtitle("Pre-enrichment taxa")

###Plot virus families according to concentration used
ggplot(filter(merged_pre_reports_full, D == "Viruses" & !startsWith(F, "unidentified")), 
       aes(x = concentration_type, y = RA, fill = f_other)) + 
  geom_bar(position = 'stack', stat = 'identity') +
  ylab("Portion of reads") +
  xlab("Nanotrap version") +
  theme_bw() +
  labs(fill = "Viral Families") + #legend title
  theme(text = element_text(size = 15)) +
  facet_wrap(~fraction)
#  facet_wrap(fraction~LIMS_ID, scales = "free")
#  theme(legend.position = "none")


###Plot virus families according to concentration used
ggplot(filter(merged_pre_reports_full, D == "Viruses" & !startsWith(F, "unidentified")), 
       aes(x = f_other, y = RA)) + 
  geom_boxplot(aes(color = fraction)) +
  ylab("Portion of reads") +
  xlab("Nanotrap version") +
  theme_bw() +
  labs(fill = "Viral Families") + #legend title
  theme(text = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5)) +
  scale_y_log10()

###Plot virus families according to fraction used
ggplot(filter(merged_pre_reports_full, D == "Viruses" & !startsWith(F, "unidentified")), 
       aes(x = fraction, y = RA, fill = f_other)) + 
  geom_bar(position = 'stack', stat = 'identity') +
  ylab("Portion of reads") +
  xlab("Fraction") +
  theme_bw() +
  labs(fill = "Viral Families") + #legend title
  theme(text = element_text(size = 15)) #+
#  facet_wrap(~LIMS_ID, scales = "free")  
#  theme(legend.position = "none")


###Same as above but without removing anything
ggplot(filter(merged_pre_reports_full, D == "Viruses"), 
       aes(x = fraction, y = RA, fill = F)) + 
  geom_bar(position = 'stack', stat = 'identity') +
  ylab("Portion of reads") +
  xlab("Fraction") +
  theme_bw() +
  labs(fill = "Viral Families") + #legend title
  theme(text = element_text(size = 15)) +
  facet_wrap(~LIMS_ID, scales = "free") +
  theme(legend.position = "none")

ggplot(filter(merged_pre_reports_full, D == "Viruses" & !startsWith(F, "unidentified")), 
       aes(x = Treatment, y = RA, fill = f_other)) + 
  geom_bar(position = 'stack', stat = 'identity') +
  ylab("Portion of reads") +
  xlab("Preparation method") +
  theme_bw() +
  labs(fill = "Viral Families") + #legend title
  theme(text = element_text(size = 15)) +
  facet_wrap(~LIMS_ID, scales = "free_y") #+
  #theme(legend.position = "none")



ggplot(filter(merged_pre_reports_full, R == "root"), 
       aes(x = Treatment, y = nodeOnly, fill = D)) + 
  geom_bar(position = 'stack', stat = 'identity') +
  ylab("Number of reads") +
  xlab("Preparation method") +
  theme_bw() +
  labs(fill = "Domains") + #legend title
  theme(text = element_text(size = 15)) +
  facet_wrap(~LIMS_ID, scales = "free_y")

#Plot all viruses that are assigned at the family level or lower (stacked bars with read counts), 
#faceted by treatment and grouped by sample
ggplot(filter(merged_pre_reports_full, D == "Viruses" & !startsWith(F, "unidentified") & RA > 0.0001), 
       aes(x = LIMS_ID, y = nodeOnly, fill = F)) + 
  geom_bar(position = 'stack', stat = 'identity') +
  ylab("Number of reads") +
  xlab("Sample No.") +
  theme_bw() +
  labs(fill = "Viral Families") + #legend title
  theme(text = element_text(size = 15)) +
  facet_grid(~Treatment) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))

ggplot(filter(merged_pre_reports_full, D == "Viruses" & Treatment != "FA" & RA > 0.000001), 
       aes(x = LIMS_ID, y = nodeOnly, fill = F)) + 
  geom_bar(position = 'stack', stat = 'identity') +
  ylab("Number of reads") +
  xlab("Sample No.") +
  theme_bw() +
  labs(fill = "Viral Families") + #legend title
  theme(text = element_text(size = 15)) +
  facet_grid(~Treatment) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))


ggplot(filter(merged_pre_reports_full, F == "Polyomaviridae"), #Plot all polyomaviruses 
       aes(x = LIMS_ID, y = nodeOnly, fill = S)) + 
  geom_bar(position = 'stack', stat = 'identity') +
  ylab("Number of reads") +
  xlab("Sample No.") +
  theme_bw() +
  labs(fill = "Polyomavirus species") + #legend title
  theme(text = element_text(size = 15)) +
  facet_grid(~Treatment) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))

ggplot(filter(merged_pre_reports_full, F == "Coronaviridae"), #Plot all coronaviruses 
       aes(x = LIMS_ID, y = nodeOnly, fill = S)) + 
  geom_bar(position = 'stack', stat = 'identity') +
  ylab("Number of reads") +
  xlab("Sample No.") +
  theme_bw() +
  labs(fill = "Polyomavirus species") + #legend title
  theme(text = element_text(size = 15)) +
  facet_grid(~Treatment) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))

ggplot(filter(merged_pre_reports_full, F == "Astroviridae"), #Plot all polyomaviruses 
       aes(x = LIMS_ID, y = nodeOnly, fill = S)) + 
  geom_bar(position = 'stack', stat = 'identity') +
  ylab("Number of reads") +
  xlab("Sample No.") +
  theme_bw() +
  labs(fill = "Polyomavirus species") + #legend title
  theme(text = element_text(size = 15)) +
  facet_grid(~Treatment) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))


virus_counts <- merged_pre_reports %>%
  filter(D == "Viruses") %>%
  group_by(SampleID) %>%
  summarise(sum = sum(nodeOnly), ra_sum = sum(RA))

