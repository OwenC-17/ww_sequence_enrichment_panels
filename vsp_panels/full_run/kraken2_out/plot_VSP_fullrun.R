library(tidyverse)
library(zoo)
library(ggpubr)
library(gridExtra)
library(stringr)


#Run script from directory where outputs are located
setwd("C:/Users/fanta/Dropbox/PC (2)/Documents/OBrienMetagenome/vsp_panels/full_run/kraken2_out")



#NOTE: THE DATA FROM THE FOLLOWING IMPORT HAS BEEN IMPORTED AND SAVED AS 
#"merged_vsp_reports_full.csv"; do not run the import loop again unless the
#original files change, as it will take hours

#This will merge several kraken2 reports into one df
rm(merged_vsp_reports_full)
for (data in list.files()){ #List files in current directory as "data"
  
  # Create the initial df if none exists yet
  if (!exists("merged_vsp_reports_full")){ 
    merged_vsp_reports_full <- generate_tax_table(data) %>% #generate_tax_table() is 
      #in functions_for_tax_analysis.R
      fill_tax_NAs() %>% #fill_tax_NAs() is also in functions_for_tax_analysis.R
      mutate(SampleID = data) 
  }
  
  # if df already exists, then append new data 
  if (exists("merged_vsp_reports_full")){
    temporary <- generate_tax_table(data) %>% #temporary table is necessary so 
      #that only the current sample is 
      #included in filling in and ID assignment
      fill_tax_NAs() %>%
      mutate(SampleID = data) #assign the file name as sample ID (it will be split later)
    merged_vsp_reports_full <- bind_rows(merged_vsp_reports_full, temporary) #dplyr row 
    #binding allows different 
    #numbers of columns
    rm(temporary) #save memory <3 
  }
}


#format imported data:
merged_vsp_reports_full <- merged_vsp_reports_full %>%
  separate(SampleID, c("QCSeqID", "LIMS_ID", "Treatment", "Other"), sep="-", 
           remove=FALSE) %>% #Make IDs into relevant parts
  mutate(site = str_replace_all( #this makes a new column even though it says "replace"
    LIMS_ID,
    c(
      "(36397|38156|41596|34970)" = "OBrien", #regex for any one of the given strings
      "(32976|32989|20040|22015)" = "OHare"
    )
  ))

write_csv(x = merged_vsp_reports_full, file = "merged_vsp_reports_full.csv")



ggplot(filter(merged_vsp_reports_full, D == "Viruses" & !startsWith(F, "unidentified")), 
       aes(x = Treatment, y = nodeOnly, fill = F)) + 
  geom_bar(position = 'stack', stat = 'identity') +
  ylab("Number of reads") +
  xlab("Preparation method") +
  theme_bw() +
  labs(fill = "Viral Families") + #legend title
  theme(text = element_text(size = 15)) +
  facet_wrap(~LIMS_ID, scales = "free_y") +
  theme(legend.position = "none")



ggplot(filter(merged_vsp_reports_full, R == "root"), 
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
ggplot(filter(merged_vsp_reports_full, D == "Viruses" & !startsWith(F, "unidentified") & RA > 0.0001), 
       aes(x = LIMS_ID, y = nodeOnly, fill = F)) + 
  geom_bar(position = 'stack', stat = 'identity') +
  ylab("Number of reads") +
  xlab("Sample No.") +
  theme_bw() +
  labs(fill = "Viral Families") + #legend title
  theme(text = element_text(size = 15)) +
  facet_grid(~Treatment) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))

ggplot(filter(merged_vsp_reports_full, D == "Viruses" & Treatment != "FA" & RA > 0.000001), 
       aes(x = LIMS_ID, y = nodeOnly, fill = F)) + 
  geom_bar(position = 'stack', stat = 'identity') +
  ylab("Number of reads") +
  xlab("Sample No.") +
  theme_bw() +
  labs(fill = "Viral Families") + #legend title
  theme(text = element_text(size = 15)) +
  facet_grid(~Treatment) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))


ggplot(filter(merged_vsp_reports_full, F == "Polyomaviridae"), #Plot all polyomaviruses 
       aes(x = LIMS_ID, y = nodeOnly, fill = S)) + 
  geom_bar(position = 'stack', stat = 'identity') +
  ylab("Number of reads") +
  xlab("Sample No.") +
  theme_bw() +
  labs(fill = "Polyomavirus species") + #legend title
  theme(text = element_text(size = 15)) +
  facet_grid(~Treatment) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))

ggplot(filter(merged_vsp_reports_full, F == "Coronaviridae"), #Plot all coronaviruses 
       aes(x = LIMS_ID, y = nodeOnly, fill = S)) + 
  geom_bar(position = 'stack', stat = 'identity') +
  ylab("Number of reads") +
  xlab("Sample No.") +
  theme_bw() +
  labs(fill = "Polyomavirus species") + #legend title
  theme(text = element_text(size = 15)) +
  facet_grid(~Treatment) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))

ggplot(filter(merged_vsp_reports_full, F == "Astroviridae"), #Plot all polyomaviruses 
       aes(x = LIMS_ID, y = nodeOnly, fill = S)) + 
  geom_bar(position = 'stack', stat = 'identity') +
  ylab("Number of reads") +
  xlab("Sample No.") +
  theme_bw() +
  labs(fill = "Polyomavirus species") + #legend title
  theme(text = element_text(size = 15)) +
  facet_grid(~Treatment) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))


virus_counts <- merged_vsp_reports %>%
  filter(D == "Viruses") %>%
  group_by(SampleID) %>%
  summarise(sum = sum(nodeOnly), ra_sum = sum(RA))
