library(tidyverse)
library(zoo)
library(ggpubr)
library(gridExtra)
library(stringr)


#Run script from directory where outputs are located

#This will merge several kraken2 reports into one df
for (data in list.files()){ #List files in current directory as "data"
  
  # Create the initial df if none exists yet
  if (!exists("merged_vsp_reports")){ 
    merged_vsp_reports <- generate_tax_table(data) %>% #generate_tax_table() is 
                                                       #in functions_for_tax_analysis.R
      fill_tax_NAs() %>% #fill_tax_NAs() is also in functions_for_tax_analysis.R
      mutate(SampleID = data) 
  }
  
  # if df already exists, then append new data 
  if (exists("merged_vsp_reports")){
    temporary <- generate_tax_table(data) %>% #temporary table is necessary so 
                                              #that only the current sample is 
                                              #included in filling in and ID assignment
      fill_tax_NAs() %>%
      mutate(SampleID = data) #assign the file name as sample ID (it will be split later)
    merged_vsp_reports <- bind_rows(merged_vsp_reports, temporary) #dplyr row 
                                                                   #binding allows different 
                                                                   #numbers of columns
    rm(temporary) #save memory <3 
  }
}

head(merged_vsp_reports) #Check that import worked

#format imported data:
merged_vsp_reports <- merged_vsp_reports %>%
  separate(SampleID, c("QCSeqID", "LIMS_ID", "Treatment", "Other"), sep="-", 
           remove=FALSE) %>% #Make IDs into relevant parts
  mutate(site = str_replace_all( #this makes a new column even though it says "replace"
    LIMS_ID,
    c(
      "(36397|38156|41596|34970)" = "OBrien", #regex for any one of the given strings
      "(32976|32989|20040|22015)" = "OHare"
    )
  ))

#not needed for now
# grouped_vsp_reports <- merged_vsp_reports %>%
#  group_by(SampleID, )


#this sample has way more viral seqs than the others, making visualization difficult
outlier_removed <- merged_vsp_reports %>%
  filter(QCSeqID != "003")


#Plot all viruses that are assigned at the family level or lower (stacked bars with read counts), 
#faceted by sample and grouped by treatment
ggplot(filter(merged_vsp_reports, D == "Viruses" & !startsWith(F, "unidentified")), 
       aes(x = Treatment, y = nodeOnly, fill = F)) + 
  geom_bar(position = 'stack', stat = 'identity') +
  ylab("Number of reads") +
  xlab("Preparation method") +
  theme_bw() +
  labs(fill = "Viral Families") + #legend title
  theme(text = element_text(size = 15)) +
  facet_wrap(~LIMS_ID, scales = "free_y")


ggplot(filter(merged_vsp_reports, R == "root"), 
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
ggplot(filter(outlier_removed, D == "Viruses" & !startsWith(F, "unidentified")), 
       aes(x = LIMS_ID, y = RA, fill = F)) + 
  geom_bar(position = 'stack', stat = 'identity') +
  ylab("Number of reads") +
  xlab("Sample No.") +
  theme_bw() +
  labs(fill = "Viral Families") + #legend title
  theme(text = element_text(size = 15)) +
  facet_grid(~Treatment) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))


virus_counts <- merged_vsp_reports %>%
  filter(D == "Viruses") %>%
  group_by(SampleID) %>%
  summarise(sum = sum(nodeOnly), ra_sum = sum(RA))
  

