library(tidyverse)
library(zoo)
library(ggpubr)
library(gridExtra)
library(stringr)


#Run script from directory where outputs are located
setwd("/projects/bios_microbe/cowen20/targeted_panels/rpip_panels/raw_fastqs/fastp_out_no_dedup/kraken2_out/k2_nt_20240530/rrna_labeled/")



#This will merge several kraken2 reports into one df
rm(merged_rpip_reports_full)
for (data in list.files(pattern = "*.tsv")){ #List files in current directory as "data"
  
  # Create the initial df if none exists yet
  if (!exists("merged_rpip_reports_full")){ 
    merged_rpip_reports_full <- generate_tax_table(data) %>% #generate_tax_table() is 
      #in functions_for_tax_analysis.R
      fill_tax_NAs() %>% #fill_tax_NAs() is also in functions_for_tax_analysis.R
      mutate(SampleID = data) 
  }
  
  # if df already exists, then append new data 
  if (exists("merged_rpip_reports_full")){
    temporary <- generate_tax_table(data) %>% #temporary table is necessary so 
      #that only the current sample is 
      #included in filling in and ID assignment
      fill_tax_NAs() %>%
      mutate(SampleID = data) #assign the file name as sample ID (it will be split later)
    merged_rpip_reports_full <- bind_rows(merged_rpip_reports_full, temporary) #dplyr row 
    #binding allows different 
    #numbers of columns
    rm(temporary) #save memory <3 
  }
}


#format imported data:
merged_rpip_reports_full <- merged_rpip_reports_full %>%
  separate(SampleID, c("QCSeqID", "LIMS_ID", "Treatment", "Other"), sep="-", 
           remove=FALSE) %>% #Make IDs into relevant parts
  mutate(site = str_replace_all( #this makes a new column even though it says "replace"
    LIMS_ID,
    c(
      "(36397|38156|41596|34970)" = "OBrien", #regex for any one of the given strings
      "(32976|32989|20040|22015)" = "OHare"
      
    )
  ))

write_csv(x = merged_rpip_reports_full, file = "merged_rpip_reports_full.csv")


merged_rpip_reports_full <- read_csv("merged_rpip_reports_full.csv")
