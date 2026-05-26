library(tidyverse)
library(zoo)

k2_rpip_main_dir <- paste0("input/link_to_raw_data/rpip_panels/raw_fastqs/",
                           "fastp_out_no_dedup/kraken2_out/k2_nt_20240530/",
                           "rrna_labeled/")

source("helper_functions.R")

#This will merge several kraken2 reports into one df
if(exists("merged_rpip_reports_full")) {rm(merged_rpip_reports_full)}
for (data in list.files(k2_rpip_main_dir,
                        pattern = "*.tsv",
                        full.names = TRUE)) { 
  
  # Create the initial df if none exists yet
  if (!exists("merged_rpip_reports_full")) { 
    merged_rpip_reports_full <- generate_tax_table(data) %>% 
      fill_tax_NAs() %>% 
      mutate(SampleID = basename(data)) 
  }
  
  # if df already exists, then append new data 
  if (exists("merged_rpip_reports_full")) {
    temporary <- generate_tax_table(data) %>% #temporary table is necessary so 
      #that only the current sample is included for filling in and ID assignment
      fill_tax_NAs() %>%
      mutate(SampleID = basename(data))
    merged_rpip_reports_full <- bind_rows(merged_rpip_reports_full, temporary)
    rm(temporary) #save memory <3 
  }
}


#format imported data:
merged_rpip_reports_full <- merged_rpip_reports_full %>%
  separate(SampleID, c("QCSeqID", "LIMS_ID", "Treatment", 
                       NA, NA, NA, NA, 
                       "ribosomal", "k2_confidence", NA),
           sep="(-|_)", 
           remove=FALSE) %>% #Make IDs into relevant parts
  mutate(site = str_replace_all(
    LIMS_ID,
    c("(36397|38156|41596|34970)" = "OBrien",
      "(32976|32989|20040|22015)" = "OHare")
    )
  )

#Save data for later:
write_csv(merged_rpip_reports_full,
          file = "input/modified/k2_rpip_reports_full.csv")
write_rds(merged_rpip_reports_full,
          file = "input/modified/k2_rpip_reports_full.rds")
