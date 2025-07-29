library(jsonlite)
library(tidyverse)

###Define filepaths
topdir <- "input/link_to_raw_data/"
vsp_fastp_no_dedup_dir <- paste0(topdir,
                                 "vsp_panels/raw_fastqs/fastp_out_no_dedup/",
                                 "reports/")

rpip_fastp_no_dedup_dir <- paste0(topdir,
                                  "rpip_panels/raw_fastqs/fastp_out_no_dedup/",
                                  "reports/")

unt_fastp_no_dedup_dir <- paste0(topdir, 
                                 "untargeted/raw_fastqs/fastp_out_no_dedup/",
                                 "reports/")


#Function to load one fastp report:
read_fastp_report <- function(fastp_json_path) {
  fastp_json <- fromJSON(fastp_json_path)
  fastp_df <- as.data.frame(fastp_json[-c(4:9)])
  fastp_df$insert_size_peak <- fastp_json[[4]][1]
  fastp_df$unknown_ins_size_count <- fastp_json[[4]][2]
  fastp_df$SampleID <- basename(fastp_json_path)
  return(fastp_df)
}


#Function to load ALL fastp reports:
import_fastp_summaries <- function(directory) {
  file_list <- list.files(directory, pattern = "*\\.json$", full.names = TRUE)
  fastp_reports <- lapply(file_list, read_fastp_report)
  fastp_reports <- bind_rows(fastp_reports)
  return(fastp_reports)
}

rpip_fastp_summaries_unmerged <- import_fastp_summaries(rpip_fastp_no_dedup_dir)

#Convert sample names (from file names) into columns of relevant information:
rpip_fastp_summaries <- separate(rpip_fastp_summaries,
         col = sampleName,
         into = c("QCSeqID", "LIMS_ID", "Treatment",
                  NA, NA, NA, NA, NA, NA, NA, NA),
         sep = "(-|_)")

#Calculate portion of read pairs removed by fastp:
rpip_fastp_summaries <- rpip_fastp_summaries %>%
  mutate(portion_pairs_removed = 1 - (after_filtering.total_reads / (before_filtering.total_reads/2)))

#Label as RPIP
rpip_fastp_summaries$Enrichment = "RPIP"

#These functions have been defined in deeparg import script
rpip_fastp_summaries <- parse_sample_treatments(rpip_fastp_summaries)
rpip_fastp_summaries <- parse_locations(rpip_fastp_summaries)

#Import untargeted fastp results:
setwd("/projects/bios_microbe/cowen20/targeted_panels/untargeted/raw_fastqs/fastp_merged_no_dedup/deeparg_out/")
rm(untargeted_fastp_summaries)
for (fastp_report in list.files("../reports/", pattern = "*.json")) {
  fp_summary <- read_fastp_report(paste0("../reports/", fastp_report))
  fp_summary$sampleName <- fastp_report
  if (!exists("untargeted_fastp_summaries")) {
    untargeted_fastp_summaries <- fp_summary
  } else {
    untargeted_fastp_summaries <- rbind(untargeted_fastp_summaries, fp_summary)
  }
}

untargeted_fastp_summaries <- separate(untargeted_fastp_summaries,
                                 col = sampleName,
                                 into = c("QCSeqID", "LIMS_ID", "Treatment",
                                          NA, NA, NA, NA, NA, NA, NA, NA),
                                 sep = "(-|_)")

untargeted_fastp_summaries <- untargeted_fastp_summaries %>%
  mutate(portion_pairs_removed = 1 - (after_filtering.total_reads / (before_filtering.total_reads/2)))

untargeted_fastp_summaries$Enrichment = "None"

untargeted_fastp_summaries <- parse_sample_treatments(untargeted_fastp_summaries)
untargeted_fastp_summaries <- parse_locations(untargeted_fastp_summaries)


#combine them
fastp_reports_rpip_and_none <- rbind(rpip_fastp_summaries, untargeted_fastp_summaries)



########################################
#########PAIRED-END#####################
########################################

#Import vsp fastp results:
setwd("/projects/bios_microbe/cowen20/targeted_panels/vsp_panels/raw_fastqs/fastp_out_no_dedup/reports/")
rm(vsp_paired_fastp_summaries)
for (fastp_report in list.files(pattern = "*.json")) {
  fp_summary <- read_fastp_report(fastp_report)
  fp_summary$sampleName <- fastp_report
  if (!exists("vsp_paired_fastp_summaries")) {
    vsp_paired_fastp_summaries <- fp_summary
  } else {
    vsp_paired_fastp_summaries <- rbind(vsp_paired_fastp_summaries, fp_summary)
  }
}
vsp_paired_fastp_summaries <- separate(vsp_paired_fastp_summaries,
                                       col = sampleName,
                                       into = c("QCSeqID", "LIMS_ID", "Treatment",
                                                NA, NA, NA, NA, NA),
                                       sep = "(-|_)")

vsp_paired_fastp_summaries <- vsp_paired_fastp_summaries %>%
  mutate(portion_pairs_removed = 1 - (after_filtering.total_reads / (before_filtering.total_reads/2)))

vsp_paired_fastp_summaries$Enrichment = "VSP"

vsp_paired_fastp_summaries <- parse_sample_treatments(vsp_paired_fastp_summaries)
vsp_paired_fastp_summaries <- parse_locations(vsp_paired_fastp_summaries)


#Import rpip fastp results:
setwd("/projects/bios_microbe/cowen20/targeted_panels/rpip_panels/raw_fastqs/fastp_out_no_dedup/reports/")
rm(rpip_paired_fastp_summaries)
for (fastp_report in list.files(pattern = "*.json")) {
  fp_summary <- read_fastp_report(fastp_report)
  fp_summary$sampleName <- fastp_report
  if (!exists("rpip_paired_fastp_summaries")) {
    rpip_paired_fastp_summaries <- fp_summary
  } else {
    rpip_paired_fastp_summaries <- rbind(rpip_paired_fastp_summaries, fp_summary)
  }
}
rpip_paired_fastp_summaries <- separate(rpip_paired_fastp_summaries,
                                       col = sampleName,
                                       into = c("QCSeqID", "LIMS_ID", "Treatment",
                                                NA, NA, NA, NA, NA),
                                       sep = "(-|_)")

rpip_paired_fastp_summaries <- rpip_paired_fastp_summaries %>%
  mutate(portion_pairs_removed = 1 - (after_filtering.total_reads / (before_filtering.total_reads/2)))

rpip_paired_fastp_summaries$Enrichment = "RPIP"

rpip_paired_fastp_summaries <- parse_sample_treatments(rpip_paired_fastp_summaries)
rpip_paired_fastp_summaries <- parse_locations(rpip_paired_fastp_summaries)

#Import untargeted fastp results:
setwd("/projects/bios_microbe/cowen20/targeted_panels/untargeted/raw_fastqs/fastp_out_no_dedup/reports/")
rm(unt_paired_fastp_summaries)
for (fastp_report in list.files(pattern = "*.json")) {
  fp_summary <- read_fastp_report(fastp_report)
  fp_summary$sampleName <- fastp_report
  if (!exists("unt_paired_fastp_summaries")) {
    unt_paired_fastp_summaries <- fp_summary
  } else {
    unt_paired_fastp_summaries <- rbind(unt_paired_fastp_summaries, fp_summary)
  }
}
unt_paired_fastp_summaries <- separate(unt_paired_fastp_summaries,
                                       col = sampleName,
                                       into = c("QCSeqID", "LIMS_ID", "Treatment",
                                                NA, NA, NA, NA, NA),
                                       sep = "(-|_)")

unt_paired_fastp_summaries <- unt_paired_fastp_summaries %>%
  mutate(portion_pairs_removed = 1 - (after_filtering.total_reads / (before_filtering.total_reads/2)))

unt_paired_fastp_summaries$Enrichment = "None"

unt_paired_fastp_summaries <- parse_sample_treatments(unt_paired_fastp_summaries)
unt_paired_fastp_summaries <- parse_locations(unt_paired_fastp_summaries)


all_paired_fastp_summaries <- bind_rows(vsp_paired_fastp_summaries,
                                        rpip_paired_fastp_summaries,
                                        unt_paired_fastp_summaries)
