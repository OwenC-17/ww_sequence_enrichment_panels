library(jsonlite)
library(tidyverse)

source("helper_functions.R")

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
  fastp_df$insert_size_peak <- fastp_json[[4]][[1]]
  fastp_df$unknown_ins_size_count <- fastp_json[[4]][[2]]
  fastp_df$filename <- basename(fastp_json_path)
  return(fastp_df)
}



#Function to load ALL fastp reports:
import_fastp_summaries <- function(directory) {
  file_list <- list.files(directory, pattern = "*\\.json$", full.names = TRUE)
  fastp_reports <- lapply(file_list, read_fastp_report)
  fastp_reports <- bind_rows(fastp_reports)
  fastp_reports <- separate(fastp_reports,
                            col = filename,
                            into = c("QCSeqID", "LIMS_ID", "Treatment",
                                     NA, NA, NA, NA, NA),
                            sep = "(-|_)",
                            remove = FALSE)
  #parse_sample_treatments is defined in import_k2_results.R
  fastp_reports <- parse_sample_treatments(fastp_reports)
  fastp_reports <- parse_locations(fastp_reports)
  fastp_reports <- mutate(fastp_reports,
                    num_reads_removed = summary.before_filtering.total_reads - 
                      summary.after_filtering.total_reads,
                    portion_reads_removed = num_reads_removed / 
                      summary.before_filtering.total_reads,
                    num_bases_removed = summary.before_filtering.total_bases -
                      summary.after_filtering.total_bases,
                    portion_bases_removed = num_bases_removed /
                      summary.before_filtering.total_bases
                    )
  return(fastp_reports)
}

###Import the fastp reports:
rpip_fastp_summaries_unmerged <- import_fastp_summaries(
  rpip_fastp_no_dedup_dir
  ) %>%
  mutate(Enrichment = "RPIP")

vsp_fastp_summaries_unmerged <- import_fastp_summaries(
  vsp_fastp_no_dedup_dir
  ) %>%
  mutate(Enrichment = "VSP")

unt_fastp_summaries_unmerged <- import_fastp_summaries(
  unt_fastp_no_dedup_dir
  ) %>%
  mutate(Enrichment = "None")

#Combine them:
all_fastp_summaries_unmerged <- rbind(rpip_fastp_summaries_unmerged,
                                      vsp_fastp_summaries_unmerged,
                                      unt_fastp_summaries_unmerged)


#Write files for later use:
write_rds(all_fastp_summaries_unmerged,
          "input/modified/fastp_summaries_unmerged_before_dedup.rds")

write_csv(all_fastp_summaries_unmerged,
          "input/modified/fastp_summaries_unmerged_before_dedup.csv")


#######################################
#####Read in deduplication reports#####
#######################################
vsp_dedup_report_dir <- paste0(topdir,
  "vsp_panels/raw_fastqs/fastp_out_no_dedup/deduped/")
rpip_dedup_report_dir <- paste0(topdir,
  "rpip_panels/raw_fastqs/fastp_out_no_dedup/deduped/")
unt_dedup_report_dir <- paste0(topdir,
  "untargeted/raw_fastqs/fastp_out_no_dedup/deduped/")

vsp_dedup_report <- import_fastp_summaries(vsp_dedup_report_dir) %>%
  mutate(Enrichment = "VSP")
rpip_dedup_report <- import_fastp_summaries(rpip_dedup_report_dir) %>%
  mutate(Enrichment = "RPIP")
unt_dedup_report <- import_fastp_summaries(unt_dedup_report_dir) %>%
  mutate(Enrichment = "None")

all_dedup_reports <- rbind(vsp_dedup_report,
                           rpip_dedup_report,
                           unt_dedup_report)


write_rds(all_dedup_reports, "input/modified/fastp_dedup_reports.rds")

write_csv(all_dedup_reports, "input/modified/fastp_dedup_reports.csv")

################################
#####Join dup/dedup reports#####
################################

#Pick columns that are not redundant:
dedup_relevant <- all_dedup_reports %>%
  select(summary.after_filtering.gc_content:summary.after_filtering.total_reads,
         insert_size_peak:Enrichment)

#Rename columns to reflect that they represent post-deduplication stats:
colnames(dedup_relevant) <- str_replace_all(colnames(dedup_relevant), 
                                        c("filtering" = "dedup", 
                                          "removed" = "removed_dedup",
                                          "insert_size_peak" = "insert_size_peak_dedup",
                                          "filename" = "filename_dedup",
                                          "unknown_ins_size_count" = "unknown_ins_size_count_dedup"))

#Join to pre-deduplicated data:
all_fastp_reports_dupdedup <- all_fastp_summaries_unmerged %>%
  left_join(dedup_relevant)

all_fastp_reports_dupdedup <- all_fastp_reports_dupdedup %>%
  mutate(num_reads_removed_qcANDdedup = summary.before_filtering.total_reads -
           summary.after_dedup.total_reads,
         portion_reads_removed_qcANDdedup = num_reads_removed_qcANDdedup /
           summary.before_filtering.total_reads,
         num_bases_removed_qcANDdedup = summary.before_filtering.total_bases -
           summary.after_dedup.total_bases,
         portion_bases_removed_qcANDdedup = num_bases_removed_qcANDdedup /
           summary.before_filtering.total_bases)

#Write files for later:
write_rds(all_fastp_reports_dupdedup, 
          "input/modified/all_fastp_reports_dupdedup.rds")

write_csv(all_fastp_reports_dupdedup,
          "input/modified/all_fastp_reports_dupdedup.csv")

####################################
#####Import merged read reports#####
####################################
#Function to load one merged fastp report:
read_merged_fastp_report <- function(fastp_json_path) {
  fastp_json <- fromJSON(fastp_json_path)
  fastp_df <- as.data.frame(fastp_json[-c(4:8)])
  fastp_df$insert_size_peak <- fastp_json[[4]][[1]]
  fastp_df$unknown_ins_size_count <- fastp_json[[4]][[2]]
  fastp_df$filename <- basename(fastp_json_path)
  return(fastp_df)
}

#Function to load ALL merged fastp reports:
import_merged_fastp_summaries <- function(directory) {
  file_list <- list.files(directory, pattern = "*\\.json$", full.names = TRUE)
  fastp_reports <- lapply(file_list, read_merged_fastp_report)
  fastp_reports <- bind_rows(fastp_reports)
  fastp_reports <- separate(fastp_reports,
                            col = filename,
                            into = c("QCSeqID", "LIMS_ID", "Treatment"),
                            sep = "(-|_)",
                            remove = FALSE,
                            extra = "drop")
  #parse_sample_treatments and parse_locationse are defined in 
  #import_k2_results.R
  fastp_reports <- parse_sample_treatments(fastp_reports)
  fastp_reports <- parse_locations(fastp_reports)
  fastp_reports <- mutate(fastp_reports,
                    num_reads_removed = summary.before_filtering.total_reads - 
                    (2 * summary.after_filtering.total_reads),
                    portion_reads_removed = num_reads_removed / 
                    summary.before_filtering.total_reads,
                    num_bases_removed = summary.before_filtering.total_bases -
                    summary.after_filtering.total_bases,
                    portion_bases_removed = num_bases_removed /
                    summary.before_filtering.total_bases
  )
  return(fastp_reports)
}


#Specify report locations:
rpip_fastp_merged_no_dedup_dir <- paste0(topdir,
                                         "rpip_panels/raw_fastqs/",
                                         "fastp_merged_no_dedup/reports/")
unt_fastp_merged_no_dedup_dir <- paste0(topdir,
                                        "untargeted/raw_fastqs/",
                                        "fastp_merged_no_dedup/reports/")

rpip_fastp_merged_then_dedup_dir <- paste0(
  topdir, "rpip_panels/raw_fastqs/fastp_merged_no_dedup/deduped"
)

unt_fastp_merged_then_dedup_dir <- paste0(
  topdir, "untargeted/raw_fastqs/fastp_merged_no_dedup/deduped"
)


#Import the reports:
rpip_fastp_summaries_merged <- import_merged_fastp_summaries(
  rpip_fastp_merged_no_dedup_dir
  ) %>%
  mutate(Enrichment = "RPIP")

rpip_fastp_summaries_merged_dedup <- import_merged_fastp_summaries(
  rpip_fastp_merged_then_dedup_dir
) %>%
  mutate(Enrichment = "RPIP")

unt_fastp_summaries_merged <- import_merged_fastp_summaries(
  unt_fastp_merged_no_dedup_dir
  ) %>%
  mutate(Enrichment = "None")

unt_fastp_summaries_merged_dedup <- import_merged_fastp_summaries(
  unt_fastp_merged_then_dedup_dir
) %>%
  mutate(Enrichment = "None")

#Combine
both_fastp_summaries_merged <- rbind(rpip_fastp_summaries_merged,
                                     unt_fastp_summaries_merged)

both_fastp_summaries_merged_dedup <- rbind(rpip_fastp_summaries_merged_dedup,
                                           unt_fastp_summaries_merged_dedup)

#Write files for later use:
write_rds(both_fastp_summaries_merged, 
          "input/modified/fastp_summaries_merged_before_dedup.rds")

write_csv(both_fastp_summaries_merged, 
          "input/modified/fastp_summaries_merged_before_dedup.csv")

write_rds(both_fastp_summaries_merged_dedup, 
          "input/modified/fastp_merged_dedup_reports.rds")

write_csv(both_fastp_summaries_merged_dedup,
          "input/modified/fastp_merged_dedup_reports.csv")
