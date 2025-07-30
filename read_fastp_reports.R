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
  fastp_df$insert_size_peak <- fastp_json[[4]][[1]]
  fastp_df$unknown_ins_size_count <- fastp_json[[4]][[2]]
  fastp_df$SampleID <- basename(fastp_json_path)
  return(fastp_df)
}



#Function to load ALL fastp reports:
import_fastp_summaries <- function(directory) {
  file_list <- list.files(directory, pattern = "*\\.json$", full.names = TRUE)
  fastp_reports <- lapply(file_list, read_fastp_report)
  fastp_reports <- bind_rows(fastp_reports)
  fastp_reports <- separate(fastp_reports,
                            col = SampleID,
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

#combine them
all_fastp_summaries_unmerged <- rbind(rpip_fastp_summaries_unmerged,
                                      vsp_fastp_summaries_unmerged,
                                      unt_fastp_summaries_unmerged)

ggplot(all_fastp_summaries_unmerged, aes(x = Enrichment, y = portion_reads_removed)) + geom_boxplot()
ggplot(all_fastp_summaries_unmerged, aes(x = Enrichment, y = portion_bases_removed)) + geom_boxplot()
ggplot(all_fastp_summaries_unmerged, aes(x = Enrichment, y = summary.after_filtering.total_reads)) + geom_boxplot()


################################
###Import merged read reports###
################################
#Function to load one fastp report:
read_merged_fastp_report <- function(fastp_json_path) {
  fastp_json <- fromJSON(fastp_json_path)
  fastp_df <- as.data.frame(fastp_json[-c(4:8)])
  fastp_df$insert_size_peak <- fastp_json[[4]][[1]]
  fastp_df$unknown_ins_size_count <- fastp_json[[4]][[2]]
  fastp_df$SampleID <- basename(fastp_json_path)
  return(fastp_df)
}

#Function to load ALL fastp reports:
import_merged_fastp_summaries <- function(directory) {
  file_list <- list.files(directory, pattern = "*\\.json$", full.names = TRUE)
  fastp_reports <- lapply(file_list, read_merged_fastp_report)
  fastp_reports <- bind_rows(fastp_reports)
  fastp_reports <- separate(fastp_reports,
                            col = SampleID,
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


rpip_fastp_merged_no_dedup_dir <- paste0(topdir,
                                         "rpip_panels/raw_fastqs/",
                                         "fastp_merged_no_dedup/reports/")

unt_fastp_merged_no_dedup_dir <- paste0(topdir,
                                         "untargeted/raw_fastqs/",
                                         "fastp_merged_no_dedup/reports/")

rpip_fastp_summaries_merged <- import_merged_fastp_summaries(
  rpip_fastp_merged_no_dedup_dir
  ) %>%
  mutate(Enrichment = "RPIP")

unt_fastp_summaries_merged <- import_merged_fastp_summaries(
  unt_fastp_merged_no_dedup_dir
  ) %>%
  mutate(Enrichment = "None")

#Combine
both_fastp_summaries_merged <- rbind(rpip_fastp_summaries_merged,
                                     unt_fastp_summaries_merged)

ggplot(both_fastp_summaries_merged, aes(x = Enrichment, y = portion_reads_removed)) + geom_boxplot()
ggplot(both_fastp_summaries_merged, aes(x = Enrichment, y = portion_bases_removed)) + geom_boxplot()
ggplot(both_fastp_summaries_merged, aes(x = Enrichment, y = summary.after_filtering.total_reads)) + geom_boxplot()
ggplot(both_fastp_summaries_merged, aes(x = Enrichment, y = summary.before_filtering.total_reads)) + geom_boxplot()
ggplot(both_fastp_summaries_merged, aes(x = site, y = summary.after_filtering.total_reads)) + geom_boxplot()