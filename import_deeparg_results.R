library(tidyverse)
source("helper_functions.R")

topdir <- "input/link_to_raw_data/"
rpip_deeparg_dir <- paste0(topdir,
                           "rpip_panels/raw_fastqs/fastp_merged_no_dedup/",
                           "deeparg_out/")

rpip_deeparg_dedup_dir <- paste0(rpip_deeparg_dir, "deduped/")
list.files(rpip_deeparg_dir, pattern = "*.ARG$")

unt_deeparg_dir <- paste0(topdir,
                          "untargeted/raw_fastqs/fastp_merged_no_dedup/",
                          "deeparg_out/")

unt_deeparg_dedup_dir <- paste0(unt_deeparg_dir, "deduped/")

import_ARG_mapfile <- function(filename) {
  imported_ARG_table <- readr::read_tsv(filename, col_types = "cddcccdddddd")
  imported_ARG_table$SampleID <- basename(filename)
  imported_ARG_table <- imported_ARG_table %>%
    rename(ARG = `#ARG`)
  colnames(imported_ARG_table) <- str_replace_all(colnames(imported_ARG_table),
                                                  pattern = "-", 
                                                  replacement = "_")
  return(imported_ARG_table)
}

import_ARG_mapfiles <- function(directory) {
  file_list <- list.files(directory, pattern = "*.ARG$", full.names = TRUE)
  imported_ARG_list <- lapply(file_list, import_ARG_mapfile)
  combined_ARG_df <- bind_rows(imported_ARG_list)
  combined_ARG_df <- separate(combined_ARG_df,
                              col = SampleID,
                              into = c("QCSeqID", "LIMS_ID", "Treatment"),
                              sep = "(-|_)",
                              remove = FALSE,
                              extra = "drop")
  #parse_sample_treatments and parse_locations are defined in 
  #import_k2_results.R
  combined_ARG_df <- parse_sample_treatments(combined_ARG_df)
  combined_ARG_df <- parse_locations(combined_ARG_df)
  return(combined_ARG_df)
}

#Import the Deeparg read mapping files:
rpip_deeparg_results_table <- import_ARG_mapfiles(rpip_deeparg_dir)
rpip_deeparg_results_table_dedup <- import_ARG_mapfiles(rpip_deeparg_dedup_dir)

unt_deeparg_results_table <- import_ARG_mapfiles(unt_deeparg_dir)
unt_deeparg_results_table_dedup <- import_ARG_mapfiles(unt_deeparg_dedup_dir)

#Separate potential from highly-probably ARG seqs **(move this to import fun)**
add_potential <- function(arg_df) {
  mutate(arg_df, potential = str_detect(SampleID, "potential"))
}

rpip_deeparg_results_table <- add_potential(rpip_deeparg_results_table)
rpip_deeparg_results_table_dedup <- add_potential(rpip_deeparg_results_table_dedup)
unt_deeparg_results_table <- add_potential(unt_deeparg_results_table)
unt_deeparg_results_table_dedup <- add_potential(unt_deeparg_results_table_dedup)

#Add enrichment column
rpip_deeparg_results_table$Enrichment <- "RPIP"
rpip_deeparg_results_table_dedup$Enrichment <- "RPIP"
unt_deeparg_results_table$Enrichment <- "None"
unt_deeparg_results_table_dedup$Enrichment <- "None"

#combine them
rpip_vs_unt_deeparg_results_table <- rbind(rpip_deeparg_results_table, 
                                           unt_deeparg_results_table)

rpip_vs_unt_deeparg_results_table_dedup <- rbind(
  rpip_deeparg_results_table_dedup,
  unt_deeparg_results_table_dedup
)
#Create an ID that contains all treatment information and is unique:
rpip_vs_unt_deeparg_results_table <- rpip_vs_unt_deeparg_results_table %>%
  mutate(UniqueID = paste(LIMS_ID, Treatment, Enrichment, sep = "-"))
rpip_vs_unt_deeparg_results_table_dedup <- rpip_vs_unt_deeparg_results_table_dedup %>%
  mutate(UniqueID = paste(LIMS_ID, Treatment, Enrichment, sep = "-"))


#####IMPORTANT: Run read_fastp_reports.R before continuing#####
source("read_fastp_reports.R")

#Make matching column to join with deeparg results:
both_fastp_summaries_merged <- both_fastp_summaries_merged %>%
  mutate(UniqueID = paste(LIMS_ID, Treatment, Enrichment, sep = "-"))

both_fastp_summaries_merged_dedup <- both_fastp_summaries_merged_dedup %>%
  mutate(UniqueID = paste(LIMS_ID, Treatment, Enrichment, sep = "-"))

#Sum the number of reads for each ARG in each sample (some of them are repeated
#due to differences in match stats from Diamond alignment)
samplewise_deeparg_results_table <- rpip_vs_unt_deeparg_results_table %>% 
  group_by(ARG, predicted_ARG_class, LIMS_ID, Treatment, Enrichment, Fraction, 
           Nanotrap_type, site, UniqueID) %>%
  summarize(NumReads = sum(counts))

samplewise_deeparg_results_table_dedup <- rpip_vs_unt_deeparg_results_table_dedup %>% 
  group_by(ARG, predicted_ARG_class, LIMS_ID, Treatment, Enrichment, Fraction, 
           Nanotrap_type, site, UniqueID) %>%
  summarize(NumReads = sum(counts))


#Calculate how many reads are classified as ARGs in each sample:
total_args_per_sample <- samplewise_deeparg_results_table %>%
  group_by(UniqueID) %>%
  summarize(AllArgs = sum(NumReads))

total_args_per_sample_dedup <- samplewise_deeparg_results_table_dedup %>%
  group_by(UniqueID) %>%
  summarize(AllArgs = sum(NumReads))

#Join total ARG counts to fastp reports, so now we have total nreads and total
#ARG-classified reads:
both_fastp_summaries_merged_with_argcounts <- both_fastp_summaries_merged %>%
  left_join(total_args_per_sample, by = "UniqueID")

fastp_summaries_merged_with_argcounts_dedup <- both_fastp_summaries_merged_dedup %>%
  left_join(total_args_per_sample_dedup, by = "UniqueID")

#Calculate total number of non ARG-classified reads
both_fastp_summaries_merged_with_argcounts <- 
  both_fastp_summaries_merged_with_argcounts %>%
  mutate(NonArgReads = summary.after_filtering.total_reads - AllArgs)

fastp_summaries_merged_with_argcounts_dedup <-
  fastp_summaries_merged_with_argcounts_dedup %>%
  mutate(NonArgReads = summary.after_filtering.total_reads - AllArgs)

#Create a df that has the same sample info as the main one, but for NON-ARGs
#in each sample:
total_non_args_per_sample <- both_fastp_summaries_merged_with_argcounts %>%
  select(NumReads = NonArgReads, LIMS_ID, Treatment, Enrichment, Fraction,
         Nanotrap_type, site, UniqueID)
total_non_args_per_sample$ARG = "NonArgReads"
total_non_args_per_sample$predicted_ARG_class <- "None"

total_non_args_per_sample_dedup <- fastp_summaries_merged_with_argcounts_dedup %>%
  select(NumReads = NonArgReads, LIMS_ID, Treatment, Enrichment, Fraction,
         Nanotrap_type, site, UniqueID)
total_non_args_per_sample_dedup$ARG = "NonArgReads"
total_non_args_per_sample_dedup$predicted_ARG_class <- "None"


#Combine the ARG and non-ARG dfs, so each sample now has a row containing
#the count of non-ARGs:
combined_args_and_nonargs_per_sample <- rbind(samplewise_deeparg_results_table,
                                              total_non_args_per_sample)


combined_args_and_nonargs_per_sample_dedup <- rbind(
  samplewise_deeparg_results_table_dedup,
  total_non_args_per_sample_dedup)

#Sanity check: barplot should show majority of each sample as non-ARG reads:
ggplot(combined_args_and_nonargs_per_sample, aes(x = UniqueID, y = NumReads, 
                                                 fill = predicted_ARG_class)) + 
  geom_bar(position = 'fill', stat = "identity")

ggplot(combined_args_and_nonargs_per_sample_dedup, aes(x = UniqueID, y = NumReads, 
                                                 fill = predicted_ARG_class)) + 
  geom_bar(position = 'fill', stat = "identity")


#Save the table for later:
write_csv(combined_args_and_nonargs_per_sample,
          "imported_deeparg_reports/imported_deeparg_results.csv")

write_csv(combined_args_and_nonargs_per_sample_dedup,
          "imported_deeparg_reports/imported_deeparg_results_dedup.csv")

###Prepare for EdgeR
processed_deeparg_table <- read_csv("imported_deeparg_reports/imported_deeparg_results.csv")
deeparg_metadata <- processed_deeparg_table %>%
  select(UniqueID, LIMS_ID, Treatment, site, Fraction, Nanotrap_type, 
         Enrichment) %>%
  distinct()

processed_deeparg_table_dedup <- read_csv("imported_deeparg_reports/imported_deeparg_results_dedup.csv")
deeparg_metadata_dedup <- processed_deeparg_table_dedup %>%
  select(UniqueID, LIMS_ID, Treatment, site, Fraction, Nanotrap_type, 
         Enrichment) %>%
  distinct()


deeparg_count_matrix <- processed_deeparg_table %>%
  mutate(ARG_w_class = paste(ARG, predicted_ARG_class, sep = ":")) %>%
  pivot_wider(id_cols = ARG_w_class,
              names_from = UniqueID,
              values_from = NumReads,
              values_fill = 0)

deeparg_count_matrix_dedup <- processed_deeparg_table_dedup %>%
  mutate(ARG_w_class = paste(ARG, predicted_ARG_class, sep = ":")) %>%
  pivot_wider(id_cols = ARG_w_class,
              names_from = UniqueID,
              values_from = NumReads,
              values_fill = 0)

write_csv(deeparg_metadata, "imported_deeparg_reports/deeparg_metadata.csv")
write_csv(deeparg_metadata_dedup, "imported_deeparg_reports/deeparg_metadata_dedup.csv")
write_csv(deeparg_count_matrix, "imported_deeparg_reports/deeparg_count_matrix.csv")
write_csv(deeparg_count_matrix_dedup, "imported_deeparg_reports/deeparg_count_matrix_dedup.csv")

