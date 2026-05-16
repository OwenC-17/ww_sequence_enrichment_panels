library(tidyverse)

#Data locations
rpip_rgi_dir <- paste0("input/link_to_raw_data/rpip_panels/raw_fastqs/",
                       "fastp_out_no_dedup/deduped/rgi_out/")
unt_rgi_dir <- paste0("input/link_to_raw_data/untargeted/raw_fastqs/",
                      "fastp_out_no_dedup/deduped/rgi_out/")


##############################
######Import Allele Data######
##############################

#Read file names:
rpip_allele_list <- list.files(rpip_rgi_dir, 
                               pattern = "*allele_mapping_data.txt",
                               full.names = TRUE)
unt_allele_list <- list.files(unt_rgi_dir,
                              pattern = "*allele_mapping_data.txt",
                              full.names = TRUE)

#Import each file:
rpip_allele_dfs <- lapply(rpip_allele_list, function(x) {
  read_tsv(x, col_types = "cccccccccnnnnnncncccnccc")
}
)

unt_allele_dfs <- lapply(unt_allele_list, function(x) {
  read_tsv(x, col_types = "cccccccccnnnnnncncccnccc")
}
)

#Give each dataframe a label matching the imported filename:
names(rpip_allele_dfs) <- lapply(rpip_allele_list, basename)
names(unt_allele_dfs) <- lapply(unt_allele_list, basename)

#Join data from all samples into a big dataframe:
rpip_allele_big_df <- bind_rows(rpip_allele_dfs, .id = "filename")
unt_allele_big_df <- bind_rows(unt_allele_dfs, .id = "filename")

#Label enrichment as RPIP or not:
rpip_allele_big_df$Enrichment <- "RPIP"
unt_allele_big_df$Enrichment <- "Non-targeted"

#Extract the sample IDs from the file names:
rpip_allele_big_df$sample_id <- substr(rpip_allele_big_df$filename, 5, 12)
unt_allele_big_df$sample_id <- substr(unt_allele_big_df$filename, 5, 12)

#Join the RPIP and unenriched samples into one dataframe:
rpipAndUnt_allele_big_df <- bind_rows(rpip_allele_big_df, unt_allele_big_df)

####An exploratory plot (skip when sourcing)####
#ggplot(rpipAndUnt_allele_big_df, aes(x = `All Mapped Reads`,
#                                     y = `Percent Coverage`,
#                                     colour = Enrichment)) +
#  geom_point() +
#  scale_x_log10() +
#  geom_smooth() +
#  facet_wrap(~`Reference Model Type`)


##############################
######Join to fastp data######
##############################
#Read in fastp reports as formatted in read_fastp_reports.R
all_fastp_reports_dupdedup <- read_rds(
  "input/modified/all_fastp_reports_dupdedup.rds"
  )

#Create specific column to join to ARG table:
all_fastp_reports_dupdedup$sample_id <- paste(
  all_fastp_reports_dupdedup$LIMS_ID,
  all_fastp_reports_dupdedup$Treatment,
  sep = "-"
)

#Select columns that will be useful for ARG analyses:
fastp_info_table <- all_fastp_reports_dupdedup %>%
  select(summary.after_dedup.total_reads,
         summary.after_dedup.total_bases, 
         filename:sample_id)


#Join to ARG table:
rpip_and_unt_allele_and_info <- rpipAndUnt_allele_big_df %>%
  left_join(fastp_info_table, by = c("sample_id", "Enrichment")) %>%
  mutate(RA_nmapped = `All Mapped Reads` / summary.after_dedup.total_reads,
         UniqueID = paste(sample_id, Enrichment, sep = "-"))


#Write file for later
write_rds(rpip_and_unt_allele_and_info,
          "input/modified/rpip_and_unt_rgi_allele_and_info.rds")


#########################
#####Gene-level data#####
#########################
#Read file names:
rpip_gene_list <- list.files(rpip_rgi_dir,
                             pattern = "*gene_mapping_data.txt",
                             full.names = TRUE)
unt_gene_list <- list.files(unt_rgi_dir,
                            pattern = "*gene_mapping_data.txt",
                            full.names = TRUE)

#Read in data from file names:
rpip_gene_dfs <- lapply(rpip_gene_list, function(x) {
  read_tsv(x, col_types = "ccccnccccnnnnnnn------nccc")
}
)

unt_gene_dfs <- lapply(unt_gene_list, function(x) {
  read_tsv(x, col_types = "ccccnccccnnnnnnn------nccc")
}
)

#Label the imported dataframes with corresponding filenames:
names(rpip_gene_dfs) <- lapply(rpip_gene_list, basename)
names(unt_gene_dfs) <- lapply(unt_gene_list, basename)

#Join the lists of dataframes into large single dataframes:
rpip_gene_big_df <- bind_rows(rpip_gene_dfs, .id = "filename")
unt_gene_big_df <- bind_rows(unt_gene_dfs, .id = "filename")

#Label as RPIP or not:
rpip_gene_big_df$Enrichment <- "RPIP"
unt_gene_big_df$Enrichment <- "Non-targeted"

#Extract sample id from filename:
rpip_gene_big_df$sample_id <- substr(rpip_gene_big_df$filename, 5, 12)
unt_gene_big_df$sample_id <- substr(unt_gene_big_df$filename, 5, 12)

#Join RPIP and unenriched dataframes together:
rpipAndUnt_gene_big_df <- bind_rows(rpip_gene_big_df, unt_gene_big_df)

####An exploratory plot (skip when sourcing)####
#ggplot(rpipAndUnt_gene_big_df, aes(x = `All Mapped Reads`,
#                                   y = `Average Percent Coverage`,
#                                   colour = Enrichment)) +
#  geom_point() +
#  scale_x_log10() +
#  geom_smooth() +
#  facet_wrap(~`Reference Model Type`)

############################
#####Join to fastp data#####
############################
#The fastp_info_table df should be the same as the one created for allele
#results above

#Join fastp table to ARG table:
rpip_and_unt_gene_and_info <- rpipAndUnt_gene_big_df %>%
  left_join(fastp_info_table, by = c("sample_id", "Enrichment")) %>%
  mutate(RA_nmapped = `All Mapped Reads` / summary.after_dedup.total_reads,
         UniqueID = paste(sample_id, Enrichment, sep = "-"))

#Write file for later:
write_rds(rpip_and_unt_gene_and_info, "input/modified/rpip_and_unt_rgi_gene_and_info.rds")
