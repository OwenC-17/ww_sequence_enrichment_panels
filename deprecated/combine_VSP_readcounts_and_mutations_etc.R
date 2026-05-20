library(tidyverse)
library(taxonomizr)
library(ggh4x)
library(paletteer)


sqlPath = "/projects/bios_microbe/cowen20/ref_db/taxonomizr/accessionTaxa.sql"
vsp_dir <- "/projects/bios_microbe/cowen20/targeted_panels/vsp_panels/raw_fastqs/fastp_out_no_dedup/deduped/bwa_out/mapped_and_low_complexity_removed/sorted/filtered095/pileup13mq/"
unt_vsp_dir <- "/projects/bios_microbe/cowen20/targeted_panels/untargeted/raw_fastqs/fastp_out_no_dedup/deduped/bwa_out/mapped_and_low_complexity_removed/sorted/filtered095/pileup13mq/"


vsp_count_list <- list.files(vsp_dir, pattern = "*.tsv", full.names = TRUE)
unt_vsp_count_list <- list.files(unt_vsp_dir, pattern = "*.tsv", full.names = TRUE)
vsp_count_list <- list.files(vsp_dir, pattern = "*.tsv", full.names = TRUE)
unt_rpip_count_list <- list.files(unt_rpip_dir, pattern = "*.tsv", full.names = TRUE)

vsp_count_dfs <- lapply(vsp_count_list, function(x) {
  read_tsv(x, col_names = c("Reference", "rlength", "nmapped", "assigned_unmapped"))
}
)


unt_vsp_count_dfs <- lapply(unt_vsp_count_list, function(x) {
  read_tsv(x, col_names = c("Reference", "rlength", "nmapped", "assigned_unmapped"))
}
)


names(vsp_count_dfs) <- lapply(vsp_count_list, basename)
names(unt_vsp_count_dfs) <- lapply(unt_vsp_count_list, basename)

vsp_counts_bigdf <- bind_rows(vsp_count_dfs, .id = "filename")
unt_vsp_counts_bigdf <- bind_rows(unt_vsp_count_dfs, .id = "filename")

vsp_counts_bigdf$Enrichment <- "VSP"
unt_vsp_counts_bigdf$Enrichment <- "None"

vsp_and_unt_counts_bigdf <- bind_rows(vsp_counts_bigdf, unt_vsp_counts_bigdf) %>%
  filter(Reference != "*")

vsp_and_unt_counts_bigdf <- vsp_and_unt_counts_bigdf %>%
  mutate(sample_id = substr(filename, 1, 12)) %>%
  mutate(rname = substr(Reference, 1, 11))

#Get taxonomy of reference IDs
vspunt_taxa <- data.frame(rname = unique(vsp_and_unt_counts_bigdf$rname))
vspunt_taxtable <- add_taxonomy(vspunt_taxa)
vspunt_taxtable[is.na(vspunt_taxtable$domain),]$domain <- "Viruses"

vsp_and_unt_counts_bigdf_tax <- vsp_and_unt_counts_bigdf %>%
  left_join(vspunt_taxtable)

#rpipANDunt_big_pileup_df is produced in diversity_of_mapped_reads.R (this is pileups by base)
vspANDunt_big_pileup_df <- read_rds("input/modified/vspANDunt_big_pileup_df.rds")

#Filter to bases with at least 1 mapped read
vspunt_covered_bases <- vspANDunt_big_pileup_df %>%
  filter(Depth >= 1) %>%
  group_by(sample_id, rname, Enrichment) %>%
  summarize(bases_covered = n())

#Only references with at least 50 covered bases
vspunt_50plus <- vspunt_covered_bases %>%
  filter(bases_covered >= 50)


#Find bases where mutations exist
vspunt_mutations_count <- vspANDunt_big_pileup_df %>%
  group_by(sample_id, rname, Enrichment) %>%
  filter(!Status %in% c("CONSERVED", "STAR")) %>%
  summarize(mutation_count = n())


#Combine coverage and mutation data
vsp_and_unt_counts_and_covstats <- vsp_and_unt_counts_bigdf_tax %>%
  left_join(vspunt_covered_bases) %>%
  left_join(vspunt_mutations_count)

vsp_and_unt_counts_and_covstats <- vsp_and_unt_counts_bigdf_tax %>%
  left_join(vspunt_covered_bases) %>%
  left_join(vspunt_mutations_count) %>%
  #The only taxon with no family ID assigned is NC_116928.1, which has no mapped
  #reads anyway, so easiest to filter it out:
  filter(!is.na(family))

vsp_and_unt_counts_and_covstats[
  is.na(vsp_and_unt_counts_and_covstats$bases_covered), ]$bases_covered <- 0

vsp_and_unt_counts_and_covstats[
  is.na(vsp_and_unt_counts_and_covstats$mutation_count), ]$mutation_count <- 0


####Some exploratory plots####
ggplot(vsp_and_unt_counts_covstats_and_info, aes(x = nmapped, y = bases_covered, colour = family)) +
  geom_point() +
  scale_x_log10() +
  geom_smooth() +
  scale_y_log10()

ggplot(vsp_and_unt_counts_and_covstats, aes(x = nmapped, y = mutation_count,
                                            colour = family)) +
  geom_smooth() +
  scale_x_log10() +
  scale_y_log10()

#Plot mutations, color palette 1:
ggplot(filter(vsp_and_unt_counts_and_covstats, bases_covered >= 50), aes(x = bases_covered, y = mutation_count,
                                            colour = family)) +
  geom_point(alpha = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~Enrichment) +
  theme_bw() +
  scale_colour_manual(name = "Family", 
        values = paletteer::paletteer_c("grDevices::Zissou 1", 11)[c(
          1, 11, 2, 10, 3, 9, 4, 8, 5, 7, 6)]) +
  xlab("Reference bases covered") +
  ylab("Mutation count")


#Alternative color palette:
ggplot(filter(vsp_and_unt_counts_and_covstats, bases_covered >= 50), aes(x = bases_covered, y = mutation_count,
                                                                         colour = family)) +
  geom_point(alpha = .75) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~Enrichment) +
  theme_bw() +
  scale_colour_manual(name = "Family", 
                      values = paletteer::paletteer_d("colorBlindness::PairedColor12Steps")[c(
                        1, 11, 2, 10, 3, 9, 4, 8, 5, 7, 6)]) +
  xlab("Reference bases covered") +
  ylab("Mutation count")


##############################
fastp_info_table <- read_rds("input/modified/all_fastp_reports_dupdedup.rds") %>%
  mutate(UniqueID = paste(LIMS_ID, Treatment, Enrichment, sep = "-"))


vsp_and_unt_counts_covstats_and_info <- vsp_and_unt_counts_and_covstats %>%
  mutate(UniqueID = paste(LIMS_ID))
  left_join(fastp_info_table, by = c("UniqueID", "Enrichment")) %>%
  mutate(RA_nmapped = nmapped / summary.after_filtering.total_reads,
         RA_bases_covered = bases_covered / summary.after_filtering.total_bases)


write_rds(vsp_and_unt_counts_covstats_and_info, "input/modified/vsp_and_unt_counts_covstats_and_info.rds")
vsp_and_unt_counts_covstats_and_info <- read_rds("input/modified/vsp_and_unt_counts_covstats_and_info.rds")

####Some exploratory plots (RA)####
ggplot(vsp_and_unt_counts_covstats_and_info, aes(y = RA_nmapped, x = genus, colour = Enrichment)) +
  geom_boxplot() +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90))


ggplot(vsp_and_unt_counts_covstats_and_info, aes(x = genus, y = RA_bases_covered,
                                            colour = Enrichment)) +
  geom_point() +
  geom_jitter() +
  geom_boxplot() +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90)) 
  

##############################
ggplot(vsp_and_unt_counts_covstats_and_info, aes(x = phylum, y = RA_nmapped, 
                                     fill = Enrichment, colour = Enrichment)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_manual(values = c("goldenrod2", "magenta3")) +
  scale_colour_manual(values = c("goldenrod4", "magenta4"))

