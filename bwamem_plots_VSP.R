library(tidyverse)
library(ggh4x)
library(paletteer)

###################
#####Load data#####
###################

vsp_and_unt_counts_covstats_and_info <- read_rds(
  "input/modified/vsp_and_unt_counts_covstats_and_info.rds"
  )

vspANDunt_big_pileup_df <- read_rds(
  "input/modified/vspANDunt_big_pileup_df.rds"
  )

########################################
#####Histogram of mapped base counts####
########################################

#Calculate how many bases of each reference are covered in each sample
vspunt_covered_bases <- vspANDunt_big_pileup_df %>%
  group_by(UniqueID, sample_id, rname, Enrichment) %>%
  filter(Depth >=1) %>%
  summarize(bases_covered = n())

#Include samples with 0 mapped bases
sample_id_df <- data.frame(sample_id = unique(vspunt_covered_bases$sample_id))
Enrichment_df <- data.frame(
  Enrichment = unique(vspunt_covered_bases$Enrichment)
  )
all_samples_df <- cross_join(sample_id_df, Enrichment_df) %>%
  mutate(UniqueID = paste(sample_id, Enrichment, sep = "-"))

#See how many bases in each sample map to any reference genome:
vspunt_basescovered_all_refs_combined <- vspunt_covered_bases %>%
  group_by(UniqueID, Enrichment) %>%
  summarize(totalRefBases = sum(bases_covered, na.rm = TRUE)) %>%
  right_join(all_samples_df, by = c("UniqueID", "Enrichment"))

#Remove NAs for plotting:
vspunt_basescovered_all_refs_combined[
  is.na(vspunt_basescovered_all_refs_combined$totalRefBases)
  , ]$totalRefBases <- 0

#Look at histogram of mapped bases
ggplot(filter(vspunt_basescovered_all_refs_combined, totalRefBases > 0), 
       aes(x = totalRefBases, fill = Enrichment)) + 
  geom_histogram(position = "stack", breaks = c(0, 2^(0:17)), colour = "NA") +
  theme_bw() +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 2),
                     breaks = c(0, 2^(0:17)),
                     labels = scales::label_parse()(c(0, paste0("2^", 0:17)))
  ) +
  
  scale_fill_manual(values = c("salmon", "goldenrod2")) +
  xlab("Total reference bases covered") +
  ylab("Sample count")


####################################
#####Barplots of mutation rates#####
####################################
#Join read count data to base count/mutation data:
#Find bases where mutations exist
vspunt_mutations_count <- vspANDunt_big_pileup_df %>%
  group_by(sample_id, rname, Enrichment) %>%
  filter(!Status %in% c("CONSERVED", "STAR")) %>%
  summarize(mutation_count = n())

vspANDunt_mutations_per_covered_base <- vspunt_covered_bases %>%
  left_join(vspunt_mutations_count) %>%
  mutate(mutations_per_covered_base = mutation_count / bases_covered)

vspANDunt_mutations_per_covered_base[
  is.na(vspANDunt_mutations_per_covered_base$mutation_count), 
  ]$mutation_count <- 0
vspANDunt_mutations_per_covered_base[
  is.na(vspANDunt_mutations_per_covered_base$mutations_per_covered_base),
]$mutations_per_covered_base <- 0

vspANDunt_mutations_per_covered_base <- vspANDunt_mutations_per_covered_base %>%
  left_join(vsp_and_unt_counts_covstats_and_info)
ggplot(vspANDunt_mutations_per_covered_base, 
       aes(x = species, y = mutations_per_covered_base, fill = Enrichment,
           colour = Enrichment)) +
  geom_boxplot() + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c("goldenrod2", "magenta3")) +
  scale_colour_manual(values = c("goldenrod4", "magenta4"))

ggplot(vspANDunt_mutations_per_covered_base, 
       aes(x = species, y = bases_covered, fill = Enrichment,
           colour = Enrichment)) +
  geom_boxplot() + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c("goldenrod2", "magenta3")) +
  scale_colour_manual(values = c("goldenrod4", "magenta4"))


########################################################
#####Bases covered vs. mutation counts: scatterplot#####
########################################################
taxa_read_counts_studywide <- vsp_and_unt_counts_covstats_and_info %>%
  ungroup() %>%
  group_by(Reference) %>%
  summarize(allReads = sum(nmapped)) %>%
  arrange(desc(allReads))

top25_refs <- taxa_read_counts_studywide$Reference[1:25]

taxa_most_common <- vsp_and_unt_counts_covstats_and_info %>%
  ungroup() %>%
  group_by(Reference) %>%
  filter(nmapped != 0) %>%
  summarize(nappearances = n()) %>%
  arrange(desc(nappearances))

top25_mostCommon <- taxa_most_common$Reference[1:25]

#Plot mutations, color palette 1:
ggplot(filter(vsp_and_unt_counts_covstats_and_info,
              Reference %in% top25_mostCommon),
       aes(x = bases_covered, y = mutation_count, colour = species)) +
  geom_point(alpha = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~Enrichment) +
  theme_bw() +
  scale_color_paletteer_d("ggsci::default_igv") +
  xlab("Reference bases covered") +
  ylab("Mutation count") +
  theme(strip.background = element_rect(fill = "white"))


#Alternative color palette:
ggplot(filter(vsp_and_unt_counts_covstats_and_info,
              Reference %in% top25_mostCommon),
       aes(x = bases_covered, y = mutation_count, colour = species)) +
  geom_point(alpha = .7) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~Enrichment) +
  theme_bw() +
  scale_color_paletteer_d("ggsci::springfield_simpsons") +
  xlab("Reference bases covered") +
  ylab("Mutation count") +
  theme(strip.background = element_rect(fill = "white"))





