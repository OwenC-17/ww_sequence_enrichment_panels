#Plots for data imported via bwamem_and_pileup_import_RPIP.R

library(tidyverse)
library(paletteer)
library(ggh4x)

rpip_and_unt_counts_covstats_and_info <- read_rds(
  "input/modified/rpip_and_unt_counts_covstats_and_info.rds"
)

##################################
#####Read counts vs. Coverage#####
##################################
rpipunt_50plus_by_taxid <- rpip_and_unt_counts_covstats_and_info %>% 
  group_by(Enrichment, UniqueID, Treatment, LIMS_ID, Fraction, Nanotrap_type,
           across(taxid:domain)) %>%
  summarize(bases_covered = sum(bases_covered),
            nmapped = sum(nmapped),
            total_mapped_bases = sum(total_mapped_bases),
            comb_rlength = sum(rlength),
            mutation_count = sum(mutation_count)) %>%
  mutate(prcov = bases_covered / comb_rlength) %>%
  filter(bases_covered >= 50)

#Number of mapped reads vs. reference coverage
ggplot(rpipunt_50plus_by_taxid, 
       aes(x = nmapped, y = prcov * 100, colour = Treatment)) +
  geom_point(alpha = 0.7, size = 1) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  scale_colour_paletteer_d("ggsci::springfield_simpsons") +
  facet_grid(rows = vars(Enrichment), cols = vars(domain), scales = "free_x") +
  theme(strip.background = element_rect(fill = "white")) +
  xlab("Read count") +
  ylab("Percent ref covered")


ggplot(rpipunt_50plus_by_taxid, 
       aes(x = nmapped, y = prcov * 100, colour = class)) +
  geom_point(alpha = 0.7, size = .7) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  scale_colour_paletteer_d("ggsci::default_igv") +
  facet_nested(Enrichment ~ domain) +
  theme(strip.background = element_rect(fill = "white")) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  xlab("Read count") +
  ylab("Percent ref covered")

##########################################
#####Show read count vs.bases covered#####
##########################################
#Calculate the weighted mean read lengths:
#(Dividing the total read count by a million prevents an integer overflow error)
wm <- weighted.mean(
 c(rpip_and_unt_counts_covstats_and_info$summary.after_dedup.read2_mean_length,
   rpip_and_unt_counts_covstats_and_info$summary.after_dedup.read1_mean_length),
 rep(rpip_and_unt_counts_covstats_and_info$summary.after_dedup.total_reads /
       1000000,
     2)
  )


#Generate a pseudo-data to draw a line corresponding to zero redundancy: 
line_df <- data.frame(
  x = seq(min(rpipunt_50plus_by_taxid$nmapped),
  max(rpipunt_50plus_by_taxid$nmapped),
  length.out = 200)
  ) %>%
  mutate(y = x * wm)
  

#Look at number of bases covered vs number of reads, facet by
#Enrichment and domain
ggplot(rpipunt_50plus_by_taxid, 
       aes(x = nmapped, y = bases_covered, colour = class)) +
  geom_point(alpha = 0.5, size = 0.7) +
  geom_line(data = line_df, aes(x = x, y = y), inherit.aes = FALSE,
            linetype = "dotted", colour = "grey30") +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_paletteer_d("ggsci::default_igv") +
  theme_bw() +
  facet_grid(Enrichment~domain) +
  theme(strip.background = element_rect(fill = "white")) +
  xlab("Read count") +
  ylab("Ref bases covered")

###Hex version:
ggplot(rpipunt_50plus_by_taxid, 
       aes(x = nmapped, y = bases_covered)) +
  geom_hex(aes(fill = after_stat(log10(count)))) +
  geom_line(data = line_df, aes(x = x, y = y), inherit.aes = FALSE,
            linetype = "21", colour = "grey25") +
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_paletteer_c("grDevices::Teal", direction = -1) +
  theme_bw() +
  facet_grid(Enrichment~domain) +
  theme(strip.background = element_rect(fill = "white"),
        legend.position = "none") +
  xlab("Read count") +
  ylab("Ref bases covered")


#########################
#####Mutation plots######
#########################
#Get most common families:
bact_countdf <- rpipunt_50plus_by_taxid %>%
  filter(domain == "Bacteria") %>%
  group_by(family) %>%
  summarize(total_reads = sum(nmapped)) %>%
  arrange(desc(total_reads))

viru_countdf <- rpipunt_50plus_by_taxid %>%
  filter(domain == "Viruses") %>%
  group_by(family) %>%
  summarize(total_reads = sum(nmapped)) %>%
  arrange(desc(total_reads))

fung_countdf <- rpipunt_50plus_by_taxid %>%
  filter(domain == "Eukaryota") %>%
  group_by(family) %>%
  summarize(total_reads = sum(nmapped)) %>%
  arrange(desc(total_reads))

#Select top 10 bacterial families and all viral/fungal families (there are <10
#total of each):
topfam <- c(bact_countdf$family[1:10], fung_countdf$family, viru_countdf$family)

#Plot mutations per covered base:
ggplot(
  filter(rpipunt_50plus_by_taxid,family %in% topfam), 
  aes(x = bases_covered, y = mutation_count, colour = genus)) +
  geom_point(size = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  facet_grid(Enrichment~domain) +
  scale_color_paletteer_d("ggsci::default_igv") +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(strip.background = element_rect(fill = "white")) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  xlab("Ref bases covered") +
  ylab("Unconserved bases")


#################################################################
#####Boxplots of RA, coverage, and  detection counts (virus)#####
#################################################################
#Make a boxplot df for convenience:
rpip_virus_boxdf <- rpipunt_50plus_by_taxid %>%
  left_join(rpip_and_unt_counts_covstats_and_info) %>%
  filter(domain == "Viruses") %>%
  filter(nmapped > 0)

#Find how many samples each species was detected in:
rpip_virus_detectcounts <- rpip_virus_boxdf %>%
  select(species, LIMS_ID, Enrichment, Nanotrap_type, Fraction) %>%
  distinct() %>%
  ungroup() %>% #Ungrouping is important!
  count(species, Enrichment, name = "n")

#Relative abundance boxplot:
rpip_virus_boxplot <- ggplot(
  rpip_virus_boxdf,
  aes(x = species, y = RA_nmapped, colour = Enrichment, fill = Enrichment)) +
  geom_boxplot(position = position_dodge(width = 1.2)) +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = c("#C3EDE7", "#E7CAA7")) +
  scale_colour_manual(values = c("#005C50", "#654500")) +
  ylab("Fraction of total reads")

rpip_virus_boxplot

#Percent ref coverage boxplot:
rpip_virus_coverage_boxplot <- ggplot(
  rpip_virus_boxdf,
  aes(x = species, y = prcov * 100, colour = Enrichment, fill = Enrichment)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Reference length covered (%)") +
  scale_fill_manual(values = c("#C3EDE7", "#E7CAA7")) +
  scale_colour_manual(values = c("#005C50", "#654500"))

rpip_virus_coverage_boxplot

#Detection count boxplot:
rpip_virus_count_col <- ggplot(
  rpip_virus_detectcounts,
  aes(x = species, y = n, colour = Enrichment, fill = Enrichment)) +
  geom_col(position = "dodge") +
  theme_bw() +
  ylim(0, 55) +
  geom_text(aes(label = n), vjust = -0.5, position = position_dodge(0.9),
            size = 3, show.legend = FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   face = "italic"),
        axis.title.x = element_blank()) +
  scale_fill_manual(values = c("#C3EDE7", "#E7CAA7")) +
  scale_colour_manual(values = c("#005C50", "#654500")) +
  ylab("Detection count")

rpip_virus_count_col

#Combine the boxplots:
cowplot::plot_grid(rpip_virus_boxplot, 
                   rpip_virus_coverage_boxplot,
                   rpip_virus_count_col,
                   nrow = 3,
                   align = "v")

#######################################
#####Blob-style heatmaps (viruses)#####
#######################################

#Get number of samples in each treatment category for calculating detection 
#ratios:
rpip_combo_count <- rpip_and_unt_counts_covstats_and_info %>%
  ungroup() %>%
  select(LIMS_ID, Enrichment, Nanotrap_type, Fraction) %>%
  distinct() %>%
  group_by(Enrichment, Nanotrap_type, Fraction) %>%
  summarize(nsamples = n())

#Get detect counts within each treatment category:
rpip_virus_detect_counts_by_category <- rpip_virus_boxdf %>%
  ungroup() %>%
  select(species, LIMS_ID, Enrichment, Nanotrap_type, Fraction) %>%
  distinct() %>%
  count(species, Enrichment, Nanotrap_type, Fraction, name = "n")

#Make a dataframe for the blob-style heatmap:
rpip_virus_heatblob_df <- rpip_virus_boxdf %>%
  group_by(LIMS_ID, Enrichment, Nanotrap_type, Fraction, species) %>%
  summarise(RA_nmapped = sum(RA_nmapped),
            prcov = sum(prcov), .groups = "drop") %>%
  group_by(Enrichment, Nanotrap_type, Fraction, species) %>%
  summarise(RA_nmapped = mean(RA_nmapped),
            prcov = mean(prcov)) %>%
  left_join(rpip_virus_detect_counts_by_category) %>%
  left_join(rpip_combo_count) %>%
  mutate(pos_rate = n / nsamples) %>%
  mutate(Nanotrap_type = str_replace_all(Nanotrap_type,
                                         c("none" = "DirEx", 
                                           "A$" = "NTM-A",
                                           "A&B" = "NTM-AB")),
         Fraction = str_replace_all(Fraction,
                                    c("unfiltered" = "Unf",
                                      "filtrate" = "Fil",
                                      "retentate" = "Ret")
  )) %>%
  filter(!is.na(RA_nmapped))



#Version 1 (filled by RA)
ggplot(rpip_virus_heatblob_df,
       aes(x = species, y = Fraction, fill = RA_nmapped)) +
  geom_point(shape = 21, stroke = 1, aes(size = pos_rate)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   face = "italic"),
        axis.title = element_blank(),
        strip.background = element_rect(fill = "white")) +
  facet_nested(Enrichment + Nanotrap_type ~ .)  +
  scale_fill_viridis_c(option = "F", trans = "log10",
                       name = "Mean relative\nabundance") +
  scale_size(name = "Detection\nrate")

#Version 2 (Filled by reference coverage)
ggplot(rpip_virus_heatblob_df,
       aes(x = species, y = Fraction, fill = prcov * 100)) +
  geom_point(shape = 21, stroke = 1, aes(size = pos_rate)) +
  #  scale_radius(transform = "log10") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   face = "italic"),
        axis.title = element_blank(),
        strip.background = element_rect(fill = "white")) +
  facet_nested(Enrichment + Nanotrap_type ~ .)  +
  scale_fill_viridis_c(option = "F", trans = "log10",
                       name = "Mean reference\ncoverage") +
  scale_size(name = "Detection\nrate")

#####UP TO HERE IS TESTED#####
#####Plot Bacteria RA
#This requires EdgeR to be run first 

top_quartile_logfc <- rpip_taxa_signif %>%
  filter(abs(logFC) > quantile(abs(logFC), 0.75))

top_half_logfc <- rpip_taxa_signif %>%
  filter(abs(logFC) > quantile(abs(logFC), 0.5))

#Make a data frame for box-plotting the bacteria species with greatest FCs:
rpip_bacteria_boxdf <- rpip_and_unt_counts_covstats_and_info_by_taxid %>%
  filter(domain == "Bacteria") %>%
  filter(nmapped > 0) %>%
  filter(species %in% top_half_logfc$Species)

#Count the number of samples each species was detected in:
rpip_bacteria_detectcounts <- rpip_bacteria_boxdf %>%
  select(genus, LIMS_ID, Enrichment, Nanotrap_type, Fraction) %>%
  distinct() %>%
  count(genus, Enrichment, name = "n")

rpip_bacteria_boxplot <- ggplot(rpip_bacteria_boxdf, aes(x = genus, y = RA_nmapped,
                                                         colour = Enrichment, fill = Enrichment)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = c("#C3EDE7", "#E7CAA7")) +
  scale_colour_manual(values = c("#005C50", "#654500")) +
  ylab("Fraction of total reads")


rpip_bacteria_boxplot

rpip_bacteria_coverage_boxplot <- ggplot(rpip_bacteria_boxdf, aes(x = genus, y = prcov * 100,
                                                                  colour = Enrichment,
                                                                  fill = Enrichment)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Reference length covered (%)") +
  scale_fill_manual(values = c("#C3EDE7", "#E7CAA7")) +
  scale_colour_manual(values = c("#005C50", "#654500"))

rpip_bacteria_coverage_boxplot

rpip_bacteria_count_col <- ggplot(rpip_bacteria_detectcounts, aes(x = genus, y = n, fill = Enrichment, colour = Enrichment)) +
  geom_col(position = "dodge") +
  theme_bw() +
  ylim(0, 55) +
  geom_text(aes(label = n), vjust = -0.5, position = position_dodge(0.9), size = 2.5,
            show.legend = FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        axis.title.x = element_blank()) +
  scale_fill_manual(values = c("#C3EDE7", "#E7CAA7")) +
  scale_colour_manual(values = c("#005C50", "#654500")) +
  ylab("Detection count")


cowplot::plot_grid(rpip_bacteria_boxplot, 
                   rpip_bacteria_coverage_boxplot, 
                   rpip_bacteria_count_col,
                   nrow = 3,
                   align = "v")


#Get number of samples in each treatment category:
rpip_combo_count <- rpip_bacteria_boxdf %>%
  select(LIMS_ID, Enrichment, Nanotrap_type, Fraction) %>%
  distinct() %>%
  group_by(Enrichment, Nanotrap_type, Fraction) %>%
  summarize(nsamples = n())

#Get detect counts within each treatment category:
rpip_bacteria_detect_counts_by_category <- rpip_bacteria_boxdf %>%
  select(species, LIMS_ID, Enrichment, Nanotrap_type, Fraction) %>%
  distinct() %>%
  count(species, Enrichment, Nanotrap_type, Fraction, name = "n")


#Make a dataframe for a blob-style heatmap:
rpip_bacteria_heatblob_df <- rpip_bacteria_boxdf %>%
  group_by(LIMS_ID, Enrichment, Nanotrap_type, Fraction, species) %>%
  summarise(RA_nmapped = sum(RA_nmapped),
            prcov = sum(prcov), .groups = "drop") %>%
  group_by(Enrichment, Nanotrap_type, Fraction, species) %>%
  summarise(RA_nmapped = mean(RA_nmapped),
            prcov = mean(prcov)) %>%
  left_join(rpip_bacteria_detect_counts_by_category) %>%
  left_join(rpip_combo_count) %>%
  mutate(pos_rate = n / nsamples) %>%
  mutate(Enrichment = str_replace(Enrichment, "None", "No enrichment"),
         Nanotrap_type = str_replace_all(Nanotrap_type, c("none" = "DE", "A$" = "NTM-A", "A&B" = "NTM-AB")),
         Fraction = str_replace_all(Fraction, c("unfiltered" = "Unf", "filtrate" = "Fil", "retentate" = "Ret")))

#Bubbles (filled by read RA)
ggplot(rpip_bacteria_heatblob_df, aes(x = species, y = Fraction, fill = RA_nmapped, size = pos_rate)) +
  geom_point(shape = 21, stroke = 1) +
  #  scale_radius(transform = "log10") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8, face = "italic"),
        axis.title = element_blank(), strip.background = element_rect(fill = "white")) +
  facet_nested(Enrichment + Nanotrap_type ~ .) +
  scale_fill_viridis_c(option = "A", trans = "log10", name = "Relative abundance\nof mapped reads") +
  scale_size(name = "Detection rate")

#Bubbles (filled by coverage)
ggplot(rpip_bacteria_heatblob_df, aes(x = species, y = Fraction, fill = prcov * 100, size = pos_rate)) +
  geom_point(shape = 21, stroke = 1) +
  #  scale_radius(transform = "log10") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8, face = "italic"),
        axis.title = element_blank(), strip.background = element_rect(fill = "white")) +
  facet_nested(Enrichment + Nanotrap_type ~ .) +
  scale_fill_viridis_c(option = "F", trans = "log10", name = "Mean % reference\ncovered") +
  scale_size(name = "Detection rate")


####
#####Plot Fungi RA
rpip_fungi_boxdf <- rpip_and_unt_counts_covstats_and_info_by_taxid %>%
  filter(domain == "Eukaryota") %>%
  filter(nmapped > 0) %>%
  filter(species %in% rpip_taxa_signif$Species)

#Get number of samples where each species detected in each treatment category:
rpip_fungi_detectcounts <- rpip_fungi_boxdf %>%
  select(species, LIMS_ID, Enrichment, Nanotrap_type, Fraction) %>%
  distinct() %>%
  count(species, Enrichment, name = "n")

#Make a boxplot:
rpip_fungi_boxplot <- ggplot(rpip_fungi_boxdf, aes(x = species, y = RA_nmapped,
                                                   colour = Enrichment,
                                                   fill = Enrichment)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_bw() +
  ylab("Fraction of total reads") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = c("#C3EDE7", "#E7CAA7")) +
  scale_colour_manual(values = c("#005C50", "#654500"))

rpip_fungi_boxplot



rpip_fungi_coverage_boxplot <- ggplot(rpip_fungi_boxdf, aes(x = species, y = prcov * 100,
                                                            colour = Enrichment,
                                                            fill = Enrichment)) +
  geom_boxplot() +
  scale_y_log10() +
  ylab("Reference length covered (%)") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = c("#C3EDE7", "#E7CAA7")) +
  scale_colour_manual(values = c("#005C50", "#654500"))

rpip_fungi_coverage_boxplot

rpip_fungi_count_col <- ggplot(rpip_fungi_detectcounts, aes(x = species, y = n, colour = Enrichment,
                                                            fill = Enrichment)) +
  geom_col(position = "dodge") +
  theme_bw() +
  ylim(0, 55) +
  geom_text(aes(label = n), vjust = -0.5, position = position_dodge(0.9), size = 3,
            show.legend = FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        axis.title.x = element_blank()) +
  scale_fill_manual(values = c("#C3EDE7", "#E7CAA7")) +
  scale_colour_manual(values = c("#005C50", "#654500")) +
  ylab("Detection count")

cowplot::plot_grid(rpip_fungi_boxplot, 
                   rpip_fungi_coverage_boxplot,
                   rpip_fungi_count_col, nrow = 3,
                   align = "v"
)


######Begin the bubbly heatmap

#Get a count of how many samples are in each category:
rpip_combo_count <-
  rpip_and_unt_counts_covstats_and_info_by_taxid %>%
  select(LIMS_ID, Enrichment, Nanotrap_type, Fraction) %>%
  distinct() %>%
  group_by(Enrichment, Nanotrap_type, Fraction) %>%
  summarize(nsamples = n())

#Detection count of each species within each treatment category:
rpip_fungi_detect_counts_by_category <- rpip_fungi_boxdf %>%
  select(species, LIMS_ID, Enrichment, Nanotrap_type, Fraction) %>%
  distinct() %>%
  count(species, Enrichment, Nanotrap_type, Fraction, name = "n")


rpip_fungi_heatblob_df <- rpip_fungi_boxdf %>%
  group_by(LIMS_ID, Enrichment, Nanotrap_type, Fraction, species) %>%
  summarise(RA_nmapped = sum(RA_nmapped), .groups = "drop",
            prcov = sum(bases_covered) / sum(rlength)) %>%
  group_by(Enrichment, Nanotrap_type, Fraction, species) %>%
  summarise(RA_nmapped = mean(RA_nmapped),
            prcov = mean(prcov)) %>%
  left_join(rpip_fungi_detect_counts_by_category) %>%
  left_join(rpip_combo_count) %>%
  mutate(pos_rate = n / nsamples) %>%
  mutate(Enrichment = str_replace(Enrichment, "None", "No enrichment"),
         Nanotrap_type = str_replace_all(Nanotrap_type, c("none" = "DE", "A$" = "NTA", "A&B" = "NTAB")),
         Fraction = str_replace_all(Fraction, c("unfiltered" = "Unf", "filtrate" = "Fil", "retentate" = "Ret")))


#Bubbles (filled by read RA)
ggplot(rpip_fungi_heatblob_df, aes(x = species, y = Fraction, fill = RA_nmapped, size = pos_rate)) +
  geom_point(shape = 21, stroke = 1) +
  #  scale_radius(transform = "log10") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8, face = "italic"),
        axis.title = element_blank()) +
  facet_nested(Enrichment + Nanotrap_type ~ .) +
  scale_fill_viridis_c(option = "A", trans = "log10", name = "Relative abundance\nof mapped reads") +
  scale_size(name = "Detection rate")

ggplot(rpip_fungi_heatblob_df, aes(x = species, y = Fraction, fill = prcov * 100, size = pos_rate)) +
  geom_point(shape = 21, stroke = 1) +
  #  scale_radius(transform = "log10") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        axis.title = element_blank(), strip.background = element_rect(fill = "white")) +
  facet_nested(Enrichment + Nanotrap_type ~ .) +
  scale_fill_viridis_c(option = "F", trans = "log10", name = "Mean % reference\ncovered") +
  scale_size(name = "Detection rate")

##############################

fungi_rpip_and_unt_counts_covstats_and_info <- rpip_and_unt_counts_covstats_and_info %>%
  filter(kingdom == "Fungi")

bacte_rpip_and_unt_counts_covstats_and_info <- rpip_and_unt_counts_covstats_and_info %>%
  filter(domain == "Bacteria")


##############################More basic barplots for reference

#Relative abundance for all fungal species:
ggplot(fungi_rpip_and_unt_counts_covstats_and_info, aes(x = species, y = RA_nmapped, 
                                                        fill = Enrichment, colour = Enrichment)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_manual(values = c("goldenrod2", "magenta3")) +
  scale_colour_manual(values = c("goldenrod4", "magenta4"))

#Relative abundance of all bacterial species:
ggplot(bacte_rpip_and_unt_counts_covstats_and_info, aes(x = family, y = RA_nmapped, 
                                                        fill = Nanotrap_type, colour = Nanotrap_type)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_manual(values = c("goldenrod2", "magenta3", "green4")) +
  scale_colour_manual(values = c("goldenrod4", "magenta4", "green"))

#Plot mutations, color palette 1:
ggplot(filter(rpip_and_unt_counts_covstats_and_info, bases_covered >= 50 & domain == "Eukaryota"), aes(x = bases_covered, y = mutation_count,
                                                                                                       colour = family)) +
  geom_point(alpha = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~Enrichment) +
  theme_bw() +
  scale_colour_manual(name = "Family", 
                      values = paletteer::paletteer_c("grDevices::Zissou 1", 23)[c(
                        1, 23, 2, 22, 3, 21, 4, 20, 5, 19, 6, 18, 7, 17, 
                        8, 16, 9, 15, 10, 14, 11, 13, 12)]) +
  xlab("Reference bases covered") +
  ylab("Mutation count")


#Plot mutations, color palette 1:
ggplot(filter(rpip_and_unt_counts_covstats_and_info, bases_covered >= 50 & domain == "Eukaryota"), aes(x = bases_covered, y = mutation_count,
                                                                                                       colour = species)) +
  geom_point(alpha = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~Enrichment) +
  theme_bw() +
  scale_colour_manual(name = "Family", 
                      values = paletteer::paletteer_c("grDevices::Zissou 1", 52)) +
  xlab("Reference bases covered") +
  ylab("Mutation count")

#Plot mutations, color palette 1 (fungi with significant FC):
ggplot(filter(rpip_fungi_boxdf, bases_covered >= 50 & domain == "Eukaryota"), aes(x = bases_covered, y = mutation_count,
                                                                                  colour = species)) +
  geom_point(alpha = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  scale_colour_manual(name = "Family", 
                      values = paletteer::paletteer_c("grDevices::Zissou 1", 23)[c(
                        1, 23, 2, 22, 3, 21, 4, 20, 5, 19, 6, 18, 7, 17, 
                        8, 16, 9, 15, 10, 14, 11, 13, 12)]) +
  xlab("Reference bases covered") +
  ylab("Mutation count") +
  facet_nested(Enrichment ~ Nanotrap_type + Fraction)


#Plot mutations, color palette 1 (all viruses):
ggplot(filter(rpip_virus_boxdf, bases_covered >= 50 & nmapped >= 10), aes(x = bases_covered, 
                                                                          y = mutation_count,
                                                                          colour = species)) +
  geom_point(alpha = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  scale_colour_manual(name = "Species", 
                      values = paletteer::paletteer_d("ggsci::default_igv")) +
  xlab("Reference bases covered") +
  ylab("Mutation count") +
  facet_wrap(~Enrichment)


ggplot(filter(vsp_boxdf, bases_covered >= 50 & nmapped >= 10), aes(x = bases_covered, 
                                                                   y = mutation_count,
                                                                   colour = genus)) +
  geom_point(alpha = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  scale_colour_manual(name = "Genus", 
                      values = paletteer::paletteer_d("ggsci::default_igv")) +
  xlab("Reference bases covered") +
  ylab("Mutation count") +
  facet_wrap(~Enrichment)

ggplot(filter(rpip_bacteria_boxdf, bases_covered >= 50 & nmapped >= 100), aes(x = bases_covered, 
                                                                              y = mutation_count,
                                                                              colour = genus)) +
  geom_point(alpha = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  scale_colour_manual(name = "Genus", 
                      values = paletteer::paletteer_d("ggsci::default_igv")) +
  xlab("Reference bases covered") +
  ylab("Mutation count") +
  facet_wrap(~Enrichment)

#########################
#####Faceted barplots#####
#########################
#Get every combination of sample and treatment
#This allows even spacing even though some combinations are not present
existing_rpip_combos <- rpip_bacteria_boxdf %>%
  select(LIMS_ID, Enrichment, Nanotrap_type, Fraction) %>%
  distinct() %>%
  mutate(exists = TRUE)


barplot_design <- expand_grid(
  LIMS_ID = unique(rpip_and_unt_counts_covstats_and_info_by_taxid$LIMS_ID),
  Fraction = unique(rpip_and_unt_counts_covstats_and_info_by_taxid$Fraction),
  Nanotrap_type = unique(rpip_and_unt_counts_covstats_and_info_by_taxid$Nanotrap_type),
  Enrichment = unique(rpip_and_unt_counts_covstats_and_info_by_taxid$Enrichment)
) %>%
  left_join(existing_rpip_combos) %>%
  mutate(missing = is.na(exists))

#################################
#####num. reads vs. coverage#####
#################################
rpip_and_unt_counts_covstats_and_info_by_taxid <- rpip_and_unt_counts_covstats_and_info_by_taxid %>%
  mutate(Enrichment = relevel(factor(Enrichment, ordered = FALSE), ref = "None")) 

ggplot(filter(rpip_and_unt_counts_covstats_and_info_by_taxid, bases_covered >= 50),
       aes(x = nmapped,
           y = bases_covered)) +
  geom_hex() +
  facet_grid(domain ~ Enrichment) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  scale_fill_viridis_c(option = "G", trans = "log10")

ggplot(filter(vsp_and_unt_counts_covstats_and_info, bases_covered >= 50),
       aes(x = nmapped,
           y = bases_covered)) +
  geom_hex() +
  facet_grid(domain ~ Enrichment) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  scale_fill_viridis_c(option = "H", trans = "log10")


ggplot(filter(rpip_and_unt_counts_covstats_and_info_by_taxid, bases_covered >= 50 & domain == "Bacteria"),
       aes(x = nmapped,
           y = bases_covered,
           colour = prcov)) +
  geom_point(size = 1, alpha = 0.5) +
  geom_abline(slope = 1, intercept = 2.0) +
  facet_grid(Fraction ~ kingdom) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  scale_color_viridis_c(option = "H", trans = "log10")


logtrans <- rpip_and_unt_counts_covstats_and_info_by_taxid %>%
  mutate(lnbases_covered = log(bases_covered), lnnmapped = log(nmapped), lnrlength = log(rlength))

covmodel1 <- glmmTMB(bases_covered / rlength  ~ nmapped, data = filter(logtrans,
                                                                       bases_covered >= 50),
                     family = beta_family()
)
sim_resid <- simulateResiduals(covmodel1, n = 100)
plot(sim_resid)
testDispersion(sim_resid)


###############################
#####Virus stacked barplot#####
###############################

#Make a dataframe for easier barplot visualization:
barplot_rpip_virus_full <- rpip_virus_boxdf %>%
  right_join(barplot_design) %>%
  mutate(Fraction = str_replace_all(Fraction, c("filtrate" = "Fil", "retentate" = "Ret", "unfiltered" = "Unfil")),
         Nanotrap_type = str_replace_all(Nanotrap_type, c("^A$" = "NTM-A", "A&B" = "NTM-AB", "none" = "DE")),
         Enrichment = str_replace_all(Enrichment, "None", "Non-targeted"))



barplot_rpip_no_bcov <- barplot_rpip_virus_full %>%
  filter(species != "Betacoronavirus gravedinis" | is.na(species))

#Make a stacked barplot of BWA-mapped reads, colored by resistance category:  
ggplot(barplot_rpip_no_bcov,
       aes(x = LIMS_ID, y = RA_nmapped, fill = species)) + 
  geom_col(data = barplot_rpip_no_bcov %>% filter(!missing), colour = NA) +
  geom_text(data = barplot_rpip_no_bcov %>% filter(missing),
            aes(label = "\u2020", y = 5e-06)) +
  facet_nested(Enrichment ~ Nanotrap_type + Fraction) +
  ylab("Relative abundance (mapped reads)") +
  xlab("Sample ID") +
  labs(fill = "Species: ") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.position = "bottom",
        legend.text = element_text(size = 8, face = "italic"),
        strip.background = element_rect(fill = "white")) +
  scale_fill_paletteer_d("ggsci::default_igv"[-2], na.translate = FALSE) +
  guides(fill = guide_legend(nrow = 4))


##################################
#####Bacteria stacked barplot#####
##################################
#Make a data frame for box-plotting the bacteria species with greatest FCs:
rpip_bacteria_bardf <- rpip_and_unt_counts_covstats_and_info_by_taxid %>%
  filter(domain == "Bacteria") %>%
  filter(nmapped > 0) %>%
  filter(species %in% rpip_taxa_signif$Species)


#Make a dataframe for easier barplot visualization:
barplot_rpip_bacteria_full <- rpip_bacteria_bardf %>%
  right_join(barplot_design) %>%
  mutate(Fraction = str_replace_all(Fraction, c("filtrate" = "Fil", "retentate" = "Ret", "unfiltered" = "Unfil")),
         Nanotrap_type = str_replace_all(Nanotrap_type, c("^A$" = "NTM-A", "A&B" = "NTM-AB", "none" = "DE")),
         Enrichment = str_replace_all(Enrichment, "None", "Non-targeted"))

#Make a stacked barplot of RGI-mapped reads, colored by family:  
ggplot(barplot_rpip_bacteria_full,
       aes(x = LIMS_ID, y = RA_nmapped, fill = family)) + 
  geom_col(data = barplot_rpip_bacteria_full %>% filter(!missing), colour = NA) +
  geom_text(data = barplot_rpip_bacteria_full %>% filter(missing),
            aes(label = "\u2020", y = 2.5e-03)) +
  facet_nested(Enrichment ~ Nanotrap_type + Fraction) +
  ylab("Relative abundance (mapped reads)") +
  xlab("Sample ID") +
  labs(fill = "Family:") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.position = "bottom",
        legend.text = element_text(size = 8),
        strip.background = element_rect(fill = "white")) +
  scale_fill_paletteer_d("ggsci::default_igv", na.translate = FALSE) +
  guides(fill = guide_legend(nrow = 5))



###############################
#####Fungi stacked barplot#####
###############################
#Make a data frame for box-plotting the bacteria species with greatest FCs:
rpip_fungi_bardf <- rpip_and_unt_counts_covstats_and_info_by_taxid %>%
  filter(domain == "Eukaryota") %>%
  filter(nmapped > 0)

#Make a dataframe for easier barplot visualization:
barplot_rpip_fungi_full <- rpip_fungi_bardf %>%
  right_join(barplot_design) %>%
  mutate(Fraction = str_replace_all(Fraction, c("filtrate" = "Fil", "retentate" = "Ret", "unfiltered" = "Unfil")),
         Nanotrap_type = str_replace_all(Nanotrap_type, c("^A$" = "NT-A", "A&B" = "NT-AB", "none" = "DirEx")),
         Enrichment = str_replace_all(Enrichment, "None", "Non-targeted"))

#Make a stacked barplot of RGI-mapped reads, colored by resistance category:  
ggplot(barplot_rpip_fungi_full,
       aes(x = LIMS_ID, y = RA_nmapped, fill = genus)) + 
  geom_col(data = barplot_rpip_fungi_full %>% filter(!missing), colour = NA) +
  geom_text(data = barplot_rpip_fungi_full %>% filter(missing),
            aes(label = "\u2020", y = 2e-04)) +
  facet_nested(Enrichment ~ Nanotrap_type + Fraction) +
  ylab("Relative abundance (mapped reads)") +
  xlab("Sample ID") +
  labs(fill = "Genus: ") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.position = "bottom",
        legend.text = element_text(size = 8, face = "italic"),
        strip.background = element_rect(fill = "white")) +
  scale_fill_paletteer_d("ggsci::default_igv", na.translate = FALSE) +
  guides(fill = guide_legend(nrow = 5))


