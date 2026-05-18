library(tidyverse)
library(ggh4x) #For nested faceting
library(paletteer)
library(cowplot)
library(Cairo) #For PDF export fidelity

############################
#####Allele level plots#####
############################
#Import filtered RGI results (created in EdgeR_rgi_allelelevel.R)
rpip_and_unt_allele_protein_homolog <- read_rds(
  "input/modified/rpip_and_unt_rgi_allele_protein_homolog.rds"
  )
  
RPIP_results_allele <- read_rds(
  "input/modified/edgeR_rgi_alleles_notgrouped_results.rds"
  )

#Import target info derived from RPIP documentation:
rpip_amr_targets <- read_csv("input/rpip_amr_targets.csv")
rpip_amr_targets <- rpip_amr_targets %>%
  mutate(matching_col = str_remove(antimicrobial_name,
                                   "[:space:]\\(.*generation\\)$")) %>%
  mutate(matching_col = tolower(str_remove(matching_col,
                                           "(s$)|([:space:]antibiotics)")))

#Make a string to match all RPIP-targeted drug categories:
drug_match_string <- paste(rpip_amr_targets$matching_col, collapse = "|")

#Format variables for display and add add RPIP targeting info:
rpip_allele_df <- rpip_and_unt_allele_protein_homolog %>% 
  mutate(Fraction = str_replace_all(Fraction,
                                    c("filtrate" = "Fil",
                                      "retentate" = "Ret",
                                      "unfiltered" = "Unfil")),
         Nanotrap_type = str_replace_all(Nanotrap_type,
                                         c("^A$" = "NT-A",
                                           "A&B" = "NT-AB",
                                           "none" = "DirEx")),
         Enrichment = str_replace_all(Enrichment,
                                      c("None" = "Non-targeted"))
         ) %>%
  mutate(drug_is_targeted = str_detect(`Drug Class`, drug_match_string)) %>%
  mutate(category = case_when(
    str_detect(tolower(`AMR Gene Family`), "beta-lactam") ~
      "beta-lactamase",
    str_detect(`Drug Class`, ";", negate = TRUE) ~
      str_remove(`Drug Class`, "antibiotic"),
    str_detect(`Drug Class`, ";") &
      str_detect(tolower(`Resistance Mechanism`), "efflux") ~
      "efflux pump (multivalent)",
    .default = "Other multivalent resistance genes"
  ))


#########################
#####Faceted barplot#####
#########################
#Get every combination of sample and treatment
#This allows even spacing even though some combinations are not present
barplot_design <- expand_grid(
  LIMS_ID = unique(rpip_allele_df$LIMS_ID),
  Fraction = unique(rpip_allele_df$Fraction),
  Nanotrap_type = unique(rpip_allele_df$Nanotrap_type),
  Enrichment = unique(rpip_allele_df$Enrichment)
)

#Join to the actual data:
barplot_rpip_allele_full <- barplot_design %>%
  left_join(rpip_allele_df) %>%
  mutate(missing = is.na(RA_nmapped))
barplot_rpip_allele_full[
  is.na(barplot_rpip_allele_full$drug_is_targeted),
    ]$drug_is_targeted <- FALSE

#Make a character label for drug targeting:
barplot_rpip_allele_full <- barplot_rpip_allele_full %>%
  mutate(c_drug_is_targeted = if_else(drug_is_targeted,
                                      "RPIP targets",
                                      "Not RPIP targets"
                                      )
         )


#Make a stacked barplot of RGI-mapped reads, colored by resistance category:  
ggplot(barplot_rpip_allele_full, aes(x = LIMS_ID, y = RA_nmapped, fill = category)) + 
  geom_col(data = barplot_rpip_allele_full %>% filter(!missing), colour = NA) +
  geom_text(data = barplot_rpip_allele_full %>% filter(missing),
            aes(label = "\u2020", y = 0.0125), size = 3) +
  facet_nested(Enrichment ~ Nanotrap_type + Fraction) +
  ylab("Relative abundance (mapped reads)") +
  xlab("Sample ID") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        strip.background = element_rect(fill = "white")) +
  scale_fill_paletteer_d("ggsci::default_igv", na.translate = FALSE) +
  guides(fill = guide_legend(nrow = 4))


#Set position = fill for easier viewing of low-abundance samples
ggplot(barplot_rpip_allele_full, aes(x = LIMS_ID, y = RA_nmapped, 
                                     fill = category)) + 
  geom_col(data = barplot_rpip_allele_full %>% filter(!missing), colour = NA,
           position = "fill") +
  geom_text(data = barplot_rpip_allele_full %>% filter(missing),
            aes(label = "\u2020", y = 0.125), size = 3) +
  facet_nested(Enrichment ~ Nanotrap_type + Fraction) +
  ylab("Relative abundance") +
  xlab("Sample ID") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        strip.background = element_rect(fill = "white")) +
  scale_fill_paletteer_d("ggsci::default_igv", na.translate = FALSE) +
  guides(fill = guide_legend(nrow = 5))


#Indicate whether drug is targeted:
ggplot(barplot_rpip_allele_full, aes(x = LIMS_ID, y = RA_nmapped, 
                                     fill = category)) + 
  geom_col(data = barplot_rpip_allele_full %>% filter(!missing), colour = NA,
           position = "fill") +
  geom_text(data = barplot_rpip_allele_full %>% filter(missing),
            aes(label = "\u2020", y = 0.25), size = 3) +
  facet_nested(Enrichment + c_drug_is_targeted ~ Nanotrap_type + Fraction) +
  ylab("Relative abundance (mapped reads)") +
  xlab("Sample ID") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        strip.background = element_rect(fill = "white")) +
  scale_fill_paletteer_d("ggsci::default_igv", na.translate = FALSE) +
  guides(fill = guide_legend(nrow = 4))


######################################
#####Volcano plot colored by Gene#####
######################################
volcano_plot_df <- RPIP_results_allele %>%
  mutate(drug_is_targeted = str_detect(`Drug Class`, drug_match_string)) %>%
  mutate(category = case_when(
    str_detect(tolower(`AMR Gene Family`), "beta-lactam") ~
      "beta-lactamase",
    str_detect(`Drug Class`, ";", negate = TRUE) ~
      str_remove(`Drug Class`, "antibiotic"),
    str_detect(`Drug Class`, ";") &
      str_detect(tolower(`Resistance Mechanism`), "efflux") ~
      "efflux pump (multivalent)",
    .default = "Other multivalent resistance genes"
  )) %>%
  mutate(Fraction = str_replace_all(Fraction,
                                    c("Unf$" = "Unfil")),
         Nanotrap_type = str_replace_all(Nanotrap_type,
                                         c("NA$" = "NT-A",
                                           "NB$" = "NT-AB",
                                           "DE$" = "DirEx")),
         Enrichment = str_replace_all(Enrichment,
                                      c("None" = "Non-targeted"))
  ) %>% 
  mutate(Treatment = paste(Fraction, Nanotrap_type, sep = "/")
         )
  
#Main plot (no legend):
volcano_plot <- ggplot(volcano_plot_df, 
                       aes(x = logFC, y = -log10(FDR), colour = category,
                           shape = drug_is_targeted)) + 
  geom_point(alpha = 0.67, size = 2) +
  geom_hline(yintercept = -log10(0.05)) +
  theme_bw() +
  scale_color_paletteer_d("ggsci::default_igv") +
  scale_shape_manual(values = c(4, 6)) +
  theme(legend.position = "none") +
  scale_x_continuous(position = "top", limits = c(-15, 15)) +
  ylim(0, 10) +
  xlab(expression(log[2](FC))) +
  ylab(expression(-log[10](FDR)))
volcano_plot

#Make a version with a legend:
volcano_plot_legend <- ggplot(volcano_plot_df, 
                          aes(x = logFC, y = -log10(FDR), colour = category,
                              shape = drug_is_targeted)) + 
  geom_point(alpha = 0.67, size = 2) +
  geom_hline(yintercept = -log10(0.05)) +
  theme_bw() +
  scale_color_paletteer_d("ggsci::default_igv") +
  scale_shape_manual(values = c(4, 6)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.key.spacing.y = unit(0, "points"),
        legend.box.spacing = unit(0, "points"),
        legend.margin = margin(0,0,0,-2)) +
  scale_x_continuous(position = "top", limits = c(-15, 15)) +
  ylim(0, 10) +
  guides(color = guide_legend(ncol = 4))

#Extract just the legend:
allele_legend <- as_grob(get_legend(volcano_plot_legend))
allele_legend


#Marginal density plot of LogFC:
FC_dens <- ggplot(volcano_plot_df, aes(x = logFC, y = Treatment, fill = Treatment)) +
  ggridges::geom_density_ridges(alpha = 0.5, scale = 1.25, rel_min_height = 0.001) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  scale_fill_paletteer_d("RColorBrewer::YlGnBu") +
  xlim(-15, 15)
FC_dens

#Marginal density plot of FDR:
sig_dens <- ggplot(volcano_plot_df, aes(x = -log10(FDR), y = Treatment, fill = Treatment)) +
  ggridges::geom_density_ridges2(alpha = 0.5, scale = -1.25, rel_min_height = 0.005) +
  xlim(0, 10) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        axis.text.y = element_blank()) +
  scale_fill_paletteer_d("RColorBrewer::YlGnBu") +
  scale_y_discrete(limits = rev, position = "right") +
  coord_flip()
sig_dens

#Volcano plot aligned with FDR density plot:
v_and_y <- plot_grid(
  volcano_plot, sig_dens,
  ncol = 2,
  rel_widths = c(3, 1),
  align = "h",
  axis = "tblr"
)

v_and_y

#Add legend (on top, full width) + 
#LogFC density plot (bottom, aligned with volcano plot axis):
v_and_y_and_x <- plot_grid(
  allele_legend, v_and_y, plot_grid(FC_dens, NULL, ncol = 2, rel_widths = c(3,1)),
  nrow = 3,
  rel_heights = c(1,3,1),
  align = "v",
  axis = "tblr"
)

#Final plot:
v_and_y_and_x

###############################
#####Bubble-style heatmaps#####
###############################
summarized_rpip_alleles <- RPIP_results_allele %>%
  filter(FDR <= 0.05) %>%
  mutate(category = case_when(
    str_detect(tolower(`AMR Gene Family`), "beta-lactam") ~
      "beta-lactamase",
    str_detect(`Drug Class`, ";", negate = TRUE) ~
      str_remove(`Drug Class`, "antibiotic"),
    str_detect(`Drug Class`, ";") &
      str_detect(tolower(`Resistance Mechanism`), "efflux") ~
      "efflux pump (multivalent)",
    .default = "Other multivalent resistance genes"
  )) %>%
  group_by(Nanotrap_type, Fraction, `ARO Term`, category) %>%
  summarize(mean_logFC = log2(mean(2^logFC)),
            numGenes = n_distinct(reference)) %>%
  filter(mean_logFC > 2)

ggplot(summarized_rpip_alleles, aes(x = `ARO Term`, y = Fraction,
                                    fill = mean_logFC, size = numGenes)) +
  geom_point(shape = 21, stroke = 1) +
  #  scale_radius(transform = "log10") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(rows = "Nanotrap_type") +
  scale_fill_viridis_c(option = "H")





################################################################################
######################
#####Gene level#######
######################

#Run EdgeR_rgi_genelevel.R first

#Import filtered RGI results (created in EdgeR_rgi_genelevel.R)
rpip_and_unt_gene_protein_homolog <- read_rds(
  "input/modified/rpip_and_unt_rgi_gene_protein_homolog.rds"
  )
RPIP_results_gene_groupedbyARO <- read_rds(
  "input/modified/edgeR_rgi_genesgroupedbyARO_results.rds"
  )


#Import target info derived from RPIP documentation:
if(!exists("rpip_amr_targets")) {
 rpip_amr_targets <- read_csv("input/rpip_amr_targets.csv")
 rpip_amr_targets <- rpip_amr_targets %>%
   mutate(matching_col = str_remove(antimicrobial_name,
                                   "[:space:]\\(.*generation\\)$")) %>%
   mutate(matching_col = tolower(str_remove(matching_col,
                                           "(s$)|([:space:]antibiotics)"))
          )

  #Make a string to match all RPIP-targeted drug categories:
  drug_match_string <- paste(rpip_amr_targets$matching_col, collapse = "|")
}

#Format variables for display and add add RPIP targeting info:
rpip_gene_df <- rpip_and_unt_gene_protein_homolog %>% 
  mutate(Fraction = str_replace_all(Fraction,
                                    c("filtrate" = "Fil",
                                      "retentate" = "Ret",
                                      "unfiltered" = "Unfil")),
         Nanotrap_type = str_replace_all(Nanotrap_type,
                                         c("^A$" = "NT-A",
                                           "A&B" = "NT-AB",
                                           "none" = "DirEx")),
         Enrichment = str_replace_all(Enrichment,
                                      c("None" = "Non-targeted"))
  ) %>%
  mutate(drug_is_targeted = str_detect(`Drug Class`, drug_match_string)) %>%
  mutate(category = case_when(
    str_detect(tolower(`AMR Gene Family`), "beta-lactam") ~
      "beta-lactamase",
    str_detect(`Drug Class`, ";", negate = TRUE) ~
      str_remove(`Drug Class`, "antibiotic"),
    str_detect(`Drug Class`, ";") &
      str_detect(tolower(`Resistance Mechanism`), "efflux") ~
      "efflux pump (multivalent)",
    .default = "Other multivalent resistance genes"
    )
  )

#########################
#####Faceted barplot#####
#########################

#Get every combination of sample and treatment
#This allows even spacing even though some combinations are not present
barplot_design <- expand_grid(
  LIMS_ID = unique(rpip_gene_df$LIMS_ID),
  Fraction = unique(rpip_gene_df$Fraction),
  Nanotrap_type = unique(rpip_gene_df$Nanotrap_type),
  Enrichment = unique(rpip_gene_df$Enrichment)
)

#Join to the actual data:
barplot_rpip_gene_full <- barplot_design %>%
  left_join(rpip_gene_df) %>%
  mutate(missing = is.na(RA_nmapped))

barplot_rpip_gene_full[
  is.na(barplot_rpip_gene_full$drug_is_targeted),
  ]$drug_is_targeted <- FALSE
  
barplot_rpip_gene_full <- barplot_rpip_gene_full %>%
  mutate(c_drug_is_targeted = if_else(drug_is_targeted,
                                      "RPIP targets",
                                      "Not RPIP targets"
  ))

#Make a stacked barplot of RGI-mapped reads, colored by resistance category:  
ggplot(barplot_rpip_gene_full, aes(x = LIMS_ID, y = RA_nmapped, fill = category)) + 
  geom_col(data = barplot_rpip_gene_full %>% filter(!missing), colour = NA) +
  geom_text(data = barplot_rpip_gene_full %>% filter(missing),
            aes(label = "\u2020", y = 0.0125), size = 3) +
  facet_nested(Enrichment ~ Nanotrap_type + Fraction) +
  ylab("Relative abundance (mapped reads)") +
  xlab("Sample ID") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        strip.background = element_rect(fill = "white")) +
  scale_fill_paletteer_d("ggsci::default_igv", na.translate = FALSE) +
  guides(fill = guide_legend(nrow = 4))


#Set position = fill for easier viewing of low-abundance samples
ggplot(barplot_rpip_gene_full, aes(x = LIMS_ID, y = RA_nmapped, 
                                     fill = category)) + 
  geom_col(data = barplot_rpip_gene_full %>% filter(!missing), colour = NA,
           position = "fill") +
  geom_text(data = barplot_rpip_gene_full %>% filter(missing),
            aes(label = "\u2020", y = 0.125), size = 3) +
  facet_nested(Enrichment ~ Nanotrap_type + Fraction) +
  ylab("Relative abundance") +
  xlab("Sample ID") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        strip.background = element_rect(fill = "white")) +
  scale_fill_paletteer_d("ggsci::default_igv", na.translate = FALSE) +
  guides(fill = guide_legend(nrow = 5))


#Indicate whether drug is targeted:
ggplot(barplot_rpip_gene_full, aes(x = LIMS_ID, y = RA_nmapped, 
                                     fill = category)) + 
  geom_col(data = barplot_rpip_gene_full %>% filter(!missing), colour = NA,
           position = "fill") +
  geom_text(data = barplot_rpip_gene_full %>% filter(missing),
            aes(label = "\u2020", y = 0.25), size = 3) +
  facet_nested(Enrichment + c_drug_is_targeted ~ Nanotrap_type + Fraction) +
  ylab("Relative abundance (mapped reads)") +
  xlab("Sample ID") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        strip.background = element_rect(fill = "white")) +
  scale_fill_paletteer_d("ggsci::default_igv", na.translate = FALSE) +
  guides(fill = guide_legend(nrow = 4))


##########################################
#####Volcano plot colored by category#####
##########################################

volcano_plot_df <- RPIP_results_gene_groupedbyARO %>%
  mutate(drug_is_targeted = str_detect(`Drug Class`, drug_match_string)) %>%
  mutate(category = case_when(
    str_detect(tolower(`AMR Gene Family`), "beta-lactam") ~
      "beta-lactamase",
    str_detect(`Drug Class`, ";", negate = TRUE) ~
      str_remove(`Drug Class`, "antibiotic"),
    str_detect(`Drug Class`, ";") &
      str_detect(tolower(`Resistance Mechanism`), "efflux") ~
      "efflux pump (multivalent)",
    .default = "Other multivalent resistance genes"
  )) %>%
  mutate(Fraction = str_replace_all(Fraction,
                                    c("Unf$" = "Unfil")),
         Nanotrap_type = str_replace_all(Nanotrap_type,
                                         c("NA$" = "NT-A",
                                           "NB$" = "NT-AB",
                                           "DE$" = "DirEx")),
         Enrichment = str_replace_all(Enrichment,
                                      c("None" = "Non-targeted"))
  ) %>% 
  mutate(Treatment = paste(Fraction, Nanotrap_type, sep = "/")
  )

#Main plot (no legend):
volcano_plot <- ggplot(volcano_plot_df, 
                       aes(x = logFC, y = -log10(FDR), colour = category,
                           shape = drug_is_targeted)) + 
  geom_point(alpha = 0.67, size = 2) +
  geom_hline(yintercept = -log10(0.05)) +
  theme_bw() +
  scale_color_paletteer_d("ggsci::default_igv") +
  scale_shape_manual(values = c(4, 6)) +
  theme(legend.position = "none") +
  scale_x_continuous(position = "top", limits = c(-15, 15)) +
  ylim(0, 10)
volcano_plot

#Make a version with a legend:
volcano_plot_legend <- ggplot(volcano_plot_df, 
                              aes(x = logFC, y = -log10(FDR), colour = category,
                                  shape = drug_is_targeted)) + 
  geom_point(alpha = 0.67, size = 2) +
  geom_hline(yintercept = -log10(0.05)) +
  theme_bw() +
  scale_color_paletteer_d("ggsci::default_igv") +
  scale_shape_manual(values = c(4, 6)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.key.spacing.y = unit(0, "points"),
        legend.box.spacing = unit(0, "points"),
        legend.margin = margin(0,0,0,-2)) +
  scale_x_continuous(position = "top", limits = c(-15, 15)) +
  ylim(0, 10) +
  guides(color = guide_legend(ncol = 4))

#Extract just the legend:
gene_legend <- as_grob(get_legend(volcano_plot_legend))
gene_legend


#Marginal density plot of LogFC:
FC_dens <- ggplot(volcano_plot_df, aes(x = logFC, y = Treatment, 
                                       fill = Treatment)) +
  ggridges::geom_density_ridges(alpha = 0.5, scale = 2, 
                                rel_min_height = 0.001) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  scale_fill_paletteer_d("RColorBrewer::YlGnBu") +
  xlim(-15, 15)
FC_dens

#Marginal density plot of FDR:
sig_dens <- ggplot(volcano_plot_df, aes(x = -log10(FDR), y = Treatment,
                                        fill = Treatment)) +
  ggridges::geom_density_ridges2(alpha = 0.5, scale = -2,
                                 rel_min_height = 0.005) +
  xlim(0, 10) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        axis.text.y = element_blank()) +
  scale_fill_paletteer_d("RColorBrewer::YlGnBu") +
  scale_y_discrete(limits = rev, position = "right") +
  coord_flip()
sig_dens

#Volcano plot aligned with FDR density plot:
v_and_y <- plot_grid(
  volcano_plot, sig_dens,
  ncol = 2,
  rel_widths = c(3, 1),
  align = "h",
  axis = "tblr"
)

v_and_y

#Add legend (on top, full width) + 
#LogFC density plot (bottom, aligned with volcano plot axis):
v_and_y_and_x <- plot_grid(
  gene_legend, v_and_y, plot_grid(FC_dens, NULL, ncol = 2, rel_widths = c(3,1)),
  nrow = 3,
  rel_heights = c(1,3,1),
  align = "hv",
  axis = "tblr"
)

v_and_y_and_x


###############################
#####Bubble-style heatmaps#####
###############################
summarized_rpip_genes <- RPIP_results_gene_groupedbyARO %>%
  filter(FDR <= 0.05) %>%
  group_by(Nanotrap_type, Fraction, `AMR Gene Family`) %>%
  summarize(mean_logFC = log2(mean(2^logFC)),
            numGenes = n_distinct(reference))

ggplot(summarized_rpip_genes, aes(x = `AMR Gene Family`, y = Fraction, fill = mean_logFC, size = numGenes)) +
  geom_point(shape = 21, stroke = 1) +
  #  scale_radius(transform = "log10") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(rows = "Nanotrap_type") +
  scale_fill_viridis_c(option = "H")


################################
#####Plot depth of coverage#####
################################
coverage_df <- rpip_and_unt_gene_protein_homolog %>%
  mutate(average_paired_length = summary.after_dedup.read1_mean_length +
           summary.after_dedup.read2_mean_length -
           insert_size_peak_dedup) %>%
  mutate(covered_depth = (`All Mapped Reads` * average_paired_length) /
           `Average Length Coverage (bp)`,
         total_depth = (`All Mapped Reads` * average_paired_length) /
           `Average Length Coverage (bp)`) %>%
  mutate(Enrichment = str_replace(Enrichment, "None", "No enrichment"))


ggplot(coverage_df, aes(colour = Enrichment, x = `Average Percent Coverage`)) +
  geom_density() +
  theme_bw()

#Log10 plot
ggplot(coverage_df, 
       aes(x = `All Mapped Reads`, y = `Average Percent Coverage`, 
           colour = Enrichment)) +
  geom_point() +
  theme_minimal() +
  scale_x_log10() +
  scale_colour_manual(values = c("deeppink3", "royalblue3"))


#Linear plot
ggplot(coverage_df, aes(x = `All Mapped Reads`, y = `Average Percent Coverage`,
                        colour = Enrichment)) +
  geom_point() +
  theme_minimal() +
  scale_colour_manual(values = c("deeppink3", "royalblue3"))

#Faceted
ggplot(coverage_df, aes(x = `All Mapped Reads`, y = `Average Percent Coverage`,
                        colour = `Drug Class`)) +
  geom_point() +
  theme_bw() +
  facet_grid(~Enrichment) +
  scale_color_paletteer_d("ggsci::default_igv") +
  theme(legend.position = "none")

#Faceted log scale
ggplot(coverage_df, aes(x = `All Mapped Reads`, y = `Average Percent Coverage`,
                        colour = `Drug Class`)) +
  geom_point() +
  theme_bw() +
  scale_x_log10() +
  facet_grid(~Enrichment) +
  scale_color_paletteer_d("ggsci::default_igv") +
  theme(legend.position = "none")


#Get a count of each ARO Term:
ARO_counts <- coverage_df %>%
  group_by(`ARO Term`) %>%
  summarize(total = sum(`All Mapped Reads`)) %>%
  arrange(desc(total))

#Filter to 20 most common for visual clarity:
top_20 <- ARO_counts$`ARO Term`[1:20]

coverage_df_top20 <- coverage_df %>%
  filter(`ARO Term` %in% top_20)



#Plot top 20 (mapped read count vs. percent coverage) with no axis 
#transformation:
ggplot(coverage_df_top20, aes(x = `All Mapped Reads`,
                              y = `Average Percent Coverage`,
                              colour = `ARO Term`)) +
  geom_point() +
  theme_bw() +
  facet_grid(~Enrichment) +
  scale_color_paletteer_d("ggsci::default_igv") +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white")) 


#Plot top 20 (mapped read count vs. percent coverage) with log x axis:
ggplot(coverage_df_top20,
       aes(x = `All Mapped Reads`, y = `Average Percent Coverage`, 
           colour = `ARO Term`)) +
  geom_point() +
  theme_bw() +
  scale_x_log10() +
  facet_grid(~Enrichment) +
  scale_color_paletteer_d("ggsci::default_igv") +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white")) 

#Plot top 20 (mapped read count vs. percent coverage) with log x AND y axes:
ggplot(coverage_df_top20, 
       aes(x = `All Mapped Reads`, y = `Average Percent Coverage`, 
           colour = `ARO Term`)) +
  geom_point() +
  theme_bw() +
  scale_x_log10() +
  scale_y_log10() +
  facet_grid(~Enrichment) +
  scale_color_paletteer_d("ggsci::default_igv") +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white")) 



ggplot(filter(coverage_df_top20, `ARO Term` == "msrE"), 
       aes(x = `All Mapped Reads`, y = `Average Percent Coverage`,
           colour = UniqueID)) +
  geom_point() +
  theme_bw() +
  facet_grid(~Enrichment) +
  scale_color_paletteer_d("ggsci::default_igv") +
  theme(legend.position = "none")

#####################################
#####Plot alleles vs. read count#####
#####################################
allele_counts <- rpip_and_unt_allele_protein_homolog %>%
  group_by(`ARO Term`, `ARO Accession`, LIMS_ID, Fraction, Nanotrap_type,
           Enrichment) %>%
  summarize(total_reads = sum(`All Mapped Reads`), num_alleles = n()) %>%
  mutate(Enrichment = str_replace(Enrichment, "None", "No enrichment"))

#Read count vs. allele count; log-scaled x & y:
ggplot(filter(allele_counts, `ARO Term` %in% top_20),
       aes(x = total_reads, y = num_alleles, colour = `ARO Term`)) +
  geom_point(size = 1) +
  scale_x_log10() + 
  scale_y_log10() +
  scale_color_paletteer_d("ggsci::default_igv") +
  facet_grid(~Enrichment) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12)) +
  ylab("Allele count") +
  xlab("Read count") +
  guides(colour = guide_legend(nrow = 2))

#Read count vs. allele count; log-scaled x and linear y:
ggplot(filter(allele_counts, `ARO Term` %in% top_20),
       aes(x = total_reads, y = num_alleles,
           colour = `ARO Term`)) +
  geom_point() +
  scale_x_log10() + 
  scale_color_paletteer_d("ggsci::default_igv") +
  facet_grid(~Enrichment) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12)) +
  ylab("Allele count") +
  xlab("Read count") +
  guides(colour = guide_legend(nrow = 2))


#Read count vs. allele count; linear axes
ggplot(filter(allele_counts, `ARO Term` %in% top_20),
       aes(x = total_reads, y = num_alleles,
           colour = `ARO Term`)) +
  geom_point() +
  scale_color_paletteer_d("ggsci::default_igv") +
  facet_grid(~Enrichment) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12)) +
  ylab("Percent ref covered") +
  xlab("Read count") +
  guides(colour = guide_legend(nrow = 2))
