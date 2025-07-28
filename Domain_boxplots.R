collapse_to_tax_level <- function(imported_table, level) {
  collapsed <- imported_table %>%
    group_by(!!sym(level), SampleID, LIMS_ID, Treatment, ribosomal, site) %>%
    summarise(readcount = sum(nodeOnly),
              RA = sum(RA))
  return(collapsed)
}

parse_sample_treatments <- function(tax_table, treatment_col = "Treatment") {
  parsed <- tax_table %>%
    mutate(Fraction = str_replace_all(
      !!sym(treatment_col),
      c("^F.$" = "filtrate",
        "^R.$" = "retentate",
        "^U.$" = "unfiltered")
    )) %>%
    mutate(Nanotrap_type = str_replace_all(
      !!sym(treatment_col),
      c("^.A$" = "A",
        "^.B$" = "A&B",
        "^.D$" = "none")
    ))
}

collapse_to_domain_no_rrna <- function(imported_report, treatment_name, conf_level) {
  domains <- collapse_to_tax_level(imported_report, "D")
  domains <- parse_sample_treatments(domains) 
  domains <- domains %>%
    mutate(Enrichment = treatment_name) %>%
    mutate(Kraken2_confidence = conf_level)
  domains <- filter(domains, ribosomal == "nonrrna")
}


collapse_to_domain_keep_rrna <- function(imported_report, treatment_name, conf_level) {
  domains <- collapse_to_tax_level(imported_report, "F")
  domains <- parse_sample_treatments(domains) 
  domains <- domains %>%
    mutate(Enrichment = treatment_name) %>%
    mutate(Kraken2_confidence = conf_level)
}


vsp_domains_no_rrna <- collapse_to_domain_no_rrna(vsp_imported, "VSP", 0)
rpip_domains_no_rrna <- collapse_to_domain_no_rrna(rpip_imported, "RPIP", 0)
untargeted_domains_no_rrna <- collapse_to_domain_no_rrna(untargeted_imported, "None", 0)

all_domains_no_rrna <- bind_rows(vsp_domains_no_rrna,
                                 rpip_domains_no_rrna,
                                 untargeted_domains_no_rrna)

all_domains_no_rrna[str_starts(all_domains_no_rrna$D, "unident"),]$D <- "unclassified" 
all_domains_no_rrna$Label_Enrich <- str_replace(all_domains_no_rrna$Enrichment, "None", "No enrichment")

viruses_only <- filter(all_domains_no_rrna, D == "Viruses")

ggplot(filter(all_domains_no_rrna, D == "Viruses"), aes(y = RA, x = Fraction, fill = Nanotrap_type)) + geom_boxplot() +
#  scale_y_log10() +
  facet_grid(cols = vars(Enrichment)) +
  theme_minimal() +
  ylim(0, 0.3)


###Total relative abundances of virus reads (pooled)
ggplot(viruses_only, aes(x = Fraction, y = RA, fill = Nanotrap_type, colour = Nanotrap_type)) + geom_boxplot() +
  scale_fill_manual(name = "Nanotrap protocol", 
                    values = c("aquamarine2", "goldenrod2", "lightpink2")) +
  scale_colour_manual(name = "Nanotrap protocol",
                      values = c("aquamarine4", "goldenrod4", "lightpink4")) +
  ylab("Relative abundance of virus reads") +
  theme_minimal() +
  facet_wrap(~Label_Enrich, strip.position = "top") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(face = "bold", size = 15),
        axis.text = element_text(size = 15, angle = 45, hjust = 1),
        axis.title = element_text(size = 20),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 20))


##Total RA of ARG reads (pooled)
ggplot(arg_fraction_df, aes(x = Fraction, y = portion_of_args, fill = Nanotrap_type, colour = Nanotrap_type)) + geom_boxplot() +
  scale_fill_manual(name = "Nanotrap protocol", 
                    values = c("aquamarine2", "goldenrod2", "lightpink2")) +
  scale_colour_manual(name = "Nanotrap protocol",
                      values = c("aquamarine4", "goldenrod4", "lightpink4")) +
  ylab("Relative abundance of all ARGs") +
  theme_minimal() +
  facet_wrap(~Label_Enrich, strip.position = "top") +
  theme(strip.placement = "outside", 
        strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(face = "bold", size = 15),
        axis.text = element_text(size = 15, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 20),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 20))


#############33
#Plot VSP with grid
ggplot(filter(family_box_df_00, `F` %in% at_least_doubled_VSP$`F` & `F` %in% vsp_targets$family_name, Enrichment != "RPIP"), aes(y = `F`, x = RA, fill = Enrichment, colour = Enrichment)) + geom_boxplot() +
  scale_x_log10() +
  theme_minimal() +
  #  facet_grid(rows = "plotting_cat", scales = "free_y", space = "free_y") +
  xlab("Relative abundance of reads") +
  ylab("") +
  theme(strip.placement = "outside", 
        strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(face = "bold", size = 15),
        strip.text.y.right = element_text(angle = -90),
        axis.text = element_text(size = 12, angle = 30, hjust = 1),
        axis.title = element_text(size = 20),
        axis.ticks.x.bottom = element_line(color = "grey40"),
        legend.position = "bottom",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 12),
        panel.grid.major.y = element_blank()) +
  geom_hline(yintercept=seq(1,nrow(boxplot_df_by_arg[,1]),1)+.5,color="black", linetype = "dotted", linewidth = .5) +
  geom_hline(yintercept=seq(1,nrow(boxplot_df_by_arg[,1]),1)-.5,color="black", linetype = "dotted", linewidth = .5) +
  scale_fill_manual(name = "Enrichment:", values = c("lightpink3", "goldenrod1")) +
  scale_colour_manual(name = "Enrichment:", values = c("lightpink4", "goldenrod3")) +
  facet_grid(rows = vars(Fraction), cols = vars(Nanotrap_type))


#############33
#Plot RPIP with grid
ggplot(filter(family_box_df_00, `F` %in% at_least_doubled_RPIP$`F` & `F` %in% rpip_targets$family, Enrichment != "VSP"), aes(y = `F`, x = RA, fill = Enrichment, colour = Enrichment)) + geom_boxplot() +
  scale_x_log10() +
  theme_minimal() +
  xlab("Relative abundance of reads") +
  ylab("") +
  theme(strip.placement = "outside", 
        strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(face = "bold", size = 15),
        strip.text.y.right = element_text(angle = -90),
        axis.text = element_text(size = 12, angle = 30, hjust = 1),
        axis.title = element_text(size = 20),
        axis.ticks.x.bottom = element_line(color = "grey40"),
        legend.position = "bottom",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        panel.grid.major.y = element_blank()) +
  geom_hline(yintercept=seq(1,nrow(boxplot_df_by_arg[,1]),1)+.5,color="black", linetype = "dotted", linewidth = .5) +
  geom_hline(yintercept=seq(1,nrow(boxplot_df_by_arg[,1]),1)-.5,color="black", linetype = "dotted", linewidth = .5) +
  scale_fill_manual(name = "Enrichment:", values = c("lightpink3", "aquamarine3")) +
  scale_colour_manual(name = "Enrichment:", values = c("lightpink4", "aquamarine4")) +
  facet_grid(rows = vars(Fraction), cols = vars(Nanotrap_type))



###########RPIP Virus only
at_least_doubled_RPIP <- filter(rpip_all_taxa_wider_enrichment_00_mean_ratios, mean_rpip_ratio >= 2)

RPIP_top_viruses <- filter(rpip_all_taxa_wider_enrichment_00_mean_ratios, `F` %in% RPIP_Viruses$family) %>%
  arrange(mean_rpip_ratio)

top_10_rpip_viruses <- 

RPIP_Viruses <- rpip_targets %>%
  filter(Domain == "Virus")
ggplot(filter(family_box_df_00, `F` %in% c('Coronaviridae', 'Picornaviridae') & `F` %in% rpip_targets$family, Enrichment != "VSP"), aes(y = `F`, x = RA, fill = Enrichment, colour = Enrichment)) + geom_boxplot() +
  scale_x_log10() +
  theme_minimal() +
  xlab("Relative abundance of reads") +
  ylab("") +
  theme(strip.placement = "outside", 
        strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(face = "bold", size = 15),
        strip.text.y.right = element_text(angle = -90, size = 10),
        axis.text = element_text(size = 12, angle = 30, hjust = 1),
        axis.title = element_text(size = 20),
        axis.ticks.x.bottom = element_line(color = "grey40"),
        legend.position = "bottom",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        panel.grid.major.y = element_blank()) +
  geom_hline(yintercept=seq(1,nrow(boxplot_df_by_arg[,1]),1)+.5,color="black", linetype = "dotted", linewidth = .5) +
  geom_hline(yintercept=seq(1,nrow(boxplot_df_by_arg[,1]),1)-.5,color="black", linetype = "dotted", linewidth = .5) +
  scale_fill_manual(name = "Enrichment:", values = c("lightpink3", "aquamarine3")) +
  scale_colour_manual(name = "Enrichment:", values = c("lightpink4", "aquamarine4")) +
  facet_grid(rows = vars(Fraction), cols = vars(Nanotrap_type))



###########Everything
ggplot(filter(family_box_df_00, (`F` %in% at_least_doubled_VSP$`F` | `F` %in% at_least_doubled_RPIP$F) & (`F` %in% vsp_targets$family_name | `F` %in% rpip_targets$family)), aes(y = `F`, x = RA, fill = Enrichment, colour = Enrichment)) + geom_boxplot() +
  scale_x_log10() +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12, angle = 45)) +
  xlab("Relative abundance of reads") +
  ylab("") +
  theme( theme(strip.placement = "outside", 
        strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(face = "bold", size = 15),
        strip.text.y.right = element_text(angle = -90),
        axis.text = element_text(size = 12, angle = 30, hjust = 1),
        axis.title = element_text(size = 20),
        axis.ticks.x.bottom = element_line(color = "grey40"),
        legend.position = "bottom",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 12),
        panel.grid.major.y = element_blank())
  theme(panel.grid.major.y = element_blank()) +
  geom_hline(yintercept=seq(1,nrow(boxplot_df_by_arg[,1]),1)+.5,color="black", linetype = "dotted", linewidth = .5) +
  geom_hline(yintercept=seq(1,nrow(boxplot_df_by_arg[,1]),1)-.5,color="black", linetype = "dotted", linewidth = .5) +
  scale_fill_manual(name = "Enrichment:", values = c("lightpink3", "aquamarine3", "goldenrod1")) +
  scale_colour_manual(name = "Enrichment:", values = c("lightpink4", "aquamarine4", "goldenrod3")) +
  facet_grid(rows = vars(Fraction), cols = vars(Nanotrap_type))

#############ARGs pooled by category
ggplot(boxplot_df_by_category, aes(y = plotting_cat, x = ra_sum, colour = Enrichment, fill = Enrichment)) + geom_boxplot() +
  scale_x_log10() +
  theme_minimal() +
  theme(strip.placement = "outside", 
        strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(face = "bold", size = 15),
        strip.text.y.right = element_text(angle = -90),
        axis.text = element_text(size = 12, angle = 30, hjust = 1),
        axis.title = element_text(size = 20),
        axis.ticks.x.bottom = element_line(color = "grey40"),
        legend.position = "bottom",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 12),
        panel.grid.major.y = element_blank()) +
  #  facet_grid(rows = "plotting_cat", scales = "free_y", space = "free_y") +
  xlab("Relative abundance of reads") +
  ylab("") +
  theme(legend.position = "bottom") +
  theme(panel.grid.major.y = element_blank()) +
  geom_hline(yintercept=seq(1,nrow(boxplot_df_by_arg[,1]),1)+.5,color="black", linetype = "dotted", linewidth = .5) +
  geom_hline(yintercept=seq(1,nrow(boxplot_df_by_arg[,1]),1)-.5,color="black", linetype = "dotted", linewidth = .5) +
  scale_color_manual(name = "Enrichment:", values = c("lightpink4", "aquamarine4")) +
  scale_fill_manual(name = "Enrichment:", values = c("lightpink2", "aquamarine2"))


#######ARGs only in RPIP
ggplot(ARGs_detected_only_w_rpip_25plus_any_sample, aes(x = RA, y = ARG)) + geom_boxplot(fill = "aquamarine2", colour = "aquamarine4") +
  scale_x_log10() +
  xlab("Relative abundance of reads") +
  ylab("") +
  theme_bw() +
  theme(strip.placement = "outside", 
        strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(face = "bold", size = 8),
        strip.text.y.right = element_text(angle = -90),
        axis.text = element_text(size = 12, angle = 30, hjust = 1),
        axis.title = element_text(size = 20),
        axis.ticks.x.bottom = element_line(color = "grey40"),
        legend.position = "bottom",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 12),
        panel.grid.major.y = element_line(color = "grey85"),
        panel.border = element_rect(linetype = "dotted", linewidth = 0.5),
        panel.grid.minor.y = element_blank()) +
  facet_grid(rows = "plotting_cat", scales = "free_y", space = "free_y")

