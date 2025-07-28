library(tidyverse)

family_box_df <- all_families
family_box_df_w_rrna <- all_families_w_rrna %>%
  group_by(LIMS_ID, Treatment, site, Fraction, Nanotrap_type, Enrichment, Kraken2_confidence, UniqueID) %>%
  mutate(RA = readcount / sum(readcount))

family_box_df_90 <- family_box_df %>%
  filter(Kraken2_confidence == "0.9")
family_box_df_90_w_rrna <- family_box_df_w_rrna %>%
  filter(Kraken2_confidence == "0.9")



targeted_families_90 <- family_box_df_90 %>%
  filter(`F` %in% vsp_targets$family_name | `F` %in% rpip_targets$family)


targeted_families_90_w_rrna <- family_box_df_90_w_rrna %>%
  filter(`F` %in% vsp_targets$family_name | `F` %in% rpip_targets$family)


ggplot(targeted_families_90, aes(x = RA, y = `F`, colour = Enrichment)) +
  geom_boxplot() +
  scale_x_log10() +
  theme_bw()

ggplot(targeted_families_90_w_rrna, aes(x = RA, y = `F`, colour = Enrichment)) +
  geom_boxplot() +
  scale_x_log10() +
  theme_bw()


mean_family_readcounts <- targeted_families_90 %>%
  group_by(`F`, Fraction, Nanotrap_type, Enrichment) %>%
  reframe(nreads = readcount) %>%
  group_by(`F`, Enrichment) %>%
  reframe(readmeans = mean(nreads), sd = sd(nreads))

families_w_over_25_mean_reads <- filter(mean_family_readcounts, readmeans > 25)

targeted_families_90_25 <- filter(targeted_families_90, `F` %in% families_w_over_25_mean_reads$`F`)

ggplot(targeted_families_90_25, aes(x = RA, y = `F`, colour = Enrichment)) +
  geom_boxplot() +
  scale_x_log10() +
  theme_bw()

vspvsnone_family_res_90conf$contrast <- "VSPvNONE"
rpipvsnone_family_res_90conf$contrast <- "RPIPvNONE"

vsp_contrast <- data.frame(vspvsnone_family_res_90conf)
vsp_contrast$Family = rownames(vsp_contrast)
rownames(vsp_contrast) <- NULL

rpip_contrast <- data.frame(rpipvsnone_family_res_90conf)
rpip_contrast$Family = rownames(rpip_contrast)
rownames(rpip_contrast) <- NULL

enrichment_contrasts <- bind_rows(vsp_contrast, rpip_contrast)

sig_enrichment_contrasts <- filter(enrichment_contrasts, padj <= 0.05)

targeted_families_90_sig_enrichment <- filter(targeted_families_90, `F` %in% sig_enrichment_contrasts$Family) %>%
  left_join(rpip_targets, by = join_by("F" == "family"))
targeted_families_90_sig_enrichment[is.na(targeted_families_90_sig_enrichment$Domain), ]$Domain <- "Virus"

ggplot(targeted_families_90_sig_enrichment, aes(x = RA, y = `F`, colour = Enrichment)) +
  geom_boxplot() +
  scale_x_log10() +
  theme_minimal() +
  facet_grid(rows = "Domain", scales = "free_y", space = "free_y") +
  xlab("Relative abundance") +
  ylab("") +  
  theme(strip.placement = "outside", strip.text.y = element_text(angle = 0, size = 10), 
                    strip.background = element_rect(fill = "gray90"),
                    strip.text = element_text(face = "bold"),
                    legend.position = "bottom") +
  geom_hline(yintercept=seq(1,nrow(boxplot_df_by_arg[,1]),1)+.5,color="black", linetype = "dotted", linewidth = .2) +
  geom_hline(yintercept=seq(1,nrow(boxplot_df_by_arg[,1]),1)-.5,color="black", linetype = "dotted", linewidth = .2)

rpip_targets$Domain


#Taxa only in enrichments
##########################

enriched_families_no_rrna <- bind_rows(rpip_families_no_rrna, vsp_families_no_rrna)
enriched_families_with_rrna <- bind_rows()

untargeted_families_no_rrna_no_zeros <- untargeted_families_no_rrna[untargeted_families_no_rrna$readcount != 0,]

enriched_only_no_rrna <- enriched_families_no_rrna %>%
  filter(!(`F` %in% untargeted_families_no_rrna_no_zeros$`F`))

table(enriched_only_90_no_rrna$`F`)

enriched_only_90_no_rrna_25 <- enriched_only_90_no_rrna %>%
  filter(`F` %in% families_w_over_25_mean_reads$`F`) %>%
  filter()

ggplot(enriched_only_90_no_rrna, aes(x = readcount, y = `F`, colour = Enrichment)) + 
  geom_boxplot() +
  scale_x_log10()

untargeted_RAs <- untargeted_families_no_rrna %>%
  select(Family = `F`, RA) %>%
  mutate(Enrichment = "None")
rpip_RAs <- rpip_families_no_rrna %>%
  select(Family = `F`, RA) %>%
  mutate(Enrichment = "RPIP")

RA_comparison <- bind_rows(untargeted_RAs, rpip_RAs)

RA_comparison_wide <- RA_comparison %>%
  pivot_wider(id_cols = c(SampleID, Family), names_from = Enrichment, values_from = RA) %>%
  mutate(RA_ratio = RPIP / None)

RA_increased_rpip <- filter(RA_comparison_wide, RA_ratio > 1)
unique(RA_increased_rpip$Family)