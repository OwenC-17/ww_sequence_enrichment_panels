#Compare presence-absence by K2 Confidence

all_families <- all_families %>%
  mutate(FamilyName = `F`)

no_fam_indices <- str_which(all_families$FamilyName, "^unidentified")

all_families[no_fam_indices, ]$FamilyName <- "Unidentified"

all_families_0 <- all_families %>%
  filter(Kraken2_confidence == 0)

all_families_90 <- all_families %>%
  filter(Kraken2_confidence == 0.9)

removed_by_conf_filtering <- all_families_0[!(all_families_0$F %in% all_families_90$F),]

removed_abund <- removed_by_conf_filtering %>%
  group_by(UniqueID, LIMS_ID, Treatment, site, Fraction, Nanotrap_type, Enrichment) %>%
  summarize(NumReadsRemoved = sum(readcount), AbundRemoved = sum(RA))

ggplot(removed_abund, aes(y = AbundRemoved, x = Nanotrap_type)) + geom_boxplot()

removed_taxa_count <- removed_by_conf_filtering %>%
  group_by(UniqueID, LIMS_ID, Treatment, site, Fraction, Nanotrap_type, Enrichment) %>%
  summarize(NumTaxaRemoved = n())

total_taxa_count <- all_families_0 %>%
  group_by(UniqueID, LIMS_ID, Treatment, site, Fraction, Nanotrap_type, Enrichment) %>%
  summarize(TotalNumTaxa = n()) %>%
  full_join(removed_taxa_count) %>%
  mutate(PropTaxaRemoved = NumTaxaRemoved / TotalNumTaxa) %>%
  full_join(removed_abund)

pivot_longer(total_taxa_count, cols = c(AbundRemoved, PropTaxaRemoved)) %>%
  ggplot(aes(x = LIMS_ID, y = value, fill = name)) + geom_boxplot() + scale_y_log10() + facet_wrap(~ Enrichment)

ggplot(total_taxa_count, aes(y = NumTaxaRemoved, fill = Enrichment)) + geom_boxplot()

range(total_taxa_count$AbundRemoved)

families_removed <- unique(removed_by_conf_filtering$F)



###reticulate
library(reticulate)
use_condaenv("cowenv", require = TRUE)


#

ggplot(pileup_example, aes(x = read_count)) + geom_histogram()

ggplot(test_, aes(y = D, x = RA)) + geom_boxplot() #+ scale_x_log10()

Viruses <- filter(test_, D == "Viruses")




##############
simulate_drift <- function(p0, generations, N) {
  freqs <- numeric(generations)
  freqs[1] <- p0
  for (i in 2:generations) {
    freqs[i] <- rbinom(1, N, freqs[i-1]) / N
  }
  return(freqs)
}


simulate_selection <- function(p0, generations, N, s) {
  freqs <- numeric(generations)
  freqs[1] <- p0
  for (i in 2:generations) {
    w_A <- 1 + s #fitness of allele 1
    w_a <- 1     #fitness of allele 2
    
  p_prev <- freqs[i - 1]
  p_sel <- (p_prev * w_A) / (p_prev * w_A + (1 - p_prev) * w_a)
  
  freqs[i] <- rbinom(1, N, p_sel) / N
  }
  return(freqs)
}

simulate_selection(0.9, 10000, 100000, -0.001) %>% plot()


simulate_selection_multiple <- function(times, p0, generations, N, s) {
  x <- logical(length = times)
  for (time in 1:times) {
    y <- simulate_selection(p0, generations, N, s)
    x[time] <- y[generations]
  }
  return(x)
}

bleble <- simulate_selection_multiple(100, 0.001, 10000, 100000, 0.001)
sum(bleble)

Calici <- allFamilies_anyConf_noRrna_dedup %>%
  filter(F == "Caliciviridae" & Kraken2_confidence == "0.9")

ggplot(Calici, aes(x = Enrichment, y = readcount, color = Nanotrap_type)) +
  geom_point() +
  geom_line(aes(group = interaction(Nanotrap_type, Fraction)))


rpip_allele_graph_df <- barplot_rpip_allele_df %>%
  mutate(full_treatment = paste(site, Fraction, Nanotrap_type, Enrichment, sep = "_")) %>%
  filter(`Percent Coverage` > 90) %>%
  select(`ARO Term`, site, Fraction, Nanotrap_type, Enrichment, full_treatment) %>%
  distinct()

rpip_allele_count_df <- rpip_allele_graph_df %>%
  group_by(site, Fraction, Nanotrap_type, Enrichment, full_treatment) %>%
  summarize(n_genes = n())

rpip_gene_lists <- rpip_allele_graph_df %>%
  group_by(site, Fraction, Nanotrap_type, Enrichment, full_treatment) %>%
  summarize(genes = list(`ARO Term`), .groups = "drop")

rpip_pairs <- t(combn(rpip_gene_lists$full_treatment, 2))

rpip_edge_df <- map_dfr(seq_len(nrow(rpip_pairs)), function(i) {
  s1 <- rpip_pairs[i, 1]
  s2 <- rpip_pairs[i, 2]
  
  g1 <- rpip_gene_lists$genes[[which(rpip_gene_lists$full_treatment == s1)]]
  g2 <- rpip_gene_lists$genes[[which(rpip_gene_lists$full_treatment == s2)]]
  
  tibble(
    from = s1,
    to = s2,
    shared = length(intersect(g1, g2))
  )
}) %>%
  filter(shared > 0)

vertex_df <- rpip_allele_count_df %>%
  select(full_treatment, site, Fraction, Nanotrap_type, Enrichment, full_treatment, n_genes)

library(igraph)
g <- graph_from_data_frame(
  d = rpip_edge_df,
  vertices = vertex_df,
  directed = FALSE
)

library(ggraph)

ggraph(g, layout = "hive", axis = Nanotrap_type) +
  geom_edge_link(aes(width = shared), alpha = 0.6, colour = "grey40") +
  geom_node_point(aes(size = n_genes, colour = site)) +
  scale_edge_width(range = c(0.2, 3)) +
  scale_size(range = c(0.2, 12)) +
  theme_bw()

rpip_allele_ohare_filt_venn_df <- barplot_rpip_allele_df %>%
  filter(site == "OHare") %>%
  mutate(full_treatment = paste(site, Fraction, sep = "_")) %>%
  filter(`Percent Coverage` > 90)

treatment_list <- as.list(unique(rpip_allele_ohare_filt_venn_df$full_treatment))
names(treatment_list) <- unlist(treatment_list)
rpip_allele_venn_list <- lapply(treatment_list, function(x) {
  filter(rpip_allele_venn_df, full_treatment == x)$`Reference Sequence` %>%
    unique()
}
)

ggVennDiagram(rpip_allele_venn_list)

rpip_ARG_readcounts <- rpip_allele_df %>%
  group_by(UniqueID, `ARO Term`) %>%
  summarize(nmapped = sum(`All Mapped Reads`)) %>%
  pivot_wider(names_from = UniqueID, values_from = nmapped, id_cols = `ARO Term`)

rpip_lib.size <- rpip_allele_group_data$summary.after_filtering.total_reads








####Some exploratory plots (RA)####
rpip_and_unt_allele_and_info %>% 
  filter(`All Mapped Reads` > 4 & `Reference Model Type` == "protein homolog model") %>%
  ggplot(aes(y = `AMR Gene Family`, x = sample_id, fill = log10(RA_nmapped))) +
  geom_tile() +
  theme(axis.text.y = element_blank()) +
  facet_grid(~ Enrichment)

allele_RA_matrix <- rpip_and_unt_allele_and_info %>%
  pivot_wider(names_from = `AMR Gene Family`, values_from = RA_nmapped,
              id_cols = UniqueID, values_fn = sum) %>%
  column_to_rownames('UniqueID') %>%
  as.matrix()

allele_RA_matrix[is.na(allele_RA_matrix)] <- 0
log10_allele_RA_matrix <- log10(allele_RA_matrix)
log10_allele_RA_matrix[is.na(log10_allele_RA_matrix) | is.infinite(log10_allele_RA_matrix)] <- 0


heatmap(log10_allele_RA_matrix)




#######Random plots for RGI gene-level results: 
ggplot(RPIP_genelevel_results, aes(x = logFC, y = -log10(FDR), colour = Fraction, shape = Nanotrap_type)) + 
  geom_point(alpha = 0.67, size = 2) +
  geom_hline(yintercept = -log10(0.05)) +
  theme_bw() +
  scale_color_manual(values = c("#A50021FF", "#3FA0FFFF", "#FFAD72FF")) +
  scale_shape_manual(values = c(1, 6, 4))

ggplot(RPIP_results, aes(x = logFC, y = Treatment, fill = Fraction)) + geom_boxplot()
ggplot(RPIP_results, aes(x = Treatment, y = -log10(FDR), fill = Fraction)) + geom_boxplot()


summarized_rpip_genes_by_family <- RPIP_results %>%
  filter(FDR <= 0.05) %>%
  select(logFC, reference, Nanotrap_type, Fraction, `AMR Gene Family`) %>%
  group_by(Nanotrap_type, Fraction, `AMR Gene Family`) %>%
  summarize(mean_logFC = mean(logFC),
            #            unique_ref = list(unique(reference)),
            numGenes = n_distinct(reference))

summarized_arg_annotations <- RPIP_results %>%
  select(reference, `AMR Gene Family`, `Drug Class`, `Resistance Mechanism`)

ggplot(summarized_rpip_genes_by_family, aes(x = `AMR Gene Family`, y = Fraction, fill = mean_logFC, size = numGenes)) +
  geom_point(shape = 21, stroke = 1) +
  #  scale_radius(transform = "log10") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(rows = "Nanotrap_type") +
  scale_fill_viridis_c(option = "H")


summarized_rpip_genes_by_target <- RPIP_results %>%
  filter(FDR <= 0.05) %>%
  select(logFC, reference, Nanotrap_type, Fraction, `Drug Class`, `Resistance Mechanism`) %>%
  group_by(Nanotrap_type, Fraction, `Drug Class`, `Resistance Mechanism`) %>%
  summarize(mean_logFC = mean(logFC),
            #            unique_ref = list(unique(reference)),
            numGenes = n_distinct(reference)) 

summarized_rpip_genes_by_target <- summarized_rpip_genes_by_target %>%
  mutate(drugClassLabel = if_else(str_count(`Drug Class`, ";") > 1,
                                  paste0("multi (", `Resistance Mechanism`, ")"),
                                  `Drug Class`) %>% str_remove_all("antibiotic"))

ggplot(summarized_rpip_genes_by_target, aes(x = `drugClassLabel`, y = Fraction, fill = mean_logFC, size = numGenes)) +
  geom_point(shape = 21, stroke = 1) +
  #  scale_radius(transform = "log10") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(rows = "Nanotrap_type") +
  scale_fill_viridis_c(option = "H")# +
theme(axis.text.x = element_blank())

summarized_rpip_genes_by_mechanism <- RPIP_results %>%
  filter(FDR <= 0.05) %>%
  select(logFC, reference, Nanotrap_type, Fraction, `Resistance Mechanism`) %>%
  group_by(Nanotrap_type, Fraction, `Resistance Mechanism`) %>%
  summarize(mean_logFC = mean(logFC),
            #            unique_ref = list(unique(reference)),
            numGenes = n_distinct(reference))

ggplot(summarized_rpip_genes_by_mechanism, aes(x = `Resistance Mechanism`, y = Fraction, fill = mean_logFC, size = numGenes)) +
  geom_point(shape = 21, stroke = 1) +
  #  scale_radius(transform = "log10") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(rows = "Nanotrap_type") +
  scale_fill_viridis_c(option = "H")


rarefaction_table <- rpip_and_unt_gene_and_info$`Resistomes & Variants: Observed Pathogen(s)` %>%
  table() %>%
  as.data.frame() %>%
  arrange(desc(Freq))
rarefaction_table <- rarefaction_table %>%
  mutate(X = seq(1:nrow(rarefaction_table))) %>%
  mutate(cums = cumsum(Freq)) %>%
  mutate(cump = 100 * (cums / sum(Freq)))

top_90perc <- rarefaction_table %>%
  filter(cump <=  80)

RPIP_results <- RPIP_results %>%
  mutate(tax_viz_var = ifelse(
    `Resistomes & Variants: Observed Pathogen(s)` %in% top_90perc$`.`, `Resistomes & Variants: Observed Pathogen(s)`, "Other"
  )) 

RPIP_results %>%
  ggplot(aes(x = logFC, y = -log10(FDR), colour = tax_viz_var, shape = Nanotrap_type)) + 
  geom_point(alpha = 0.67, size = 2) +
  geom_hline(yintercept = -log10(0.05)) +
  theme_bw() +
  #    scale_color_manual(values = c("#A50021FF", "#3FA0FFFF", "#FFAD72FF")) +
  scale_shape_manual(values = c(1, 6, 4))

signif <- RPIP_results %>%
  filter(FDR <= 0.05)

most_important <- signif %>%
  group_by(reference) %>%
  summarize(mean_logfc = mean(logFC)) %>%
  arrange(desc(mean_logfc))

RPIP_results_important <- RPIP_results %>%
  filter(reference %in% most_important$reference)

ggplot(RPIP_results_important, aes(x = Nanotrap_type, 
                                   y = Fraction, 
                                   fill = logFC)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkred") +
  facet_wrap(~ reference, ncol = 10) +
  theme_minimal()

RPIP_results_important$stars <- cut(
  RPIP_results_important$FDR,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", "-")
)



ggplot(RPIP_results_important, aes(x = Nanotrap_type, 
                                   y = Fraction, 
                                   fill = logFC)) +
  geom_tile(color = "black", linewidth = .5) +
  scale_fill_gradient(low = "white", high = "darkred") + 
  facet_wrap(~ Species, ncol = 10) +
  theme_minimal() +
  theme(strip.text = element_text(size = 12), 
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 10))+
  geom_text(aes(label = stars), size = 5, vjust = 0.5, fontface = "bold") +
  xlab("Concentration method")

rpipunt_fungi <- rpipunt_taxtable %>%
  filter(kingdom == "Fungi")

rpipunt_virus <- rpipunt_taxtable %>%
  filter(domain == "Viruses")

rpipunt_bacte <- rpipunt_taxtable %>%
  filter(domain == "Bacteria")

ggplot(filter(RPIP_results_important, Species %in% rpipunt_virus$species), aes(x = Nanotrap_type, 
                                                                               y = Fraction, 
                                                                               fill = logFC)) +
  geom_tile() +
  scale_fill_gradient2(low = "navyblue", mid = "white", high = "darkred", midpoint = 0) + 
  facet_wrap(~ Species, ncol = 3) +
  theme_minimal()

RPIP_results_important <- RPIP_results_important %>% 
  mutate(FDR = p.adjust(PValue, method = "BH"))

RPIP_results_important$stars <- cut(
  RPIP_results_important$FDR,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", "-")
)



ggplot(RPIP_results_important, aes(x = Nanotrap_type, 
                                   y = Fraction, 
                                   fill = logFC)) +
  geom_tile(color = "black", linewidth = .5) +
  scale_fill_gradient2(low = "deepskyblue4", mid = "white", high = "darkred", midpoint = 0) + 
  facet_wrap(~ Family, ncol = 3) +
  theme_minimal() +
  theme(strip.text = element_text(size = 12), 
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 10))+
  geom_text(aes(label = stars), size = 5, vjust = 0.5, fontface = "bold") +
  xlab("Concentration method")





rpip_and_unt_counts_covstats_and_info <- read_rds("input/modified/rpip_and_unt_counts_covstats_and_info.rds")

rpip_and_unt_counts_covstats_and_info %>%
  select(Reference, rlength) %>%
  distinct() %>%
  select(Reference) %>%
  table() %>%
  View()




rpip_and_unt_counts_covstats_and_info %>%
  select(taxid, species:name) %>%
  distinct() %>%
  select(taxid) %>%
  table() %>%
  View()









vsp_group_data <- vsp_group_data %>%
  mutate(scaled_shannon = shannon / log(libsize))

vsp_group_data <- vsp_group_data %>%
  mutate(scaled_richness = richness / libsize)

descdist(vsp_group_data$scaled_richness, boot = 1000)

shannon_beta_model <- glmmTMB(scaled_shannon ~ Nanotrap_type*Enrichment + Fraction*Enrichment +
                                Nanotrap_type * Fraction + (1 | LIMS_ID) + (1 | LIMS_ID:Nanotrap_type:Fraction),
                              ziformula = ~Enrichment,
                              family = beta_family(link = "logit"),
                              data = filter(vsp_group_data, scaled_shannon > 0))
AIC(shannon_beta_model)

sim_resid <- simulateResiduals(shannon_beta_model, n = 1000)
plot(sim_resid)
testDispersion(sim_resid)

emmeans(shannon_beta_model, pairwise ~ Nanotrap_type, adjust = "fdr")
emmeans(shannon_beta_model, pairwise ~ Fraction, adjust = "fdr")
emmeans(shannon_beta_model, pairwise ~ Nanotrap_type | Enrichment, adjust = "fdr")
emmeans(shannon_beta_model, pairwise ~ Fraction | Enrichment, adjust = "fdr")
emmeans(shannon_beta_model, pairwise ~ Fraction | Nanotrap_type | Enrichment, adjust = "fdr")


shannon_beta_model_2 <- glmmTMB(scaled_shannon ~ Nanotrap_type*Enrichment + Fraction*Enrichment +
                                  Nanotrap_type * Fraction + (1 | LIMS_ID),
                                ziformula = ~Enrichment,
                                family = beta_family(link = "logit"),
                                data = vsp_group_data)
AIC(shannon_beta_model_2)

sim_resid <- simulateResiduals(shannon_beta_model_2, n = 1000)
plot(sim_resid)
testDispersion(sim_resid)

emmeans(shannon_beta_model_2, pairwise ~ Nanotrap_type, adjust = "fdr")
emmeans(shannon_beta_model_2, pairwise ~ Fraction, adjust = "fdr")
emmeans(shannon_beta_model_2, pairwise ~ Nanotrap_type | Enrichment, adjust = "fdr")
emmeans(shannon_beta_model_2, pairwise ~ Fraction | Enrichment, adjust = "fdr")
emmeans(shannon_beta_model_2, pairwise ~ Fraction | Nanotrap_type | Enrichment, adjust = "fdr")



richness_nbinom_model <- glmmTMB(richness ~ Nanotrap_type*Enrichment + Fraction*Enrichment +
                                   Nanotrap_type * Fraction + (1 | LIMS_ID) + (1 | LIMS_ID:Nanotrap_type:Fraction),
                                 family = nbinom2,
                                 data = vsp_group_data)
sim_resid <- simulateResiduals(richness_nbinom_model, n = 1000)
plot(sim_resid)
testDispersion(sim_resid)
testZeroInflation(sim_resid)

richness_nbinom_model_zi <- glmmTMB(richness ~ Nanotrap_type*Enrichment + Fraction*Enrichment +
                                      Nanotrap_type * Fraction + (1 | LIMS_ID) + (1 | LIMS_ID:Nanotrap_type:Fraction),
                                    ziformula = ~Enrichment,
                                    family = nbinom2,
                                    data = vsp_group_data)
sim_resid <- simulateResiduals(richness_nbinom_model_zi, n = 1000)
plot(sim_resid)
testDispersion(sim_resid)
testZeroInflation(sim_resid)

emmeans(richness_nbinom_model_zi, pairwise ~ Nanotrap_type, adjust = "fdr")
emmeans(richness_nbinom_model_zi, pairwise ~ Fraction, adjust = "fdr")
emmeans(richness_nbinom_model_zi, pairwise ~ Nanotrap_type | Enrichment, adjust = "fdr")
emmeans(richness_nbinom_model_zi, pairwise ~ Fraction | Enrichment, adjust = "fdr")

descdist(vsp_group_data$simpson, boot = 1000)



richness_gamma_model_zi <- glmmTMB(richness ~ Nanotrap_type*Enrichment + Fraction*Enrichment +
                                     Nanotrap_type * Fraction + (1 | LIMS_ID),
                                   ziformula = ~Enrichment,
                                   family = ziGamma(link = "inverse"),
                                   data = vsp_group_data)
sim_resid <- simulateResiduals(richness_gamma_model_zi, n = 1000)
plot(sim_resid)
testDispersion(sim_resid)
testZeroInflation(sim_resid)

emmeans(richness_nbinom_model_zi, pairwise ~ Nanotrap_type, adjust = "fdr")
emmeans(richness_nbinom_model_zi, pairwise ~ Fraction, adjust = "fdr")
emmeans(richness_nbinom_model_zi, pairwise ~ Nanotrap_type | Enrichment, adjust = "fdr")
emmeans(richness_nbinom_model_zi, pairwise ~ Fraction | Enrichment, adjust = "fdr")

descdist(vsp_group_data$scaled_richness, boot = 1000)
diagnose(richness_beta_model_zi)


#Duplication rate GLMM:
library(glmmTMB)
glmm_dupRate_Enrichment_model <- glmmTMB(rate ~ Enrichment + Fraction +
                                           (1|LIMS_ID),
                                         data = all_dedup_reports, 
                                         family = beta_family())

glmm_dupRate_Enrichment_model %>% summary()

glmmB_dupRate_Enrichment_model <- glmmTMB(cbind(num_reads_removed,
                                                summary.after_filtering.total_reads) ~ Enrichment +
                                            (1|LIMS_ID),
                                          data = all_dedup_reports, 
                                          family = binomial)

glmmB_dupRate_Enrichment_model %>% summary()

#Diagnostic plots normality tests
qqnorm(residuals(glmmB_dupRate_Enrichment_model))
qqline(residuals(glmmB_dupRate_Enrichment_model))

shapiro.test(residuals(glmmB_dupRate_Enrichment_model)) 
#Normal-ish

nortest::ad.test(residuals(glmmB_dupRate_Enrichment_model))
#Non-normal

car::leveneTest(rate ~ Enrichment, data = all_fastp_summaries_unmerged)
#Still non-normal


##############Non-parametric tests
#Friedman (Rank-based) test; sample_id controls for original sample ID and
#filtration/nanotrap pretreatments:
fr_dupRate_Enrichment_model <- friedman.test(rate ~ Enrichment | sample_id,
                                             data = all_dedup_reports)
#Show model output:
fr_dupRate_Enrichment_model
#Print full p-value:
fr_dupRate_Enrichment_model$p.value

#Post-hoc Wilcoxon Signed-Rank Tests
all_dedup_reports <- all_dedup_reports %>%
  group_by(Enrichment) %>%
  arrange(sample_id)

wilcox_dupRate_Enrichment_post <- pairwise.wilcox.test(
  all_dedup_reports$rate, 
  all_dedup_reports$Enrichment,
  paired = TRUE,
  p.adjust.method = "none")

wilcox_dupRate_Enrichment_post$p.value

################################################################################
###Duplication in rRNA vs. non-rRNA
dupdedup_withRrna <- bind_rows(allFamilies_anyConf_withRrna, allFamilies_anyConf_withRrna_dedup)

dupdedup_noRrna <- bind_rows(allFamilies_anyConf_noRrna, allFamilies_anyConf_noRrna_dedup)

dupdedup_corona <- filter(dupdedup_noRrna, `F` == "Coronaviridae")



ggplot(dupdedup_corona, aes(x = dedup, y = readcount, group = UniqueID)) + geom_point() + geom_line() + scale_y_log10()


dupdedup_onlyRrna <- dupdedup_withRrna %>%
  filter()

collapse_to_family_keep_rrna_separated <- function(imported_k2_report) {
  families <- collapse_to_tax_level(imported_k2_report, "F")
  families <- families %>%
    group_by(across(-c(readcount, RA))) %>%
    summarise(readcount = sum(readcount), .groups = "drop") %>%
    group_by(across(-c(`F`, readcount))) %>%
    mutate(RA = readcount / sum(readcount))
}

allFamilies_anyConf_rrnaSep <- collapse_to_family_keep_rrna_separated(
  all_k2_reports_anyConf
)

allFamilies_anyConf_rrnaSep_dedup <- collapse_to_family_keep_rrna_separated(
  all_k2_reports_deduped
)

dupdedup_rrna_separated <- bind_rows(allFamilies_anyConf_rrnaSep,
                                     allFamilies_anyConf_rrnaSep_dedup)

dupdedup_rrnaOnly <- dupdedup_rrna_separated %>%
  filter(ribosomal == "rrnaOnly") %>%
  filter(Kraken2_confidence == "0.9")

rrna_counts_dupdedup <- dupdedup_rrnaOnly %>%
  mutate(UniqueID = paste(LIMS_ID, Treatment, Enrichment, dedup, sep = "_")) %>%
  group_by(across(-c(`F`, readcount, RA, Kraken2_confidence))) %>%
  summarize(total_rrna_reads = sum(readcount))

nonRrna_counts_dupdedup <- dupdedup_noRrna %>%
  filter(Kraken2_confidence == "0.9") %>%
  mutate(UniqueID = paste(LIMS_ID, Treatment, Enrichment, dedup, sep = "_")) %>%
  group_by(across(-c(`F`, readcount, RA, Kraken2_confidence))) %>%
  summarize(total_nonRrna_reads = sum(readcount))


rrna_counts_dupdedup_wide <- rrna_counts_dupdedup %>%
  ungroup() %>%
  select(-UniqueID) %>%
  pivot_wider(names_from = dedup, values_from = total_rrna_reads) %>%
  mutate(duprate = (`FALSE` - `TRUE`) / `FALSE`)  %>%
  arrange(duprate) %>%
  mutate(x = 1:144)


nonRrna_counts_dupdedup_wide <- nonRrna_counts_dupdedup %>%
  ungroup() %>%
  select(-UniqueID) %>%
  pivot_wider(names_from = dedup, values_from = total_nonRrna_reads) %>%
  mutate(duprate = (`FALSE` - `TRUE`) / `FALSE`) %>%
  arrange(duprate) %>%
  mutate(x = 1:144)

rrna_vs_nonRrna_duprates <- bind_rows(rrna_counts_dupdedup_wide, 
                                      nonRrna_counts_dupdedup_wide)

ggplot(rrna_vs_nonRrna_duprates, aes(x = x, y = duprate, 
                                     color = ribosomal)) + 
  geom_point()

