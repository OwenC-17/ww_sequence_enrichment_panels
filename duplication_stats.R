library(lme4)
library(tidyverse)

#source("read_fastp_reports.R")

all_dedup_reports$PreEnrichmentID <- paste(
  all_dedup_reports$LIMS_ID,
  all_dedup_reports$Fraction,
  all_dedup_reports$Nanotrap_type,
  sep = "_"
)

all_dedup_reports$Nanotrap_type <- relevel(factor(all_dedup_reports$Nanotrap_type, ordered = FALSE), ref = "none")
all_dedup_reports$Fraction <- relevel(factor(all_dedup_reports$Fraction, ordered = FALSE), ref = "unfiltered")
all_dedup_reports$Enrichment <- relevel(factor(all_dedup_reports$Enrichment, ordered = FALSE), ref = "None")



library(glmmTMB)
glmm_dupRate_Enrichment_model <- glmmTMB(rate ~ Enrichment + Fraction +
                                  (1|LIMS_ID),
                                data = all_dedup_reports, 
                                family = beta_family())

glmm_dupRate_Enrichment_model %>% summary()

glmmB_dupRate_Enrichment_model <- glmmTMB(cbind(num_reads_removed,
                                              summary.after_filtering.total_reads) ~ Enrichment + Fraction +
                                           (1|LIMS_ID),
                                         data = all_dedup_reports, 
                                         family = binomial)

glmmB_dupRate_Enrichment_model %>% summary()

#Diagnostic plots normality tests
qqnorm(residuals(lmer_dupRate_Enrichment_model))
qqline(residuals(lmer_dupRate_Enrichment_model))

shapiro.test(residuals(lmer_dupRate_Enrichment_model)) 
#Non-normal

nortest::ad.test(residuals(lmer_dupRate_Enrichment_model))
#Still non-normal

car::leveneTest(rate ~ Enrichment, data = all_fastp_summaries_unmerged)
#Still non-normal


##############Non-parametric tests
#Friedman (Rank-based) test:
fr_dupRate_Enrichment_model <- friedman.test(rate ~ Enrichment | PreEnrichmentID,
                                             data = all_dedup_reports)
#Show model output:
fr_dupRate_Enrichment_model
#Print full p-value:
fr_dupRate_Enrichment_model$p.value

#Post-hoc Wilcoxon Signed-Rank Tests
all_dedup_reports <- all_dedup_reports %>%
  group_by(Enrichment) %>%
  arrange(PreEnrichmentID)

wilcox_dupRate_Enrichment_post <- pairwise.wilcox.test(
  all_dedup_reports$rate, 
  all_dedup_reports$Enrichment,
  paired = TRUE,
  p.adjust.method = "none")

wilcox_dupRate_Enrichment_post$p.value

####Plot dup comparison:
library(ggpubr)
library(ggplot2)

ggplot(all_dedup_reports, aes(x = Enrichment, y = rate)) +
  geom_boxplot() + 
  theme_minimal() +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("RPIP", "None"),
                                        c("VSP", "RPIP"),
                                        c("VSP", "None")),
                     paired = TRUE,
                     label = "p.format",
                     digits = 10,
                     size = 5) +
  ylab("Read duplication rate") +
  xlab("Enrichment panel") +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 20))


ggplot(all_dedup_reports, aes(x = Enrichment, y = rate, fill = Fraction)) +
  geom_boxplot() + 
  theme_minimal() +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("RPIP", "None"),
                                        c("VSP", "RPIP"),
                                        c("VSP", "None")),
                     paired = TRUE,
                     label = "p.format",
                     digits = 10,
                     size = 5) +
  ylab("Read duplication rate") +
  xlab("Enrichment panel") +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 20))





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