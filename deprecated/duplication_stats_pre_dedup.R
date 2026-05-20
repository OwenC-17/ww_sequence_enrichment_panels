all_fastp_summaries_unmerged$PreEnrichmentID <- paste(
  all_fastp_summaries_unmerged$LIMS_ID,
  all_fastp_summaries_unmerged$Fraction,
  all_fastp_summaries_unmerged$Nanotrap_type,
  sep = "_"
)

library(lme4)


####Based on pre-deduplication reports:
lmer_dupRate_Enrichment_model <- lmer(rate ~ Enrichment + 
                                  (1|PreEnrichmentID),
                                data = all_fastp_summaries_unmerged)

qqnorm(residuals(lmer_dupRate_Enrichment_model))
qqline(residuals(lmer_dupRate_Enrichment_model))

shapiro.test(residuals(lmer_dupRate_Enrichment_model)) 
#Sensitive to sample size? Probably overly sensitive.

nortest::ad.test(residuals(lmer_dupRate_Enrichment_model))

car::leveneTest(rate ~ Enrichment, data = all_fastp_summaries_unmerged)

plot(lmer_dupRate_Enrichment_model, which = 1)

######Arcsin-sqrt of duplication rate:
all_fastp_summaries_unmerged$arcsin_sqrt_dup <- asin(sqrt(all_fastp_summaries_unmerged$rate))

trans_lmer_dupRate_Enrichment_model <- lmer(arcsin_sqrt_dup ~ Enrichment +
                                              (1|PreEnrichmentID),
                                            data = all_fastp_summaries_unmerged)

qqnorm(residuals(trans_lmer_dupRate_Enrichment_model))
qqline(residuals(trans_lmer_dupRate_Enrichment_model))

shapiro.test(residuals(trans_lmer_dupRate_Enrichment_model))
nortest::ad.test(residuals(trans_lmer_dupRate_Enrichment_model))
car::leveneTest(arcsin_sqrt_dup ~ Enrichment, data = all_fastp_summaries_unmerged)

##############Non-parametric tests
#Friedman (Rank-based) test:
fr_dupRate_Enrichment_model <- friedman.test(rate ~ Enrichment | PreEnrichmentID,
                                             data = all_fastp_summaries_unmerged)
fr_dupRate_Enrichment_model

#Post-hoc Wilcoxon Signed-Rank Tests
all_fastp_summaries_unmerged_sorted <- all_fastp_summaries_unmerged %>%
  group_by(Enrichment) %>%
  arrange(PreEnrichmentID)

wilcox_dupRate_Enrichment <- pairwise.wilcox.test(
  all_fastp_summaries_unmerged_sorted$rate, 
  all_fastp_summaries_unmerged_sorted$Enrichment,
  paired = TRUE,
  p.adjust.method = "none")

wilcox_dupRate_Enrichment$p.value


##Based on post-duplication reports:
all_dedup_reports$PreEnrichmentID <- paste(
  all_dedup_reports$LIMS_ID,
  all_dedup_reports$Fraction,
  all_dedup_reports$Nanotrap_type,
  sep = "_"
)

fr_dupRate_Enrichment_model_post <- friedman.test(rate ~ Enrichment | PreEnrichmentID,
                                                  data = all_dedup_reports)
fr_dupRate_Enrichment_model_post$p.value

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
library(tidyverse)

custom_label <- function(p) {
  vapply(p, function(x) formatC(x, format = "f", digits = 4), character(1))
}

ggplot(all_fastp_summaries_unmerged, aes(x = Enrichment, y = rate)) +
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
