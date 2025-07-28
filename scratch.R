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
