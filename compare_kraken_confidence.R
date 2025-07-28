conf_comparison_df <- all_families %>%
  mutate(Kraken2_confidence = as.character(Kraken2_confidence)) %>%
  group_by(Kraken2_confidence) %>% 
  select(!UniqueID) %>%
  pivot_wider(names_from = Kraken2_confidence, values_from = c(readcount, RA)) %>%
  mutate(readcount_0.9 = replace_na(readcount_0.9, 0)) %>%
  mutate(RA_0.9 = replace_na(RA_0.9, 0)) %>%
  mutate(readcount_0 = replace_na(readcount_0, 0)) %>%
  mutate(RA_0 = replace_na(RA_0, 0)) %>%
  mutate(Percent_reads_removed = 100 * (1 - readcount_0.9/readcount_0)) %>%
  mutate(Num_reads_removed = readcount_0 - readcount_0.9)




ggplot(filter(conf_comparison_df, `F` == "Coronaviridae"), aes(x = Num_reads_removed, fill = Enrichment)) + geom_density(alpha = 0.5) +
  scale_x_log10()
ggplot(filter(conf_comparison_df, `F` == "Coronaviridae"), aes(x = Percent_reads_removed, fill = Enrichment)) + geom_histogram(position = "dodge", binwidth = 10)


ggplot(filter(conf_comparison_df, `F` == "Caliciviridae"), aes(x = Num_reads_removed, fill = Enrichment)) + geom_density(alpha = 0.5) +
  scale_x_log10()
ggplot(filter(conf_comparison_df, `F` == "Caliciviridae"), aes(x = Percent_reads_removed, fill = Enrichment)) + geom_histogram(position = "dodge", binwidth = 10)

ggplot(filter(conf_comparison_df, `F` == "Caliciviridae"), aes(y = Percent_reads_removed, fill = Enrichment)) + geom_boxplot() + scale_y_log10()


View(filter(conf_comparison_df, `F` == "Coronaviridae"))



setwd("../../../../../../../vsp_panels/raw_fastqs/fastp_out_no_dedup/kraken2_out/k2_nt_20240530/")
nomin <- read_tsv("rrna_labeled/003-38156-FA-RP_S115_L007_k2out_nonrrna.out", col_names = c("is_classified", "read_id", "tax_id", "length", "matchloc"))
min90 <- read_tsv("90conf/rrna_labeled/003-38156-FA-RP_S115_L007_90confk2out_nonrrna.out", col_names = c("is_classified", "read_id", "tax_id", "length", "matchloc"))

coronaviridae_nomin <- nomin %>%
  filter(tax_id == 11118)
coronaviridae_min90 <- min90 %>%
  filter(tax_id == 11118)

reads_identified_only_at_90conf <- coronaviridae_min90 %>%
  filter(!(read_id %in% coronaviridae_nomin$read_id))

nomin_but_its_reads_ided_as_corona_in_90min <- nomin %>%
  filter(read_id %in% reads_identified_only_at_90conf$read_id)

betacoronavirus <- vsp_imported %>%
  filter(G == "Betacoronavirus")
