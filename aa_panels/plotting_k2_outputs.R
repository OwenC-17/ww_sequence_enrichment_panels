k36397_FA <- generate_tax_table("165-36397-FA_S1_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "36397_FA")
k36397_RA <- generate_tax_table("166-36397-RA_S2_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "36397_RA")

k36397_FB <- generate_tax_table("167-36397-FB_S3_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "36397_FB")

k36397_RB <- generate_tax_table("168-36397-RB_S4_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "36397_RB")

k36397_UA <- generate_tax_table("169-36397-UA_S5_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "36397_UA")

k36397_all <- bind_rows(k36397_FA, k36397_FB, k36397_RA, k36397_RB, k36397_UA)


k38156_UA <- generate_tax_table("170-38156-UA_S6_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "38156_UA")

k32976_UA <- generate_tax_table("171-32976-UA_S7_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "32976_UA")

k32989_UA <- generate_tax_table("172-32989-UA_S8_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "32989_UA")

k36397_UB <- generate_tax_table("173-36397-UB_S9_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "36397_UB")

k36397_UB <- generate_tax_table("174-38156-UB_S10_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "36397_UB")

k32976_UB <- generate_tax_table("175-32976-UB_S11_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "32976_UB")

k34970_FA <- generate_tax_table("177-34970-FA_S13_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "34970_FA")

k34970_RA <- generate_tax_table("178-34970-RA_S14_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "34970_RA")

k22015_UA <- generate_tax_table("179-22015-UA_S15_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "22015_UA")

k22015_UA <- generate_tax_table("179-22015-UA_S15_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "22015_UA")

k20040_UA <- generate_tax_table("180-20040-UA_S16_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "20040_UA")

k41596_UA <- generate_tax_table("181-41596-UA_S17_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "41596_UA")

k34970_UA <- generate_tax_table("182-34970-UA_S18_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "34970_UA")

k34970_UD <- generate_tax_table("183-34970-UD_S19_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "34970_UD")

k41596_UD <- generate_tax_table("184-41596-UD_S20_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "41596_UD")

k20040_UD <- generate_tax_table("185-20040-UD_S21_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "20040_UD")

k22015_UD <- generate_tax_table("186-22015-UD_S22_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "22015_UD")

k34970_RD <- generate_tax_table("187-34970-RD_S23_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "34970_RD")


k34970_FD <- generate_tax_table("188-34970-FD_S24_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "34970_FD")

all_vs_aa_db <- bind_rows(k20040_UA,k20040_UD,k22015_UA,k22015_UD,
                          k32976_UA,k32976_UB,k32989_UA,k34970_FA,
                          k34970_FD,k34970_RA,k34970_RD,k34970_UA,
                          k34970_UD,k36397_FA,k36397_FB,k36397_RA,
                          k36397_RB,k36397_UA,k36397_UB,k38156_UA,
                          k41596_UA,k41596_UD)


plot_subset(all_vs_aa_db, "nodeOnly", "R", "root", "S") +
  labs(fill = "Species") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ggplot(filter(k36397_all, R == "root"), aes(x = SampleID, y = RA, fill = S)) + geom_bar(position = 'stack', stat = 'identity') +
  ylab("Relative abundance of reads") +
  xlab("sampleID") + 
  theme_bw() +
  labs(fill = "Species") + 
  theme(text = element_text(size = 15)) +
  guides(fill = guide_legend(ncol = 1))


ggplot(filter(k36397_all, R == "root"), aes(x = "", y = RA, fill = S)) + geom_bar(position = 'stack', stat = 'identity') +
  coord_polar("y", start = 0)
  ylab("Relative abundance of reads") +
  xlab("sampleID") + 
  theme_bw() +
  labs(fill = "Species") + 
  theme(text = element_text(size = 15)) +
  guides(fill = guide_legend(ncol = 1))
