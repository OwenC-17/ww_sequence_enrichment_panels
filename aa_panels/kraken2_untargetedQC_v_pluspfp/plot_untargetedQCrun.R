#Run script from directory where outputs are located

k36397_FA_QC <- generate_tax_table("165-36397-FA_S165_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "36397_FA")

k36397_RA_QC <- generate_tax_table("166-36397-RA_S166_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "36397_RA")

k36397_FB_QC <- generate_tax_table("167-36397-FB_S167_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "36397_FB")

k36397_RB_QC <- generate_tax_table("168-36397-RB_S168_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "36397_RB")

k36397_UA_QC <- generate_tax_table("169-36397-UA_S169_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "36397_UA")

k36397_all_QC <- bind_rows(k36397_FA_QC, k36397_FB_QC, k36397_RA_QC, k36397_RB_QC, k36397_UA_QC)


k38156_UA_QC <- generate_tax_table("170-38156-UA_S170_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "38156_UA")

k32976_UA_QC <- generate_tax_table("171-32976-UA_S171_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "32976_UA")

k32989_UA_QC <- generate_tax_table("172-32989-UA_S172_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "32989_UA")

k36397_UB_QC <- generate_tax_table("173-36397-UB_S173_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "36397_UB")

k36397_UB_QC <- generate_tax_table("174-38156-UB_S174_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "36397_UB")

k32976_UB_QC <- generate_tax_table("175-32976-UB_S175_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "32976_UB")

k34970_FA_QC <- generate_tax_table("177-34970-FA_S177_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "34970_FA")

k34970_RA_QC <- generate_tax_table("178-34970-RA_S178_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "34970_RA")

k22015_UA_QC <- generate_tax_table("179-22015-UA_S179_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "22015_UA")

k20040_UA_QC <- generate_tax_table("180-20040-UA_S180_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "20040_UA")

k41596_UA_QC <- generate_tax_table("181-41596-UA_S181_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "41596_UA")

k34970_UA_QC <- generate_tax_table("182-34970-UA_S182_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "34970_UA")

k34970_UD_QC <- generate_tax_table("183-34970-UD_S183_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "34970_UD")

k41596_UD_QC <- generate_tax_table("184-41596-UD_S184_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "41596_UD")

k20040_UD_QC <- generate_tax_table("185-20040-UD_S185_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "20040_UD")

k22015_UD_QC <- generate_tax_table("186-22015-UD_S186_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "22015_UD")

k34970_RD_QC <- generate_tax_table("187-34970-RD_S187_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "34970_RD")


k34970_FD_QC <- generate_tax_table("188-34970-FD_S188_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "34970_FD")

all_vs_pluspfp_QC <- bind_rows(k20040_UA_QC,k20040_UD_QC,k22015_UA_QC,k22015_UD_QC,
                               k32976_UA_QC,k32976_UB_QC,k32989_UA_QC,k34970_FA_QC,
                               k34970_FD_QC,k34970_RA_QC,k34970_RD_QC,k34970_UA_QC,
                               k34970_UD_QC,k36397_FA_QC,k36397_FB_QC,k36397_RA_QC,
                               k36397_RB_QC,k36397_UA_QC,k36397_UB_QC,k38156_UA_QC,
                               k41596_UA_QC,k41596_UD_QC)

all_vs_pluspfp_QC$rrna_removed <- "rRNA not removed"

plot_subset(all_vs_pluspfp_QC, "RA", "U", "unclassified", "D") +
  labs(fill = "Domain") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

plot_subset(all_vs_pluspfp_QC, "RA", "D", "Viruses", "P") +
  labs(fill = "Phylum") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Plot all viruses that are assigned at the family level or lower (stacked bars with read counts), 
#faceted by treatment and grouped by sample
ggplot(filter(all_vs_pluspfp_QC, D == "Viruses" & !startsWith(F, "unidentified")), 
       aes(x = SampleID, y = nodeOnly, fill = F)) + 
  geom_bar(position = 'stack', stat = 'identity') +
  ylab("Number of reads") +
  theme_bw() +
  labs(fill = "Viral Families") + #legend title
  theme(text = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))
