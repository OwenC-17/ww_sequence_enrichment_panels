#Run script from directory where outputs are located

k36397_FA <- generate_tax_table("165-36397-FA_S224_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "36397_FA")
k36397_RA <- generate_tax_table("166-36397-RA_S225_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "36397_RA")

k36397_FB <- generate_tax_table("167-36397-FB_S226_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "36397_FB")

k36397_RB <- generate_tax_table("168-36397-RB_S227_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "36397_RB")

k36397_UA <- generate_tax_table("169-36397-UA_S228_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "36397_UA")

k36397_all <- bind_rows(k36397_FA, k36397_FB, k36397_RA, k36397_RB, k36397_UA)


k38156_UA <- generate_tax_table("170-38156-UA_S229_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "38156_UA")

k32976_UA <- generate_tax_table("171-32976-UA_S230_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "32976_UA")

k32989_UA <- generate_tax_table("172-32989-UA_S231_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "32989_UA")

k36397_UB <- generate_tax_table("173-36397-UB_S232_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "36397_UB")

k36397_UB <- generate_tax_table("174-38156-UB_S233_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "36397_UB")

k32976_UB <- generate_tax_table("175-32976-UB_S234_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "32976_UB")

k34970_FA <- generate_tax_table("177-34970-FA_S236_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "34970_FA")

k34970_RA <- generate_tax_table("178-34970-RA_S237_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "34970_RA")

k22015_UA <- generate_tax_table("179-22015-UA_S238_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "22015_UA")

k20040_UA <- generate_tax_table("180-20040-UA_S239_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "20040_UA")

k41596_UA <- generate_tax_table("181-41596-UA_S240_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "41596_UA")

k34970_UA <- generate_tax_table("182-34970-UA_S241_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "34970_UA")

k34970_UD <- generate_tax_table("183-34970-UD_S242_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "34970_UD")

k41596_UD <- generate_tax_table("184-41596-UD_S243_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "41596_UD")

k20040_UD <- generate_tax_table("185-20040-UD_S244_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "20040_UD")

k22015_UD <- generate_tax_table("186-22015-UD_S245_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "22015_UD")

k34970_RD <- generate_tax_table("187-34970-RD_S246_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "34970_RD")


k34970_FD <- generate_tax_table("188-34970-FD_S247_L001_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "34970_FD")

all_vs_pluspfp <- bind_rows(k20040_UA,k20040_UD,k22015_UA,k22015_UD,
                          k32976_UA,k32976_UB,k32989_UA,k34970_FA,
                          k34970_FD,k34970_RA,k34970_RD,k34970_UA,
                          k34970_UD,k36397_FA,k36397_FB,k36397_RA,
                          k36397_RB,k36397_UA,k36397_UB,k38156_UA,
                          k41596_UA,k41596_UD)
all_vs_pluspfp$rrna_removed <- "rRNA not removed"


plot_subset(all_vs_aa_db, "RA", "D", "Viruses", "P") +
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
  
  
setwd("C:/Users/fanta/Dropbox/PC/Documents/OBrienMetagenome/aa_panels/kraken2_jun6seqs_v_pluspfp/rRNA_removed/")

k36397_FA_no_rRNA <- generate_tax_table("36397-FA_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "36397_FA_no_rRNA")

k36397_RA_no_rRNA <- generate_tax_table("36397-RA_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "36397_RA_no_rRNA")

k36397_FB_no_rRNA <- generate_tax_table("36397-FB_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "36397_FB_no_rRNA")

k36397_RB_no_rRNA <- generate_tax_table("36397-RB_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "36397_RB_no_rRNA")

k36397_UA_no_rRNA <- generate_tax_table("36397-UA_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "36397_UA_no_rRNA")

k36397_all_no_rRNA <- bind_rows(k36397_FA_no_rRNA, k36397_FB_no_rRNA, k36397_RA_no_rRNA, k36397_RB_no_rRNA, k36397_UA_no_rRNA)


k38156_UA_no_rRNA <- generate_tax_table("38156-UA_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "38156_UA_no_rRNA")

k32976_UA_no_rRNA <- generate_tax_table("32976-UA_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "32976_UA_no_rRNA")

k32989_UA_no_rRNA <- generate_tax_table("32989-UA_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "32989_UA_no_rRNA")

k36397_UB_no_rRNA <- generate_tax_table("36397-UB_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "36397_UB_no_rRNA")

k36397_UB_no_rRNA <- generate_tax_table("38156-UB_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "36397_UB_no_rRNA")

k32976_UB_no_rRNA <- generate_tax_table("32976-UB_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "32976_UB_no_rRNA")

k34970_FA_no_rRNA <- generate_tax_table("34970-FA_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "34970_FA_no_rRNA")

k34970_RA_no_rRNA <- generate_tax_table("34970-RA_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "34970_RA_no_rRNA")

k22015_UA_no_rRNA <- generate_tax_table("22015-UA_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "22015_UA_no_rRNA")

k20040_UA_no_rRNA <- generate_tax_table("20040-UA_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "20040_UA_no_rRNA")

k41596_UA_no_rRNA <- generate_tax_table("41596-UA_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "41596_UA_no_rRNA")

k34970_UA_no_rRNA <- generate_tax_table("34970-UA_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "34970_UA_no_rRNA")

k34970_UD_no_rRNA <- generate_tax_table("34970-UD_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "34970_UD_no_rRNA")

k41596_UD_no_rRNA <- generate_tax_table("41596-UD_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "41596_UD_no_rRNA")

k20040_UD_no_rRNA <- generate_tax_table("20040-UD_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "20040_UD_no_rRNA")

k22015_UD_no_rRNA <- generate_tax_table("22015-UD_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "22015_UD_no_rRNA")

k34970_RD_no_rRNA <- generate_tax_table("34970-RD_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "34970_RD_no_rRNA")

k34970_FD_no_rRNA <- generate_tax_table("34970-FD_k2report.tsv") %>%
  fill_tax_NAs() %>%
  mutate(SampleID = "34970_FD_no_rRNA")

all_vs_pluspfp_db_no_rRNA <- bind_rows(k20040_UA_no_rRNA,k20040_UD_no_rRNA,k22015_UA_no_rRNA,k22015_UD_no_rRNA,
                          k32976_UA_no_rRNA,k32976_UB_no_rRNA,k32989_UA_no_rRNA,k34970_FA_no_rRNA,
                          k34970_FD_no_rRNA,k34970_RA_no_rRNA,k34970_RD_no_rRNA,k34970_UA_no_rRNA,
                          k34970_UD_no_rRNA,k36397_FA_no_rRNA,k36397_FB_no_rRNA,k36397_RA_no_rRNA,
                          k36397_RB_no_rRNA,k36397_UA_no_rRNA,k36397_UB_no_rRNA,k38156_UA_no_rRNA,
                          k41596_UA_no_rRNA,k41596_UD_no_rRNA)

all_vs_pluspfp_db_no_rRNA$rrna_removed <- "rRNA removed"

plot_subset(all_vs_pluspfp_db_no_rRNA, "RA", "U", "unclassified", "D") +
  labs(fill = "Domain") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

plot_subset(filter(all_vs_pluspfp_db_no_rRNA, G != "Homo"), "RA", "D", "Eukaryota", "F") +
  labs(fill = "Family") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


rRNA_YesNo <- bind_rows(all_vs_pluspfp, all_vs_pluspfp_db_no_rRNA)


ggplot(filter(rRNA_YesNo, R == "root"), aes(x = rrna_removed, y = RA, fill = D)) + geom_bar(position = 'stack', stat = 'identity') +
  ylab("Relative abundance of reads") +
  xlab("sampleID") + 
  theme_bw() +
  labs(fill = "Domain") + 
  theme(text = element_text(size = 15)) +
  guides(fill = guide_legend(ncol = 1)) +
  facet_grid(~SampleID, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

