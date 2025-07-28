library(tidyverse)

nrow(family_box_df)
nrow(family_box_df_90)

family_box_df_00 <- family_box_df %>%
  filter(Kraken2_confidence == "0")


wider_enrichment_00 <- family_box_df_00 %>%
  select(`F`, readcount, Nanotrap_type, ribosomal, RA, Enrichment, LIMS_ID, site, Fraction) %>%
  pivot_wider(names_from = Enrichment, values_from = c("RA", "readcount"))

wider_enrichment_00[is.na(wider_enrichment_00$RA_RPIP),]$RA_RPIP <- 0
wider_enrichment_00[is.na(wider_enrichment_00$RA_VSP),]$RA_VSP <- 0
wider_enrichment_00[is.na(wider_enrichment_00$RA_None),]$RA_None <- 0

wider_enrichment_00 <- wider_enrichment_00 %>%
  mutate(RPIP_untarg_RA_ratio = RA_RPIP/RA_None) %>%
  mutate(VSP_untarg_RA_ratio = RA_VSP/RA_None)



widen_and_calculate_ratios <- function(df, names_from) {
  wider <- df %>%
    select(`F`, readcount, Nanotrap_type, ribosomal, RA, Enrichment, LIMS_ID, 
           site, Fraction) %>%
    pivot_wider(names_from = names_from, values_from = c("RA", "readcount"))
  return(wider)
}


#Enrichment ratio RPIP

rpip_all_taxa_wider_enrichment_00_100plus <- wider_enrichment_00 %>%
  filter(readcount_RPIP >= 100)

rpip_all_taxa_wider_enrichment_00_mean_ratios <- rpip_all_taxa_wider_enrichment_00_100plus %>%
  group_by(`F`) %>%
  summarize(mean_rpip_ratio = mean(RPIP_untarg_RA_ratio, na.rm = TRUE), mean_rpip_reads = mean(readcount_RPIP))

at_least_doubled_RPIP <- filter(rpip_all_taxa_wider_enrichment_00_mean_ratios, mean_rpip_ratio >= 2)


ggplot(filter(rpip_all_taxa_wider_enrichment_00_100plus, `F` %in% at_least_doubled_RPIP$`F` & `F` %in% rpip_targets$family), aes(y = `F`, x = RPIP_untarg_RA_ratio, fill = Fraction, colour = Fraction)) + geom_boxplot() +
  scale_x_log10() +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12)) +
  #  facet_grid(rows = "plotting_cat", scales = "free_y", space = "free_y") +
  xlab("Ratio of relative abundances") +
  ylab("") +
  theme(strip.placement = "outside", strip.text.y = element_text(angle = 0, size = 10), 
        strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom") +
  theme(panel.grid.major.y = element_blank()) +
  geom_hline(yintercept=seq(1,nrow(boxplot_df_by_arg[,1]),1)+.5,color="black", linetype = "dotted", linewidth = .5) +
  geom_hline(yintercept=seq(1,nrow(boxplot_df_by_arg[,1]),1)-.5,color="black", linetype = "dotted", linewidth = .5) +
  scale_fill_manual(name = "Fraction:", values = c("lightpink3", "aquamarine3", "goldenrod1")) +
  scale_colour_manual(name = "Fraction:", values = c("lightpink4", "aquamarine4", "goldenrod3")) +
  facet_grid(rows = "Nanotrap_type")




#Enrichment ratio VSP

vsp_all_taxa_wider_enrichment_00_100plus <- wider_enrichment_00 %>%
  filter(readcount_VSP >= 100)

vsp_all_taxa_wider_enrichment_00_mean_ratios <- vsp_all_taxa_wider_enrichment_00_100plus %>%
  group_by(`F`) %>%
  summarize(mean_vsp_ratio = mean(VSP_untarg_RA_ratio, na.rm = TRUE), mean_vsp_reads = mean(readcount_VSP))

at_least_doubled_VSP <- filter(vsp_all_taxa_wider_enrichment_00_mean_ratios, mean_vsp_ratio >= 2)


ggplot(filter(vsp_all_taxa_wider_enrichment_00_100plus, `F` %in% at_least_doubled_VSP$`F` & `F` %in% vsp_targets$family_name), aes(y = `F`, x = VSP_untarg_RA_ratio, fill = Fraction, colour = Fraction)) + geom_boxplot() +
  scale_x_log10() +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12)) +
  #  facet_grid(rows = "plotting_cat", scales = "free_y", space = "free_y") +
  xlab("Ratio of relative abundances") +
  ylab("") +
  theme(strip.placement = "outside", strip.text.y = element_text(angle = 0, size = 10), 
        strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom") +
  theme(panel.grid.major.y = element_blank()) +
  geom_hline(yintercept=seq(1,nrow(boxplot_df_by_arg[,1]),1)+.5,color="black", linetype = "dotted", linewidth = .5) +
  geom_hline(yintercept=seq(1,nrow(boxplot_df_by_arg[,1]),1)-.5,color="black", linetype = "dotted", linewidth = .5) +
  scale_fill_manual(name = "Fraction:", values = c("lightpink3", "aquamarine3", "goldenrod1")) +
  scale_colour_manual(name = "Fraction:", values = c("lightpink4", "aquamarine4", "goldenrod3")) +
  facet_grid(rows = "Nanotrap_type")


ggplot(filter(vsp_all_taxa_wider_enrichment_00_100plus, `F` %in% at_least_doubled_VSP$`F` & `F` %in% vsp_targets$family_name), aes(y = `F`, x = RA_VSP - RA_None, fill = Fraction, colour = Fraction)) + geom_boxplot() +
  scale_x_log10() +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12)) +
  #  facet_grid(rows = "plotting_cat", scales = "free_y", space = "free_y") +
  xlab("Ratio of relative abundances") +
  ylab("") +
  theme(strip.placement = "outside", strip.text.y = element_text(angle = 0, size = 10), 
        strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom") +
  theme(panel.grid.major.y = element_blank()) +
  geom_hline(yintercept=seq(1,nrow(boxplot_df_by_arg[,1]),1)+.5,color="black", linetype = "dotted", linewidth = .5) +
  geom_hline(yintercept=seq(1,nrow(boxplot_df_by_arg[,1]),1)-.5,color="black", linetype = "dotted", linewidth = .5) +
  scale_fill_manual(name = "Fraction:", values = c("lightpink3", "aquamarine3", "goldenrod1")) +
  scale_colour_manual(name = "Fraction:", values = c("lightpink4", "aquamarine4", "goldenrod3")) +
  facet_grid(rows = "Nanotrap_type")


#############33
#Plot VSP with grid
ggplot(filter(family_box_df_00, `F` %in% at_least_doubled_VSP$`F` & `F` %in% vsp_targets$family_name, Enrichment != "RPIP"), aes(y = `F`, x = RA, fill = Enrichment, colour = Enrichment)) + geom_boxplot() +
  scale_x_log10() +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12)) +
  #  facet_grid(rows = "plotting_cat", scales = "free_y", space = "free_y") +
  xlab("Relative abundance of reads") +
  ylab("") +
  theme(strip.placement = "outside", strip.text.y = element_text(angle = 0, size = 10), 
        strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(face = "bold", size = 15),
        legend.position = "bottom") +
  theme(panel.grid.major.y = element_blank()) +
  geom_hline(yintercept=seq(1,nrow(boxplot_df_by_arg[,1]),1)+.5,color="black", linetype = "dotted", linewidth = .5) +
  geom_hline(yintercept=seq(1,nrow(boxplot_df_by_arg[,1]),1)-.5,color="black", linetype = "dotted", linewidth = .5) +
  scale_fill_manual(name = "Enrichment:", values = c("lightpink3", "aquamarine3", "goldenrod1")) +
  scale_colour_manual(name = "Enrichment:", values = c("lightpink4", "aquamarine4", "goldenrod3")) +
  facet_grid(rows = vars(Fraction), cols = vars(Nanotrap_type))

#Plot RPIP with grid
ggplot(filter(family_box_df_00, `F` %in% at_least_doubled_RPIP$`F` & `F` %in% rpip_targets$family, Enrichment != "VSP"), aes(y = `F`, x = RA, fill = Enrichment, colour = Enrichment)) + geom_boxplot() +
  scale_x_log10() +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12)) +
  #  facet_grid(rows = "plotting_cat", scales = "free_y", space = "free_y") +
  xlab("Relative abundance of reads") +
  ylab("") +
  theme(strip.placement = "outside", strip.text.y = element_text(angle = 0, size = 10), 
        strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom") +
  theme(panel.grid.major.y = element_blank()) +
  geom_hline(yintercept=seq(1,nrow(boxplot_df_by_arg[,1]),1)+.5,color="black", linetype = "dotted", linewidth = .5) +
  geom_hline(yintercept=seq(1,nrow(boxplot_df_by_arg[,1]),1)-.5,color="black", linetype = "dotted", linewidth = .5) +
  scale_fill_manual(name = "Enrichment:", values = c("lightpink3","goldenrod1")) +
  scale_colour_manual(name = "Enrichment:", values = c("lightpink4", "goldenrod3")) +
  facet_grid(rows = vars(Fraction), cols = vars(Nanotrap_type))

################################################################################
#Taxa only in RPIP
rpip_only <- filter(wider_enrichment_00, RA_RPIP != 0 & RA_None == 0) %>%
  filter(readcount_RPIP >= 100)

vsp_only <- filter(wider_enrichment_00, RA_VSP != 0 & RA_None == 0) %>%
  filter(readcount_VSP >= 100) %>%
  filter(`F` %in% vsp_targets$family_name)

ggplot(vsp_only, aes(x = RA_VSP, y = `F`, fill = Fraction)) + geom_boxplot() +
  scale_x_log10() +
  facet_grid(rows = "Nanotrap_type")
  


######All targeted taxa with filtration and Nano
all_targeted_00 <- family_box_df_00 %>%
  filter(`F` %in% rpip_targets$family | `F` %in% vsp_targets$family_name) %>%
  ungroup() %>%
  group_by(SampleID, LIMS_ID, Treatment, Fraction, Nanotrap_type, Enrichment) %>%
  summarize(total_targeted_RA = sum(RA))

ggplot(all_targeted_00, aes(x = Enrichment, y = total_targeted_RA, fill = Treatment)) + geom_boxplot()



#############
#Same for nano
wider_nano_00 <- family_box_df_00 %>%
  ungroup() %>%
  select(`F`, readcount, Nanotrap_type, ribosomal, RA, Enrichment, LIMS_ID, site, Fraction) %>%
  pivot_wider(names_from = Nanotrap_type, values_from = c("RA", "readcount"))

wider_nano_00[is.na(wider_nano_00$RA_A),]$RA_A <- 0
wider_nano_00[is.na(wider_nano_00$`RA_A&B`),]$`RA_A&B` <- 0
wider_nano_00[is.na(wider_nano_00$RA_none),]$RA_none <- 0

wider_nano_00 <- wider_nano_00 %>%
  mutate(A_untarg_RA_ratio = RA_A/RA_none) %>%
  mutate(BA_untarg_RA_ratio = `RA_A&B`/RA_none)


wider_nano_00$RPIP_Targeted <- wider_nano_00$`F` %in% rpip_targets$family
wider_nano_00$VSP_Targeted <- wider_nano_00$`F` %in% vsp_targets$family_name
wider_nano_00 <- wider_nano_00 %>%
  mutate(Targeted_by = case_when(RPIP_Targeted & VSP_Targeted ~ "Both",
                                 RPIP_Targeted & !VSP_Targeted ~ "RPIP",
                                 !RPIP_Targeted & VSP_Targeted ~ "VSP",
                                 !(RPIP_Targeted | VSP_Targeted) ~ "Neither"))
#nano ratio

ntA_all_taxa_wider_nano_00_100plus <- wider_nano_00 %>%
  filter(readcount_A >= 100)

ntA_all_taxa_wider_nano_00_mean_ratios <- ntA_all_taxa_wider_nano_00_100plus %>%
  group_by(`F`) %>%
  summarize(mean_ntA_ratio = mean(A_untarg_RA_ratio, na.rm = TRUE), mean_ntA_reads = mean(readcount_A))

at_least_doubled_ntA <- filter(ntA_all_taxa_wider_nano_00_mean_ratios, mean_ntA_ratio >= 2)
biggest_diff_ntA <- ntA_all_taxa_wider_nano_00_mean_ratios %>%
  filter(!str_starts(`F`, "unid")) %>%
  arrange(desc(abs(log(mean_ntA_ratio)))) %>%
  `[`(1:10,)



ggplot(filter(ntA_all_taxa_wider_nano_00_100plus, `F` %in% at_least_doubled_ntA$`F` & (`F` %in% vsp_targets$family_name | `F` %in% rpip_targets$family)), aes(y = `F`, x = A_untarg_RA_ratio, fill = Targeted_by)) + geom_boxplot() +
  scale_x_log10() +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12)) +
  #  facet_grid(rows = "plotting_cat", scales = "free_y", space = "free_y") +
  xlab("Ratio of relative abundances") +
  ylab("") +
  theme(strip.placement = "outside", strip.text.y = element_text(angle = 0, size = 10), 
        strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom") +
  theme(panel.grid.major.y = element_blank()) +
  geom_hline(yintercept=seq(1,nrow(boxplot_df_by_arg[,1]),1)+.5,color="black", linetype = "dotted", linewidth = .5) +
  geom_hline(yintercept=seq(1,nrow(boxplot_df_by_arg[,1]),1)-.5,color="black", linetype = "dotted", linewidth = .5) +
  scale_fill_manual(name = "Targeted by:", values = c("lightpink3", "aquamarine3", "goldenrod1")) +
  scale_colour_manual(name = "Targeted by:", values = c("lightpink4", "aquamarine4", "goldenrod3"))




#Top 10 differences Nanotrap
ggplot(filter(family_box_df_00, `F` %in% biggest_diff_ntA$`F`), aes(y = `F`, x = RA, fill = Nanotrap_type, colour = Nanotrap_type)) + geom_boxplot() +
  scale_x_log10() +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12)) +
  #  facet_grid(rows = "plotting_cat", scales = "free_y", space = "free_y") +
  xlab("Relative abundance of reads") +
  ylab("") +
  theme(strip.placement = "outside", strip.text.y = element_text(angle = 0, size = 10), 
        strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom") +
  theme(panel.grid.major.y = element_blank()) +
  geom_hline(yintercept=seq(1,nrow(boxplot_df_by_arg[,1]),1)+.5,color="black", linetype = "dotted", linewidth = .5) +
  geom_hline(yintercept=seq(1,nrow(boxplot_df_by_arg[,1]),1)-.5,color="black", linetype = "dotted", linewidth = .5) +
  scale_fill_manual(name = "Nanotrap type:", values = c("lightpink3", "aquamarine3", "goldenrod1")) +
  scale_colour_manual(name = "Nanotrap type:", values = c("lightpink4", "aquamarine4", "goldenrod3"))

#Enrichment ratio VSP

vsp_all_taxa_wider_enrichment_00_100plus <- wider_enrichment_00 %>%
  filter(readcount_VSP >= 100)

vsp_all_taxa_wider_enrichment_00_mean_ratios <- vsp_all_taxa_wider_enrichment_00_100plus %>%
  group_by(`F`) %>%
  summarize(mean_vsp_ratio = mean(VSP_untarg_RA_ratio, na.rm = TRUE), mean_vsp_reads = mean(readcount_VSP))

at_least_doubled_VSP <- filter(vsp_all_taxa_wider_enrichment_00_mean_ratios, mean_vsp_ratio >= 2)


ggplot(filter(vsp_all_taxa_wider_enrichment_00_100plus, `F` %in% at_least_doubled_VSP$`F` & `F` %in% vsp_targets$family_name), aes(y = `F`, x = VSP_untarg_RA_ratio, fill = Fraction, colour = Fraction)) + geom_boxplot() +
  scale_x_log10() +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12)) +
  #  facet_grid(rows = "plotting_cat", scales = "free_y", space = "free_y") +
  xlab("Ratio of relative abundances") +
  ylab("") +
  theme(strip.placement = "outside", strip.text.y = element_text(angle = 0, size = 10), 
        strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom") +
  theme(panel.grid.major.y = element_blank()) +
  geom_hline(yintercept=seq(1,nrow(boxplot_df_by_arg[,1]),1)+.5,color="black", linetype = "dotted", linewidth = .5) +
  geom_hline(yintercept=seq(1,nrow(boxplot_df_by_arg[,1]),1)-.5,color="black", linetype = "dotted", linewidth = .5) +
  scale_fill_manual(name = "Fraction:", values = c("lightpink3", "aquamarine3", "goldenrod1")) +
  scale_colour_manual(name = "Fraction:", values = c("lightpink4", "aquamarine4", "goldenrod3")) +
  facet_grid(rows = "Nanotrap_type")


ggplot(filter(vsp_all_taxa_wider_enrichment_00_100plus, `F` %in% at_least_doubled_VSP$`F` & `F` %in% vsp_targets$family_name), aes(y = `F`, x = RA_VSP - RA_None, fill = Fraction, colour = Fraction)) + geom_boxplot() +
  scale_x_log10() +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12)) +
  #  facet_grid(rows = "plotting_cat", scales = "free_y", space = "free_y") +
  xlab("Ratio of relative abundances") +
  ylab("") +
  theme(strip.placement = "outside", strip.text.y = element_text(angle = 0, size = 10), 
        strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom") +
  theme(panel.grid.major.y = element_blank()) +
  geom_hline(yintercept=seq(1,nrow(boxplot_df_by_arg[,1]),1)+.5,color="black", linetype = "dotted", linewidth = .5) +
  geom_hline(yintercept=seq(1,nrow(boxplot_df_by_arg[,1]),1)-.5,color="black", linetype = "dotted", linewidth = .5) +
  scale_fill_manual(name = "Fraction:", values = c("lightpink3", "aquamarine3", "goldenrod1")) +
  scale_colour_manual(name = "Fraction:", values = c("lightpink4", "aquamarine4", "goldenrod3")) +
  facet_grid(rows = "Nanotrap_type")


#############33
#Plot VSP with grid
ggplot(filter(family_box_df_00, `F` %in% at_least_doubled_VSP$`F` & `F` %in% vsp_targets$family_name, Enrichment != "RPIP"), aes(y = `F`, x = RA, fill = Enrichment, colour = Enrichment)) + geom_boxplot() +
  scale_x_log10() +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12)) +
  #  facet_grid(rows = "plotting_cat", scales = "free_y", space = "free_y") +
  xlab("Relative abundance of reads") +
  ylab("") +
  theme(strip.placement = "outside", strip.text.y = element_text(angle = 0, size = 10), 
        strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom") +
  theme(panel.grid.major.y = element_blank()) +
  geom_hline(yintercept=seq(1,nrow(boxplot_df_by_arg[,1]),1)+.5,color="black", linetype = "dotted", linewidth = .5) +
  geom_hline(yintercept=seq(1,nrow(boxplot_df_by_arg[,1]),1)-.5,color="black", linetype = "dotted", linewidth = .5) +
  scale_fill_manual(name = "Enrichment:", values = c("lightpink3", "aquamarine3", "goldenrod1")) +
  scale_colour_manual(name = "Enrichment:", values = c("lightpink4", "aquamarine4", "goldenrod3")) +
  facet_grid(rows = vars(Fraction), cols = vars(Nanotrap_type))



#############################################3
#Taxa barplot
ggplot(filter(family_box_df_00, `F` %in% vsp_targets$family_name | `F` %in% rpip_targets$family), aes(x = Fraction, y = RA, fill = `F`)) + geom_bar(stat = "identity", position = "fill") +
  theme(legend.position = "none")


