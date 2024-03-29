library(tidyverse)

setwd("C:/Users/fanta/Dropbox/PC (2)/Documents/OBrienMetagenome/rpip_panels/deeparg_predictions/")

deeparg001 <- read.csv("001-36397-FA-RP_S265_L005_R1_merged__merged_deeparg_ss.out.mapping.ARG", header = TRUE, sep = "\t")

rm(deeparg_reports_rpip) #(If it is still in the workspace, re-importing will just duplicate the data)

for (data in list.files()){ #List files in current directory as "data"
  
  # Create the initial df if none exists yet
  if (!exists("deeparg_reports_rpip")){ 
    deeparg_reports_rpip <- read.csv(data, header = TRUE, sep = "\t") %>%
      mutate(SampleID = substr(data, 1, 12)) %>%
      mutate(SampleID = substr(data, 1, 12)) %>%
      mutate(treatment = substr(data, 11, 12))  
  }
  
  # if df already exists, then append new data 
  if (exists("deeparg_reports_rpip")){
    temporary <- read.csv(data, header = TRUE, sep = "\t") %>%
      mutate(SampleID = substr(data, 1, 12)) %>%
      mutate(treatment = substr(data, 11, 12))
    deeparg_reports_rpip <- bind_rows(deeparg_reports_rpip, temporary) 
    #dplyr row binding allows different numbers of columns
    rm(temporary) #save memory <3 
  }
}







setwd("C:/Users/fanta/Dropbox/PC (2)/Documents/OBrienMetagenome/pre_enrichment/deeparg_predictions")

rm(deeparg_reports_pre_enrichment) #Avoid duplicating
for (data in list.files()){ #List files in current directory as "data"
  
  # Create the initial df if none exists yet
  if (!exists("deeparg_reports_pre_enrichment")){ 
    deeparg_reports_pre_enrichment <- read.csv(data, header = TRUE, sep = "\t") %>%
      mutate(SampleID = substr(data, 1, 12)) %>% #This is 1-indexed, not 0-indexed
      mutate(treatment = substr(data, 11, 12)) %>%
      mutate(X.ARG = as.character(X.ARG)) %>%
      mutate(read_id = as.character(read_id)) %>%
      mutate(predicted_ARG.class = as.character(predicted_ARG.class)) %>%
      mutate(best.hit = as.character(best.hit))  
  }
  
  # if df already exists, then append new data 
  if (exists("deeparg_reports_pre_enrichment")){
    temporary <- read.csv(data, header = TRUE, sep = "\t") %>%
      mutate(SampleID = substr(data, 1, 12)) %>%
      mutate(treatment = substr(data, 11, 12)) %>%
      mutate(X.ARG = as.character(X.ARG)) %>%
      mutate(read_id = as.character(read_id)) %>%
      mutate(predicted_ARG.class = as.character(predicted_ARG.class)) %>%
      mutate(best.hit = as.character(best.hit))
    deeparg_reports_pre_enrichment <- bind_rows(deeparg_reports_pre_enrichment, temporary) 
    #dplyr row binding allows different numbers of columns
    rm(temporary) #save memory <3 
  }
}

ggplot(deeparg_reports_rpip, aes(x = SampleID, fill = predicted_ARG.class)) + geom_bar(position = "stack", stat = "count") +
  facet_wrap(~treatment)
ggplot(deeparg_reports_pre_enrichment, aes(x = SampleID, fill = predicted_ARG.class)) + geom_bar(position = "stack", stat = "count") +
  facet_wrap(~treatment)

################################################################################
###Add total read counts 
setwd("C:/Users/fanta/Dropbox/PC (2)/Documents/OBrienMetagenome/")
rpip_read_counts <- read.csv("rpip_panels/rrna_filtered_read_counts.csv", header = FALSE, col.names = c("id", "count")) %>%
  mutate(SampleID = substr(id, 1, 12))

pre_enrichment_read_counts <- read.csv("pre_enrichment/rrna_filtered_read_counts.csv", header = FALSE, col.names = c("id", "count")) %>%
  mutate(SampleID = substr(id, 1, 12))



################################################################################
###Summarize counts/sample

pre_enrichment_arg_counts <- deeparg_reports_pre_enrichment %>%
  group_by(SampleID, X.ARG) %>%
  summarize(sampleCount = n()) %>%
  full_join(pre_enrichment_read_counts, by = "SampleID") %>%
  mutate(RA = sampleCount/count)

pre_enrichment_target_counts <- deeparg_reports_pre_enrichment %>%
  group_by(SampleID, predicted_ARG.class) %>%
  summarize(sampleCount = n()) %>%
  full_join(pre_enrichment_read_counts, by = "SampleID") %>%
  mutate(RA = sampleCount/count)

rpip_arg_counts <- deeparg_reports_rpip %>%
  group_by(SampleID, X.ARG) %>%
  summarize(sampleCount = n()) %>%
  full_join(rpip_read_counts, by = "SampleID") %>%
  mutate(RA = sampleCount/count)

rpip_target_counts <- deeparg_reports_rpip %>%
  group_by(SampleID, predicted_ARG.class) %>%
  summarize(sampleCount = n()) %>%
  full_join(rpip_read_counts, by = "SampleID") %>%
  mutate(RA = sampleCount/count)



################################################################################
###Comparing enriched to non-enriched (targets)
rpip_targets = c("aminoglycoside", "beta-lactam", "carbapenem", "cephalosporin", 
                 "diaminopyrimidine", "fluoroquinolone", "fosfomycin", "glycopeptide",
                 "lincosamides", "macrolides", "MLS", "oxazolidinone", "penicillin", 
                 "polymyxin", "sulfonamide", "tetracycline", "isoniazid", 
                 "polyamine:peptide", "pyrazinamide", "rifamycin", "ethionamide", 
                 "para-aminosalisylic acid")

pre_enrichment_target_summary <- pre_enrichment_target_counts %>%
  group_by(predicted_ARG.class) %>%
  summarize(mean = mean(RA), sd = sd(RA), iqr = IQR(RA)) %>%
  mutate(treatment = "No enrichment", rpip_target = (predicted_ARG.class %in% rpip_targets)) %>%
  mutate(rpip_target = factor(rpip_target, levels = c("TRUE", "FALSE"))) %>%
  mutate(rpip_target = recode(rpip_target, "TRUE" = "RPIP targets", "FALSE" = "Non-targets"))

rpip_target_summary <- rpip_target_counts %>%
  group_by(predicted_ARG.class) %>%
  summarize(mean = mean(RA), sd = sd(RA), iqr = IQR(RA, na.rm = TRUE)) %>%
  mutate(treatment = "RPIP", rpip_target = (predicted_ARG.class %in% rpip_targets)) %>%
  mutate(rpip_target = factor(rpip_target, levels = c("TRUE", "FALSE"))) %>%
  mutate(rpip_target = recode(rpip_target, "TRUE" = "RPIP targets", "FALSE" = "Non-targets"))

complete_target_summary <- rbind(pre_enrichment_target_summary, rpip_target_summary)

ggplot(complete_target_summary, aes(x = predicted_ARG.class, y = mean, col = treatment)) + 
  geom_point() + 
  #geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd))+
  facet_grid(~rpip_target, scales = "free_x", space = "free") +
  scale_y_log10() +
  scale_color_manual(values = c("blue", "orange")) +
  ylab("mean relative abundance") +
  xlab("") +
  theme(axis.text.x = element_text(angle = -60, hjust = 0, vjust = 0.5))
################################################################################
###Comparing enriched to non-enriched (genes)
pre_enrichment_summary <- pre_enrichment_arg_counts %>%
  group_by(X.ARG) %>%
  summarize(mean = mean(RA), sd = sd(RA), iqr = IQR(RA)) %>%
  mutate(treatment = "Pre-Enrichment")

rpip_summary <- rpip_arg_counts %>%
  group_by(X.ARG) %>%
  summarize(mean = mean(RA), sd = sd(RA), iqr = IQR(RA, na.rm = TRUE)) %>%
  mutate(treatment = "RPIP")

complete_summary <- rbind(pre_enrichment_summary, rpip_summary)

comparison_df <- complete_summary %>%
  pivot_wider(id_cols = X.ARG, names_from = treatment, values_from = c(mean, sd, iqr)) %>%
  mutate(difference = mean_RPIP - `mean_Pre-Enrichment`, na.rm = TRUE) %>%
  arrange(difference)

ggplot(complete_summary, aes(x = X.ARG, y = mean, col = treatment)) + geom_point() +
  scale_y_log10() +
  scale_color_manual(values = c("blue", "orange")) +
  theme(axis.text.x = element_blank()) +
  xlab("") +
  ylab("relative abundances")

################################################################################

sulfonamides <- filter(complete_summary, predicted_ARG.class == "sulfonamide")
aminoglycosides <- filter(complete_summary, predicted_ARG.class == "aminoglycoside")

ggplot(aminoglycosides, aes(x = SampleID, fill = X.ARG)) + geom_bar(position = "stack", stat = "count")


head(complete_summary$X.ARG)
