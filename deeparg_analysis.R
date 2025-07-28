#library(stringdist)

setwd("/projects/bios_microbe/cowen20/targeted_panels/rpip_panels/raw_fastqs/fastp_merged_no_dedup/deeparg_out/")

#Start from a clear slate:
rm(rpip_deeparg_results_table)

#Import RPIP ARG results:
for (ARG in list.files(pattern = "*mapping.ARG")) {
  print(ARG)
  arg_report <- read_tsv(ARG)
  arg_report$sampleName <- ARG #(The file name)
  if (!exists("rpip_deeparg_results_table")) {
    rpip_deeparg_results_table <- arg_report
  } else {
    rpip_deeparg_results_table <- rbind(rpip_deeparg_results_table, arg_report)
  }
}

#Add enrichment column
rpip_deeparg_results_table$Enrichment <- "RPIP"

#Remove the weird character
colnames(rpip_deeparg_results_table)[colnames(rpip_deeparg_results_table) == "#ARG"] <- "ARG"

#Convert "sampleName" into separate columns
rpip_deeparg_results_table <- separate(rpip_deeparg_results_table,
                                 col = sampleName,
                                 into = c("QCSeqID", "LIMS_ID", "Treatment",
                                          NA, NA, NA, NA, NA, NA, NA, NA),
                                 sep = "(-|_)")

#Function to convert 2-letter codes into readable treatment info
parse_sample_treatments <- function(tax_table, treatment_col = "Treatment") {
  parsed <- tax_table %>%
    mutate(Fraction = str_replace_all(
      !!sym(treatment_col),
      c("^F.$" = "filtrate",
        "^R.$" = "retentate",
        "^U.$" = "unfiltered")
    )) %>%
    mutate(Nanotrap_type = str_replace_all(
      !!sym(treatment_col),
      c("^.A$" = "A",
        "^.B$" = "A&B",
        "^.D$" = "none")
    ))
}

#Convert the 2-letter codes
rpip_deeparg_results_table <- parse_sample_treatments(rpip_deeparg_results_table)

#Indicate which sample site each sample came from
parse_locations <- function(taxtable) {
  mutate(taxtable, site = str_replace_all(
    LIMS_ID,
    c(
      "(36397|38156|41596|34970)" = "OBrien", #regex for any one of the given strings
      "(32976|32989|20040|22015)" = "OHare"
    )
  )
  )
}

#Add location info to the df
rpip_deeparg_results_table <- parse_locations(rpip_deeparg_results_table)


#For untargeted sequencing:
setwd("/projects/bios_microbe/cowen20/targeted_panels/untargeted/raw_fastqs/fastp_merged_no_dedup/deeparg_out/")

#Start from clear slate:
rm(untargeted_deeparg_results_table)

#Import the unargeted deeparg reports:
for (ARG in list.files(pattern = "*mapping.ARG")) {
  print(ARG)
  arg_report <- read_tsv(ARG) #(The file name)
  arg_report$sampleName <- ARG #(The file name is the sample name)
  if (!exists("untargeted_deeparg_results_table")) {
    untargeted_deeparg_results_table <- arg_report
  } else {
    untargeted_deeparg_results_table <- rbind(untargeted_deeparg_results_table, arg_report)
  }
}

#Label untargeted as no enrichment
untargeted_deeparg_results_table$Enrichment <- "None"

#Remove weird character:
colnames(untargeted_deeparg_results_table)[colnames(untargeted_deeparg_results_table) == "#ARG"] <- "ARG"

#Convert sampleName column into separate columns for relevant info:
untargeted_deeparg_results_table <- separate(untargeted_deeparg_results_table,
                                       col = sampleName,
                                       into = c("QCSeqID", "LIMS_ID", "Treatment",
                                                NA, NA, NA, NA, NA, NA, NA, NA),
                                       sep = "(-|_)")

#Function to convert 2-letter codes into readable treatments (don't need to 
#re-define if already defined during RPIP import)
parse_sample_treatments <- function(tax_table, treatment_col = "Treatment") {
  parsed <- tax_table %>%
    mutate(Fraction = str_replace_all(
      !!sym(treatment_col),
      c("^F.$" = "filtrate",
        "^R.$" = "retentate",
        "^U.$" = "unfiltered")
    )) %>%
    mutate(Nanotrap_type = str_replace_all(
      !!sym(treatment_col),
      c("^.A$" = "A",
        "^.B$" = "A&B",
        "^.D$" = "none")
    ))
}

#Convert the 2-letter codes:
untargeted_deeparg_results_table <- parse_sample_treatments(untargeted_deeparg_results_table)

#Funtion to label sample sites based on IDs (don't need to re-define if already
#run. Just leaving here in case I need to run through importing untargeted dataset again)
parse_locations <- function(taxtable) {
  mutate(taxtable, site = str_replace_all(
    LIMS_ID,
    c(
      "(36397|38156|41596|34970)" = "OBrien", #regex for any one of the given strings
      "(32976|32989|20040|22015)" = "OHare"
    )
  )
  )
}

untargeted_deeparg_results_table <- parse_locations(untargeted_deeparg_results_table)

#combine them
rpip_vs_none_deeparg_results_table <- rbind(rpip_deeparg_results_table, 
                                            untargeted_deeparg_results_table)

#Rename that one weird column 
colnames(rpip_vs_none_deeparg_results_table)[colnames(rpip_vs_none_deeparg_results_table) == "predicted_ARG-class"] <- "predicted_ARG_class"

#Create an ID that contains all treatment information and is unique:
rpip_vs_none_deeparg_results_table <- rpip_vs_none_deeparg_results_table %>%
  mutate(UniqueID = paste0(LIMS_ID, Treatment, Enrichment))

#####IMPORTANT: Run read_fastp_reports.R before continuing#####

#Make matching column to join with deeparg results:
fastp_reports_rpip_and_none <- fastp_reports_rpip_and_none %>%
  mutate(UniqueID = paste0(LIMS_ID, Treatment, Enrichment))

#Sum the number of reads for each ARG in each sample (some of them are repeated
#due to differences in match stats from Diamond alignment)
samplewise_deeparg_results_table <- rpip_vs_none_deeparg_results_table %>% 
  group_by(ARG, predicted_ARG_class, LIMS_ID, Treatment, Enrichment, Fraction, Nanotrap_type, site, UniqueID) %>%
  summarize(NumReads = sum(counts))


#Calculate how many reads are classified as ARGs in each sample:
all_args_per_sample <- samplewise_deeparg_results_table %>%
  group_by(UniqueID) %>%
  summarize(AllArgs = sum(NumReads))

#Join TOTAL ARG counts to fastp reports, so now we have total nreads and total
#ARG-classified reads:
fastp_reports_rpip_and_none <- fastp_reports_rpip_and_none %>%
  left_join(all_args_per_sample, by = "UniqueID")

#Calculate total number of non ARG-classified reads
fastp_reports_rpip_and_none$NonArgReads = fastp_reports_rpip_and_none$after_filtering.total_reads - fastp_reports_rpip_and_none$AllArgs

#Create a df that has the same sample info as the main one, but for NON-ARGs
#in each sample:
non_args_per_sample <- fastp_reports_rpip_and_none %>%
  select(NumReads = NonArgReads, LIMS_ID, Treatment, Enrichment, Fraction,
         Nanotrap_type, site, UniqueID)
non_args_per_sample$ARG = "NonArgReads"
non_args_per_sample$predicted_ARG_class <- "None"

#Combine the ARG and non-ARG dfs, so each sample now has a row containing
#the count of non-ARGs:
all_args_and_nonargs_per_sample <- rbind(samplewise_deeparg_results_table,
                                         non_args_per_sample)
  
#Sanity check: barplot should show majority of each sample as non-ARG reads:
ggplot(all_args_and_nonargs_per_sample, aes(x = UniqueID, y = NumReads, 
                                            fill = predicted_ARG_class)) + 
  geom_bar(position = 'fill', stat = "identity")

#Save the table for later:
write_csv(all_args_and_nonargs_per_sample, "/projects/bios_microbe/cowen20/rprojects/targeted_panels/imported_deeparg_results_df.csv")

################################################################################
#Look at fraction of all ARGs 
arg_fraction_df <- full_join(all_args_per_sample, non_args_per_sample, by = "UniqueID") %>%
  dplyr::rename(NonArg_Reads = NumReads, Arg_Reads = AllArgs)
arg_fraction_df <- arg_fraction_df %>%
  mutate(portion_of_args = Arg_Reads / (Arg_Reads + NonArg_Reads))

arg_fraction_df_wide <- arg_fraction_df %>%
  pivot_wider(names_from = Enrichment, values_from = c(UniqueID, portion_of_args, Arg_Reads, NonArg_Reads))

arg_fraction_df_wide <- arg_fraction_df_wide %>%
  mutate(RPIP_Unta_RA_ratio = portion_of_args_RPIP/portion_of_args_None)

ggplot(arg_fraction_df_wide, aes(x = paste(Fraction, Nanotrap_type, sep = "_"), y = RPIP_Unta_RA_ratio)) + geom_boxplot()
ggplot(arg_fraction_df_wide, aes(x = paste(Fraction, Nanotrap_type, sep = "_"), y = portion_of_args_RPIP)) + geom_boxplot()

arg_fraction_df <- arg_fraction_df %>%
  mutate(LabelCol = paste(str_to_sentence(Fraction), "+", Nanotrap_type)) %>%
  mutate(Label_Enrich = str_replace(Enrichment, "None", "No enrichment"))

ggplot(arg_fraction_df, aes(x = Fraction, y = portion_of_args, fill = Nanotrap_type, colour = Nanotrap_type)) + geom_boxplot() +
  scale_fill_manual(name = "Nanotrap protocol", 
                     values = c("aquamarine2", "goldenrod2", "lightpink2")) +
  scale_colour_manual(name = "Nanotrap protocol",
                     values = c("aquamarine4", "goldenrod4", "lightpink4")) +
  ylab("Relative abundance of all ARGs") +
  theme_minimal() +
  facet_wrap(~Label_Enrich, strip.position = "top") +
  theme(strip.placement = "outside", strip.text.y = element_text(size = 10), 
        strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(size = 10, angle = 30, hjust = 1),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom")
  
####Prepare for DESeq2
deeparg_sample_metadata <- all_args_and_nonargs_per_sample %>%
  ungroup() %>%
  select(UniqueID, LIMS_ID, Treatment, site, Fraction, Nanotrap_type, Enrichment) %>%
  distinct()

deeparg_count_matrix <- all_args_and_nonargs_per_sample %>%
  ungroup() %>%
  mutate(ARG_w_class = paste(ARG, predicted_ARG_class, sep = ":")) %>%
  pivot_wider(id_cols = ARG_w_class,
              names_from = UniqueID,
              values_from = NumReads,
              values_fill = 0)

write_csv(deeparg_sample_metadata, "/projects/bios_microbe/cowen20/rprojects/targeted_panels/deeparg_DS2_sample_metadata.csv")
write_csv(deeparg_count_matrix, "/projects/bios_microbe/cowen20/rprojects/targeted_panels/deeparg_count_matrix_DS2.csv")


###DESeq2 analysis
library(tidyverse)
library(DESeq2)
library(pheatmap)
setwd("/projects/bios_microbe/cowen20/rprojects/targeted_panels/")

deeparg_count_table <- read_csv("deeparg_count_matrix_DS2.csv", guess_max = Inf)

#Convert to matrix
deeparg_count_matrix <- deeparg_count_table %>%
  select(!ARG_w_class) %>%
  as.matrix()
row.names(deeparg_count_matrix) <- deeparg_count_table$ARG_w_class

#Load metadata  
deeparg_sample_metadata <- read.csv("deeparg_DS2_sample_metadata.csv", row.names = 1, stringsAsFactors = TRUE)
deeparg_sample_metadata$LIMS_ID <- factor(deeparg_sample_metadata$LIMS_ID)

#Check sample names in same order
sum(colnames(deeparg_count_matrix) != rownames(deeparg_sample_metadata))

#Perform DESeq2 analysis
deeparg_dds <- DESeqDataSetFromMatrix(countData = deeparg_count_matrix,
                                     colData = deeparg_sample_metadata,
                                     design = ~ LIMS_ID + Enrichment + Fraction + Nanotrap_type)



deeparg_dds <- DESeq(deeparg_dds)

deeparg_result_table <- data.frame(results(deeparg_dds))


#Set reference levels
deeparg_dds$Enrichment <- relevel(deeparg_dds$Enrichment, ref = "None")
deeparg_dds$Fraction <- relevel(deeparg_dds$Fraction, ref = "unfiltered")
deeparg_dds$Nanotrap_type <- relevel(deeparg_dds$Nanotrap_type, ref = "none")

deeparg_dds <- DESeq(deeparg_dds)

View(data.frame(mcols(deeparg_dds)))

#heatmap

dsheatmap <- function(ds2_obj, an_source, an_name) {
  vsd_dds <- varianceStabilizingTransformation(ds2_obj, blind = TRUE)
  vsd_mat_dds <- assay(vsd_dds)
  vsd_cor_dds <- cor(vsd_mat_dds)
  pheatmap(vsd_cor_dds, annotation = select(an_source, an_name))
}

dsheatmap(deeparg_dds, deeparg_sample_metadata, "Enrichment")
dsheatmap(deeparg_dds, deeparg_sample_metadata, "Nanotrap_type")
dsheatmap(deeparg_dds, deeparg_sample_metadata, "Fraction")
dsheatmap(deeparg_dds, deeparg_sample_metadata, "LIMS_ID")


vsd_deeparg_dds <- varianceStabilizingTransformation(deeparg_dds, blind = TRUE)
vsd_mat_deeparg_dds <- assay(vsd_deeparg_dds)
vsd_cor_deeparg_dds <- cor(vsd_mat_deeparg_dds)
pheatmap(vsd_cor_deeparg_dds, annotation = select(deeparg_sample_metadata, Enrichment))

#PCA
plotPCA(vsd_deeparg_dds, intgroup = "LIMS_ID", pcsToUse = c(1,2))
plotPCA(vsd_deeparg_dds, intgroup = "Enrichment")
#plotPCA(vsd_deeparg_dds, intgroup = "site")
plotPCA(vsd_deeparg_dds, intgroup = "Nanotrap_type")
plotPCA(vsd_deeparg_dds, intgroup = "Fraction")

#Dispersion
plotDispEsts(deeparg_dds)


#Contrast RPIP
rpipvsnone_deeparg_res <- results(deeparg_dds,
                                       contrast = c("Enrichment", 
                                                    "RPIP",
                                                    "None"),
                                       alpha = 0.05)
rpipvsnone_deeparg_res <- lfcShrink(deeparg_dds, contrast = c("Enrichment", "RPIP", "None"),
                                         res = rpipvsnone_deeparg_res, type = "ashr")

resultsNames(deeparg_dds)

summary(rpipvsnone_deeparg_res)
View(data.frame(rpipvsnone_deeparg_res))

plotMA(rpipvsnone_deeparg_res)
mcols(rpipvsnone_deeparg_res)


#Contrast filtration
filtratevsunfiltered_deeparg_res <- results(deeparg_dds,
                                  contrast = c("Fraction", 
                                               "filtrate",
                                               "unfiltered"),
                                  alpha = 0.05)
filtratevsunfiltered_deeparg_res <- lfcShrink(deeparg_dds, contrast = c("Fraction", "filtrate", "unfiltered"),
                                    res = filtratevsunfiltered_deeparg_res, type = "ashr")

summary(filtratevsunfiltered_deeparg_res)
View(data.frame(filtratevsunfiltered_deeparg_res))

plotMA(filtratevsunfiltered_deeparg_res)

retentatevsunfiltered_deeparg_res <- results(deeparg_dds,
                                             contrast = c("Fraction",
                                                          "retentate",
                                                          "unfiltered"),
                                             alpha = 0.05)

retentatevsunfiltered_deeparg_res <- lfcShrink(deeparg_dds, contrast = c("Fraction", "retentate", "unfiltered"),
                                               res = retentatevsunfiltered_deeparg_res, type = "ashr")

summary(retentatevsunfiltered_deeparg_res)
View(data.frame(retentatevsunfiltered_deeparg_res))

plotMA(retentatevsunfiltered_deeparg_res)


#Contrast Nanotraps
nanotrapAvsnone_deeparg_res <- results(deeparg_dds,
                                       contrast = c("Nanotrap_type",
                                                    "A",
                                                    "none"),
                                       alpha = 0.05)

nanotrapAvsnone_deeparg_res <- lfcShrink(deeparg_dds, contrast = c("Nanotrap_type", "A", "none"),
                                         res = nanotrapAvsnone_deeparg_res, type = "ashr")


summary(nanotrapAvsnone_deeparg_res)
View(data.frame(nanotrapAvsnone_deeparg_res))
plotMA(nanotrapAvsnone_deeparg_res)


nanotrapBAvsnone_deeparg_res <- results(deeparg_dds,
                                       contrast = c("Nanotrap_type",
                                                    "A&B",
                                                    "none"),
                                       alpha = 0.05)

nanotrapBAvsnone_deeparg_res <- lfcShrink(deeparg_dds, contrast = c("Nanotrap_type", "A&B", "none"),
                                         res = nanotrapBAvsnone_deeparg_res, type = "ashr")


summary(nanotrapBAvsnone_deeparg_res)
View(data.frame(nanotrapBAvsnone_deeparg_res))
plotMA(nanotrapBAvsnone_deeparg_res)




nanotrapBAvsA_deeparg_res <- results(deeparg_dds,
                                        contrast = c("Nanotrap_type",
                                                     "A&B",
                                                     "A"),
                                        alpha = 0.05)

nanotrapBAvsA_deeparg_res <- lfcShrink(deeparg_dds, contrast = c("Nanotrap_type", "A&B", "A"),
                                          res = nanotrapBAvsA_deeparg_res, type = "ashr")


summary(nanotrapBAvsA_deeparg_res)
View(data.frame(nanotrapBAvsA_deeparg_res))
plotMA(nanotrapBAvsA_deeparg_res)



###Look at targeted taxa (RPIP)

#Create a df of DESeq results:
rpip_deeparg_res_df <- data.frame(rpipvsnone_deeparg_res)
rpip_deeparg_res_df$ARG_w_class <- rownames(rpipvsnone_deeparg_res)
rpip_deeparg_res_df$padj[is.na(rpip_deeparg_res_df$padj)] <- 1
rpip_deeparg_res_df <- rpip_deeparg_res_df %>%
  separate(ARG_w_class, into = c("ARG", "class"), remove = FALSE, sep = ":")

#Import ambiguous/unclear category translations:
mls_categories <- read_csv("input/rpip_panels/mls_labeled_categories.csv", col_types = "ccc-")
multidrug_categories <- read_csv("input/rpip_panels/multidrug_labeled_categories.csv")
unclassified_categories <- read_csv("input/rpip_panels/unclassified_labeled_categories.csv")

#Join them together:
all_weird_categories <- rbind(mls_categories, multidrug_categories, unclassified_categories)

#Join them to the main data frame:
rpip_deeparg_res_df <- rpip_deeparg_res_df %>%
  left_join(all_weird_categories, by = join_by(ARG == Gene_name))

#Sanity check
sum(is.na(rpip_deeparg_res_df$class)) == 0

#Copy the class column in case we're about to mess it up:
rpip_deeparg_res_df$original_class <- rpip_deeparg_res_df$class

#Rename the ambiguous/unclear values in the class column:    
rpip_deeparg_res_df[rpip_deeparg_res_df$class == "MLS",]$class <- rpip_deeparg_res_df[rpip_deeparg_res_df$class == "MLS",]$Parent_drug_category
rpip_deeparg_res_df[rpip_deeparg_res_df$class == "multidrug",]$class <- rpip_deeparg_res_df[rpip_deeparg_res_df$class == "multidrug",]$Parent_drug_category
rpip_deeparg_res_df[rpip_deeparg_res_df$class == "unclassified",]$class <- rpip_deeparg_res_df[rpip_deeparg_res_df$class == "unclassified",]$Parent_drug_category


#Open the table that translates Deeparg labels to RPIP labels:
rpip_deeparg_equivalents <- read_csv("input/rpip_panels/rpip_targets_vs_deeparg_detects.csv", col_types = "-ccc")

rpip_deeparg_res_df <- rpip_deeparg_res_df %>%
  left_join(rpip_deeparg_equivalents, by = join_by(class == Deeparg_CARD_category))
rpip_deeparg_res_df$Targeted <- rpip_deeparg_res_df$RPIP_documentation_category != "NONE"
rpip_deeparg_res_df$UpperARG <- toupper(rpip_deeparg_res_df$ARG)

rpip_deeparg_res_df$MultipleTargets <- str_detect(rpip_deeparg_res_df$Parent_drug_category, ";")
rpip_deeparg_res_df$MultipleTargets[is.na(rpip_deeparg_res_df$MultipleTargets)] <- FALSE


ggplot(rpip_deeparg_res_df, 
       aes(x = log2FoldChange, y = -log10(padj), col = Targeted)) + 
  geom_point(aes(shape = MultipleTargets), size = 3, alpha = 0.5) +
  geom_text(aes(label = ifelse(abs(log2FoldChange) > 10 | -log10(padj) > 100, ARG_w_class, ""), vjust = 1.5, angle = -30), size = 3.5, show.legend = FALSE) +
  theme_minimal() +
  xlim(-12, 30) +
  ylim(0, 300) +
  ggtitle("RPIP vs. No Enrichment") +
  xlab(expression(paste("log"[2], " fold change"))) +
  ylab(expression(paste("-log"[10], " adjusted p-value"))) +
  geom_hline(yintercept = 1.30103) +
  annotate("text", x = -8, y = 10, label = "adjusted p-value = 0.05") +
  scale_color_manual(values = c("grey", "#0072B2")) 

#########################Boxplots###############################################
#Make a new df so we don't mess up the main one:
boxplot_df <- all_args_and_nonargs_per_sample

setwd("/projects/bios_microbe/cowen20/rprojects/targeted_panels/")

#This table indicates which targets from RPIP documentation correspond to which
#resistance categories from deeparg results:
rpip_deeparg_equivalents <- read_csv("input/rpip_panels/rpip_targets_vs_deeparg_detects.csv", col_types = "-cc-")

#Create a df for ambiguous/unclear categories from deeparg output:
mls_categories <- read_csv("input/rpip_panels/mls_labeled_categories.csv", 
                           col_types = "ccc-")
multidrug_categories <- read_csv("input/rpip_panels/multidrug_labeled_categories.csv")
unclassified_categories <- read_csv("input/rpip_panels/unclassified_labeled_categories.csv")
weird_categories <- bind_rows(mls_categories, multidrug_categories, unclassified_categories)

#And connect it to the main df:
boxplot_df <- boxplot_df %>%
  left_join(weird_categories, by = join_by(ARG == Gene_name))


#Create a column where ambiguous categories are replaced by looked-up CARD
#categories:
boxplot_df$lookup_ARG_class <- boxplot_df$predicted_ARG_class
boxplot_df[boxplot_df$lookup_ARG_class %in% c("MLS",
                                              "multidrug",
                                              "unclassified"),
           ]$lookup_ARG_class <- boxplot_df[boxplot_df$lookup_ARG_class %in% 
                                              c("MLS", 
                                                "multidrug", 
                                                "unclassified"),
                                            ]$Parent_drug_category

#Now add corresponding RPIP target names:
boxplot_df <- boxplot_df %>%
  left_join(rpip_deeparg_equivalents, by = join_by(lookup_ARG_class == Deeparg_CARD_category))
boxplot_df$Targeted <- boxplot_df$RPIP_documentation_category != "NONE"
boxplot_df <- boxplot_df %>%
  mutate(ARG_w_class = paste(ARG, predicted_ARG_class, sep = ":"))

#Join to DESeq results:
boxplot_df <- boxplot_df %>%
  inner_join(rpip_deeparg_res_df, by = "ARG_w_class")

# Define significance labels based on thresholds
boxplot_df$significance_stars <- ifelse(boxplot_df$padj < 0.001, "***",
                          ifelse(boxplot_df$padj < 0.01, "**",
                                 ifelse(boxplot_df$padj < 0.05, "*", "ns")))

#Make ARG/significance combined (for labeling plots):
boxplot_df <- boxplot_df %>%
  mutate(ARG_sig = paste0(ARG.x, " (", significance_stars, ")")) %>%
  group_by(LIMS_ID, Fraction, Nanotrap_type, Enrichment) %>%
  mutate(RA = NumReads / sum(NumReads))

#Boxplot of CATEGORIES
boxplot_df %>%
  group_by(LIMS_ID, Fraction, Nanotrap_type, Enrichment, RPIP_documentation_category.x) %>%
  filter(Targeted.x == TRUE) %>%
  mutate(plotting_cat = case_when(MultipleTargets == TRUE ~ "Multiple",
                                  .default = RPIP_documentation_category.x)) %>%
  filter(padj < 0.05) %>%
  group_by(LIMS_ID, Fraction, Nanotrap_type, Enrichment, plotting_cat) %>%
  summarize(nreads = NumReads) %>%
  ggplot(aes(x = plotting_cat, y = nreads, colour = Enrichment)) + geom_boxplot() +
  scale_y_log10() +
  theme_minimal() +
  xlab("") +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 0)) 

#Boxplot of ARGs
mean_readcounts <- boxplot_df %>%
  filter(Targeted.x == TRUE) %>%
  filter(padj < 0.05) %>%
  mutate(plotting_cat = case_when(MultipleTargets == TRUE ~ "Multiple",
                                  .default = RPIP_documentation_category.x)) %>%
  mutate(plottng_cat = str_replace(plotting_cat, "lactam", "lactams")) %>%
  group_by(ARG.x, Fraction, Nanotrap_type, Enrichment, plotting_cat) %>%
  reframe(nreads = NumReads) %>%
  group_by(ARG.x, Enrichment) %>%
  reframe(readmeans = mean(nreads), sd = sd(nreads))

ARGs_w_over_25_mean_reads <- filter(mean_readcounts, readmeans > 25)
  
boxplot_df_by_arg <- boxplot_df %>%
  mutate(RPIP_documentation_category.x = str_replace(RPIP_documentation_category.x, "Beta-lactam \\+ beta-lactamase inhibitor", "Beta-lactam")) %>%
  mutate(plotting_cat = case_when(MultipleTargets == TRUE ~ "Multiple",
                                  .default = RPIP_documentation_category.x)) %>%
  mutate(plotting_cat = str_replace_all(plotting_cat, "lactam", "lactams")) %>%
  filter(Targeted.x == TRUE) %>%
  filter(ARG.x %in% ARGs_w_over_25_mean_reads$ARG.x) %>% 
  filter(!ARG.x %in% ARGs_detected_only_w_rpip_25plus_any_sample$ARG) %>%
  group_by(LIMS_ID, ARG_sig, Fraction, Nanotrap_type, Enrichment, plotting_cat) %>%
  summarize(nreads = sum(NumReads), ra_sum = sum(RA))

ggplot(filter(boxplot_df_by_arg, plotting_cat != "Multiple"), aes(y = ARG_sig, x = ra_sum, colour = Enrichment)) + geom_boxplot() +
  scale_x_log10() +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12)) +
  facet_grid(rows = "plotting_cat", scales = "free_y", space = "free_y") +
  xlab("Relative abundance of reads") +
  ylab("") +
  theme(strip.placement = "outside", strip.text.y = element_text(angle = 0, size = 10), 
        strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom") +
  theme(panel.grid.major.y = element_blank()) +
  geom_hline(yintercept=seq(1,nrow(boxplot_df_by_arg[,1]),1)+.5,color="black", linetype = "dotted", linewidth = .2) +
  geom_hline(yintercept=seq(1,nrow(boxplot_df_by_arg[,1]),1)-.5,color="black", linetype = "dotted", linewidth = .2)



#####
boxplot_df_by_category <- boxplot_df %>%
  mutate(RPIP_documentation_category.x = str_replace(RPIP_documentation_category.x, "Beta-lactam \\+ beta-lactamase inhibitor", "Beta-lactam")) %>%
  mutate(plotting_cat = case_when(MultipleTargets == TRUE ~ "Multiple",
                                  .default = RPIP_documentation_category.x)) %>%
  mutate(plotting_cat = str_replace_all(plotting_cat, "lactam", "lactams")) %>%
  filter(Targeted.x == TRUE) %>%
  filter(ARG.x %in% ARGs_w_over_25_mean_reads$ARG.x) %>% 
  filter(!ARG.x %in% ARGs_detected_only_w_rpip_25plus_any_sample$ARG) %>%
  group_by(LIMS_ID, plotting_cat, Fraction, Nanotrap_type, Enrichment) %>%
  summarize(nreads = sum(NumReads), ra_sum = sum(RA))

ggplot(boxplot_df_by_category, aes(y = plotting_cat, x = ra_sum, colour = Enrichment, fill = Enrichment)) + geom_boxplot() +
  scale_x_log10() +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 18)) +
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
  scale_color_manual(name = "Enrichment:", values = c("lightpink4", "aquamarine4")) +
  scale_fill_manual(name = "Enrichment:", values = c("lightpink2", "aquamarine2"))


#Look for genes that are ONLY in rpip:
only_rpip <- all_args_and_nonargs_per_sample %>%
  group_by(LIMS_ID, Fraction, Nanotrap_type, Enrichment) %>%
  mutate(RA = NumReads / sum(NumReads)) %>%
  filter(Enrichment == "RPIP") %>%
  left_join(mls_categories, by = join_by(ARG == Gene_name))

only_rpip[only_rpip$predicted_ARG_class == "MLS",]$predicted_ARG_class <- only_rpip[only_rpip$predicted_ARG_class == "MLS",]$Parent_drug_category


only_rpip <- only_rpip %>%
  left_join(rpip_deeparg_equivalents, by = join_by(predicted_ARG_class == Deeparg_CARD_category)) %>%
  mutate(RPIP_documentation_category = str_replace(RPIP_documentation_category, "Beta-lactam \\+ beta-lactamase inhibitor", "Beta-lactams"))

only_rpip$MultipleTargets <- str_detect(only_rpip$Parent_drug_category, ";")
only_rpip$MultipleTargets[is.na(only_rpip$MultipleTargets)] <- FALSE
  
  




only_untargeted <- all_args_and_nonargs_per_sample %>%
  filter(Enrichment == "None")

ARGs_detected_only_w_rpip <- only_rpip %>%
  filter(!(ARG %in% only_untargeted$ARG))

ARGs_detected_only_untargeted <- only_untargeted %>%
  filter(!(ARG %in% only_rpip$ARG))

ARGs_detected_only_w_rpip_25plus <- ARGs_detected_only_w_rpip %>%
  filter(NumReads >= 25)

ARGs_detected_only_w_rpip_25plus_any_sample <- ARGs_detected_only_w_rpip %>%
  filter(ARG %in% ARGs_detected_only_w_rpip_25plus$ARG) %>%
  mutate(plotting_cat = case_when(MultipleTargets == TRUE ~ "Multiple",
                                  .default = RPIP_documentation_category))
  



ggplot(ARGs_detected_only_w_rpip_25plus_any_sample, aes(x = RA, y = ARG)) + geom_boxplot(fill = "aquamarine2", colour = "aquamarine4") +
  scale_x_log10() +
  facet_grid(rows = "plotting_cat", scales = "free_y", space = "free_y") +
  xlab("Relative abundance of reads") +
  ylab("") +
  theme_bw() +
  theme(strip.placement = "outside", strip.text.y = element_text(angle = 270, size = 10), 
        strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(face = "bold")) +
  theme(axis.text.y = element_text(size = 12))

ggplot()
colnames(boxplot_df_by_arg)


all_ARGs_summed <- boxplot_df_by_arg %>%
  group_by(LIMS_ID, Fraction, Nanotrap_type, Enrichment) %>%
  summarize(total_targeted_reads = sum(ra_sum))

ggplot(all_ARGs_summed, aes(x = Fraction, y = total_targeted_reads, fill = Enrichment)) + geom_boxplot() +
  scale_y_log10() +
  facet_wrap(~Nanotrap_type)
