###Run import_deeparg_results.R before this script
library(tidyverse)
library(edgeR)

#Read in the imported/formatted table
rpip_v_unt_deeparg_results <- read_csv(
  "imported_deeparg_reports/imported_deeparg_results.csv"
  )

###EdgeR analysis
edger_arg_count_table <- read_csv(
  "imported_deeparg_reports/deeparg_count_matrix.csv", guess_max = Inf
  )

#Load metadata  
edger_arg_metadata <- read_csv(
  "imported_deeparg_reports/deeparg_metadata.csv",
  guess_max = Inf
  ) %>%
  arrange(UniqueID)

#edger_arg_metadata$LIMS_ID <- factor(edger_arg_metadata$LIMS_ID)


prepare_arg_count_table_for_edgeR <- function(count_table, group) {
  #This will ensure columns and their order match in counts/metadata:
  acounts <- count_table[, c("ARG_w_class", group$UniqueID)]
  
  #Sanity check (make sure sample names and order match in counts and group 
  #data):
  mismatches <- sum(colnames(acounts[,2:ncol(acounts)]) != group$UniqueID)
  if (mismatches != 0) {
    stop("Something is wrong; the samples in the count table don't match the \
         samples in the group data.")
  }
  return(acounts)
}

edger_arg_count_table <- prepare_arg_count_table_for_edgeR(
  edger_arg_count_table,
  edger_arg_metadata
  )

###Generate a design matrix for EdgeR:
arg_generate_levels <- function(group_df) {
  #Use LIMS_ID to control for which sample we started with:
  LIMS_ID <- factor(group_df$LIMS_ID)
  
  #Define explanatory variables as factors and order them to have an appropriate
  #baseline:  
  Fraction <- factor(group_df$Fraction, levels = c("unfiltered", 
                                                   "retentate", 
                                                   "filtrate"))
  
  Nanotrap_type <- factor(group_df$Nanotrap_type, levels = c("none", 
                                                             "A", 
                                                             "A&B"))
  
  Enrichment <- factor(group_df$Enrichment, levels = c("None", 
                                                       "RPIP"))
  
  #Create a tibble containing all combinations of treatment variables:
  treat_tb <- expand.grid(Fraction = levels(Fraction),
                          Nanotrap_type = levels(Nanotrap_type),
                          Enrichment = levels(Enrichment),
                          stringsAsFactors = FALSE)
  
  #Create an empty list, then add every possible combination of treatment
  #variables:
  treat_list <- vector("list", length = nrow(treat_tb) - 1)
  
  for(row in 2:nrow(treat_tb)) {
    treat_vec <- (Fraction == treat_tb[row, 1] & 
                    Nanotrap_type == treat_tb[row, 2] &
                    Enrichment == treat_tb[row, 3])
    treat_list[[row - 1]] <- treat_vec
    names(treat_list) <- paste(treat_tb[, 1], 
                               treat_tb[, 2], 
                               treat_tb[, 3], 
                               sep = ".")[-1]
  }
  
  treat_mat <- lapply(treat_list, as.numeric) %>% as_tibble() %>% as.matrix()
  
  #Start model with only LIMS_ID as explanatory variable:
  design <- model.matrix(~LIMS_ID)
  
  #Append all of the boolean treatment combinations to the model matrix:
  design <- cbind(design, treat_mat)
  
  #And there we have it, the model matrix! 
  return(design)
}


#Now run the design generator to get a model matrix:
design_arg <- arg_generate_levels(edger_arg_metadata)


###Fit the model
#Create a DGEList object (what EdgeR works with) from the count matrix:
edger_arg_dge <- DGEList(
  counts = edger_arg_count_table
)

#Find low-frequency args that don't give us enough information to be useful but
#mess with the analysis:
edger_arg_dge_lf_remover <- filterByExpr(
  edger_arg_dge, 
  design = design_arg)

edger_arg_dge_lfRemoved <- edger_arg_dge[
  edger_arg_dge_lf_remover, , keep.lib.sizes = FALSE
]

#Fit the model:
edger_arg_disp_lfRemoved <- estimateDisp(
  y = edger_arg_dge_lfRemoved,
  design = design_arg
)

edger_arg_fit_lfRemoved <- glmQLFit(
  edger_arg_disp_lfRemoved, 
  design_arg
)


#Create Boolean vectors indicating which columns of the design matrix correspond
#to individual treatments:
index_main_effects <- function(design) {
  is.RPIP <<- str_detect(colnames(design), "RPIP")
  is.VSP <<- str_detect(colnames(design), "VSP")
  is.Untargeted <<- str_detect(colnames(design), "None")
  is.Filtrate <<- str_detect(colnames(design), "filtrate")
  is.Retentate <<- str_detect(colnames(design), "retentate")
  is.Unfiltered <<- str_detect(colnames(design), "unfiltered")
  is.NanoA <<- str_detect(colnames(design), "A\\.")
  is.NanoB <<- str_detect(colnames(design), "A&B")
  is.DirectExt <<- str_detect(colnames(design), "none")
}

index_main_effects(design_arg)

top25 <- function(result) {
  topTags(result, n = 25, sort.by = "PValue")
}

###RPIP - Untargeted:
#Main effect:
RPIP_all_treatments_args <- glmQLFTest(edger_arg_fit_lfRemoved, 
                                       contrast = is.RPIP - is.Untargeted)
top25(RPIP_all_treatments_args)
#No NT, no filtration:
RPIP_vs_Unt_DirExt_Unfiltered_arg <- glmQLFTest(edger_arg_fit_lfRemoved, 
 contrast = (is.RPIP & is.Unfiltered & is.DirectExt) - 
   (is.Untargeted & is.Unfiltered & is.DirectExt))
top25(RPIP_vs_Unt_DirExt_Unfiltered_arg)

#NTA, no filtration:
RPIP_vs_Unt_NanoA_Unfiltered_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.RPIP & is.Unfiltered & is.NanoA) - (is.Untargeted & is.Unfiltered & is.NanoA))
top25(RPIP_vs_Unt_NanoA_Unfiltered_90conf_no_rrna)

#NTB, no filtration:
RPIP_vs_Unt_NanoBA_Unfiltered_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.RPIP & is.NanoB & is.Unfiltered) - (is.Untargeted & is.NanoB & is.Unfiltered))
top25(RPIP_vs_Unt_NanoBA_Unfiltered_90conf_no_rrna)

#No NT, filtrate:
RPIP_vs_Unt_DirExt_Filtrate_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.RPIP & is.Filtrate & is.DirectExt) - (is.Untargeted & is.Filtrate & is.DirectExt))
top25(RPIP_vs_Unt_DirExt_Filtrate_90conf_no_rrna)

#NTA, filtrate:
RPIP_vs_Unt_NanoA_Filtrate_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.RPIP & is.NanoA & is.Filtrate) - (is.Untargeted & is.NanoA & is.Filtrate))
top25(RPIP_vs_Unt_NanoA_Filtrate_90conf_no_rrna)

#NTB, filtrate:
RPIP_vs_Unt_NanoBA_Filtrate_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.RPIP & is.NanoB & is.Filtrate) - (is.Untargeted & is.NanoB & is.Filtrate))
top25(RPIP_vs_Unt_NanoBA_Filtrate_90conf_no_rrna)

#No NT, retentate:
RPIP_vs_Unt_DirExt_Retentate_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.RPIP & is.Retentate & is.DirectExt) - (is.Untargeted & is.Retentate & is.DirectExt))
top25(RPIP_vs_Unt_DirExt_Retentate_90conf_no_rrna)

#NTA, retentate:
RPIP_vs_Unt_NanoA_Retentate_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.RPIP & is.NanoA & is.Retentate) - (is.Untargeted & is.NanoA & is.Retentate))
top25(RPIP_vs_Unt_NanoA_Retentate_90conf_no_rrna)

#NTB, retentate:
RPIP_vs_Unt_NanoBA_Retentate_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.RPIP & is.NanoB & is.Retentate) - (is.Untargeted & is.NanoB & is.Retentate))
top25(RPIP_vs_Unt_NanoBA_Retentate_90conf_no_rrna)

########################################3
##########################################
########################################
######################################3
##########################################
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
