library(tidyverse)
library(DESeq2)
library(pheatmap)
setwd("/projects/bios_microbe/cowen20/rprojects/targeted_panels/")
family_count_table <- read_csv("DS2_family_count_matrix_no_rrna.csv", guess_max = Inf)

#Find rows not identified to family level
no_fam_indices <- str_which(family_count_table$`F`, ("^unidentified|^unclassified"))
#Name them "Unidentified"
family_count_table[no_fam_indices, ]$`F` <- "Unidentified"

#Merge unidentified taxa into 1 row
family_count_table <- family_count_table %>%
  group_by(`F`) %>%
  summarise(across(everything(), sum))

#Convert to matrix
family_count_matrix <- family_count_table %>%
  select(!F) %>%
  as.matrix()
row.names(family_count_matrix) <- family_count_table$F

#Load metadata  
sample_metadata <- read.csv("DS2_sample_metadata.csv", row.names = 1, stringsAsFactors = TRUE)
sample_metadata$LIMS_ID <- factor(sample_metadata$LIMS_ID)

#Filter metadata for samples where Kraken2 is confident
sample_metadata_90conf <- filter(sample_metadata, Kraken2_confidence == 0.9)

#Filter family counts to samples where Kraken2 is confident
family_count_matrix_90conf <- family_count_matrix[, rownames(sample_metadata_90conf)]
family_count_matrix_90conf <- family_count_matrix_90conf[rowSums(family_count_matrix_90conf) > 0, ]

#Check sample names in same order
sum(colnames(family_count_matrix) != rownames(sample_metadata))
sum(colnames(family_count_matrix_90conf) != rownames(sample_metadata_90conf))


#Perform DESeq2 analysis
family_dds <- DESeqDataSetFromMatrix(countData = family_count_matrix,
                                     colData = sample_metadata,
                                     design = ~ 0 + LIMS_ID + Fraction + Nanotrap_type + Enrichment)

family_dds$Enrichment <- relevel(family_dds$Enrichment, ref = "None")
family_dds$Fraction <- relevel(family_dds$Fraction, ref = "unfiltered")
family_dds$Nanotrap_type <- relevel(family_dds$Nanotrap_type, ref = "none")

family_dds <- DESeq(family_dds)

family_result_table <- data.frame(results(family_dds))


#Perform DESeq2 analysis on 90conf set:
family_dds_90conf <- DESeqDataSetFromMatrix(countData = family_count_matrix_90conf,
                                            colData = sample_metadata_90conf,
                                            design = ~ 0 + 
                                              LIMS_ID + 
                                              Enrichment + 
                                              Fraction + 
                                              Nanotrap_type
                                            )


#Set reference levels
family_dds_90conf$Enrichment <- relevel(family_dds_90conf$Enrichment, ref = "None")
family_dds_90conf$Fraction <- relevel(family_dds_90conf$Fraction, ref = "unfiltered")
family_dds_90conf$Nanotrap_type <- relevel(family_dds_90conf$Nanotrap_type, ref = "none")

family_dds_90conf <- DESeq(family_dds_90conf)

View(data.frame(mcols(family_dds_90conf)))

#Try with low-abundance removed
abund_family_count_matrix_90conf <- family_count_matrix_90conf[rowSums(family_count_matrix_90conf) >= 10,]

abund_family_dds_90conf <- DESeqDataSetFromMatrix(countData = abund_family_count_matrix_90conf,
                                           colData = sample_metadata_90conf,
                                           design = ~ 0 + 
                                             LIMS_ID + 
                                             Enrichment + 
                                             Fraction + 
                                             Nanotrap_type
                                           )

abund_family_dds_90conf <- DESeq(abund_family_dds_90conf)



#heatmap

dsheatmap <- function(ds2_obj, an_source, an_name) {
  vsd_dds <- varianceStabilizingTransformation(ds2_obj, blind = TRUE)
  vsd_mat_dds <- assay(vsd_dds)
  vsd_cor_dds <- cor(vsd_mat_dds)
  pheatmap(vsd_cor_dds, annotation = select(an_source, an_name))
}

dsheatmap(family_dds_90conf, sample_metadata_90conf, "LIMS_ID")
vsd_family_dds_90 <- varianceStabilizingTransformation(family_dds_90conf, blind = TRUE)
vsd_mat_family_dds <- assay(vsd_family_dds_90)
vsd_cor_family_dds <- cor(vsd_mat_family_dds)
pheatmap(vsd_cor_family_dds, annotation = select(sample_metadata, LIMS_ID))

#PCA
plotPCA(vsd_family_dds_90, intgroup = "LIMS_ID", pcsToUse = c(1,2))
plotPCA(vsd_family_dds_90, intgroup = "Enrichment")
#plotPCA(vsd_family_dds_90, intgroup = "site")
plotPCA(vsd_family_dds_90, intgroup = "Nanotrap_type")
plotPCA(vsd_family_dds_90, intgroup = "Fraction")

#Dispersion
plotDispEsts(family_dds_90conf)


#Contrast VSP
vspvsnone_family_res_90conf <- results(family_dds_90conf,
                        contrast = c("Enrichment", 
                                     "VSP",
                                     "None"),
                        alpha = 0.05)
vspvsnone_family_res_90conf <- lfcShrink(family_dds_90conf, contrast = c("Enrichment", "VSP", "None"),
                        res = vspvsnone_family_res_90conf, type = "ashr")

resultsNames(family_dds_90conf)

summary(vspvsnone_family_res_90conf)
View(data.frame(vspvsnone_family_res_90conf))

plotMA(vspvsnone_family_res_90conf, ylim = c(-8,8))
mcols(vspvsnone_family_res_90conf)


#Contrast RPIP
rpipvsnone_family_res_90conf <- results(family_dds_90conf,
                             contrast = c("Enrichment", 
                                          "RPIP",
                                          "None"),
                             alpha = 0.05)
rpipvsnone_family_res_90conf <- lfcShrink(family_dds_90conf, contrast = c("Enrichment", "RPIP", "None"),
                               res = rpipvsnone_family_res_90conf, type = "ashr")

summary(rpipvsnone_family_res_90conf)
View(data.frame(rpipvsnone_family_res_90conf))


###Look at targeted taxa (VSP)
vsp_targets <- read_csv("input/vsp_panels/vsp_target_info.csv")
vsp_family_res_90conf <- data.frame(vspvsnone_family_res_90conf)
vsp_family_res_90conf$Family_name <- rownames(vsp_family_res_90conf)
vsp_family_res_90conf$Targeted <- vsp_family_res_90conf$Family_name %in% vsp_targets$family_name
vsp_family_res_90conf$padj[is.na(vsp_family_res_90conf$padj)] <- 1

#Volcano colored by targeted
ggplot(vsp_family_res_90conf, 
       aes(x = log2FoldChange, y = -log10(padj), col = Targeted)) + 
  geom_point(aes(color = Targeted)) +
  geom_text(aes(
    label = ifelse((-log10(padj) > 1 & abs(log2FoldChange) > 1), Family_name, ""),
    vjust = -0.5, angle = -30), 
    size = 3.5, show.legend = FALSE) +
  theme_minimal() +
  xlim(-5, 15) +
  ylim(0, 150) +
  ggtitle("VSP vs. No Enrichment") +
  xlab(expression(paste("log"[2], " fold change"))) +
  ylab(expression(paste("-log"[10], " adjusted p-value"))) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  geom_hline(yintercept = 1.30103) +
  annotate("text", x = 10, y = 6, label = "adjusted p-value = 0.05")


###Look at targeted taxa (RPIP)
rpip_targets <- read_csv("input/rpip_panels/rpip_target_info.csv")
rpip_family_res_90conf <- data.frame(rpipvsnone_family_res_90conf)
rpip_family_res_90conf$Family_name <- rownames(rpip_family_res_90conf)
rpip_family_res_90conf$Targeted <- rpip_family_res_90conf$Family_name %in% rpip_targets$family
rpip_family_res_90conf$padj[is.na(rpip_family_res_90conf$padj)] <- 1

rpip_family_res_90conf %>%
  filter(Targeted == TRUE) %>%
  View()

ggplot(rpip_family_res_90conf, 
       aes(x = log2FoldChange, y = -log10(padj), col = Targeted)) + 
  geom_point() +
  geom_text(aes(label = ifelse((-log10(padj) > 5 & abs(log2FoldChange) > 1), Family_name, ""), vjust = -0.5, angle = -10), size = 3.5, show.legend = FALSE) +
  theme_minimal() +
  xlim(-25, 10) +
  ggtitle("RPIP vs. No Enrichment") +
  xlab(expression(paste("log"[2], " fold change"))) +
  ylab(expression(paste("-log"[10], " adjusted p-value"))) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  geom_hline(yintercept = 1.30103) +
  annotate("text", x = -22, y = 2, label = "adjusted p-value = 0.05")





#Contrast Nanotrap A
nanoAvsnone_family_res_90conf <- results(family_dds_90conf,
                                        contrast = c("Nanotrap_type", 
                                                     "A",
                                                     "none"),
                                        alpha = 0.05)
nanoAvsnone_family_res_90conf <- lfcShrink(family_dds_90conf, contrast = c("Nanotrap_type", "A", "none"),
                                          res = nanoAvsnone_family_res_90conf, type = "ashr")

summary(nanoAvsnone_family_res_90conf)
View(data.frame(nanoAvsnone_family_res_90conf))


rpip_targets <- read_csv("input/rpip_panels/rpip_target_info.csv")
nanoA_by_targeted_family_res_90conf <- data.frame(nanoAvsnone_family_res_90conf)
nanoA_by_targeted_family_res_90conf$Family_name <- rownames(nanoA_by_targeted_family_res_90conf)
nanoA_by_targeted_family_res_90conf$RPIP_Targeted <- nanoA_by_targeted_family_res_90conf$Family_name %in% rpip_targets$family
nanoA_by_targeted_family_res_90conf$VSP_Targeted <- nanoA_by_targeted_family_res_90conf$Family_name %in% vsp_targets$family_name
nanoA_by_targeted_family_res_90conf$padj[is.na(nanoA_by_targeted_family_res_90conf$padj)] <- 1
nanoA_by_targeted_family_res_90conf <- nanoA_by_targeted_family_res_90conf %>%
  mutate(Targeted_by = case_when(RPIP_Targeted & VSP_Targeted ~ "Both",
                                 RPIP_Targeted & !VSP_Targeted ~ "RPIP",
                                 !RPIP_Targeted & VSP_Targeted ~ "VSP",
                                 !(RPIP_Targeted | VSP_Targeted) ~ "Neither"))

ggplot(nanoA_by_targeted_family_res_90conf, 
       aes(x = log2FoldChange, y = -log10(padj), col = Targeted_by)) + 
  geom_point(alpha = 0.8, size = 0.5) +
  #geom_text(aes(label = ifelse((-log10(padj) > 1.30103), Family_name, ""), vjust = -0.5, angle = -30), size = 1, show.legend = FALSE) +
  theme_minimal() +
  xlim(-5, 3) +
  #ylim(0, 6) +
  ggtitle("Nanotrap A vs. No Concentration") +
  xlab(expression(paste("log"[2], " fold change"))) +
  ylab(expression(paste("-log"[10], " adjusted p-value"))) +
  scale_color_manual(values = c("magenta", "grey", "forestgreen", "#D55E00")) +
  geom_hline(yintercept = 1.30103) +
  annotate("text", x = -4, y = 2, label = "adjusted p-value = 0.05")



#Contrast Nanotrap A
nanoAvsnone_family_res_90conf <- results(family_dds_90conf,
                                         contrast = c("Nanotrap_type", 
                                                      "A",
                                                      "none"),
                                         alpha = 0.05)
nanoAvsnone_family_res_90conf <- lfcShrink(family_dds_90conf, contrast = c("Nanotrap_type", "A", "none"),
                                           res = nanoAvsnone_family_res_90conf, type = "ashr")

summary(nanoAvsnone_family_res_90conf)
View(data.frame(nanoAvsnone_family_res_90conf))


rpip_targets <- read_csv("input/rpip_panels/rpip_target_info.csv")
nanoA_by_targeted_family_res_90conf <- data.frame(nanoAvsnone_family_res_90conf)
nanoA_by_targeted_family_res_90conf$Family_name <- rownames(nanoA_by_targeted_family_res_90conf)
nanoA_by_targeted_family_res_90conf$RPIP_Targeted <- nanoA_by_targeted_family_res_90conf$Family_name %in% rpip_targets$family
nanoA_by_targeted_family_res_90conf$VSP_Targeted <- nanoA_by_targeted_family_res_90conf$Family_name %in% vsp_targets$family_name
nanoA_by_targeted_family_res_90conf$padj[is.na(nanoA_by_targeted_family_res_90conf$padj)] <- 1
nanoA_by_targeted_family_res_90conf <- nanoA_by_targeted_family_res_90conf %>%
  mutate(Targeted_by = case_when(RPIP_Targeted & VSP_Targeted ~ "Both",
                                 RPIP_Targeted & !VSP_Targeted ~ "RPIP",
                                 !RPIP_Targeted & VSP_Targeted ~ "VSP",
                                 !(RPIP_Targeted | VSP_Targeted) ~ "Neither"))

ggplot(nanoA_by_targeted_family_res_90conf, 
       aes(x = log2FoldChange, y = -log10(padj), col = Targeted_by)) + 
  geom_point(alpha = 0.8, size = 1.5) +
#  geom_text(aes(label = ifelse((-log10(padj) > 1.30103), Family_name, ""), vjust = -0.5, angle = -30), size = 3, show.legend = FALSE) +
  theme_minimal() +
  xlim(-5, 3) +
  #ylim(0, 6) +
  ggtitle("Nanotrap A vs. No Concentration") +
  xlab(expression(paste("log"[2], " fold change"))) +
  ylab(expression(paste("-log"[10], " adjusted p-value"))) +
  scale_color_manual(values = c("magenta", "grey", "forestgreen", "#D55E00")) +
  geom_hline(yintercept = 1.30103) +
  annotate("text", x = -4, y = 2, label = "adjusted p-value = 0.05")


make_df_for_plotting_from_dds <- function(dds_object, contrast_col, 
                                          contrast_lvl1, contrast_lvl2) {
  res_df <- results(dds_object, contrast = c(contrast_col,
                                             contrast_lvl1,
                                             contrast_lvl2),
                    alpha = 0.05)
  res_df <- lfcShrink(dds_object, contrast = c(contrast_col,
                                               contrast_lvl1,
                                               contrast_lvl2),
                      res = res_df, type = "ashr")
  res_df_by_targeted <- data.frame(res_df)
  res_df_by_targeted$Family_name <- rownames(res_df_by_targeted)
  res_df_by_targeted$RPIP_Targeted <- res_df_by_targeted$Family_name %in% rpip_targets$family
  res_df_by_targeted$VSP_Targeted <- res_df_by_targeted$Family_name %in% vsp_targets$family_name
  res_df_by_targeted$padj[is.na(res_df_by_targeted$padj)] <- 1
  res_df_by_targeted <- res_df_by_targeted %>%
    mutate(Targeted_by = case_when(RPIP_Targeted & VSP_Targeted ~ "Both",
                                   RPIP_Targeted & !VSP_Targeted ~ "RPIP",
                                   !RPIP_Targeted & VSP_Targeted ~ "VSP",
                                   !(RPIP_Targeted | VSP_Targeted) ~ "Neither"))
  return(res_df_by_targeted)
}


nanoAvsBA_by_targeted_family_res_90conf <- make_df_for_plotting_from_dds(family_dds_90conf, "Nanotrap_type", "A", "A&B")

ggplot(nanoAvsBA_by_targeted_family_res_90conf, 
       aes(x = log2FoldChange, y = -log10(padj), col = Targeted_by)) + 
  geom_point(alpha = 0.8, size = 0.5) +
  #geom_text(aes(label = ifelse((-log10(padj) > 1.30103), Family_name, ""), vjust = -0.5, angle = -30), size = 1, show.legend = FALSE) +
  theme_minimal() +
  #xlim(-5, 3) +
  #ylim(0, 6) +
  ggtitle("Nanotrap A vs. No Concentration") +
  xlab(expression(paste("log"[2], " fold change"))) +
  ylab(expression(paste("-log"[10], " adjusted p-value"))) +
  scale_color_manual(values = c("magenta", "grey", "forestgreen", "#D55E00")) +
  geom_hline(yintercept = 1.30103) +
  annotate("text", x = -1, y = 1.4, label = "adjusted p-value = 0.05")


ggplot(nanoBA_by_targeted_family_res_90conf, 
       aes(x = log2FoldChange, y = -log10(padj), col = Targeted_by)) + 
  geom_point(alpha = 0.8, size = 1.5) +
  #geom_text(aes(label = ifelse((-log10(padj) > 1.30103), Family_name, ""), vjust = -0.5, angle = -30), size = 1, show.legend = FALSE) +
  theme_minimal() +
  #xlim(-5, 3) +
  #ylim(0, 6) +
  ggtitle("Nanotrap B vs. No Concentration") +
  xlab(expression(paste("log"[2], " fold change"))) +
  ylab(expression(paste("-log"[10], " adjusted p-value"))) +
  scale_color_manual(values = c("magenta", "grey", "forestgreen", "#D55E00")) +
  geom_hline(yintercept = 1.30103) +
  annotate("text", x = -1.5, y = 1.5, label = "adjusted p-value = 0.05")



#Contrast Filtrate##############################################################
FractionFiltratevsUnfiltered_family_res_90conf <- results(family_dds_90conf,
                                         contrast = c("Fraction", 
                                                      "filtrate",
                                                      "unfiltered"),
                                         alpha = 0.05)
FractionFiltratevsUnfiltered_family_res_90conf <- lfcShrink(family_dds_90conf, contrast = c("Fraction", "filtrate", "unfiltered"),
                                           res = FractionFiltratevsUnfiltered_family_res_90conf, type = "ashr")

summary(FractionFiltratevsUnfiltered_family_res_90conf)
View(data.frame(FractionFiltratevsUnfiltered_family_res_90conf))


rpip_targets <- read_csv("input/rpip_panels/rpip_target_info.csv")
FractionFiltered_by_targeted_family_res_90conf <- data.frame(FractionFiltratevsUnfiltered_family_res_90conf)
FractionFiltered_by_targeted_family_res_90conf$Family_name <- rownames(FractionFiltered_by_targeted_family_res_90conf)
FractionFiltered_by_targeted_family_res_90conf$RPIP_Targeted <- FractionFiltered_by_targeted_family_res_90conf$Family_name %in% rpip_targets$family
FractionFiltered_by_targeted_family_res_90conf$VSP_Targeted <- FractionFiltered_by_targeted_family_res_90conf$Family_name %in% vsp_targets$family_name
FractionFiltered_by_targeted_family_res_90conf$padj[is.na(FractionFiltered_by_targeted_family_res_90conf$padj)] <- 1
FractionFiltered_by_targeted_family_res_90conf <- FractionFiltered_by_targeted_family_res_90conf %>%
  mutate(Targeted_by = case_when(RPIP_Targeted & VSP_Targeted ~ "Both",
                                 RPIP_Targeted & !VSP_Targeted ~ "RPIP",
                                 !RPIP_Targeted & VSP_Targeted ~ "VSP",
                                 !(RPIP_Targeted | VSP_Targeted) ~ "Neither"))

ggplot(FractionFiltered_by_targeted_family_res_90conf, 
       aes(x = log2FoldChange, y = -log10(padj), col = Targeted_by)) + 
  geom_point(alpha = 0.8, size = 1.5) +
  #  geom_text(aes(label = ifelse((-log10(padj) > 1.30103), Family_name, ""), vjust = -0.5, angle = -30), size = 3, show.legend = FALSE) +
  theme_minimal() +
#  xlim(-5, 3) +
  #ylim(0, 6) +
  ggtitle("Filtrate vs. Unfiltered") +
  xlab(expression(paste("log"[2], " fold change"))) +
  ylab(expression(paste("-log"[10], " adjusted p-value"))) +
  scale_color_manual(values = c("magenta", "grey", "forestgreen", "#D55E00")) +
  geom_hline(yintercept = 1.30103) +
  annotate("text", x = -4, y = 2.5, label = "adjusted p-value = 0.05")


make_df_for_plotting_from_dds <- function(dds_object, contrast_col, 
                                          contrast_lvl1, contrast_lvl2) {
  res_df <- results(dds_object, contrast = c(contrast_col,
                                             contrast_lvl1,
                                             contrast_lvl2),
                    alpha = 0.05)
  res_df <- lfcShrink(dds_object, contrast = c(contrast_col,
                                               contrast_lvl1,
                                               contrast_lvl2),
                      res = res_df, type = "ashr")
  res_df_by_targeted <- data.frame(res_df)
  res_df_by_targeted$Family_name <- rownames(res_df_by_targeted)
  res_df_by_targeted$RPIP_Targeted <- res_df_by_targeted$Family_name %in% rpip_targets$family
  res_df_by_targeted$VSP_Targeted <- res_df_by_targeted$Family_name %in% vsp_targets$family_name
  res_df_by_targeted$padj[is.na(res_df_by_targeted$padj)] <- 1
  res_df_by_targeted <- res_df_by_targeted %>%
    mutate(Targeted_by = case_when(RPIP_Targeted & VSP_Targeted ~ "Both",
                                   RPIP_Targeted & !VSP_Targeted ~ "RPIP",
                                   !RPIP_Targeted & VSP_Targeted ~ "VSP",
                                   !(RPIP_Targeted | VSP_Targeted) ~ "Neither"))
  return(res_df_by_targeted)
}


FractionRetentateVsUnfiltered <- make_df_for_plotting_from_dds(family_dds_90conf, "Fraction", "retentate", "unfiltered")
ggplot(FractionRetentateVsUnfiltered, 
       aes(x = log2FoldChange, y = -log10(padj), col = Targeted_by)) + 
  geom_point(alpha = 0.8, size = 1.5) +
  #  geom_text(aes(label = ifelse((-log10(padj) > 1.30103), Family_name, ""), vjust = -0.5, angle = -30), size = 3, show.legend = FALSE) +
  theme_minimal() +
  xlim(-10, 3) +
  #ylim(0, 6) +
  ggtitle("Retentate vs. Unfiltered") +
  xlab(expression(paste("log"[2], " fold change"))) +
  ylab(expression(paste("-log"[10], " adjusted p-value"))) +
  scale_color_manual(values = c("magenta", "grey", "forestgreen", "#D55E00")) +
  geom_hline(yintercept = 1.30103) +
  annotate("text", x = -4, y = 1.5, label = "adjusted p-value = 0.05")
