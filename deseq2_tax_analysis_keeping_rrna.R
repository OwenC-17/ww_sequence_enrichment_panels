library(tidyverse)
library(DESeq2)
library(pheatmap)
setwd("/projects/bios_microbe/cowen20/rprojects/targeted_panels/")
family_count_table_w_rrna <- read_csv("DS2_family_count_matrix_w_rrna.csv", guess_max = Inf)

#Find rows not identified to family level
no_fam_indices_w_rrna <- str_which(family_count_table_w_rrna$`F`, ("^unidentified|^unclassified"))
#Name them "Unidentified"
family_count_table_w_rrna[no_fam_indices_w_rrna, ]$`F` <- "Unidentified"

#Merge unidentified taxa into 1 row
family_count_table_w_rrna <- family_count_table_w_rrna %>%
  group_by(`F`) %>%
  summarise(across(everything(), sum))

#Convert to matrix
family_count_matrix_w_rrna <- family_count_table_w_rrna %>%
  select(!F) %>%
  as.matrix()
row.names(family_count_matrix_w_rrna) <- family_count_table_w_rrna$F

#Load metadata  
sample_metadata_w_rrna <- read.csv("DS2_sample_metadata_w_rrna.csv", row.names = 1, stringsAsFactors = TRUE)
sample_metadata_w_rrna$LIMS_ID <- factor(sample_metadata_w_rrna$LIMS_ID)

#Filter metadata for samples where Kraken2 is confident
sample_metadata_90conf_w_rrna <- filter(sample_metadata_w_rrna, Kraken2_confidence == 0.9)

#Filter family counts to samples where Kraken2 is confident
family_count_matrix_90conf_w_rrna <- family_count_matrix_w_rrna[, rownames(sample_metadata_90conf_w_rrna)]
family_count_matrix_90conf_w_rrna <- family_count_matrix_90conf_w_rrna[rowSums(family_count_matrix_90conf_w_rrna) > 0, ]

#Check sample names in same order
sum(colnames(family_count_matrix_w_rrna) != rownames(sample_metadata_w_rrna))
sum(colnames(family_count_matrix_90conf_w_rrna) != rownames(sample_metadata_90conf_w_rrna))


#Perform DESeq2 analysis
family_dds_incl_rrna <- DESeqDataSetFromMatrix(countData = family_count_matrix_w_rrna,
                                     colData = sample_metadata_w_rrna,
                                     design = ~ 0 + LIMS_ID + Fraction + Nanotrap_type + Enrichment)

family_dds_incl_rrna$Enrichment <- relevel(family_dds_incl_rrna$Enrichment, ref = "None")
family_dds_incl_rrna$Fraction <- relevel(family_dds_incl_rrna$Fraction, ref = "unfiltered")
family_dds_incl_rrna$Nanotrap_type <- relevel(family_dds_incl_rrna$Nanotrap_type, ref = "none")

family_dds_incl_rrna <- DESeq(family_dds_incl_rrna)

family_result_table_incl_rrna <- data.frame(results(family_dds_incl_rrna))


#Perform DESeq2 analysis on 90conf set:
family_dds_90conf_incl_rrna <- DESeqDataSetFromMatrix(countData = family_count_matrix_90conf_w_rrna,
                                            colData = sample_metadata_90conf_w_rrna,
                                            design = ~ 0 + 
                                              LIMS_ID + 
                                              Fraction + 
                                              Nanotrap_type +
                                              Enrichment
)


#Set reference levels
family_dds_90conf_incl_rrna$Enrichment <- relevel(family_dds_90conf_incl_rrna$Enrichment, ref = "None")
family_dds_90conf_incl_rrna$Fraction <- relevel(family_dds_90conf_incl_rrna$Fraction, ref = "unfiltered")
family_dds_90conf_incl_rrna$Nanotrap_type <- relevel(family_dds_90conf_incl_rrna$Nanotrap_type, ref = "none")

family_dds_90conf_incl_rrna <- DESeq(family_dds_90conf_incl_rrna)

View(data.frame(mcols(family_dds_90conf_incl_rrna)))

#Try with low-abundance removed
#abund_family_count_matrix_90conf <- family_count_matrix_90conf[rowSums(family_count_matrix_90conf) >= 10,]

#abund_family_dds_90conf <- DESeqDataSetFromMatrix(countData = abund_family_count_matrix_90conf,
#                                                  colData = sample_metadata_90conf,
#                                                  design = ~ 0 + 
#                                                    LIMS_ID + 
#                                                    Enrichment + 
#                                                    Fraction + 
#                                                    Nanotrap_type
#)

#abund_family_dds_90conf <- DESeq(abund_family_dds_90conf)



#heatmap

dsheatmap <- function(ds2_obj, an_source, an_name) {
  vsd_dds <- varianceStabilizingTransformation(ds2_obj, blind = TRUE)
  vsd_mat_dds <- assay(vsd_dds)
  vsd_cor_dds <- cor(vsd_mat_dds)
  pheatmap(vsd_cor_dds, annotation = select(an_source, an_name))
}

dsheatmap(family_dds_90conf_incl_rrna, sample_metadata_90conf_w_rrna, "LIMS_ID")
vsd_family_dds_90_incl_rrna <- varianceStabilizingTransformation(family_dds_90conf_incl_rrna, blind = TRUE)
vsd_mat_family_dds_incl_rrna <- assay(vsd_family_dds_90_incl_rrna)
vsd_cor_family_dds_incl_rrna <- cor(vsd_mat_family_dds_incl_rrna)
pheatmap(vsd_cor_family_dds_incl_rrna, annotation = select(sample_metadata_90conf_w_rrna, LIMS_ID))

#PCA
plotPCA(vsd_family_dds_90_incl_rrna, intgroup = "LIMS_ID", pcsToUse = c(1,2))
plotPCA(vsd_family_dds_90_incl_rrna, intgroup = "Enrichment")
#plotPCA(vsd_family_dds_90, intgroup = "site")
plotPCA(vsd_family_dds_90_incl_rrna, intgroup = "Nanotrap_type")
plotPCA(vsd_family_dds_90_incl_rrna, intgroup = "Fraction")

#Dispersion
plotDispEsts(family_dds_90conf_incl_rrna)


#Contrast VSP
vspvsnone_family_res_90conf_incl_rrna <- results(family_dds_90conf_incl_rrna,
                                       contrast = c("Enrichment", 
                                                    "VSP",
                                                    "None"),
                                       alpha = 0.05)
vspvsnone_family_res_90conf_incl_rrna <- lfcShrink(family_dds_90conf_incl_rrna, contrast = c("Enrichment", "VSP", "None"),
                                         res = vspvsnone_family_res_90conf_incl_rrna, type = "ashr")

resultsNames(family_dds_90conf_incl_rrna)

summary(vspvsnone_family_res_90conf_incl_rrna)
View(data.frame(vspvsnone_family_res_90conf_incl_rrna))

plotMA(vspvsnone_family_res_90conf_incl_rrna, ylim = c(-8,8))
mcols(vspvsnone_family_res_90conf_incl_rrna)


#Contrast RPIP
rpipvsnone_family_res_90conf_incl_rrna <- results(family_dds_90conf_incl_rrna,
                                        contrast = c("Enrichment", 
                                                     "RPIP",
                                                     "None"),
                                        alpha = 0.05)
rpipvsnone_family_res_90conf_incl_rrna <- lfcShrink(family_dds_90conf_incl_rrna, contrast = c("Enrichment", "RPIP", "None"),
                                          res = rpipvsnone_family_res_90conf_incl_rrna, type = "ashr")

summary(rpipvsnone_family_res_90conf_incl_rrna)
View(data.frame(rpipvsnone_family_res_90conf_incl_rrna))


###Look at targeted taxa (VSP)
vsp_targets <- read_csv("input/vsp_panels/vsp_target_info.csv")
vsp_family_res_90conf_incl_rrna <- data.frame(vspvsnone_family_res_90conf_incl_rrna)
vsp_family_res_90conf_incl_rrna$Family_name <- rownames(vsp_family_res_90conf_incl_rrna)
vsp_family_res_90conf_incl_rrna$Targeted <- vsp_family_res_90conf_incl_rrna$Family_name %in% vsp_targets$family_name
vsp_family_res_90conf_incl_rrna$padj[is.na(vsp_family_res_90conf_incl_rrna$padj)] <- 1

#Volcano colored by targeted
ggplot(vsp_family_res_90conf_incl_rrna, 
       aes(x = log2FoldChange, y = -log10(padj), col = Targeted)) + 
  geom_point(aes(color = Targeted)) +
  geom_text(aes(
    label = ifelse((-log10(padj) > 1 & abs(log2FoldChange) > 1), Family_name, ""),
    vjust = -0.5, angle = -30), 
    size = 3.5, show.legend = FALSE) +
  theme_minimal() +
  xlim(-1, 15) +
  ylim(0, 30) +
  ggtitle("VSP vs. No Enrichment") +
  xlab(expression(paste("log"[2], " fold change"))) +
  ylab(expression(paste("-log"[10], " adjusted p-value"))) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  geom_hline(yintercept = 1.30103) +
  annotate("text", x = 0.5, y = 2, label = "adjusted p-value = 0.05")


###Look at targeted taxa (RPIP)
rpip_targets <- read_csv("input/rpip_panels/rpip_target_info.csv")
rpip_family_res_90conf_incl_rrna <- data.frame(rpipvsnone_family_res_90conf_incl_rrna)
rpip_family_res_90conf_incl_rrna$Family_name <- rownames(rpip_family_res_90conf_incl_rrna)
rpip_family_res_90conf_incl_rrna$Targeted <- rpip_family_res_90conf_incl_rrna$Family_name %in% rpip_targets$family
rpip_family_res_90conf_incl_rrna$padj[is.na(rpip_family_res_90conf_incl_rrna$padj)] <- 1

rpip_family_res_90conf_incl_rrna %>%
  filter(Targeted == TRUE) %>%
  View()

ggplot(rpip_family_res_90conf_incl_rrna, 
       aes(x = log2FoldChange, y = -log10(padj), col = Targeted)) + 
  geom_point() +
  geom_text(aes(label = ifelse((-log10(padj) > 5 & abs(log2FoldChange) > .1), Family_name, ""), vjust = -0.5, angle = -10), size = 3.5, show.legend = FALSE) +
  theme_minimal() +
  xlim(-25, 10) +
  ggtitle("RPIP vs. No Enrichment") +
  xlab(expression(paste("log"[2], " fold change"))) +
  ylab(expression(paste("-log"[10], " adjusted p-value"))) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  geom_hline(yintercept = 1.30103) +
  annotate("text", x = -22, y = 2, label = "adjusted p-value = 0.05")





#Contrast Nanotrap A
nanoAvsnone_family_res_90conf_incl_rrna <- results(family_dds_90conf_incl_rrna,
                                         contrast = c("Nanotrap_type", 
                                                      "A",
                                                      "none"),
                                         alpha = 0.05)
nanoAvsnone_family_res_90conf_incl_rrna <- lfcShrink(family_dds_90conf_incl_rrna, contrast = c("Nanotrap_type", "A", "none"),
                                           res = nanoAvsnone_family_res_90conf_incl_rrna, type = "ashr")

summary(nanoAvsnone_family_res_90conf_incl_rrna)
View(data.frame(nanoAvsnone_family_res_90conf_incl_rrna))


rpip_targets <- read_csv("input/rpip_panels/rpip_target_info.csv")
nanoA_by_targeted_family_res_90conf_incl_rrna <- data.frame(nanoAvsnone_family_res_90conf_incl_rrna)
nanoA_by_targeted_family_res_90conf_incl_rrna$Family_name <- rownames(nanoA_by_targeted_family_res_90conf_incl_rrna)
nanoA_by_targeted_family_res_90conf_incl_rrna$RPIP_Targeted <- nanoA_by_targeted_family_res_90conf_incl_rrna$Family_name %in% rpip_targets$family
nanoA_by_targeted_family_res_90conf_incl_rrna$VSP_Targeted <- nanoA_by_targeted_family_res_90conf_incl_rrna$Family_name %in% vsp_targets$family_name
nanoA_by_targeted_family_res_90conf_incl_rrna$padj[is.na(nanoA_by_targeted_family_res_90conf_incl_rrna$padj)] <- 1
nanoA_by_targeted_family_res_90conf_incl_rrna <- nanoA_by_targeted_family_res_90conf_incl_rrna %>%
  mutate(Targeted_by = case_when(RPIP_Targeted & VSP_Targeted ~ "Both",
                                 RPIP_Targeted & !VSP_Targeted ~ "RPIP",
                                 !RPIP_Targeted & VSP_Targeted ~ "VSP",
                                 !(RPIP_Targeted | VSP_Targeted) ~ "Neither"))

ggplot(nanoA_by_targeted_family_res_90conf_incl_rrna, 
       aes(x = log2FoldChange, y = -log10(padj), col = Targeted_by)) + 
  geom_point() +
#  geom_text(aes(label = ifelse((-log10(padj) > 1.30103), Family_name, ""), vjust = -0.5, angle = -30), size = 1, show.legend = FALSE) +
  theme_minimal() +
  xlim(-5, 3) +
  #ylim(0, 6) +
  ggtitle("Nanotrap A vs. No Concentration") +
  xlab(expression(paste("log"[2], " fold change"))) +
  ylab(expression(paste("-log"[10], " adjusted p-value"))) +
  scale_color_manual(values = c("#0072B2", "grey", "forestgreen", "#D55E00")) +
  geom_hline(yintercept = 1.30103) +
  annotate("text", x = -4, y = 2, label = "adjusted p-value = 0.05")


#Contrast Nanotrap AB v none
nanoBAvsnone_family_res_90conf_incl_rrna <- results(family_dds_90conf_incl_rrna,
                                                   contrast = c("Nanotrap_type", 
                                                                "A&B",
                                                                "none"),
                                                   alpha = 0.05)
nanoBAvsnone_family_res_90conf_incl_rrna <- lfcShrink(family_dds_90conf_incl_rrna, contrast = c("Nanotrap_type", "A&B", "none"),
                                                     res = nanoBAvsnone_family_res_90conf_incl_rrna, type = "ashr")

summary(nanoBAvsnone_family_res_90conf_incl_rrna)
View(data.frame(nanoBAvsnone_family_res_90conf_incl_rrna))

nanoBA_by_targeted_family_res_90conf_incl_rrna <- data.frame(nanoBAvsnone_family_res_90conf_incl_rrna)
nanoBA_by_targeted_family_res_90conf_incl_rrna$Family_name <- rownames(nanoBA_by_targeted_family_res_90conf_incl_rrna)
nanoBA_by_targeted_family_res_90conf_incl_rrna$RPIP_Targeted <- nanoBA_by_targeted_family_res_90conf_incl_rrna$Family_name %in% rpip_targets$family
nanoBA_by_targeted_family_res_90conf_incl_rrna$VSP_Targeted <- nanoBA_by_targeted_family_res_90conf_incl_rrna$Family_name %in% vsp_targets$family_name
nanoBA_by_targeted_family_res_90conf_incl_rrna$padj[is.na(nanoBA_by_targeted_family_res_90conf_incl_rrna$padj)] <- 1
nanoBA_by_targeted_family_res_90conf_incl_rrna <- nanoBA_by_targeted_family_res_90conf_incl_rrna %>%
  mutate(Targeted_by = case_when(RPIP_Targeted & VSP_Targeted ~ "Both",
                                 RPIP_Targeted & !VSP_Targeted ~ "RPIP",
                                 !RPIP_Targeted & VSP_Targeted ~ "VSP",
                                 !(RPIP_Targeted | VSP_Targeted) ~ "Neither"))

ggplot(nanoBA_by_targeted_family_res_90conf_incl_rrna, 
       aes(x = log2FoldChange, y = -log10(padj), col = Targeted_by)) + 
  geom_point() +
  #  geom_text(aes(label = ifelse((-log10(padj) > 1.30103), Family_name, ""), vjust = -0.5, angle = -30), size = 1, show.legend = FALSE) +
  theme_minimal() +
  xlim(-5, 3) +
  #ylim(0, 6) +
  ggtitle("Nanotrap A + B vs. No Concentration") +
  xlab(expression(paste("log"[2], " fold change"))) +
  ylab(expression(paste("-log"[10], " adjusted p-value"))) +
  scale_color_manual(values = c("#0072B2", "grey", "forestgreen", "#D55E00")) +
  geom_hline(yintercept = 1.30103) +
  annotate("text", x = -4, y = 2, label = "adjusted p-value = 0.05")


#Contrast Nanotrap AB vs A only
nanoBAvsA_family_res_90conf_incl_rrna <- results(family_dds_90conf_incl_rrna,
                                                    contrast = c("Nanotrap_type", 
                                                                 "A&B",
                                                                 "A"),
                                                    alpha = 0.05)
nanoBAvsA_family_res_90conf_incl_rrna <- lfcShrink(family_dds_90conf_incl_rrna, contrast = c("Nanotrap_type", "A&B", "A"),
                                                      res = nanoBAvsA_family_res_90conf_incl_rrna, type = "ashr")

summary(nanoBAvsA_family_res_90conf_incl_rrna)
View(data.frame(nanoBAvsA_family_res_90conf_incl_rrna))

nanoBAvsA_by_targeted_family_res_90conf_incl_rrna <- data.frame(nanoBAvsA_family_res_90conf_incl_rrna)
nanoBAvsA_by_targeted_family_res_90conf_incl_rrna$Family_name <- rownames(nanoBAvsA_by_targeted_family_res_90conf_incl_rrna)
nanoBAvsA_by_targeted_family_res_90conf_incl_rrna$RPIP_Targeted <- nanoBAvsA_by_targeted_family_res_90conf_incl_rrna$Family_name %in% rpip_targets$family
nanoBAvsA_by_targeted_family_res_90conf_incl_rrna$VSP_Targeted <- nanoBAvsA_by_targeted_family_res_90conf_incl_rrna$Family_name %in% vsp_targets$family_name
nanoBAvsA_by_targeted_family_res_90conf_incl_rrna$padj[is.na(nanoBAvsA_by_targeted_family_res_90conf_incl_rrna$padj)] <- 1
nanoBAvsA_by_targeted_family_res_90conf_incl_rrna <- nanoBAvsA_by_targeted_family_res_90conf_incl_rrna %>%
  mutate(Targeted_by = case_when(RPIP_Targeted & VSP_Targeted ~ "Both",
                                 RPIP_Targeted & !VSP_Targeted ~ "RPIP",
                                 !RPIP_Targeted & VSP_Targeted ~ "VSP",
                                 !(RPIP_Targeted | VSP_Targeted) ~ "Neither"))

ggplot(nanoBAvsA_by_targeted_family_res_90conf_incl_rrna, 
       aes(x = log2FoldChange, y = -log10(padj), col = Targeted_by)) + 
  geom_point() +
  geom_text(aes(label = ifelse((-log10(padj) > 1.30103), Family_name, ""), vjust = -0.5, angle = -30), size = 4, show.legend = FALSE) +
  theme_minimal() +
  xlim(-5, 3) +
  #ylim(0, 6) +
  ggtitle("Nanotrap A + B vs. A only") +
  xlab(expression(paste("log"[2], " fold change"))) +
  ylab(expression(paste("-log"[10], " adjusted p-value"))) +
  scale_color_manual(values = c("#0072B2", "grey", "forestgreen", "#D55E00")) +
  geom_hline(yintercept = 1.30103) +
  annotate("text", x = -4, y = 2, label = "adjusted p-value = 0.05")


#Contrast Nanotrap AB v none
nanoBAvsnone_family_res_90conf_incl_rrna <- results(family_dds_90conf_incl_rrna,
                                                    contrast = c("Nanotrap_type", 
                                                                 "A&B",
                                                                 "none"),
                                                    alpha = 0.05)
nanoBAvsnone_family_res_90conf_incl_rrna <- lfcShrink(family_dds_90conf_incl_rrna, contrast = c("Nanotrap_type", "A&B", "none"),
                                                      res = nanoBAvsnone_family_res_90conf_incl_rrna, type = "ashr")

summary(nanoBAvsnone_family_res_90conf_incl_rrna)
View(data.frame(nanoBAvsnone_family_res_90conf_incl_rrna))

nanoBA_by_targeted_family_res_90conf_incl_rrna <- data.frame(nanoBAvsnone_family_res_90conf_incl_rrna)
nanoBA_by_targeted_family_res_90conf_incl_rrna$Family_name <- rownames(nanoBA_by_targeted_family_res_90conf_incl_rrna)
nanoBA_by_targeted_family_res_90conf_incl_rrna$RPIP_Targeted <- nanoBA_by_targeted_family_res_90conf_incl_rrna$Family_name %in% rpip_targets$family
nanoBA_by_targeted_family_res_90conf_incl_rrna$VSP_Targeted <- nanoBA_by_targeted_family_res_90conf_incl_rrna$Family_name %in% vsp_targets$family_name
nanoBA_by_targeted_family_res_90conf_incl_rrna$padj[is.na(nanoBA_by_targeted_family_res_90conf_incl_rrna$padj)] <- 1
nanoBA_by_targeted_family_res_90conf_incl_rrna <- nanoBA_by_targeted_family_res_90conf_incl_rrna %>%
  mutate(Targeted_by = case_when(RPIP_Targeted & VSP_Targeted ~ "Both",
                                 RPIP_Targeted & !VSP_Targeted ~ "RPIP",
                                 !RPIP_Targeted & VSP_Targeted ~ "VSP",
                                 !(RPIP_Targeted | VSP_Targeted) ~ "Neither"))

ggplot(nanoBA_by_targeted_family_res_90conf_incl_rrna, 
       aes(x = log2FoldChange, y = -log10(padj), col = Targeted_by)) + 
  geom_point() +
  #  geom_text(aes(label = ifelse((-log10(padj) > 1.30103), Family_name, ""), vjust = -0.5, angle = -30), size = 1, show.legend = FALSE) +
  theme_minimal() +
  xlim(-5, 3) +
  #ylim(0, 6) +
  ggtitle("Nanotrap A + B vs. No Concentration") +
  xlab(expression(paste("log"[2], " fold change"))) +
  ylab(expression(paste("-log"[10], " adjusted p-value"))) +
  scale_color_manual(values = c("#0072B2", "grey", "forestgreen", "#D55E00")) +
  geom_hline(yintercept = 1.30103) +
  annotate("text", x = -4, y = 2, label = "adjusted p-value = 0.05")


#Contrast filtrate vs unfiltered
filtratevsUnfiltered_family_res_90conf_incl_rrna <- results(family_dds_90conf_incl_rrna,
                                                 contrast = c("Fraction", 
                                                              "filtrate",
                                                              "unfiltered"),
                                                 alpha = 0.05)
filtratevsUnfiltered_family_res_90conf_incl_rrna <- lfcShrink(family_dds_90conf_incl_rrna, contrast = c("Fraction", "filtrate", "unfiltered"),
                                                   res = filtratevsUnfiltered_family_res_90conf_incl_rrna, type = "ashr")

summary(filtratevsUnfiltered_family_res_90conf_incl_rrna)
View(data.frame(filtratevsUnfiltered_family_res_90conf_incl_rrna))

filtratevsUnfiltered_by_targeted_family_res_90conf_incl_rrna <- data.frame(filtratevsUnfiltered_family_res_90conf_incl_rrna)
filtratevsUnfiltered_by_targeted_family_res_90conf_incl_rrna$Family_name <- rownames(filtratevsUnfiltered_by_targeted_family_res_90conf_incl_rrna)
filtratevsUnfiltered_by_targeted_family_res_90conf_incl_rrna$RPIP_Targeted <- filtratevsUnfiltered_by_targeted_family_res_90conf_incl_rrna$Family_name %in% rpip_targets$family
filtratevsUnfiltered_by_targeted_family_res_90conf_incl_rrna$VSP_Targeted <- filtratevsUnfiltered_by_targeted_family_res_90conf_incl_rrna$Family_name %in% vsp_targets$family_name
filtratevsUnfiltered_by_targeted_family_res_90conf_incl_rrna$padj[is.na(filtratevsUnfiltered_by_targeted_family_res_90conf_incl_rrna$padj)] <- 1
filtratevsUnfiltered_by_targeted_family_res_90conf_incl_rrna <- filtratevsUnfiltered_by_targeted_family_res_90conf_incl_rrna %>%
  mutate(Targeted_by = case_when(RPIP_Targeted & VSP_Targeted ~ "Both",
                                 RPIP_Targeted & !VSP_Targeted ~ "RPIP",
                                 !RPIP_Targeted & VSP_Targeted ~ "VSP",
                                 !(RPIP_Targeted | VSP_Targeted) ~ "Neither"))

filtratevsUnfiltered_by_targeted_family_res_90conf_incl_rrna$Targeted_by = relevel(as.factor(filtratevsUnfiltered_by_targeted_family_res_90conf_incl_rrna$Targeted_by), ref = "Neither")

ggplot(arrange(filtratevsUnfiltered_by_targeted_family_res_90conf_incl_rrna, Targeted_by), 
       aes(x = log2FoldChange, y = -log10(padj), col = Targeted_by)) + 
  geom_point() +
  geom_text(aes(label = ifelse((-log10(padj) > 1.30103 & abs(log2FoldChange) > 2 & Targeted_by != "Neither"), Family_name, ""), vjust = -0.5, angle = -30), size = 4, show.legend = FALSE) +
  theme_minimal() +
#  xlim(-5, 3) +
  #ylim(0, 6) +
  ggtitle("Filtrate vs. unfiltered samples") +
  xlab(expression(paste("log"[2], " fold change"))) +
  ylab(expression(paste("-log"[10], " adjusted p-value"))) +
  scale_color_manual(values = c("grey", "#0072B2", "forestgreen", "#D55E00")) +
  geom_hline(yintercept = 1.30103) +
  annotate("text", x = -4, y = 2, label = "adjusted p-value = 0.05")


#Contrast retentate vs unfiltered
retentatevsUnfiltered_family_res_90conf_incl_rrna <- results(family_dds_90conf_incl_rrna,
                                                            contrast = c("Fraction", 
                                                                         "retentate",
                                                                         "unfiltered"),
                                                            alpha = 0.05)
retentatevsUnfiltered_family_res_90conf_incl_rrna <- lfcShrink(family_dds_90conf_incl_rrna, contrast = c("Fraction", "retentate", "unfiltered"),
                                                              res = retentatevsUnfiltered_family_res_90conf_incl_rrna, type = "ashr")

summary(retentatevsUnfiltered_family_res_90conf_incl_rrna)
View(data.frame(retentatevsUnfiltered_family_res_90conf_incl_rrna))

retentatevsUnfiltered_by_targeted_family_res_90conf_incl_rrna <- data.frame(retentatevsUnfiltered_family_res_90conf_incl_rrna)
retentatevsUnfiltered_by_targeted_family_res_90conf_incl_rrna$Family_name <- rownames(retentatevsUnfiltered_by_targeted_family_res_90conf_incl_rrna)
retentatevsUnfiltered_by_targeted_family_res_90conf_incl_rrna$RPIP_Targeted <- retentatevsUnfiltered_by_targeted_family_res_90conf_incl_rrna$Family_name %in% rpip_targets$family
retentatevsUnfiltered_by_targeted_family_res_90conf_incl_rrna$VSP_Targeted <- retentatevsUnfiltered_by_targeted_family_res_90conf_incl_rrna$Family_name %in% vsp_targets$family_name
retentatevsUnfiltered_by_targeted_family_res_90conf_incl_rrna$padj[is.na(retentatevsUnfiltered_by_targeted_family_res_90conf_incl_rrna$padj)] <- 1
retentatevsUnfiltered_by_targeted_family_res_90conf_incl_rrna <- retentatevsUnfiltered_by_targeted_family_res_90conf_incl_rrna %>%
  mutate(Targeted_by = case_when(RPIP_Targeted & VSP_Targeted ~ "Both",
                                 RPIP_Targeted & !VSP_Targeted ~ "RPIP",
                                 !RPIP_Targeted & VSP_Targeted ~ "VSP",
                                 !(RPIP_Targeted | VSP_Targeted) ~ "Neither"))

retentatevsUnfiltered_by_targeted_family_res_90conf_incl_rrna$Targeted_by = relevel(as.factor(retentatevsUnfiltered_by_targeted_family_res_90conf_incl_rrna$Targeted_by), ref = "Neither")


ggplot(arrange(retentatevsUnfiltered_by_targeted_family_res_90conf_incl_rrna, Targeted_by), 
       aes(x = log2FoldChange, y = -log10(padj), col = Targeted_by)) + 
  geom_point(aes(alpha = 0.5)) +
  #geom_text(aes(label = ifelse((-log10(padj) > 1.30103 | abs(log2FoldChange) > 5), Family_name, ""), vjust = -0.5, angle = -30), size = 4, show.legend = FALSE) +
  theme_minimal() +
  #  xlim(-5, 3) +
  #ylim(0, 6) +
  ggtitle("Retentate vs. unfiltered samples") +
  xlab(expression(paste("log"[2], " fold change"))) +
  ylab(expression(paste("-log"[10], " adjusted p-value"))) +
  scale_color_manual(values = c("grey", "#0072B2", "forestgreen", "#D55E00")) +
  geom_hline(yintercept = 1.30103) +
  annotate("text", x = -4, y = 2, label = "adjusted p-value = 0.05")

#Boxplots
