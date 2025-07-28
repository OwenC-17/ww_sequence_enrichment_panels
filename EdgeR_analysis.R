library(edgeR)
library(tidyverse)
library(ordinal)

edger_family_count_table_w_rrna <- read_csv("DS2_family_count_matrix_w_rrna.csv")

edger_family_count_table_w_rrna <- edger_family_count_table_w_rrna %>%
  dplyr::rename(Family = `F`)

colnames(edger_family_count_table_w_rrna)

group_data <- read_csv("DS2_sample_metadata_w_rrna.csv") %>%
  arrange(by = UniqueID)

######K2 0 confidence EdgeR
group_data_0conf <- group_data %>%
  filter(Kraken2_confidence == 0)

edger_family_count_table_w_rrna_0conf <- edger_family_count_table_w_rrna[,c("Family", group_data_0conf$UniqueID)]

##Confirm orders are correct
colnames(edger_family_count_table_w_rrna_0conf[,2:ncol(edger_family_count_table_w_rrna_0conf)]) == group_data_0conf$UniqueID

family_0conf <- DGEList(counts = edger_family_count_table_w_rrna_0conf)

colnames(group_data_0conf)

LIMS_ID <- factor(group_data_0conf$LIMS_ID)
Fraction <- factor(group_data_0conf$Fraction, levels = c("unfiltered", "retentate", "filtrate"))
Nanotrap_type <- factor(group_data_0conf$Nanotrap_type, levels = c("none", "A", "A&B"))
Enrichment <- factor(group_data_0conf$Enrichment, levels = c("None", "RPIP", "VSP"))


treat_list <- list()

for (x in levels(Fraction)) {
  for (y in levels(Nanotrap_type)) {
    for (z in levels(Enrichment)) {
      tx_name <- paste(x, y, z, sep = ".")
      assign(tx_name, (Fraction == x & Nanotrap_type == y & Enrichment == z))
      treat_list <- c(treat_list, tx_name)
    }
  }
}


design <- model.matrix(~LIMS_ID)
additional_groups <- mget(unlist(treat_list))

for (group in 1:27) {
  design <- cbind(design, unlist(additional_groups[group]))
}

colnames(design)[9:ncol(design)] <- names(additional_groups)

design <- design[,c(1:8,10:ncol(design))]

design <- drop.coef(design)
View(design)
disp_families_0conf <- estimateDisp(y = family_0conf, design = design)
fit_families_0conf <- glmQLFit(disp_families_0conf, design)


is.RPIP <- str_detect(colnames(design), "RPIP")
is.VSP <- str_detect(colnames(design), "VSP")
is.Untargeted <- str_detect(colnames(design), "None")
is.Filtrate <- str_detect(colnames(design), "filtrate")
is.Retentate <- str_detect(colnames(design), "retentate")
is.Unfiltered <- str_detect(colnames(design), "unfiltered")
is.NanoA <- str_detect(colnames(design), "A\\.")
is.NanoB <- str_detect(colnames(design), "A&B")
is.DirectExt <- str_detect(colnames(design), "none")


colnames(design)
RPIP_all_treatments_with_rrna <- glmQLFTest(fit_families_0conf, contrast = is.RPIP - is.Untargeted)
topTags(RPIP_all_treatments_with_rrna2)

VSP_all_treatments_with_rrna <- glmQLFTest(fit_families_0conf, contrast = is.VSP - is.Untargeted)
topTags(VSP_all_treatments_with_rrna)

NanoA_all_other_treatments_with_rrna <- glmQLFTest(fit_families_0conf, contrast = is.NanoA - is.DirectExt)
topTags(NanoA_all_other_treatments_with_rrna)

NanoBA_all_other_treatments_with_rrna <- glmQLFTest(fit_families_0conf, contrast = is.NanoB - is.DirectExt)
topTags(NanoBA_all_other_treatments_with_rrna)

Filtrate_all_other_treatments_with_rrna <- glmQLFTest()

################################################################################
######K2 90 confidence EdgeR


group_data_90conf <- group_data %>%
  filter(Kraken2_confidence == 0.9)

edger_family_count_table_w_rrna_90conf <- edger_family_count_table_w_rrna[,c("Family", group_data_90conf$UniqueID)]

##Confirm orders are correct
sum(colnames(edger_family_count_table_w_rrna_90conf[,2:ncol(edger_family_count_table_w_rrna_90conf)]) != group_data_90conf$UniqueID)

family_90conf <- DGEList(counts = edger_family_count_table_w_rrna_90conf)

colnames(group_data_90conf)

generate_levels <- function(group_df) {
  LIMS_ID <- factor(group_df$LIMS_ID)
  Fraction <- factor(group_df$Fraction, levels = c("unfiltered", "retentate", "filtrate"))
  Nanotrap_type <- factor(group_df$Nanotrap_type, levels = c("none", "A", "A&B"))
  Enrichment <- factor(group_df$Enrichment, levels = c("None", "RPIP", "VSP"))
  treat_list <- list()
  for (x in levels(Fraction)) {
    for (y in levels(Nanotrap_type)) {
      for (z in levels(Enrichment)) {
        tx_name <- paste(x, y, z, sep = ".")
        assign(tx_name, (Fraction == x & Nanotrap_type == y & Enrichment == z))
        treat_list <- c(treat_list, tx_name)
      }
    }
  }
  additional_groups <- mget(unlist(treat_list))
  design <- model.matrix(~LIMS_ID)
  for (group in 1:length(additional_groups)) {
    design <- cbind(design, unlist(additional_groups[group]))
  }
  colnames(design)[(length(levels(LIMS_ID)) + 1):ncol(design)] <- names(additional_groups)
  design <- design[, colnames(design) != "unfiltered.none.None"]
  return(design)
}


design <- generate_levels(group_data_90conf)
#Assuming group_data is sorted, this design matrix will be identical to the one 
#produced from 0conf data.

disp_families_90conf <- estimateDisp(y = family_90conf, design = design)
fit_families_90conf <- glmQLFit(disp_families_90conf, design)


is.RPIP <- str_detect(colnames(design), "RPIP")
is.VSP <- str_detect(colnames(design), "VSP")
is.Untargeted <- str_detect(colnames(design), "None")
is.Filtrate <- str_detect(colnames(design), "filtrate")
is.Retentate <- str_detect(colnames(design), "retentate")
is.Unfiltered <- str_detect(colnames(design), "unfiltered")
is.NanoA <- str_detect(colnames(design), "A\\.")
is.NanoB <- str_detect(colnames(design), "A&B")
is.DirectExt <- str_detect(colnames(design), "none")


colnames(design)
RPIP_all_treatments_90conf_with_rrna <- glmQLFTest(fit_families_90conf, contrast = is.RPIP - is.Untargeted)
topTags(RPIP_all_treatments_90conf_with_rrna)

VSP_all_treatments_90conf_with_rrna <- glmQLFTest(fit_families_90conf, contrast = is.VSP - is.Untargeted)
topTags(VSP_all_treatments_90conf_with_rrna)

VSP_filtrate_allNT_90conf_with_rrna <- glmQLFTest(fit_families_90conf, contrast = (is.Filtrate & is.VSP) - (is.Untargeted & is.VSP))
topTags(VSP_filtrate_allNT_90conf_with_rrna)

NanoA_all_other_treatments_90conf_with_rrna <- glmQLFTest(fit_families_90conf, contrast = is.NanoA - is.DirectExt)
topTags(NanoA_all_other_treatments_90conf_with_rrna)

NanoBA_all_other_treatments_90conf_with_rrna <- glmQLFTest(fit_families_90conf, contrast = is.NanoB - is.DirectExt)
topTags(NanoBA_all_other_treatments_90conf_with_rrna)


Filtrate_all_other_treatments_90conf_with_rrna <- glmQLFTest(fit_families_90conf, contrast = is.Filtrate - is.Unfiltered)
topTags(Filtrate_all_other_treatments_90conf_with_rrna)

################################################################################
######K2 90 confidence no rrna EdgeR
edger_family_count_table_no_rrna <- read_csv("DS2_family_count_matrix_no_rrna.csv")
edger_family_count_table_no_rrna <- edger_family_count_table_no_rrna %>%
  dplyr::rename(Family = `F`)
edger_family_count_table_no_rrna[str_starts(edger_family_count_table_no_rrna$Family, "unid"), "Family"] <- "unclassified"


colnames(edger_family_count_table_no_rrna)

edger_family_count_table_no_rrna <- edger_family_count_table_no_rrna %>%
  group_by(Family) %>%
  summarize(across(everything(), sum), .groups = "drop")

group_data_no_rrna <- read_csv("DS2_sample_metadata.csv") %>%
  arrange(by = UniqueID)
group_data_90conf_no_rrna <- group_data_no_rrna %>%
  filter(Kraken2_confidence == 0.9)

edger_family_count_table_no_rrna_90conf <- edger_family_count_table_no_rrna[,c("Family", group_data_90conf_no_rrna$UniqueID)]

##Confirm orders are correct
sum(colnames(edger_family_count_table_no_rrna_90conf[,2:ncol(edger_family_count_table_no_rrna_90conf)]) != group_data_90conf$UniqueID)

family_90conf_no_rrna <- DGEList(counts = edger_family_count_table_no_rrna_90conf)
expressed <- filterByExpr(family_90conf_no_rrna, design = design)
family_90conf_no_rrna_expressed <- family_90conf_no_rrna[expressed, , keep.lib.sizes = FALSE]

colnames(group_data_90conf_no_rrna)

generate_levels <- function(group_df) {
  LIMS_ID <- factor(group_df$LIMS_ID)
  Fraction <- factor(group_df$Fraction, levels = c("unfiltered", "retentate", "filtrate"))
  Nanotrap_type <- factor(group_df$Nanotrap_type, levels = c("none", "A", "A&B"))
  Enrichment <- factor(group_df$Enrichment, levels = c("None", "RPIP", "VSP"))
  treat_list <- list()
  for (x in levels(Fraction)) {
    for (y in levels(Nanotrap_type)) {
      for (z in levels(Enrichment)) {
        tx_name <- paste(x, y, z, sep = ".")
        assign(tx_name, (Fraction == x & Nanotrap_type == y & Enrichment == z))
        treat_list <- c(treat_list, tx_name)
      }
    }
  }
  additional_groups <- mget(unlist(treat_list))
  design <- model.matrix(~LIMS_ID)
  for (group in 1:length(additional_groups)) {
    design <- cbind(design, unlist(additional_groups[group]))
  }
  colnames(design)[(length(levels(LIMS_ID)) + 1):ncol(design)] <- names(additional_groups)
  design <- design[, colnames(design) != "unfiltered.none.None"]
  return(design)
}


design <- generate_levels(group_data_90conf_no_rrna)
#Assuming group_data is sorted, this design matrix will be identical to the one 
#produced previously (i.e. with rrna or with 0conf).

disp_families_90conf_no_rrna <- estimateDisp(y = family_90conf_no_rrna, design = design)
fit_families_90conf_no_rrna <- glmQLFit(disp_families_90conf_no_rrna, design)

disp_families_90conf_no_rrna_expressed <- estimateDisp(y = family_90conf_no_rrna_expressed, design = design)
fit_families_90conf_no_rrna_expressed <- glmQLFit(disp_families_90conf_no_rrna_expressed, design)


is.RPIP <- str_detect(colnames(design), "RPIP")
is.VSP <- str_detect(colnames(design), "VSP")
is.Untargeted <- str_detect(colnames(design), "None")
is.Filtrate <- str_detect(colnames(design), "filtrate")
is.Retentate <- str_detect(colnames(design), "retentate")
is.Unfiltered <- str_detect(colnames(design), "unfiltered")
is.NanoA <- str_detect(colnames(design), "A\\.")
is.NanoB <- str_detect(colnames(design), "A&B")
is.DirectExt <- str_detect(colnames(design), "none")


colnames(design)
#############################
#########CONTRASTS###########
#############################

top25 <- function(result) {
  topTags(result, n = 25, sort.by = "logFC")
}

###RPIP - Untargeted:
#Main effect:
RPIP_all_treatments_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = is.RPIP - is.Untargeted)
top25(RPIP_all_treatments_90conf_no_rrna)
#No NT, no filtration:
RPIP_vs_Unt_DirExt_Unfiltered_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.RPIP & is.Unfiltered & is.DirectExt) - (is.Untargeted & is.Unfiltered & is.DirectExt))
top25(RPIP_vs_Unt_DirExt_Unfiltered_90conf_no_rrna)

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




#####VSP - Untargeted
#Main effect:
VSP_all_treatments_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = is.VSP - is.Untargeted)
top25(VSP_all_treatments_90conf_no_rrna)

#No NT, no filtration:
VSP_vs_Unt_DirExt_Unfiltered_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.VSP & is.Unfiltered & is.DirectExt) - (is.Untargeted & is.Unfiltered & is.DirectExt))
top25(VSP_vs_Unt_DirExt_Unfiltered_90conf_no_rrna)

#NTA, no filtration:
VSP_vs_Unt_NanoA_Unfiltered_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.VSP & is.Unfiltered & is.NanoA) - (is.Untargeted & is.Unfiltered & is.NanoA))
top25(VSP_vs_Unt_NanoA_Unfiltered_90conf_no_rrna)

#NTB, no filtration:
VSP_vs_Unt_NanoBA_Unfiltered_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.VSP & is.NanoB & is.Unfiltered) - (is.Untargeted & is.NanoB & is.Unfiltered))
top25(VSP_vs_Unt_NanoBA_Unfiltered_90conf_no_rrna)

#No NT, filtrate:
VSP_vs_Unt_DirExt_Filtrate_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.VSP & is.Filtrate & is.DirectExt) - (is.Untargeted & is.Filtrate & is.DirectExt))
top25(VSP_vs_Unt_DirExt_Filtrate_90conf_no_rrna)

#NTA, filtrate:
VSP_vs_Unt_NanoA_Filtrate_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.VSP & is.NanoA & is.Filtrate) - (is.Untargeted & is.NanoA & is.Filtrate))
top25(VSP_vs_Unt_NanoA_Filtrate_90conf_no_rrna)

#NTB, filtrate:
VSP_vs_Unt_NanoBA_Filtrate_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.VSP & is.NanoB & is.Filtrate) - (is.Untargeted & is.NanoB & is.Filtrate))
top25(VSP_vs_Unt_NanoBA_Filtrate_90conf_no_rrna)

#No NT, retentate:
VSP_vs_Unt_DirExt_Retentate_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.VSP & is.Retentate & is.DirectExt) - (is.Untargeted & is.Retentate & is.DirectExt))
top25(VSP_vs_Unt_DirExt_Retentate_90conf_no_rrna)

#NTA, retentate:
VSP_vs_Unt_NanoA_Retentate_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.VSP & is.NanoA & is.Retentate) - (is.Untargeted & is.NanoA & is.Retentate))
top25(VSP_vs_Unt_NanoA_Retentate_90conf_no_rrna)

#NTB, retentate:
VSP_vs_Unt_NanoBA_Retentate_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.VSP & is.NanoB & is.Retentate) - (is.Untargeted & is.NanoB & is.Retentate))
top25(VSP_vs_Unt_NanoBA_Retentate_90conf_no_rrna)


########################
##########Nanotrap contrasts
########################
#NanoA main effect:
NanoA_vs_DirExt_main_effect_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = is.NanoA - is.DirectExt)
top25(NanoA_vs_DirExt_main_effect_90conf_no_rrna)

#NanoA, unfiltered, untargeted:
NanoA_vs_DirExt_Unfiltered_Untargeted_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.NanoA & is.Unfiltered & is.Untargeted) - (is.DirectExt & is.Unfiltered & is.Untargeted))
top25(NanoA_vs_DirExt_Unfiltered_Untargeted_90conf_no_rrna)

#NanoA, unfiltered, VSP:
NanoA_vs_DirExt_Unfiltered_VSP_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.NanoA & is.Unfiltered & is.VSP) - (is.DirectExt & is.Unfiltered & is.VSP))
top25(NanoA_vs_DirExt_Unfiltered_VSP_90conf_no_rrna)

#NanoA, unfiltered, RPIP:
NanoA_vs_DirExt_Unfiltered_RPIP_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.NanoA & is.Unfiltered & is.RPIP) - (is.DirectExt & is.Unfiltered & is.RPIP))
top25(NanoA_vs_DirExt_Unfiltered_RPIP_90conf_no_rrna)

#NanoA, filtrate, untargeted:
NanoA_vs_DirExt_Filtrate_Untargeted_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.NanoA & is.Filtrate & is.Untargeted) - (is.DirectExt & is.Filtrate & is.Untargeted))
top25(NanoA_vs_DirExt_Filtrate_Untargeted_90conf_no_rrna)

#NanoA, filtrate, VSP:
NanoA_vs_DirExt_Filtrate_VSP_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.NanoA & is.Filtrate & is.VSP) - (is.DirectExt & is.Filtrate & is.VSP))
top25(NanoA_vs_DirExt_Filtrate_VSP_90conf_no_rrna)

#NanoA, filtrate, RPIP:
NanoA_vs_DirExt_Filtrate_RPIP_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.NanoA & is.Filtrate & is.RPIP) - (is.DirectExt & is.Filtrate & is.RPIP))
top25(NanoA_vs_DirExt_Filtrate_RPIP_90conf_no_rrna)

#NanoA, Retentate, untargeted:
NanoA_vs_DirExt_Retentate_Untargeted_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.NanoA & is.Retentate & is.Untargeted) - (is.DirectExt & is.Retentate & is.Untargeted))
top25(NanoA_vs_DirExt_Retentate_Untargeted_90conf_no_rrna)

#NanoA, Retentate, VSP:
NanoA_vs_DirExt_Retentate_VSP_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.NanoA & is.Retentate & is.VSP) - (is.DirectExt & is.Retentate & is.VSP))
top25(NanoA_vs_DirExt_Retentate_VSP_90conf_no_rrna)

#NanoA, Retentate, RPIP:
NanoA_vs_DirExt_Retentate_RPIP_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.NanoA & is.Retentate & is.RPIP) - (is.DirectExt & is.Retentate & is.RPIP))
top25(NanoA_vs_DirExt_Retentate_RPIP_90conf_no_rrna)

#NanoB main effect:
NanoB_vs_DirExt_main_effect_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = is.NanoB - is.DirectExt)
top25(NanoB_vs_DirExt_main_effect_90conf_no_rrna)

#NanoB, unfiltered, untargeted:
NanoB_vs_DirExt_Unfiltered_Untargeted_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.NanoB & is.Unfiltered & is.Untargeted) - (is.DirectExt & is.Unfiltered & is.Untargeted))
top25(NanoB_vs_DirExt_Unfiltered_Untargeted_90conf_no_rrna)

#NanoB, unfiltered, VSP:
NanoB_vs_DirExt_Unfiltered_VSP_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.NanoB & is.Unfiltered & is.VSP) - (is.DirectExt & is.Unfiltered & is.VSP))
top25(NanoB_vs_DirExt_Unfiltered_VSP_90conf_no_rrna)

#NanoB, unfiltered, RPIP:
NanoB_vs_DirExt_Unfiltered_RPIP_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.NanoB & is.Unfiltered & is.RPIP) - (is.DirectExt & is.Unfiltered & is.RPIP))
top25(NanoB_vs_DirExt_Unfiltered_RPIP_90conf_no_rrna)

#NanoB, filtrate, untargeted:
NanoB_vs_DirExt_Filtrate_Untargeted_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.NanoB & is.Filtrate & is.Untargeted) - (is.DirectExt & is.Filtrate & is.Untargeted))
top25(NanoB_vs_DirExt_Filtrate_Untargeted_90conf_no_rrna)

#NanoB, filtrate, VSP:
NanoB_vs_DirExt_Filtrate_VSP_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.NanoB & is.Filtrate & is.VSP) - (is.DirectExt & is.Filtrate & is.VSP))
top25(NanoB_vs_DirExt_Filtrate_VSP_90conf_no_rrna)

#NanoB, filtrate, RPIP:
NanoB_vs_DirExt_Filtrate_RPIP_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.NanoB & is.Filtrate & is.RPIP) - (is.DirectExt & is.Filtrate & is.RPIP))
top25(NanoB_vs_DirExt_Filtrate_RPIP_90conf_no_rrna)

#NanoB, Retentate, untargeted:
NanoB_vs_DirExt_Retentate_Untargeted_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.NanoB & is.Retentate & is.Untargeted) - (is.DirectExt & is.Retentate & is.Untargeted))
top25(NanoB_vs_DirExt_Retentate_Untargeted_90conf_no_rrna)

#NanoB, Retentate, VSP:
NanoB_vs_DirExt_Retentate_VSP_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.NanoB & is.Retentate & is.VSP) - (is.DirectExt & is.Retentate & is.VSP))
top25(NanoB_vs_DirExt_Retentate_VSP_90conf_no_rrna)

#NanoB, Retentate, RPIP:
NanoB_vs_DirExt_Retentate_RPIP_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.NanoB & is.Retentate & is.RPIP) - (is.DirectExt & is.Retentate & is.RPIP))
top25(NanoB_vs_DirExt_Retentate_RPIP_90conf_no_rrna)

##################
#Compare Nanotrap Types
##################

#NanoB vs NanoA main effect:
NanoB_vs_NanoA_main_effect_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = is.NanoB - is.NanoA)
top25(NanoB_vs_NanoA_main_effect_90conf_no_rrna)

#NanoB vs NanoA, unfiltered, untargeted:
NanoB_vs_NanoA_Unfiltered_Untargeted_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.NanoB & is.Unfiltered & is.Untargeted) - (is.NanoA & is.Unfiltered & is.Untargeted))
top25(NanoB_vs_NanoA_Unfiltered_Untargeted_90conf_no_rrna)

#NanoB vs NanoA, unfiltered, VSP:
NanoB_vs_NanoA_Unfiltered_VSP_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.NanoB & is.Unfiltered & is.VSP) - (is.NanoA & is.Unfiltered & is.VSP))
top25(NanoB_vs_NanoA_Unfiltered_VSP_90conf_no_rrna)

#NanoB vs NanoA, unfiltered, RPIP:
NanoB_vs_NanoA_Unfiltered_RPIP_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.NanoB & is.Unfiltered & is.RPIP) - (is.NanoA & is.Unfiltered & is.RPIP))
top25(NanoB_vs_NanoA_Unfiltered_RPIP_90conf_no_rrna)

#NanoB vs NanoA, filtrate, untargeted:
NanoB_vs_NanoA_Filtrate_Untargeted_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.NanoB & is.Filtrate & is.Untargeted) - (is.NanoA & is.Filtrate & is.Untargeted))
top25(NanoB_vs_NanoA_Filtrate_Untargeted_90conf_no_rrna)

#NanoB vs NanoA, filtrate, VSP:
NanoB_vs_NanoA_Filtrate_VSP_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.NanoB & is.Filtrate & is.VSP) - (is.NanoA & is.Filtrate & is.VSP))
top25(NanoB_vs_NanoA_Filtrate_VSP_90conf_no_rrna)

#NanoB vs NanoA, filtrate, RPIP:
NanoB_vs_NanoA_Filtrate_RPIP_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.NanoB & is.Filtrate & is.RPIP) - (is.NanoA & is.Filtrate & is.RPIP))
top25(NanoB_vs_NanoA_Filtrate_RPIP_90conf_no_rrna)

#NanoB vs NanoA, Retentate, untargeted:
NanoB_vs_NanoA_Retentate_Untargeted_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.NanoB & is.Retentate & is.Untargeted) - (is.NanoA & is.Retentate & is.Untargeted))
top25(NanoB_vs_NanoA_Retentate_Untargeted_90conf_no_rrna)

#NanoB vs NanoA, Retentate, VSP:
NanoB_vs_NanoA_Retentate_VSP_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.NanoB & is.Retentate & is.VSP) - (is.NanoA & is.Retentate & is.VSP))
top25(NanoB_vs_NanoA_Retentate_VSP_90conf_no_rrna)

#NanoB vs NanoA, Retentate, RPIP:
NanoB_vs_NanoA_Retentate_RPIP_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.NanoB & is.Retentate & is.RPIP) - (is.NanoA & is.Retentate & is.RPIP))
top25(NanoB_vs_NanoA_Retentate_RPIP_90conf_no_rrna)

##################
###Effect of Filtration
##################
#Filtrate vs. unfiltered main effect:
Filtrate_vs_Unfiltered_main_effect_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = is.Filtrate - is.Unfiltered)
top25(Filtrate_vs_Unfiltered_main_effect_90conf_no_rrna)

#Filtrate vs. unfiltered, no Nanotrap, Untargeted
Filtrate_vs_Unfiltered_DirExt_Untargeted <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.Filtrate & is.DirectExt & is.Untargeted) - (is.Unfiltered & is.DirectExt & is.Untargeted))
top25(Filtrate_vs_Unfiltered_DirExt_Untargeted)

#Filtrate vs. unfiltered, no Nanotrap, VSP
Filtrate_vs_Unfiltered_DirExt_VSP <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.Filtrate & is.DirectExt & is.VSP) - (is.Unfiltered & is.DirectExt & is.VSP))
top25(Filtrate_vs_Unfiltered_DirExt_VSP)

#Filtrate vs. unfiltered, no Nanotrap, RPIP
Filtrate_vs_Unfiltered_DirExt_RPIP <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.Filtrate & is.DirectExt & is.RPIP) - (is.Unfiltered & is.DirectExt & is.RPIP))
top25(Filtrate_vs_Unfiltered_DirExt_RPIP)

#Filtrate vs. unfiltered, NanoA, Untargeted
Filtrate_vs_Unfiltered_NanoA_Untargeted <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.Filtrate & is.NanoA & is.Untargeted) - (is.Unfiltered & is.NanoA & is.Untargeted))
top25(Filtrate_vs_Unfiltered_NanoA_Untargeted)

#Filtrate vs. unfiltered, NanoA, VSP
Filtrate_vs_Unfiltered_NanoA_VSP <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.Filtrate & is.NanoA & is.VSP) - (is.Unfiltered & is.NanoA & is.VSP))
top25(Filtrate_vs_Unfiltered_NanoA_VSP)

#Filtrate vs. unfiltered, NanoA, RPIP
Filtrate_vs_Unfiltered_NanoA_RPIP <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.Filtrate & is.NanoA & is.RPIP) - (is.Unfiltered & is.NanoA & is.RPIP))
top25(Filtrate_vs_Unfiltered_NanoA_RPIP)

#Filtrate vs. unfiltered, NanoB, Untargeted
Filtrate_vs_Unfiltered_NanoB_Untargeted <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.Filtrate & is.NanoB & is.Untargeted) - (is.Unfiltered & is.NanoB & is.Untargeted))
top25(Filtrate_vs_Unfiltered_NanoB_Untargeted)

#Filtrate vs. unfiltered, NanoB, VSP
Filtrate_vs_Unfiltered_NanoB_VSP <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.Filtrate & is.NanoB & is.VSP) - (is.Unfiltered & is.NanoB & is.VSP))
top25(Filtrate_vs_Unfiltered_NanoB_VSP)

#Filtrate vs. unfiltered, NanoB, RPIP
Filtrate_vs_Unfiltered_NanoB_RPIP <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.Filtrate & is.NanoB & is.RPIP) - (is.Unfiltered & is.NanoB & is.RPIP))
top25(Filtrate_vs_Unfiltered_NanoB_RPIP)

##################
#Retentate vs. unfiltered main effect:
Retentate_vs_Unfiltered_main_effect_90conf_no_rrna <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = is.Retentate - is.Unfiltered)
top25(Retentate_vs_Unfiltered_main_effect_90conf_no_rrna)

#Retentate vs. unfiltered, no Nanotrap, Untargeted
Retentate_vs_Unfiltered_DirExt_Untargeted <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.Retentate & is.DirectExt & is.Untargeted) - (is.Unfiltered & is.DirectExt & is.Untargeted))
top25(Retentate_vs_Unfiltered_DirExt_Untargeted)

#Retentate vs. unfiltered, no Nanotrap, VSP
Retentate_vs_Unfiltered_DirExt_VSP <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.Retentate & is.DirectExt & is.VSP) - (is.Unfiltered & is.DirectExt & is.VSP))
top25(Retentate_vs_Unfiltered_DirExt_VSP)

#Retentate vs. unfiltered, no Nanotrap, RPIP
Retentate_vs_Unfiltered_DirExt_RPIP <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.Retentate & is.DirectExt & is.RPIP) - (is.Unfiltered & is.DirectExt & is.RPIP))
top25(Retentate_vs_Unfiltered_DirExt_RPIP)

#Retentate vs. unfiltered, NanoA, Untargeted
Retentate_vs_Unfiltered_NanoA_Untargeted <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.Retentate & is.NanoA & is.Untargeted) - (is.Unfiltered & is.NanoA & is.Untargeted))
top25(Retentate_vs_Unfiltered_NanoA_Untargeted)

#Retentate vs. unfiltered, NanoA, VSP
Retentate_vs_Unfiltered_NanoA_VSP <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.Retentate & is.NanoA & is.VSP) - (is.Unfiltered & is.NanoA & is.VSP))
top25(Retentate_vs_Unfiltered_NanoA_VSP)

#Retentate vs. unfiltered, NanoA, RPIP
Retentate_vs_Unfiltered_NanoA_RPIP <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.Retentate & is.NanoA & is.RPIP) - (is.Unfiltered & is.NanoA & is.RPIP))
top25(Retentate_vs_Unfiltered_NanoA_RPIP)

#Retentate vs. unfiltered, NanoB, Untargeted
Retentate_vs_Unfiltered_NanoB_Untargeted <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.Retentate & is.NanoB & is.Untargeted) - (is.Unfiltered & is.NanoB & is.Untargeted))
top25(Retentate_vs_Unfiltered_NanoB_Untargeted)

#Retentate vs. unfiltered, NanoB, VSP
Retentate_vs_Unfiltered_NanoB_VSP <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.Retentate & is.NanoB & is.VSP) - (is.Unfiltered & is.NanoB & is.VSP))
top25(Retentate_vs_Unfiltered_NanoB_VSP)

#Retentate vs. unfiltered, NanoB, RPIP
Retentate_vs_Unfiltered_NanoB_RPIP <- glmQLFTest(fit_families_90conf_no_rrna_expressed, contrast = (is.Retentate & is.NanoB & is.RPIP) - (is.Unfiltered & is.NanoB & is.RPIP))
top25(Retentate_vs_Unfiltered_NanoB_RPIP)

