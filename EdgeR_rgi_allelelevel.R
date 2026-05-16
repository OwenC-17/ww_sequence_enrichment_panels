library(tidyverse)
library(edgeR)
library(ggh4x)
library(paletteer)

##########################################################
#####Prepare allele-level RGI data for EdgeR analysis#####
##########################################################
#Load the explanatory variables (group data) from Kraken analysis
#and modify for current data
group_data <- read_csv("DS2_sample_metadata_w_rrna.csv") %>%
  arrange(by = UniqueID) %>%
  select(LIMS_ID, Treatment, site, Fraction, Nanotrap_type, Enrichment) %>%
  distinct() %>%
  mutate(UniqueID = paste(LIMS_ID, Treatment, Enrichment, sep=("-")),
         LIMS_ID = as.character(LIMS_ID),
         Fraction = str_replace_all(Fraction,
                                    c("filtrate" = "Fil",
                                      "retentate" = "Ret",
                                      "unfiltered" = "Unfil")),
         Nanotrap_type = str_replace_all(Nanotrap_type,
                                         c("^A$" = "NT-A",
                                           "A&B" = "NT-AB",
                                           "none" = "DirEx")),
         Enrichment = str_replace_all(Enrichment,
                                      c("None" = "Non-targeted"))
  )

#Read in RGI allele results: (from rgi_ARG_import.R)
rpip_and_unt_allele_and_info <- read_rds(
  "input/modified/rpip_and_unt_rgi_allele_and_info.rds"
  )

#Filter to only protein homolog model with >= 13 MAPQ (95% mapping accuracy)
#AND over 50 bp covered:
rpip_and_unt_allele_protein_homolog <- rpip_and_unt_allele_and_info %>%
  filter(`Reference Model Type` == "protein homolog model")# %>%
#  filter(`Average MAPQ (Completely Mapped Reads)` >= 13) %>%
#  filter(`Length Coverage (bp)` >= 50)

#Get read counts
rpip_libsizes = rpip_and_unt_allele_protein_homolog %>%
  select(UniqueID, summary.after_dedup.total_reads) %>%
  distinct()

#Remove VSP samples (irrelevant here)
rpip_allele_group_data <- group_data %>%
  filter(Enrichment %in% c("RPIP", "Non-targeted")) %>%
  left_join(rpip_libsizes)

rpip_allele_df <- rpip_and_unt_allele_protein_homolog %>%
  filter(
    `Length Coverage (bp)` >= 50.0 &
    `Average MAPQ (Completely Mapped Reads)` >= 13
    ) %>%
  mutate(
    AroTerm_Pathogen = paste(
      `ARO Term`, 
      `Resistomes & Variants: Observed Pathogen(s)`
      )
    ) %>%
  left_join(group_data)

#Write the filtered ARG results to a file for plotting:
write_rds(rpip_and_unt_allele_protein_homolog,
          "input/modified/rpip_and_unt_rgi_allele_protein_homolog.rds"
          )

###############################################
#####Some exploratory stuff (okay to skip)#####
###############################################

#Investigate which genes are detected or not in RPIP vs unenriched samples:
detect_nondetect_comparison <- rpip_and_unt_allele_protein_homolog %>%
  select(Enrichment, `Reference Sequence`, `All Mapped Reads`) %>%
  group_by(Enrichment, `Reference Sequence`) %>%
  summarize(total_reads = sum(`All Mapped Reads`)) %>%
  pivot_wider(names_from = Enrichment, values_from = total_reads, values_fill = 0)

#Number of genes detected only in RPIP-enriched samples:
nrow(filter(detect_nondetect_comparison, RPIP > 0 & None == 0))

#Number of genes detected only in unenriched samples:
nrow(filter(detect_nondetect_comparison, RPIP == 0 & None > 0))

#Number of genes detected in both:
nrow(filter(detect_nondetect_comparison, RPIP > 0 & None > 0))

###########################################
#####Continue preparing data for EdgeR#####
###########################################
#Make count matrix for EdgeR (every row is an allele, every column is a sample):
rpip_allele_readcounts <- rpip_allele_df %>%
  group_by(UniqueID, `Reference Sequence`) %>%
  summarize(nmapped = sum(`All Mapped Reads`)) %>%
  pivot_wider(names_from = UniqueID,
              values_from = nmapped,
              id_cols = `Reference Sequence`)

#Formatting function
prepare_rgi_allele_count_table_for_edgeR <- function(count_table, group) {
  
  #Rename `Reference Sequence` column to "reference" for easier EdgeR use:
  fcounts <- count_table %>%
    dplyr::rename(reference = `Reference Sequence`)
  
  fcounts <- fcounts[, c("reference", group$UniqueID)]
  
  #Sanity check (make sure sample names/order match in counts and group data):
  mismatches <- sum(colnames(fcounts[,2:ncol(fcounts)]) != group$UniqueID)
  if (mismatches != 0) {
    stop("Something is wrong; the samples in the count table don't match the \
         samples in the group data.")
  }
  return(fcounts)
}


#Get the corresponding family count data:
rpip_allele_count_table <- prepare_rgi_allele_count_table_for_edgeR(
  rpip_allele_readcounts, rpip_allele_group_data
)

#Recode missing observations as zero counts:
rpip_allele_count_table[is.na(rpip_allele_count_table)] <- 0


############################################
#####Generate a design matrix for EdgeR#####
############################################
generate_levels <- function(group_df) {
  #Use LIMS_ID to control for which sample we started with:
  LIMS_ID <- factor(group_df$LIMS_ID)
  
  #Define explanatory variables as factors:
  Fraction <- factor(group_df$Fraction,
                     levels = c("Unfil", "Ret", "Fil"))
  
  Nanotrap_type <- factor(group_df$Nanotrap_type,
                          levels = c("DirEx", "NT-A", "NT-AB"))
  
  Enrichment <- factor(group_df$Enrichment) %>%
    relevel(ref = "Non-targeted")
  
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
  
  #Append all of the Boolean treatment combinations to the model matrix:
  design <- cbind(design, treat_mat)
  
  #And there we have it, the model matrix! 
  return(design)
}

#Now run the design generator to get a model matrix:
rpip_allele_model_design <- generate_levels(rpip_allele_group_data)

#######################################################
#####Build EdgeR DGElist object and fit main model#####
#######################################################
#Get library size vector from group data:
rpip_lib.size <- rpip_allele_group_data$summary.after_dedup.total_reads

#Create a DGEList object (what EdgeR works with) from the count matrix:
rpip_allele_DGElist <- DGEList(
  counts = rpip_allele_count_table, lib.size = rpip_lib.size
)

#Find low-frequency taxa that don't give us enough information to be useful but 
#mess with the analysis:
rpip_allele_lfRemover <- filterByExpr(
  rpip_allele_DGElist,
  design = rpip_allele_model_design)

rpip_allele_DGElist_lfremoved <- rpip_allele_DGElist[
  rpip_allele_lfRemover, , keep.lib.sizes = TRUE
]

#See which genes were removed (okay to skip)
which.removed = bind_cols(rpip_allele_DGElist$genes,
                          rpip_allele_lfRemover) %>%
  left_join(detect_nondetect_comparison,
            by = c("reference" ="Reference Sequence")) %>%
  mutate(removed = !(`...2`)) %>%
  select(-`...2`)

#Fit the model:
rpip_allele_disp_lfremoved <- estimateDisp(
  y = rpip_allele_DGElist_lfremoved,
  design = rpip_allele_model_design
)

rpip_allele_fit_lfremoved <- glmQLFit(
  rpip_allele_disp_lfremoved,
  rpip_allele_model_design
)

####################################################
#####Get useful contrasts from the fitted model#####
####################################################
#Create Boolean vectors indicating which columns of the design matrix correspond
#to individual treatments:
is.RPIP <- str_detect(colnames(rpip_allele_model_design), "RPIP")
is.Untargeted <- str_detect(colnames(rpip_allele_model_design), "Non-")
is.Filtrate <- str_detect(colnames(rpip_allele_model_design), "Fil")
is.Retentate <- str_detect(colnames(rpip_allele_model_design), "Ret")
is.Unfiltered <- str_detect(colnames(rpip_allele_model_design), "Unf")
is.NanoA <- str_detect(colnames(rpip_allele_model_design), "A\\.")
is.NanoB <- str_detect(colnames(rpip_allele_model_design), "AB")
is.DirectExt <- str_detect(colnames(rpip_allele_model_design), "DirEx")

#Useful function for looking at a few genes with large FC for fast exploration:
top25 <- function(result) {
  topTags(result, n = 25, sort.by = "logFC")
}

#Function to get results for all alleles from a given contrast: 
get_contrast_table <- function(edger_glm) {
  bind_cols(edger_glm$genes, edger_glm$table)
}

#########################
###Perform each contrast#
#########################
###RPIP vs. Untargeted###
#########################
#within direct extraction from unfiltered sample:
RPIP_DE_Unf <- glmQLFTest(
  rpip_allele_fit_lfremoved,
  contrast = (is.RPIP & is.DirectExt & is.Unfiltered) - 
    (is.Untargeted & is.DirectExt & is.Unfiltered))

#within direct extraction from filtrate:
RPIP_DE_Fil <- glmQLFTest(
  rpip_allele_fit_lfremoved,
  contrast = (is.RPIP & is.DirectExt & is.Filtrate) - 
    (is.Untargeted & is.DirectExt & is.Filtrate))

#within direct extraction from retentate:
RPIP_DE_Ret <- glmQLFTest(
  rpip_allele_fit_lfremoved,
  contrast = (is.RPIP & is.DirectExt & is.Retentate) - 
    (is.Untargeted & is.DirectExt & is.Retentate))


#within Nanotrap A concentrated from unfiltered sample:
RPIP_NA_Unf <- glmQLFTest(
  rpip_allele_fit_lfremoved,
  contrast = (is.RPIP & is.NanoA & is.Unfiltered) - 
    (is.Untargeted & is.NanoA & is.Unfiltered))


#within Nanotrap A concentrated from filtrate:
RPIP_NA_Fil <- glmQLFTest(
  rpip_allele_fit_lfremoved,
  contrast = (is.RPIP & is.NanoA & is.Filtrate) - 
    (is.Untargeted & is.NanoA & is.Filtrate))


#within Nanotrap A concentrated from retentate:
RPIP_NA_Ret <- glmQLFTest(
  rpip_allele_fit_lfremoved,
  contrast = (is.RPIP & is.NanoA & is.Retentate) - 
    (is.Untargeted & is.NanoA & is.Retentate))



#within Nanotrap A+B concentrated from unfiltered sample:
RPIP_NB_Unf <- glmQLFTest(
  rpip_allele_fit_lfremoved,
  contrast = (is.RPIP & is.NanoB & is.Unfiltered) - 
    (is.Untargeted & is.NanoB & is.Unfiltered))

#within Nanotrap A+B concentrated from filtrate:
RPIP_NB_Fil <- glmQLFTest(
  rpip_allele_fit_lfremoved,
  contrast = (is.RPIP & is.NanoB & is.Filtrate) - 
    (is.Untargeted & is.NanoB & is.Filtrate))

#within Nanotrap A+B concentrated from retentate:
RPIP_NB_Ret <- glmQLFTest(
  rpip_allele_fit_lfremoved,
  contrast = (is.RPIP & is.NanoB & is.Retentate) - 
    (is.Untargeted & is.NanoB & is.Retentate))

#Join results for each category
RPIP_allele_results <- list(
  RPIP_DE_Unf = RPIP_DE_Unf,
  RPIP_DE_Fil = RPIP_DE_Fil,
  RPIP_DE_Ret = RPIP_DE_Ret,
  RPIP_NA_Unf = RPIP_NA_Unf,
  RPIP_NA_Fil = RPIP_NA_Fil,
  RPIP_NA_Ret = RPIP_NA_Ret,
  RPIP_NB_Unf = RPIP_NB_Unf,
  RPIP_NB_Fil = RPIP_NB_Fil,
  RPIP_NB_Ret = RPIP_NB_Ret
)

#Get the contrast tables for all categories
RPIP_allele_results <- lapply(RPIP_allele_results, get_contrast_table) %>%
  imap(~ mutate(.x, ID = .y)) %>%
  bind_rows() %>%
  separate(ID, into = c("Enrichment", "Nanotrap_type", "Fraction"),
           sep = "_")  %>% 
  mutate(FDR = p.adjust(PValue, method = "BH")) %>%
  mutate(Treatment = paste(Fraction, Nanotrap_type, sep = "/"))

#Get relevant CARD annotations from the RGI output:
arg_allele_annotations <- rpip_and_unt_allele_and_info %>%
  filter(`Reference Model Type` == "protein homolog model") %>%
  select(`Reference Sequence`, `ARO Term`, `AMR Gene Family`, `Drug Class`, `Resistance Mechanism`) %>%
  rename("reference" = "Reference Sequence") %>%
  distinct()

arg_allele_annotations_cat <- arg_allele_annotations %>%
  mutate(category = case_when(
    str_detect(tolower(`AMR Gene Family`), "beta-lactam") ~
      "beta-lactamase",
    str_detect(`Drug Class`, ";", negate = TRUE) ~
      str_remove(`Drug Class`, "antibiotic"),
    str_detect(`Drug Class`, ";") &
      str_detect(tolower(`Resistance Mechanism`), "efflux") ~
      "efflux pump (multivalent)",
    .default = "Other multivalent resistance genes"
  ))

#Join the annotations to the results table
RPIP_allele_results <- RPIP_allele_results %>%
  left_join(arg_allele_annotations_cat)

#Write file for later:
write_rds(RPIP_allele_results,
          "input/modified/edgeR_rgi_alleles_notgrouped_results.rds")
