library(tidyverse)
library(edgeR)

rpip_and_unt_counts_covstats_and_info_by_taxid <- read_rds(
  "input/modified/rpip_and_unt_counts_covstats_and_info_by_taxid.rds"
  )

#Load the explanatory variables (group data) from Kraken analysis
#and modify for current data
group_data <- read_csv("DS2_sample_metadata_w_rrna.csv") %>%
  arrange(by = UniqueID) %>%
  select(LIMS_ID, Treatment, site, Fraction, Nanotrap_type, Enrichment) %>%
  distinct() %>%
  mutate(UniqueID = paste(LIMS_ID, Treatment, Enrichment, sep=("-")),
         LIMS_ID = as.character(LIMS_ID))

#Get library sizes
rpip_libsizes = rpip_and_unt_counts_covstats_and_info_by_taxid %>%
  mutate(UniqueID = paste(LIMS_ID, Treatment, Enrichment, sep = "-")) %>%
  select(UniqueID, summary.after_dedup.total_reads) %>%
  distinct()

#Filter group data to RPIP and unenriched samples
rpip_group_data <- group_data %>%
  filter(Enrichment %in% c("RPIP", "None")) %>%
  left_join(rpip_libsizes)

#Get library size vector for EdgeR (it will be in the same order as group data)
rpip_lib.size <- rpip_group_data$summary.after_dedup.total_reads

#Get read counts by species
rpip_species_readcounts <- rpip_and_unt_counts_covstats_and_info_by_taxid %>%
  group_by(UniqueID, species) %>%
  summarize(nmapped = sum(nmapped)) %>%
  pivot_wider(names_from = UniqueID, values_from = nmapped, id_cols = species)


#Formatting function
prepare_count_table_for_edgeR <- function(count_table, group) {
  #Rename the `species` column to "Species":
  fcounts <- count_table %>%
    dplyr::rename(Species = species)
  
  #Filter to match group data (in case group data has already been filtered, 
  #e.g. by kraken2 confidence level). This will also order the columns in 
  #fcounts to match the group data sample order:
  fcounts <- fcounts[, c("Species", group$UniqueID)]
  
  
  
  #Sanity check (make sure sample names and order match in counts and group data):
  mismatches <- sum(colnames(fcounts[,2:ncol(fcounts)]) != group$UniqueID)
  if (mismatches != 0) {
    stop("Something is wrong; the samples in the count table don't match the \
         samples in the group data.")
  }
  return(fcounts)
}


#Get the corresponding species count data:
rpip_species_count_table <- prepare_count_table_for_edgeR(
  rpip_species_readcounts, rpip_group_data
)


###Generate a design matrix for EdgeR:
generate_levels <- function(group_df) {
  #Use LIMS_ID to control for which sample we started with:
  LIMS_ID <- factor(group_df$LIMS_ID)
  
  #Define explanatory variables as factors and order them to have an appropriate
  #baseline:  
  Fraction <- factor(group_df$Fraction, levels = c("unfiltered", 
                                                   "retentate", 
                                                   "filtrate"))
  
  Nanotrap_type <- factor(group_df$Nanotrap_type, levels = c("none", "A", "A&B"))
  
  Enrichment <- factor(group_df$Enrichment)
  Enrichment <- relevel(Enrichment, ref = "None")

    
  #Create a tibble containing all combinations of treatment variables:
  treat_tb <- expand.grid(Fraction = levels(Fraction),
                          Nanotrap_type = levels(Nanotrap_type),
                          Enrichment = levels(Enrichment),
                          stringsAsFactors = FALSE)
  
  #Create an empty list, then add every possible combination of treatment variables:
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
rpip_model_design <- generate_levels(rpip_group_data)

###Fit the model
#Create a DGEList object (what EdgeR works with) from the count matrix:
rpip_DGElist <- DGEList(
  counts = rpip_species_count_table, lib.size = rpip_lib.size
)

#Save this object for the diversity analyses
write_rds(rpip_DGElist, "input/modified/edgeR_bwamem_rpip_DGE.rds")

#Find low-frequency taxa that don't give us enough information to be useful but 
#mess with the analysis:
rpip_lfRemover <- filterByExpr(
  rpip_DGElist,
  design = rpip_model_design)

rpip_DGElist_lfremoved <- rpip_DGElist[
   rpip_lfRemover, , keep.lib.sizes = TRUE
]

#Fit the model:
rpip_disp_lfremoved <- estimateDisp(
  y = rpip_DGElist_lfremoved,
  design = rpip_model_design
)

rpip_fit_lfremoved <- glmQLFit(
  rpip_disp_lfremoved,
  rpip_model_design
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

index_main_effects(rpip_model_design)

#Useful for looking at the taxa with the largest fold change:
top25 <- function(result) {
  topTags(result, n = 25, sort.by = "logFC")
}

#Function to get EdgeR results for all taxa from a specific contrast:
get_contrast_table <- function(edger_glm) {
  bind_cols(edger_glm$genes, edger_glm$table)
}


#####Perform the contrasts#####

###Direct Extraction
#Unfiltered
RPIP_DE_Unf <- glmQLFTest(
  rpip_fit_lfremoved,
  contrast = (is.RPIP & is.DirectExt & is.Unfiltered) - 
    (is.Untargeted & is.DirectExt & is.Unfiltered))

#Filtrate
RPIP_DE_Fil <- glmQLFTest(
  rpip_fit_lfremoved,
  contrast = (is.RPIP & is.DirectExt & is.Filtrate) - 
    (is.Untargeted & is.DirectExt & is.Filtrate))

#Retentate
RPIP_DE_Ret <- glmQLFTest(
  rpip_fit_lfremoved,
  contrast = (is.RPIP & is.DirectExt & is.Retentate) - 
    (is.Untargeted & is.DirectExt & is.Retentate))


###RPIP, Nanotrap A
#Unfiltered
RPIP_NA_Unf <- glmQLFTest(
  rpip_fit_lfremoved,
  contrast = (is.RPIP & is.NanoA & is.Unfiltered) - 
    (is.Untargeted & is.NanoA & is.Unfiltered))

#Filtrate
RPIP_NA_Fil <- glmQLFTest(
  rpip_fit_lfremoved,
  contrast = (is.RPIP & is.NanoA & is.Filtrate) - 
    (is.Untargeted & is.NanoA & is.Filtrate))

#Retentate
RPIP_NA_Ret <- glmQLFTest(
  rpip_fit_lfremoved,
  contrast = (is.RPIP & is.NanoA & is.Retentate) - 
    (is.Untargeted & is.NanoA & is.Retentate))


###RPIP, Nanotrap A and B
#Unfiltered
RPIP_NB_Unf <- glmQLFTest(
  rpip_fit_lfremoved,
  contrast = (is.RPIP & is.NanoB & is.Unfiltered) - 
    (is.Untargeted & is.NanoB & is.Unfiltered))

#Filtrate
RPIP_NB_Fil <- glmQLFTest(
  rpip_fit_lfremoved,
  contrast = (is.RPIP & is.NanoB & is.Filtrate) - 
    (is.Untargeted & is.NanoB & is.Filtrate))

#Retentate
RPIP_NB_Ret <- glmQLFTest(
  rpip_fit_lfremoved,
  contrast = (is.RPIP & is.NanoB & is.Retentate) - 
    (is.Untargeted & is.NanoB & is.Retentate))


#Make a list of the results from all contrasts (they are DGELRT objects):
RPIP_results <- list(
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

#Get the results tables and join them into a data frame:
RPIP_results <- lapply(RPIP_results, get_contrast_table) %>%
  imap(~ mutate(.x, ID = .y)) %>%
  bind_rows() %>%
  separate(ID, into = c("Enrichment", "Nanotrap_type", "Fraction"),
           sep = "_")  %>% 
  mutate(FDR = p.adjust(PValue, method = "BH")) %>%
  mutate(Treatment = paste(Fraction, Nanotrap_type, sep = "/"))


#Read in the full taxonomy table:
rpipunt_taxtable <- read_rds("input/modified/rpipunt_taxtable.rds")

#Just use standard classification scheme columns:
tax_to_join <- rpipunt_taxtable %>%
  select(domain, kingdom, phylum, class, order, family, genus, species) %>%
  distinct()

#Join to the EdgeR results: 
RPIP_results_with_tax <- RPIP_results %>%
  left_join(tax_to_join, by = c("Species" = "species"))

#Save ths results table:
write_rds(RPIP_results_with_tax, "input/modified/rpip_bwamem_edger_results.rds")
