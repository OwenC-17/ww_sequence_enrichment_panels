library(edgeR)
library(tidyverse)

#Family count tables:
edger_family_count_table_withRrna <- read_csv("edger_tables/edger_family_count_matrix_withRrna.csv", 
                                            show_col_types = FALSE)
edger_family_count_table_noRrna <- read_csv("edger_tables/edger_family_count_matrix_noRrna.csv",
                                             show_col_types = FALSE)

#Metadata (or group data in EdgeR lingo):
group_data_withRrna <- read_csv("edger_tables/edger_sample_metadata_withRrna.csv", show_col_types = FALSE) %>%
  arrange(by = UniqueID)

group_data_noRrna <- read_csv("edger_tables/edger_sample_metadata_noRrna.csv", show_col_types = FALSE) %>%
  arrange(by = UniqueID)

filter_group_by_k2conf <- function(group, k2conf) {
  group <- filter(group, Kraken2_confidence == k2conf)
  return(group)
}

#Formatting function
prepare_count_table_for_edgeR <- function(count_table, group) {
  #Rename the `F` column to "Family" to avoid confusion with Boolean F:
  fcounts <- count_table %>%
    dplyr::rename(Family = `F`)
  
  #Designate everything that has no family-level classification as "unclassified":
  fcounts[str_starts(fcounts$Family, "unid"), "Family"] <- "unclassified"
  
  #Sum the "unclassified" counts within each sample:
  fcounts <- fcounts %>%
    group_by(Family) %>%
    summarize(across(everything(), sum), .groups = "drop")
  
  #Filter to match group data (in case group data has already been filtered, 
  #e.g. by kraken2 confidence level). This will also order the columns in 
  #fcounts to match the group data sample order:
  fcounts <- fcounts[, c("Family", group$UniqueID)]
  
  #Sanity check (make sure sample names and order match in counts and group data):
  mismatches <- sum(colnames(fcounts[,2:ncol(fcounts)]) != group$UniqueID)
  if (mismatches != 0) {
    stop("Something is wrong; the samples in the count table don't match the \
         samples in the group data.")
  }
  return(fcounts)
}

#Get the group data with K2 confidence of 0.9:
group_data_90conf_noRrna <- filter_group_by_k2conf(group_data_noRrna, 0.9)

#Get the corresponding family count data:
edger_family_count_table_90conf_noRrna <- prepare_count_table_for_edgeR(
  edger_family_count_table_noRrna, group_data_90conf_noRrna
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
  
  Enrichment <- factor(group_df$Enrichment, levels = c("None", "RPIP", "VSP"))
  
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
    names(treat_list) <- paste(treat_tb[, 1], treat_tb[, 2], treat_tb[, 3], sep = ".")[-1]
  }
  
  treat_mat <- lapply(treat_list, as.numeric) %>% as.tibble() %>% as.matrix()
  
  #Start model with only LIMS_ID as explanatory variable:
  design <- model.matrix(~LIMS_ID)
  
  #Append all of the boolean treatment combinations to the model matrix:
  design <- cbind(design, treat_mat)

  #And there we have it, the model matrix! 
  return(design)
}


#Now run the design generator to get a model matrix:
design_90conf_noRrna <- generate_levels2(group_data_90conf_noRrna)


###Fit the model
#Create a DGEList object (what EdgeR works with) from the count matrix:
edger_family_dge_90conf_noRrna <- DGEList(
  counts = edger_family_count_table_90conf_noRrna
  )

#Find low-frequency taxa that don't give us enough information to be useful but mess 
#with the analysis:
edger_family_dge_90conf_noRrna_lfRemover <- filterByExpr(
  edger_family_dge_90conf_noRrna, 
  design = design_90conf_noRrna)

edger_family_dge_lfRemoved_90conf_noRrna <- family_90conf_no_rrna[
  edger_family_dge_90conf_noRrna_lfRemover, , keep.lib.sizes = FALSE
  ]

#Fit the model:
edger_family_disp_lfRemoved_90conf_noRrna <- estimateDisp(
  y = edger_family_dge_lfRemoved_90conf_noRrna,
  design = design_90conf_noRrna
)

edger_family_fit_lfRemoved_90conf_noRrna <- glmQLFit(
  edger_family_disp_lfRemoved_90conf_noRrna, 
  design_90conf_noRrna
)
