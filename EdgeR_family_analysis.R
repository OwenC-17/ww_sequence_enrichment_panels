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

