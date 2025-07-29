library(tidyverse)
library(data.table)
library(edgeR)

#Main dir that contains all targeted panel results
saved_reports_path <- "imported_k2_reports/all_k2_reports_anyConf.csv"

#Function to apply after importing (there are too many columns to specify types
#during import)
assert_k2_report_coltypes <- function(k2_report_df) {
  fixed <- k2_report_df %>%
    mutate(taxID = as_factor(as.character(taxID)),
           nodeAndChildren = as.integer(nodeAndChildren),
           nodeOnly = as.integer(nodeOnly),
           LIMS_ID = as_factor(LIMS_ID),
           Treatment = as_factor(Treatment),
           ribosomal = as_factor(ribosomal),
           site = as_factor(site),
           Enrichment = as_factor(Enrichment),
           #sprintf is used to preserve correct digit length in type conversion
           Kraken2_confidence = as_factor(sprintf("%.1f", Kraken2_confidence)))
  return(fixed)
}

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

#######################
#####Load the data#####
#######################
all_k2_reports_anyConf <- read_csv(saved_reports_path, guess_max = Inf)
all_k2_reports_anyConf <- assert_k2_report_coltypes(all_k2_reports_anyConf)
all_k2_reports_anyConf <- parse_sample_treatments(all_k2_reports_anyConf)
#######################

#To summarize the imported reports at any particular tax level
collapse_to_tax_level <- function(imported_k2_reports, level) {
  collapsed <- imported_k2_reports
  setDT(collapsed)
  collapsed[, .(readcount = sum(nodeOnly), RA = sum(RA)),
            by = c(level, "LIMS_ID", "site", "Treatment", 
                   "Enrichment", "ribosomal", "Kraken2_confidence", "Fraction",
                   "Nanotrap_type")]
}

collapse_to_family_no_rrna <- function(imported_k2_report) {
  families <- collapse_to_tax_level(imported_k2_report, "F")
  families <- filter(families, ribosomal == "nonrrna")
}

allFamilies_anyConf_noRrna <- collapse_to_family_no_rrna(all_k2_reports_anyConf)

collapse_to_family_keep_rrna <- function(imported_k2_report) {
  families <- collapse_to_tax_level(imported_k2_report, "F")
  families <- families %>%
    group_by(across(-c(readcount, RA, ribosomal))) %>%
    summarise(readcount = sum(readcount), .groups = "drop") %>%
    group_by(across(-c(`F`, readcount))) %>%
    mutate(RA = readcount / sum(readcount))
}

allFamilies_anyConf_withRrna <- collapse_to_family_keep_rrna(
  all_k2_reports_anyConf
  )


################################################################################
allFamilies_anyConf_noRrna <- allFamilies_anyConf_noRrna %>%
  mutate(UniqueID = paste(LIMS_ID, Treatment, Enrichment, Kraken2_confidence, sep = "-"))

allFamilies_anyConf_withRrna <- allFamilies_anyConf_withRrna %>%
  mutate(UniqueID = paste(LIMS_ID, Treatment, Enrichment, Kraken2_confidence, sep = "-"))


sample_metadata_noRrna <- allFamilies_anyConf_noRrna %>%
  ungroup() %>%
  select(UniqueID, LIMS_ID, Treatment, site, Fraction, Nanotrap_type,
         Enrichment, Kraken2_confidence) %>%
  distinct()


sample_metadata_withRrna <- allFamilies_anyConf_withRrna %>%
  ungroup() %>%
  select(UniqueID, LIMS_ID, Treatment, site, Fraction, Nanotrap_type,
         Enrichment, Kraken2_confidence) %>%
  distinct()

family_count_matrix_noRrna <- allFamilies_anyConf_noRrna %>%
  pivot_wider(id_cols = !!sym("F"),
              names_from = UniqueID,
              values_from = readcount,
              values_fill = 0)


family_count_matrix_withRrna <- allFamilies_anyConf_withRrna %>%
  pivot_wider(id_cols = !!sym("F"),
              names_from = UniqueID,
              values_from = readcount,
              values_fill = 0)


write_csv(sample_metadata_noRrna, "edger_tables/edger_sample_metadata_noRrna.csv")
write_csv(sample_metadata_withRrna, "edger_tables/edger_sample_metadata_withRrna.csv")

write_csv(family_count_matrix_noRrna, "edger_tables/edger_family_count_matrix_noRrna.csv")
write_csv(family_count_matrix_withRrna, "edger_tables/edger_family_count_matrix_withRrna.csv")