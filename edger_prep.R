library(tidyverse)
library(data.table)
library(edgeR)

#Main dir that contains all targeted panel results
saved_reports_path <- "imported_k2_reports/all_k2_reports_anyConf.csv"

#Locations of preformatted tax tables relative to topdir
#vsp_loc <- "vsp_panels/vsp_all_taxa_all_samples_rrna_separated.tsv"
#rpip_loc <- "rpip_panels/rpip_all_taxa_all_samples_rrna_separated.tsv"
#untargeted_loc <- "untargeted/untargeted_all_taxa_all_samples_rrna_separated.tsv"
#vsp_90_loc <- "vsp_panels/vsp_all_taxa_all_samples_rrna_separated_90conf.tsv"
#rpip_90_loc <- "rpip_panels/rpip_all_taxa_all_samples_rrna_separated_90conf.tsv"
#untargeted_90_loc <- "untargeted/untargeted_all_taxa_all_samples_rrna_separated_90conf.tsv"


load_imported <- function(tsv_location,
                   topdir = "/projects/bios_microbe/cowen20/targeted_panels/") {
  imported_df <- read_tsv(paste0(topdir, tsv_location), guess_max = Inf)
  imported_df <- imported_df %>%
    mutate(taxID = as.character(taxID),
           LIMS_ID = as_factor(LIMS_ID),
           ribosomal = as_factor(ribosomal),
           site = as_factor(site))
  return(imported_df)
}

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
           Kraken2_confidence = as_factor(sprintf("%.1f", Kraken2_confidence)))
  return(fixed)
}



all_k2_reports_anyConf <- read_csv(saved_reports_path, guess_max = Inf)
all_k2_reports_anyConf <- assert_k2_report_coltypes(all_k2_reports_anyConf)

collapse_to_tax_level <- function(imported_table, level) {
    collapsed <- imported_table %>%
    group_by(!!sym(level), SampleID, LIMS_ID, Treatment, ribosomal, 
             Kraken2_confidence) %>%
      summarise(readcount = sum(nodeOnly),
                RA = sum(RA))
    return(collapsed)
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

collapse_to_family_no_rrna <- function(imported_report, treatment_name, conf_level) {
  families <- collapse_to_tax_level(imported_report, "F")
  families <- parse_sample_treatments(families) 
  families <- families %>%
    mutate(Enrichment = treatment_name) %>%
    mutate(Kraken2_confidence = conf_level)
  families <- filter(families, ribosomal == "nonrrna")
}


collapse_to_family_keep_rrna <- function(imported_report, treatment_name, conf_level) {
  families <- collapse_to_tax_level(imported_report, "F")
  families <- parse_sample_treatments(families) 
  families <- families %>%
    mutate(Enrichment = treatment_name) %>%
    mutate(Kraken2_confidence = conf_level)
}

vsp_families_no_rrna <- collapse_to_family_no_rrna(vsp_imported, "VSP", 0)
vsp_families_w_rrna <- collapse_to_family_keep_rrna(vsp_imported, "VSP", 0)

vsp_families_90_no_rrna <- collapse_to_family_no_rrna(vsp_imported_90, "VSP", 0.9)
vsp_families_90_w_rrna <- collapse_to_family_keep_rrna(vsp_imported_90, "VSP", 0.9)

rpip_families_no_rrna <- collapse_to_family_no_rrna(rpip_imported, "RPIP", 0)
rpip_families_w_rrna <- collapse_to_family_keep_rrna(rpip_imported, "RPIP", 0)

rpip_families_90_no_rrna <- collapse_to_family_no_rrna(rpip_imported_90, "RPIP", 0.9)
rpip_families_90_w_rrna <- collapse_to_family_keep_rrna(rpip_imported_90, "RPIP", 0.9)


untargeted_families_no_rrna <- collapse_to_family_no_rrna(untargeted_imported, "None", 0)
untargeted_families_w_rrna <- collapse_to_family_keep_rrna(untargeted_imported, "None", 0)

untargeted_families_90_no_rrna <- collapse_to_family_no_rrna(untargeted_imported_90, "None", 0.9)
untargeted_families_90_w_rrna <- collapse_to_family_keep_rrna(untargeted_imported_90, "None", 0.9)

all_families <- bind_rows(vsp_families_no_rrna,
                          rpip_families_no_rrna,
                          untargeted_families_no_rrna,
                          vsp_families_90_no_rrna,
                          rpip_families_90_no_rrna,
                          untargeted_families_90_no_rrna)

all_families_w_rrna <- bind_rows(vsp_families_w_rrna,
                          rpip_families_w_rrna,
                          untargeted_families_w_rrna,
                          vsp_families_90_w_rrna,
                          rpip_families_90_w_rrna,
                          untargeted_families_90_w_rrna)

all_families$Kraken2_confidence <- as.character(all_families$Kraken2_confidence)
all_families_w_rrna$Kraken2_confidence <- as.character(all_families_w_rrna$Kraken2_confidence)

all_families_w_rrna <- all_families_w_rrna %>%
  group_by(`F`, LIMS_ID, Treatment, site, Fraction, Nanotrap_type, Enrichment, Kraken2_confidence) %>%
  summarize(readcount = sum(readcount))
  

all_families <- all_families %>%
  mutate(UniqueID = paste(LIMS_ID, Treatment, Enrichment, as.character(Kraken2_confidence), sep = "-"))

all_families_w_rrna <- all_families_w_rrna %>%
  mutate(UniqueID = paste(LIMS_ID, Treatment, Enrichment, as.character(Kraken2_confidence), sep = "-"))


sample_metadata <- all_families %>%
  ungroup() %>%
  select(UniqueID, SampleID, LIMS_ID, Treatment, site, Fraction, Nanotrap_type,
         Enrichment, Kraken2_confidence) %>%
  distinct()


sample_metadata_w_rrna <- all_families_w_rrna %>%
  ungroup() %>%
  select(UniqueID, LIMS_ID, Treatment, site, Fraction, Nanotrap_type,
         Enrichment, Kraken2_confidence) %>%
  distinct()

family_count_matrix <- all_families %>%
  pivot_wider(id_cols = !!sym("F"),
              names_from = UniqueID,
              values_from = readcount,
              values_fill = 0)


family_count_matrix_w_rrna <- all_families_w_rrna %>%
  pivot_wider(id_cols = !!sym("F"),
              names_from = UniqueID,
              values_from = readcount,
              values_fill = 0)


write_csv(sample_metadata, "/projects/bios_microbe/cowen20/rprojects/targeted_panels/DS2_sample_metadata.csv")
write_csv(sample_metadata_w_rrna, "/projects/bios_microbe/cowen20/rprojects/targeted_panels/DS2_sample_metadata_w_rrna.csv")

write_csv(family_count_matrix, "/projects/bios_microbe/cowen20/rprojects/targeted_panels/DS2_family_count_matrix_no_rrna.csv")
write_csv(family_count_matrix_w_rrna, "/projects/bios_microbe/cowen20/rprojects/targeted_panels/DS2_family_count_matrix_w_rrna.csv")
