library(tidyverse)
library(zoo)

#Main dir that contains all targeted panel results
topdir <- "input/link_to_raw_data/"

#Subdir with k2 reports from each panel
vsp_dir_00conf <- paste0(topdir, "vsp_panels/raw_fastqs/fastp_out_no_dedup/",
                  "kraken2_out/k2_nt_20240530/rrna_labeled/")
rpip_dir_00conf <- paste0(topdir, "rpip_panels/raw_fastqs/fastp_out_no_dedup/",
                   "kraken2_out/k2_nt_20240530/rrna_labeled/")
unt_dir_00conf <- paste0(topdir, "untargeted/raw_fastqs/fastp_out_no_dedup/",
                  "kraken2_out/k2_nt_20240530/rrna_labeled/")

vsp_dir_90conf <- paste0(topdir, "vsp_panels/raw_fastqs/fastp_out_no_dedup/",
                  "kraken2_out/k2_nt_20240530/90conf/rrna_labeled/")
rpip_dir_90conf <- paste0(topdir, "rpip_panels/raw_fastqs/fastp_out_no_dedup/",
                   "kraken2_out/k2_nt_20240530/90conf/rrna_labeled/")
unt_dir_90conf <- paste0(topdir, "untargeted/raw_fastqs/fastp_out_no_dedup/",
                  "kraken2_out/k2_nt_20240530/90conf/rrna_labeled/")

generate_tax_table = function(tax_report_tsv) {
  ###create a table from a kraken2 taxonomy report whose location is specified
  ###at tax_report_tsv
  
  #Specify columns in k2 report format
  colnames = c("percentOfReads", "nodeAndChildren", "nodeOnly", "taxLevel", 
               "taxID", "name")
  coltypes = "diicff"
  
  #Import the k2 report
  tax_report = read_tsv(tax_report_tsv, col_names = colnames, 
                        col_types = coltypes)
  
  #This is needed to maintain the correct order of taxonomic levels when 
  #widening the k2 report data (also removes anything deeper than 9 nested
  #sublevels)
  generate_tax_order = function(prefixes, max_depth = 9) {
    expanded <- lapply(prefixes, 
                       function(p) c(p, paste0(p, seq_len(max_depth)))) %>%
      unlist()
    return(expanded)
  }
  
  #Get the order for all desired main levels:
  all_tax_levels = generate_tax_order(c("U", "R", "D", "P", "C", "O", "F", "G",
                                        "S"))
  
  #Filter the ordered vector to only include levels and sublevels present in the
  #k2 report
  selected_tax_levels = intersect(all_tax_levels, tax_report$taxLevel)
  
  #Widen the k2 report so that the represent each taxonomic level in order:
  tax_table = pivot_wider(tax_report, names_from = taxLevel,
                          values_from = name) %>%
    select("percentOfReads", "nodeAndChildren", "nodeOnly", "taxID", 
           selected_tax_levels)
  
  return(tax_table)
}



parse_sample_ids <- function(taxtable) {
  ####Split content of SampleID column into more useful info
  separate(taxtable,
           col = SampleID,
           into = c("QCSeqID", "LIMS_ID", "Treatment",
                    NA, NA, NA, NA, NA, 
                    "ribosomal", NA ), 
           sep="(-|_)")
}

parse_locations <- function(taxtable) {
  ###Add sample site names based on numeric IDs
  mutate(taxtable, site = str_replace_all(
    LIMS_ID,
    c(
      "(36397|38156|41596|34970)" = "OBrien",
      "(32976|32989|20040|22015)" = "OHare"
    )))
}

fill_tax_NAs = function(tax_table) {
  ###Fill in blank spaces with lowest applicable taxa IDs.
  
  #Pull out columns containing taxa labels and assert that they are strings:
  taxa = tax_table %>%
    select("U":last_col()) %>%
    mutate(across(everything(), as.character))
  
  
  col_names <- colnames(taxa)
  
  fill_missing <- function(row) {
    for (i in seq_len(length(row) - 1)) {
      #Identify most specific ID available (first missing level will be NA):
      if (!is.na(row[i]) && is.na(row[i + 1])) {
        #Fill all cells to the right of the most specific ID, indicating that
        #the taxon is an unidentified member of that ID:
        row[(i + 1):length(row)] <- paste0("unidentified (", col_names[i], ": ",
                                           row[i], ")")
        break #There will be only one most specific ID per row
      }
    }
    return(row)
  }
  
  #Apply fill_missing to all rows in the taxa df
  taxa <- t(apply(taxa, 1, function(row) fill_missing(row)))
  
  #Top row is unclassified (at all levels)
  taxa[1, ] <- "unclassified"
  
  #Fill in IDs (downward) 
  taxa <- na.locf(taxa)
  
  #Reconnect to the other data from the input tax_table
  out <- bind_cols(select(tax_table, "percentOfReads":"taxID"), taxa)
  
  #Calculate relative abundances
  out <- mutate(out, RA = nodeOnly / sum(nodeOnly))
  
  return(out)
}

import_kraken2_summary <- function(directory, pattern) {
  file_list <- list.files(path = directory, pattern = pattern)
  placeholder_list <- vector("list", length = length(file_list))
  for (k2t in seq_len(length(file_list))) {
    report <- generate_tax_table(paste0(directory, file_list[k2t])) %>%
      fill_tax_NAs() %>%
      mutate(SampleID = file_list[k2t])
    placeholder_list[[k2t]] <- report
    print(file_list[k2t])
  }
  imported <- bind_rows(placeholder_list) %>%
    parse_sample_ids() %>%
    parse_locations()
  return(imported)
}


#Import k2 reports with k2 confidence of 0.0 (default)
vsp_k2_reports_00conf <- import_kraken2_summary(vsp_dir_00conf, "report.tsv")
rpip_k2_reports_00conf <- import_kraken2_summary(rpip_dir_00conf, "report.tsv")
unt_k2_reports_00conf <- import_kraken2_summary(unt_dir_00conf, "report.tsv")

#Import k2 reports with k2 confidence of 0.9
vsp_k2_reports_90conf <- import_kraken2_summary(vsp_dir_90conf, "report.tsv")
rpip_k2_reports_90conf <- import_kraken2_summary(rpip_dir_90conf, "report.tsv")
unt_k2_reports_90conf <- import_kraken2_summary(unt_dir_90conf, "report.tsv")


vsp_k2_reports_00conf <- vsp_k2_reports_00conf %>% 
  mutate(Enrichment = "VSP",
         Kraken2_confidence = "0.0" )
rpip_k2_reports_00conf <- rpip_k2_reports_00conf %>% 
  mutate(Enrichment = "RPIP",
         Kraken2_confidence = "0.0" )
unt_k2_reports_00conf <- unt_k2_reports_00conf %>% 
  mutate(Enrichment = "None",
         Kraken2_confidence = "0.0" )



vsp_k2_reports_90conf <- vsp_k2_reports_90conf %>% 
  mutate(Enrichment = "VSP",
         Kraken2_confidence = "0.9" )
rpip_k2_reports_90conf <- rpip_k2_reports_90conf %>% 
  mutate(Enrichment = "RPIP",
         Kraken2_confidence = "0.9" )
unt_k2_reports_90conf <- unt_k2_reports_90conf %>% 
  mutate(Enrichment = "None",
         Kraken2_confidence = "0.9" )


all_k2_reports_anyConf <- bind_rows(vsp_k2_reports_00conf,
                                    rpip_k2_reports_00conf,
                                    unt_k2_reports_00conf,
                                    vsp_k2_reports_90conf,
                                    rpip_k2_reports_90conf,
                                    unt_k2_reports_90conf)

rm(vsp_k2_reports_00conf, rpip_k2_reports_00conf, unt_k2_reports_00conf,
   vsp_k2_reports_90conf, rpip_k2_reports_90conf, unt_k2_reports_90conf)

dir.create("imported_k2_reports")
write_csv(all_k2_reports_anyConf, "imported_k2_reports/all_k2_reports_anyConf.csv")