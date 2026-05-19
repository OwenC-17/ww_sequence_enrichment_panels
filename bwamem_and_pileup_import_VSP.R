#Import pileup results based on BWA-MEM read mapping to VSP reference genomes

library(tidyverse)
library(taxonomizr)


############################
#####Set data locations#####
############################

#Taxonomic information for taxonomizr:
sqlPath = "/projects/bios_microbe/cowen20/ref_db/taxonomizr/accessionTaxa.sql"

#Pileup files based on VSP reads vs. VSP references:
vsp_dir <- paste0("/projects/bios_microbe/cowen20/targeted_panels/vsp_panels/",
    "raw_fastqs/fastp_out_no_dedup/deduped/bwa_out/",
    "mapped_and_low_complexity_removed/sorted/filtered095/pileup13mq/")

#Pileup files based on untargeted reads vs. VSP references:
unt_vsp_dir <- paste0("/projects/bios_microbe/cowen20/targeted_panels/",
    "untargeted/raw_fastqs/fastp_out_no_dedup/deduped/bwa_out/",
    "mapped_and_low_complexity_removed/sorted/filtered095/pileup13mq/")


#List all the pileup files:
vsp_count_list <- list.files(
  vsp_dir, pattern = "*.tsv", full.names = TRUE
  )
unt_vsp_count_list <- list.files(
  unt_vsp_dir, pattern = "*.tsv", full.names = TRUE
  )


#################################
#####Import the pileup files#####
#################################
#From lists of files, generate lists of DFs:
vsp_count_dfs <- lapply(vsp_count_list, function(x) {
  read_tsv(x, 
           col_names = c("Reference",
                         "rlength",
                         "nmapped",
                         "assigned_unmapped"))
  })

unt_vsp_count_dfs <- lapply(unt_vsp_count_list, function(x) {
  read_tsv(x,
           col_names = c("Reference",
                         "rlength",
                         "nmapped",
                         "assigned_unmapped"))
  })

#Name the DF lists based on their file of origin (which contains sample IDs)
names(vsp_count_dfs) <- lapply(vsp_count_list, basename)
names(unt_vsp_count_dfs) <- lapply(unt_vsp_count_list, basename)

#Merge each list of DFs into one large DF with a column called "filename"
#containing the file of origin:
vsp_counts_bigdf <- bind_rows(vsp_count_dfs, .id = "filename")
unt_vsp_counts_bigdf <- bind_rows(unt_vsp_count_dfs, .id = "filename")

#Label the enrichment used for each big DF:
vsp_counts_bigdf$Enrichment <- "VSP"
unt_vsp_counts_bigdf$Enrichment <- "Non-targeted"

#Combine VSP and untargeted DFs into one DF:
vsp_and_unt_counts_bigdf <- bind_rows(vsp_counts_bigdf,
                                      unt_vsp_counts_bigdf) %>%
  #Each sample has a contains reference genome with 0 mapped reads that is
  #labeled as an asterisk. Currently this command treats the asterisk as plain
  #text rather than a regular expression so it works:
  filter(Reference != "*")

#Extract sample_id (LIMS_ID and treatment code) from filename,
#and extract reference genome accession from Reference column:
vsp_and_unt_counts_bigdf <- vsp_and_unt_counts_bigdf %>%
  mutate(sample_id = substr(filename, 5, 12)) %>%
  mutate(rname = substr(Reference, 1, 11)) %>%
  mutate(UniqueID = paste(sample_id, Enrichment, sep = "-"))


###################################
#####Add taxonomic information#####
###################################

#Function to convert taxids to full  classification
add_taxonomy <- function(coverage_table) {
  unique_accessions <- unique(coverage_table$rname)
  
  taxonomyDf <- data.frame(rname = unique_accessions)
  taxonomyDf$taxid <- accessionToTaxa(taxonomyDf$rname, sqlPath)
  taxonomyDf$taxid <- as.character(taxonomyDf$taxid)
  
  uniqueTaxId <- unique(taxonomyDf$taxid)
  
  rawTaxa <- getRawTaxonomy(uniqueTaxId, sqlPath)
  rawTaxaDf <- bind_rows(rawTaxa, .id = "taxid")
  rawTaxaDf$taxid <- str_remove_all(rawTaxaDf$taxid, "\\s")
  weird_accessions_annotated <- read_csv("/projects/bios_microbe/cowen20/rprojects/targeted_panels/input/weird_accessions_fixed.csv")
  weird_accessions_annotated <- rename(weird_accessions_annotated, rname = Accession)
  
  accessions2taxids <- taxonomyDf %>%
    filter(!is.na(taxid))
  
  rawTaxaDf <- rawTaxaDf %>%
    full_join(accessions2taxids, by = "taxid")
  
  rawTaxaDf <- bind_rows(rawTaxaDf, weird_accessions_annotated)
  
  annotated_table <- coverage_table %>%
    left_join(rawTaxaDf, by = c("rname"))
  
  return(annotated_table)
  
}

#Get taxonomy of reference IDs
vspunt_taxa <- data.frame(rname = unique(vsp_and_unt_counts_bigdf$rname))
vspunt_taxtable <- add_taxonomy(vspunt_taxa)
vspunt_taxtable[is.na(vspunt_taxtable$domain),]$domain <- "Viruses"

vsp_and_unt_counts_bigdf_tax <- vsp_and_unt_counts_bigdf %>%
  left_join(vspunt_taxtable)


##############################
#####Import basewise data#####
##############################
#Read file names:
vsp_pileup_list <- list.files(vsp_dir, pattern = "*.tab", full.names = TRUE)
unt_vsp_pileup_list <- list.files(
  unt_vsp_dir, pattern = "*.tab", full.names = TRUE
  )

#Load tables into list of dataframes:
vsp_pileup_dfs <- lapply(vsp_pileup_list, function(x) {
  read_tsv(x,
           col_types = "cdccdcdddddddddddd",
           na = c(".", "#", ""))
  })


unt_vsp_pileup_dfs <- lapply(unt_vsp_pileup_list, function(x) {
  read_tsv(x,
           col_types = "cdccdcdddddddddddd",
           na = c(".", "#", ""))
  })

#Name DFs in the lists after their file of origin:
names(vsp_pileup_dfs) <- lapply(vsp_pileup_list, basename)
names(unt_vsp_pileup_dfs) <- lapply(unt_vsp_pileup_list, basename)

#Bind the DF lists into single big DFs:
vsp_big_pileup_df <- bind_rows(vsp_pileup_dfs, .id = "filename")
unt_vsp_big_pileup_df <- bind_rows(unt_vsp_pileup_dfs, .id = "filename")

#Add Enrichment label:
vsp_big_pileup_df$Enrichment <- "VSP"
unt_vsp_big_pileup_df$Enrichment <- "Non-targeted"


#Combine VSP and untargeted DFs into one DF:
vspANDunt_big_pileup_df <- bind_rows(vsp_big_pileup_df,
                                     unt_vsp_big_pileup_df) %>%
  mutate(sample_id = substr(filename, 5, 12)) %>%
  mutate(UniqueID = paste(sample_id, Enrichment, sep = "-")) %>%
  mutate(rname = substr(Reference, 1, 11)) %>%
  left_join(vspunt_taxtable)

#Save imported file:
write_rds(vspANDunt_big_pileup_df, "input/modified/vspANDunt_big_pileup_df.rds")


###############################
#####Combine imported data#####
###############################

#Combine counts, coverage and mutation data
vsp_and_unt_counts_and_covstats <- vsp_and_unt_counts_bigdf_tax %>%
  left_join(vspunt_covered_bases) %>%
  left_join(vspunt_mutations_count) %>%
  #The only taxon with no family ID assigned is NC_116928.1, which has no mapped
  #reads anyway, so easiest to filter it out:
  filter(!is.na(family))

#Set NAs to 0
vsp_and_unt_counts_and_covstats[
  is.na(vsp_and_unt_counts_and_covstats$bases_covered), ]$bases_covered <- 0
vsp_and_unt_counts_and_covstats[
  is.na(vsp_and_unt_counts_and_covstats$mutation_count), ]$mutation_count <- 0

fastp_info_table <- read_rds(
  "input/modified/all_fastp_reports_dupdedup.rds"
  ) %>%
  mutate(UniqueID = paste(LIMS_ID, Treatment, Enrichment, sep = "-"))


vsp_and_unt_counts_covstats_and_info <- vsp_and_unt_counts_and_covstats %>%
  left_join(fastp_info_table, by = c("UniqueID", "Enrichment")) %>%
  mutate(RA_nmapped = nmapped / summary.after_filtering.total_reads,
         RA_bases_covered = bases_covered / summary.after_filtering.total_bases)

write_rds(vsp_and_unt_counts_covstats_and_info,
          "input/modified/vsp_and_unt_counts_covstats_and_info.rds")
