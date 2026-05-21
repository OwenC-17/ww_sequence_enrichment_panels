#Import pileup results based on BWA-MEM read mapping to RPIP reference genomes
#IMPORTANT: Make sure at least available RAM is >400 GB or R will crash

library(tidyverse)
library(data.table)
library(taxonomizr)


############################
#####Set data locations#####
############################

#Taxonomic information for taxonomizr:
sqlPath = "/projects/bios_microbe/cowen20/ref_db/taxonomizr/accessionTaxa.sql"

#Pileup files based on RPIP reads vs. RPIP references:
rpip_dir <- paste0("/projects/bios_microbe/cowen20/targeted_panels/",
                   "rpip_panels/raw_fastqs/fastp_out_no_dedup/deduped/",
                   "bwa_out/mapped_and_low_complexity_removed/sorted/",
                   "deduped_bams/filtered095/pileup13mq/")

#Pileup files based on untargeted reads vs. RPIP references:
unt_rpip_dir <- paste0("/projects/bios_microbe/cowen20/targeted_panels/",
                   "untargeted/raw_fastqs/fastp_out_no_dedup/deduped/bwa_out/",
                   "mapped_and_low_complexity_removed/sorted/deduped_bams/",
                   "filtered095/pileup13mq/")

#List all the pileup files:
rpip_count_list <- list.files(
  rpip_dir, pattern = "*.tsv", full.names = TRUE
)
unt_rpip_count_list <- list.files(
  unt_rpip_dir, pattern = "*.tsv", full.names = TRUE
)

#################################
#####Import the pileup files#####
#################################
#From lists of files, generate lists of DFs:
rpip_count_dfs <- lapply(rpip_count_list, function(x) {
  read_tsv(x, 
           col_names = c("Reference",
                         "rlength",
                         "nmapped",
                         "assigned_unmapped"))
})

unt_rpip_count_dfs <- lapply(unt_rpip_count_list, function(x) {
  read_tsv(x,
           col_names = c("Reference",
                         "rlength",
                         "nmapped",
                         "assigned_unmapped"))
})

#Name the DF lists based on their file of origin (which contains sample IDs)
names(rpip_count_dfs) <- lapply(rpip_count_list, basename)
names(unt_rpip_count_dfs) <- lapply(unt_rpip_count_list, basename)

#Merge each list of DFs into one large DF with a column called "filename"
#containing the file of origin:
rpip_counts_bigdf <- bind_rows(rpip_count_dfs, .id = "filename")
unt_rpip_counts_bigdf <- bind_rows(unt_rpip_count_dfs, .id = "filename")


#Label the enrichment used for each big DF:
rpip_counts_bigdf$Enrichment <- "RPIP"
unt_rpip_counts_bigdf$Enrichment <- "Non-targeted"

#Combine VSP and untargeted DFs into one DF:
rpip_and_unt_counts_bigdf <- bind_rows(rpip_counts_bigdf,
                                      unt_rpip_counts_bigdf) %>%
  #Each sample has a contains reference genome with 0 mapped reads that is
  #labeled as an asterisk. Currently this command treats the asterisk as plain
  #text rather than a regular expression so it works:
  filter(Reference != "*")

#Extract sample_id (LIMS_ID and treatment code) from filename,
#and extract reference genome accession from Reference column:
rpip_and_unt_counts_bigdf <- rpip_and_unt_counts_bigdf %>%
  mutate(sample_id = substr(filename, 5, 12)) %>%
  mutate(rname = Reference) %>%
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
  weird_accessions_annotated <- read_csv(paste0("/projects/bios_microbe/",
                                                "cowen20/rprojects/targeted_panels/input/weird_accessions_fixed.csv"))
  weird_accessions_annotated <- rename(weird_accessions_annotated,
                                       rname = Accession)
  
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
rpipunt_taxa <- data.frame(rname = unique(rpip_and_unt_counts_bigdf$rname))
rpipunt_taxtable <- add_taxonomy(rpipunt_taxa)

#Viruses won't have a domain assigned automatically:
rpipunt_taxtable[is.na(rpipunt_taxtable$domain),]$domain <- "Viruses"

#[Eubacterium] is not a stable genus name; it is the only bacterial genus in 
#the reference genome set that isn't assigned automatically: 
rpipunt_taxtable[rpipunt_taxtable$domain == "Bacteria" & 
                   is.na(rpipunt_taxtable$genus), ]$genus <- "[Eubacterium]"


rpip_and_unt_counts_bigdf_tax <- rpip_and_unt_counts_bigdf %>%
  left_join(rpipunt_taxtable)


##############################
#####Import basewise data#####
##############################
#Read file names:
rpip_pileup_list <- list.files(
  rpip_dir, pattern = "*.tab", full.names = TRUE
  )
unt_rpip_pileup_list <- list.files(
  unt_rpip_dir, pattern = "*.tab", full.names = TRUE
)

#Load tables into list of dataframes:
rpip_pileup_dfs <- lapply(rpip_pileup_list, function(x) {
  read_tsv(x,
           col_types = "cdccdcdddddddddddd",
           na = c(".", "#", ""))
})

unt_rpip_pileup_dfs <- lapply(unt_rpip_pileup_list, function(x) {
  read_tsv(x,
           col_types = "cdccdcdddddddddddd",
           na = c(".", "#", ""))
})

#Name DFs in the lists after their file of origin:
names(rpip_pileup_dfs) <- lapply(rpip_pileup_list, basename)
names(unt_rpip_pileup_dfs) <- lapply(unt_rpip_pileup_list, basename)

#Bind the DF lists into single big DFs:
rpip_big_pileup_df <- bind_rows(rpip_pileup_dfs, .id = "filename")
unt_rpip_big_pileup_df <- bind_rows(unt_rpip_pileup_dfs, .id = "filename")

#Save some memory:
rm(rpip_pileup_dfs, unt_rpip_pileup_dfs)
gc()

#Add Enrichment label:
rpip_big_pileup_df$Enrichment <- "RPIP"
unt_rpip_big_pileup_df$Enrichment <- "Non-targeted"


#Combine VSP and untargeted DFs into one DF:
rpipANDunt_big_pileup_df <- bind_rows(rpip_big_pileup_df,
                                      unt_rpip_big_pileup_df) %>%
  mutate(sample_id = substr(filename, 5, 12)) %>%
  mutate(UniqueID = paste(sample_id, Enrichment, sep = "-")) %>%
  mutate(rname = Reference) %>%
  left_join(rpipunt_taxtable)

#Save more memory:
rm(rpip_big_pileup_df, unt_rpip_big_pileup_df)
gc()

#Save imported file:
write_rds(rpipANDunt_big_pileup_df,
          "input/modified/rpipANDunt_big_pileup_df.rds")
gc()

rpipunt_big_dt <- as.data.table(rpipANDunt_big_pileup_df)
write_rds(rpipunt_big_dt, "input/modified/rpipANDunt_big_pileup_dt.rds")

#Sum number of bases covered per reference:
rpipunt_covered_bases <- rpipunt_big_dt[
  Depth >= 1,
  .(bases_covered = .N, total_mapped_bases = sum(Depth)),
  by = .(UniqueID, rname, Enrichment)
]

#Sum number of unconserved bases per reference:
rpipunt_mutations_count <- rpipunt_big_dt[
  !(Status %in% c("CONSERVED", "STAR")),
  .(mutation_count = .N),
  by = .(UniqueID, rname, Enrichment)
]

###############################
#####Combine imported data#####
###############################

#Combine counts, coverage and mutation data
rpip_and_unt_counts_and_covstats <- rpip_and_unt_counts_bigdf_tax %>%
  left_join(as_tibble(rpipunt_covered_bases)) %>%
  left_join(as_tibble(rpipunt_mutations_count)) %>%
  mutate(UniqueID = str_replace_all(UniqueID, "Non-targeted", "None"))

#Set NAs to 0
rpip_and_unt_counts_and_covstats[
  is.na(rpip_and_unt_counts_and_covstats$bases_covered), ]$bases_covered <- 0
rpip_and_unt_counts_and_covstats[
  is.na(rpip_and_unt_counts_and_covstats$mutation_count), ]$mutation_count <- 0

#Import QC info:
fastp_info_table <- read_rds(
  "input/modified/all_fastp_reports_dupdedup.rds"
) %>%
  mutate(UniqueID = paste(LIMS_ID, Treatment, Enrichment, sep = "-"),
         Enrichment = str_replace_all(Enrichment, "None", "Non-targeted"))

#Join the fastp data to the pileup data:
rpip_and_unt_counts_covstats_and_info <- rpip_and_unt_counts_and_covstats %>%
  left_join(fastp_info_table, by = c("UniqueID", "Enrichment")) %>%
  mutate(RA_nmapped = nmapped / summary.after_filtering.total_reads,
         RA_bases = total_mapped_bases / summary.after_filtering.total_bases,
         prcov = bases_covered / rlength)

#Recode E. brachy (where taxid is NA):
rpip_and_unt_counts_covstats_and_info <- 
  rpip_and_unt_counts_covstats_and_info %>%
  mutate(
    taxid = case_when(
      str_starts(
        name, "Eubacterium brachy ATCC 33089"
      ) ~ "Eubacterium_brachy_ATCC_33089",
      str_starts(
        name, "MAG\\: \\[Eubacterium",
      ) ~ "MAG:[Eubacterium]_brachy_isolate_JCVI_32_bin.35",
      .default = taxid)
  )

#Save the file:
write_rds(rpip_and_unt_counts_covstats_and_info,
          "input/modified/rpip_and_unt_counts_covstats_and_info.rds")
