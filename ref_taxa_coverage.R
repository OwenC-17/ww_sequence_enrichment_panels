#install.packages("taxonomizr")
library(taxonomizr)
#prepareDatabase('accessionTaxa.sql')

sqlPath = "/projects/bios_microbe/cowen20/ref_db/taxonomizr/accessionTaxa.sql"
setwd("/projects/bios_microbe/cowen20/targeted_panels/vsp_panels/raw_fastqs/fastp_out_no_dedup/bwa_out/mapped_and_low_complexity_removed/coverage_tables/")
unt_vs_vsp_dir <- "/projects/bios_microbe/cowen20/targeted_panels/untargeted/raw_fastqs/fastp_out_no_dedup/bwa_out/vs_vsp/mapped_and_low_complexity_removed/readcoverage_tagged/coverage_tables/"
unt_vs_rpip_dir <- "/projects/bios_microbe/cowen20/targeted_panels/untargeted/raw_fastqs/fastp_out_no_dedup/bwa_out/vs_rpip/mapped_and_low_complexity_removed/readcoverage_tagged/coverage_tables/"
vsp_dir <- "/projects/bios_microbe/cowen20/targeted_panels/vsp_panels/raw_fastqs/fastp_out_no_dedup/bwa_out/mapped_and_low_complexity_removed/readcoverage_tagged/coverage_tables/"
rpip_dir <- "/projects/bios_microbe/cowen20/targeted_panels/rpip_panels/raw_fastqs/fastp_out_no_dedup/bwa_out/mapped_and_low_complexity_removed/readcoverage_tagged/coverage_tables/"



import_coverage_of_vsp_ref <- function(coveragePath, taxPath = sqlPath) {
  
  coverage <- read_tsv(coveragePath, col_types = "ciiiidddd") %>%
    rename(rname = `#rname`)
  
  coverage$taxid <- accessionToTaxa(coverage$rname, taxPath)
  
  coverage$taxid <- as.character(coverage$taxid)
  
  rawTaxa <- getRawTaxonomy(unique(coverage$taxid), taxPath)
  
  rawTaxaDf <- bind_rows(rawTaxa, .id = "taxid")
  
  rawTaxaDf$taxid <- str_remove_all(rawTaxaDf$taxid, "\\s")
  
  rawTaxaDf$tname <- apply(rawTaxaDf[c("species", "no rank", "no rank.1", 
                                       "no rank.2", "serotype", "clade", 
                                       "clade.1", "genotype", "serotype")], 
                           1, function(x) {
                             
                             paste(na.omit(x), collapse = " / ")
                             
                           })
  
  coverage <- coverage %>%
    left_join(rawTaxaDf, by = "taxid")
  
  sampleInfo <- strsplit(coveragePath, "-")[[1]]
  
  coverage$LIMS_ID <- sampleInfo[2]
  coverage$Treatment <- sampleInfo[3]
  
  return(coverage)
  
}

cov1 <- import_coverage_of_vsp_ref("001-36397-FA-RP_S113_L007_bwa_vs_vsp_ref_mapped_and_filtered_coverage.tsv")
cov2 <- import_coverage_of_vsp_ref("002-36397-RA-RP_S114_L007_bwa_vs_vsp_ref_mapped_and_filtered_coverage.tsv")


rm(coverage_vsp_vs_vsp)

for (cvt in list.files(pattern = "*coverage.tsv$")) {
  print(cvt)
  }
for (cvt in list.files(pattern = "*coverage.tsv$")) {
  x = import_coverage_of_vsp_ref(cvt)
  if(exists("coverage_vsp_vs_vsp")) {
    coverage_vsp_vs_vsp <- bind_rows(coverage_vsp_vs_vsp, x)
  } else {
    coverage_vsp_vs_vsp <- x
  }
}

coverage_vsp_vs_vsp$Enrichment <- "VSP"

rm(coverage_unt_vs_vsp)
for (cvt in list.files(unt_vs_vsp_dir, pattern = "*coverage.tsv$")) {
  x = import_coverage_of_vsp_ref(paste0(unt_vs_vsp_dir, cvt))
  if(exists("coverage_unt_vs_vsp")) {
    coverage_unt_vs_vsp <- bind_rows(coverage_unt_vs_vsp, x)
  } else {
    coverage_unt_vs_vsp <- x
  }
}
coverage_unt_vs_vsp$Enrichment <- "None"

coverage_unt_and_vsp_vs_vsp <- bind_rows(coverage_unt_vs_vsp, coverage_vsp_vs_vsp)

coverage_unt_and_vsp_vs_vsp_by_tname <- coverage_unt_and_vsp_vs_vsp %>%
  group_by(Enrichment, LIMS_ID, Treatment, tname) %>%
  summarize(total_length = sum(endpos), numreads = sum(numreads), covbases = sum(covbases)) %>%
  mutate(coverage = (covbases / total_length) * 100)

ggplot(coverage_unt_and_vsp_vs_vsp_by_tname, aes(y = tname, x = coverage, fill = Enrichment)) + 
  geom_boxplot()



####################################################
###Impact of read portion aligned on ref coverage###
####################################################

setwd("../readcoverage_tagged/coverage_tables/")

import_readportion_coverage_of_vsp_ref <- function(coveragePath, taxPath = sqlPath) {
  
  coverage <- read_tsv(coveragePath, col_types = "ciiiidddd") %>%
    rename(rname = `#rname`)
  
  coverage$taxid <- accessionToTaxa(coverage$rname, taxPath)
  
  coverage$taxid <- as.character(coverage$taxid)
  
  rawTaxa <- getRawTaxonomy(unique(coverage$taxid), taxPath)
  
  rawTaxaDf <- bind_rows(rawTaxa, .id = "taxid")
  
  rawTaxaDf$taxid <- str_remove_all(rawTaxaDf$taxid, "\\s")
  
  rawTaxaDf$tname <- apply(rawTaxaDf[c("species", "no rank", "no rank.1", 
                                       "no rank.2", "serotype", "clade", 
                                       "clade.1", "genotype", "serotype")], 
                           1, function(x) {
                             
                             paste(na.omit(x), collapse = " / ")
                             
                           })
  
  coverage <- coverage %>%
    left_join(rawTaxaDf, by = "taxid")
  
  sampleInfo <- strsplit(coveragePath, "-")[[1]]
  
  coverage$LIMS_ID <- sampleInfo[2]
  coverage$Treatment <- sampleInfo[3]
  coverage$Threshold_readCov <- as.double(str_extract(sampleInfo[4], "(?<=PR)[:digit:](\\.?[:digit:]){0,2}"))
  
  return(coverage)
  
}

rm(read_coverage_vs_ref_coverage_vsp_vs_vsp)

for (cvt in list.files(path = vsp_dir, pattern = "*.tsv$")) {
  x = import_readportion_coverage_of_vsp_ref(paste0(vsp_dir, cvt))
  if(exists("read_coverage_vs_ref_coverage_vsp_vs_vsp")) {
    read_coverage_vs_ref_coverage_vsp_vs_vsp <- bind_rows(read_coverage_vs_ref_coverage_vsp_vs_vsp, x)
  } else {
    read_coverage_vs_ref_coverage_vsp_vs_vsp <- x
  }
}

read_coverage_vs_ref_coverage_vsp_vs_vsp$Enrichment <- "VSP"

grouped_read_coverage_vs_ref_coverage_vsp_vs_vsp <- read_coverage_vs_ref_coverage_vsp_vs_vsp %>%
  mutate(UniqueID = paste(LIMS_ID, Treatment, Enrichment, sep = "_")) %>%
  group_by(UniqueID, species, tname, Threshold_readCov) %>%
  summarize(total_length = sum(endpos), numreads = sum(numreads), covbases = sum(covbases)) %>%
  mutate(coverage = (covbases / total_length) * 100)


mean_across_samples_vsp_vs_vsp <- grouped_read_coverage_vs_ref_coverage_vsp_vs_vsp %>%
  ungroup() %>%
  group_by(Threshold_readCov, species, tname) %>%
  summarise(numreads = mean(numreads), covbases = mean(covbases), coverage = mean(coverage)) %>%
  mutate(Enrichment = "VSP")


ggplot(mean_across_samples_vsp_vs_vsp, aes(x = Threshold_readCov, y = numreads, colour = tname)) + geom_point() + geom_line() +
  theme(legend.position = "none")

ggplot(mean_across_samples_vsp_vs_vsp, aes(x = Threshold_readCov, y = coverage, colour = tname)) + geom_point() + geom_line() +
  theme(legend.position = "none")

ggplot(mean_across_samples_vsp_vs_vsp, aes(x = Threshold_readCov, y = covbases, colour = tname)) + geom_point() + geom_line() +
  geom_text(aes(label = species)) +
#  guides(color = guide_legend(ncol = 1))

  theme(legend.position = "none")

View(mean_across_samples_vsp_vs_vsp)




rm(read_coverage_vs_ref_coverage_unt_vs_vsp)

for (cvt in list.files(path = unt_vs_vsp_dir, pattern = "*.tsv$")) {
  x = import_readportion_coverage_of_vsp_ref(paste0(unt_vs_vsp_dir, cvt))
  if(exists("read_coverage_vs_ref_coverage_unt_vs_vsp")) {
    read_coverage_vs_ref_coverage_unt_vs_vsp <- bind_rows(read_coverage_vs_ref_coverage_unt_vs_vsp, x)
  } else {
    read_coverage_vs_ref_coverage_unt_vs_vsp <- x
  }
}

read_coverage_vs_ref_coverage_unt_vs_vsp$Enrichment <- "None"

grouped_read_coverage_vs_ref_coverage_unt_vs_vsp <- read_coverage_vs_ref_coverage_unt_vs_vsp %>%
  mutate(UniqueID = paste(LIMS_ID, Treatment, Enrichment, sep = "_")) %>%
  group_by(UniqueID, species, tname, Threshold_readCov) %>%
  summarize(total_length = sum(endpos), numreads = sum(numreads), covbases = sum(covbases)) %>%
  mutate(coverage = (covbases / total_length) * 100)


mean_across_samples_unt_vs_vsp <- grouped_read_coverage_vs_ref_coverage_unt_vs_vsp %>%
  ungroup() %>%
  group_by(Threshold_readCov, species, tname) %>%
  summarise(numreads = mean(numreads), covbases = mean(covbases), coverage = mean(coverage)) %>%
  mutate(Enrichment = "None")


ggplot(mean_across_samples, aes(x = Threshold_readCov, y = numreads, colour = tname)) + geom_point() + geom_line() +
  theme(legend.position = "none") + scale_y_log10()

ggplot(mean_across_samples, aes(x = Threshold_readCov, y = coverage, colour = tname)) + geom_point() + geom_line() +
  theme(legend.position = "none") + scale_y_log10()

ggplot(mean_across_samples, aes(x = Threshold_readCov, y = covbases, colour = tname)) + geom_point() + geom_line() +
  theme(legend.position = "none")



mean_across_samples_both_vs_vsp <- bind_rows(mean_across_samples_unt_vs_vsp, mean_across_samples_vsp_vs_vsp)

ggplot(mean_across_samples_both_vs_vsp, aes(x = Threshold_readCov, y = numreads, colour = tname)) + geom_point() + geom_line() +
  theme(legend.position = "none") +
  facet_wrap(~Enrichment, ncol = 1)

ggplot(mean_across_samples_both_vs_vsp, aes(x = Threshold_readCov, y = coverage, colour = tname)) + geom_point() + geom_line() +
  theme(legend.position = "none") +
  facet_wrap(~Enrichment, ncol = 1)

ggplot(mean_across_samples_both_vs_vsp, aes(x = Threshold_readCov, y = covbases, colour = tname)) + geom_point() + geom_line() +
  theme(legend.position = "none") +
  facet_wrap(~Enrichment, ncol = 1)


ggplot(mean_across_samples_vsp_vs_vsp, aes(x = numreads, y = coverage, colour = species, fill = species)) + geom_point() +
#  geom_smooth() + 
  theme(legend.position = "none")

ggplot(mean_across_samples_both_vs_vsp, aes(x = numreads, y = coverage, colour = species)) + geom_point() +
  theme(legend.position = "none") +
  facet_wrap(~Enrichment, ncol = 1, scales = "free") +
  scale_x_log10() +
  scale_y_log10()

ggplot(mean_across_samples_both_vs_vsp, aes(x = covbases, y = coverage, colour = species)) + geom_point() +
  theme(legend.position = "none") +
  facet_wrap(~Enrichment, ncol = 1, scales = "free")

###############
###RPIP EXAMPLE
###############

setwd("/projects/bios_microbe/cowen20/targeted_panels/rpip_panels/raw_fastqs/fastp_out_no_dedup/bwa_out/mapped_and_low_complexity_removed/coverage_tables/")



import_coverage_of_rpip_ref <- function(coveragePath, taxPath = sqlPath) {
  require(data.table)
  print(paste("Working on", coveragePath))
  coverage <- fread(coveragePath)
  setnames(coverage, "#rname", "rname")
  sampleInfo <- strsplit(coveragePath, "-")[[1]]
  coverage$LIMS_ID <- sampleInfo[2]
  coverage$Treatment <- sampleInfo[3]
  coverage$Threshold_readCov <- as.double(str_extract(sampleInfo[4], "(?<=PR)[:digit:](\\.?[:digit:]){0,2}"))
  return(coverage)
  
}

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

  rawTaxaDf <- rawTaxaDf %>%
    unite("tname", c("species", "subspecies", "no rank", "no rank.1", "no rank.2",
                   "serotype", "name", "varietas"), sep = " / ", na.rm = TRUE, 
        remove = FALSE) 

  annotated_table <- coverage_table %>%
  left_join(rawTaxaDf, by = c("rname"))
  
  return(annotated_table)

}

rm(read_coverage_vs_ref_coverage_rpip_vs_rpip)
file_list <- list.files(path = rpip_dir, pattern = "*.tsv")
placeholder_list <- vector("list", length = length(file_list))
for (cvt in seq_len(length(file_list))) {
  x = import_coverage_of_rpip_ref(paste0(rpip_dir, file_list[cvt]))
  placeholder_list[[cvt]] <- x
}

read_coverage_vs_ref_coverage_rpip_vs_rpip <- bind_rows(placeholder_list)

read_coverage_vs_ref_coverage_rpip_vs_rpip <- add_taxonomy(read_coverage_vs_ref_coverage_rpip_vs_rpip)
read_coverage_vs_ref_coverage_rpip_vs_rpip$Enrichment <- "RPIP"
viruses_only_read_coverage_vs_ref_coverage_rpip_vs_rpip <- read_coverage_vs_ref_coverage_rpip_vs_rpip %>%
  filter(`acellular root` == "Viruses")


rm(read_coverage_vs_ref_coverage_unt_vs_rpip)
file_list <- list.files(path = unt_vs_rpip_dir, pattern = "*.tsv")
placeholder_list <- vector("list", length = length(file_list))
for (cvt in seq_len(length(file_list))) {
  x = import_coverage_of_rpip_ref(paste0(unt_vs_rpip_dir, file_list[cvt]))
  placeholder_list[[cvt]] <- x
}

read_coverage_vs_ref_coverage_unt_vs_rpip <- bind_rows(placeholder_list)

read_coverage_vs_ref_coverage_unt_vs_rpip <- add_taxonomy(read_coverage_vs_ref_coverage_unt_vs_rpip)
read_coverage_vs_ref_coverage_unt_vs_rpip$Enrichment <- "None"
viruses_only_read_coverage_vs_ref_coverage_unt_vs_rpip <- read_coverage_vs_ref_coverage_unt_vs_rpip %>%
  filter(`acellular root` == "Viruses")


viruses_only_read_coverage_vs_ref_coverage_both_vs_rpip <- bind_rows(
  viruses_only_read_coverage_vs_ref_coverage_rpip_vs_rpip,
  viruses_only_read_coverage_vs_ref_coverage_unt_vs_rpip
)

viruses_only_90rc_both_vs_rpip <- viruses_only_read_coverage_vs_ref_coverage_both_vs_rpip %>%
  filter(Threshold_readCov == 0.9)


grouped_viruses_read_coverage_vs_ref_coverage_both_vs_rpip <- viruses_only_read_coverage_vs_ref_coverage_both_vs_rpip %>%
  mutate(UniqueID = paste(LIMS_ID, Treatment, Enrichment, sep = "_")) %>%
  group_by(Enrichment, UniqueID, species, tname, Threshold_readCov) %>%
  summarize(total_length = sum(endpos), numreads = sum(numreads), covbases = sum(covbases)) %>%
  mutate(coverage = (covbases / total_length) * 100)


mean_viruses_across_samples_both_vs_rpip <- grouped_viruses_read_coverage_vs_ref_coverage_both_vs_rpip %>%
  ungroup() %>%
  group_by(Enrichment, Threshold_readCov, species, tname) %>%
  summarise(numreads = mean(numreads), covbases = mean(covbases), coverage = mean(coverage))

ggplot(mean_viruses_across_samples_both_vs_rpip, aes(x = Threshold_readCov, y = numreads, colour = tname)) + geom_point() + geom_line() +
  theme(legend.position = "none") +
  facet_wrap(~Enrichment, ncol = 1)

ggplot(mean_viruses_across_samples_both_vs_rpip, aes(x = Threshold_readCov, y = coverage, colour = tname)) + geom_point() + geom_line() +
#  geom_text(aes(label = species)) +
  theme(legend.position = "none") +
  facet_wrap(~Enrichment, ncol = 1)

ggplot(mean_viruses_across_samples_both_vs_rpip, aes(x = Threshold_readCov, y = covbases, colour = tname)) + geom_point() + geom_line() +
  theme(legend.position = "none") +
  facet_wrap(~Enrichment, ncol = 1)


ggplot(read_coverage_90_rpip_vs_rpip, aes(x = numreads, y = covbases, colour = species, fill = species)) + geom_point() +
  #  geom_smooth() + 
  theme(legend.position = "none") +
  scale_x_log10() +
  scale_y_log10()

ggplot(read_coverage_90_unt_vs_rpip, aes(x = numreads, y = covbases, colour = species, fill = species)) + geom_point() +
  #  geom_smooth() + 
  theme(legend.position = "none") +
  scale_x_log10() +
  scale_y_log10()

ggplot(viruses_only_90rc_both_vs_rpip, aes(x = numreads, y = covbases, colour = Enrichment)) + geom_point(alpha = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  theme_minimal() +
  scale_color_manual(values = c("magenta", "yellow"))

ggplot(viruses_only_90rc_both_vs_rpip, aes(x = numreads, y = coverage, colour = Enrichment)) + geom_point(alpha = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  theme_minimal() +
  scale_color_manual(values = c("magenta", "yellow"))

ggplot(mean_across_samples_both_vs_vsp, aes(x = numreads, y = coverage, colour = species)) + geom_point() +
  theme(legend.position = "none") +
  facet_wrap(~Enrichment, ncol = 1, scales = "free") +
  scale_x_log10() +
  scale_y_log10()

ggplot(mean_across_samples_both_vs_vsp, aes(x = covbases, y = coverage, colour = species)) + geom_point() +
  theme(legend.position = "none") +
  facet_wrap(~Enrichment, ncol = 1, scales = "free")



read_coverage_90_rpip_vs_rpip <- read_coverage_vs_ref_coverage_rpip_vs_rpip %>%
  filter(numreads != 0 & Threshold_readCov == 0.9)

read_coverage_90_unt_vs_rpip <- read_coverage_vs_ref_coverage_unt_vs_rpip %>%
  filter(numreads != 0 & Threshold_readCov == 0.9)

read_coverage_90_both_vs_rpip <- bind_rows(read_coverage_90_rpip_vs_rpip, 
                                           read_coverage_90_unt_vs_rpip)

read_coverage_90_both_vs_rpip_betacorona <- filter(read_coverage_90_both_vs_rpip, genus == "Betacoronavirus")
