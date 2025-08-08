library(Rsamtools)
library(taxonomizr)
library(tidyverse)        
sqlPath = "/projects/bios_microbe/cowen20/ref_db/taxonomizr/accessionTaxa.sql"

pparam <- PileupParam(max_depth = 999999, include_insertions = TRUE, distinguish_strands = FALSE)

makePileup <- function(bamloc) {
  pileup(bamloc, index = paste0(bamloc, ".bai"), pileupParam = pparam) %>%
    mutate(SampleID = basename(bamloc))
}

freqs_from_pileup <- function(pileup) {
  pileup %>%
    as.data.frame() %>%
    group_by(SampleID, seqnames, pos, nucleotide) %>%
    summarize(countfr = sum(count)) %>%
    mutate(freq = countfr / sum(countfr))
}

calculate_h <- function(freqdf) {
  hdf <- freqdf %>%
    pivot_wider(names_from = nucleotide, values_from = c(freq, countfr), 
                values_fill = 0) %>%
    dplyr::rename(freq_Ins = `freq_+`, 
                  freq_Del = `freq_-`, 
                  countfr_Ins = `countfr_+`, 
                  countfr_Del = `countfr_-`)
    num_site_vars <- rowSums(select(ungroup(hdf), starts_with("freq_")) > 0)
    depth <- rowSums(select(ungroup(hdf), starts_with("countfr_")))
    hdf$total_nts <- num_site_vars
    hdf$n <- depth
    hdf$h <- ifelse(hdf$n > 1, ((hdf$n / (hdf$n - 1))*(1 - rowSums(select(ungroup(hdf),
                                                                          starts_with("freq"))
                                                                   ^2))),
                    0)
    return(hdf)
}

calculate_pi_from_h <- function(hdf) {
  hdf %>%
    group_by(SampleID, Enrichment, seqnames) %>%
    summarise(num_covsites = n(),
              num_segsites = sum(h != 0),
              totalpi = sum(h),
              persitepi = totalpi/n())
}

addTax <- function(df, taxPath = sqlPath) {
  df$taxid <- accessionToTaxa(df$seqnames, taxPath)
  df$taxid <- as.character(df$taxid)
  rawTaxa <- getRawTaxonomy(unique(df$taxid), taxPath)
  rawTaxaDf <- bind_rows(rawTaxa, .id = "taxid")
  rawTaxaDf$taxid <- str_remove_all(rawTaxaDf$taxid, "\\s")
  df <- df %>%
    left_join(rawTaxaDf, by = "taxid")
  return(df)
}

import_readportion_coverage_of_vsp_ref <- function(coveragePath, taxPath = sqlPath) {
  
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

vsp_bam_location <- paste0("input/link_to_raw_data/vsp_panels/raw_fastqs/",
                           "fastp_out_no_dedup/bwa_out/",
                           "mapped_and_low_complexity_removed")
vsp_bams <- list.files(vsp_bam_location, pattern = "*sorted.bam$", 
                       full.names = TRUE)
vsp_pileups <- lapply(vsp_bams, makePileup)
vsp_freqs <- lapply(vsp_pileups, freqs_from_pileup)
vsp_freqs <- bind_rows(vsp_freqs)
vsp_freqs$Enrichment <- "VSP"
vsp_h <- calculate_h(vsp_freqs)
vsp_pi <- calculate_pi_from_h(vsp_h)

unt_vs_vsp_bam_location <- paste0("input/link_to_raw_data/untargeted/",
                                  "raw_fastqs/fastp_out_no_dedup/bwa_out/",
                                  "vs_vsp/mapped_and_low_complexity_removed")
unt_vs_vsp_bams <- list.files(unt_vs_vsp_bam_location, pattern = "*sorted.bam$", 
                       full.names = TRUE)
unt_vs_vsp_pileups <- lapply(unt_vs_vsp_bams, makePileup)
unt_vs_vsp_freqs <- lapply(unt_vs_vsp_pileups, freqs_from_pileup)
unt_vs_vsp_freqs <- bind_rows(unt_vs_vsp_freqs)
unt_vs_vsp_freqs$Enrichment <- "None"
unt_vs_vsp_h <- calculate_h(unt_vs_vsp_freqs)
unt_vs_vsp_pi <- calculate_pi_from_h(unt_vs_vsp_h)

unt_and_vsp_pi <- bind_rows(unt_vs_vsp_pi, vsp_pi) %>%
  parse_sample_ids() %>%
  parse_locations()

unt_and_vsp_pi <- addTax(unt_and_vsp_pi)

ggplot(unt_and_vsp_pi, aes(x = species, y = persitepi, colour = Enrichment)) + 
  geom_boxplot() + 
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



################################################################################
rpip_bam_location <- paste0("input/link_to_raw_data/rpip_panels/raw_fastqs/",
                           "fastp_out_no_dedup/bwa_out/",
                           "mapped_and_low_complexity_removed")
rpip_bams <- list.files(rpip_bam_location, pattern = "*sorted.bam$", 
                       full.names = TRUE)
rpip_pileups <- lapply(rpip_bams, makePileup)

library(doParallel)
library(foreach)

cl <- makeCluster(48)
registerDoParallel(cl)

rpip_freqs <- foreach(pu = rpip_pileups, .packages = c("dplyr")) %dopar% {
  message(pu$sampleID[1])
  freqs_from_pileup(pu)
}

stopCluster(cl)

rpip_freqs <- bind_rows(rpip_freqs)
rpip_freqs$Enrichment <- "RPIP"
rpip_h <- calculate_h(rpip_freqs)
vsp_pi <- calculate_pi_from_h(vsp_h)

unt_vs_vsp_bam_location <- paste0("input/link_to_raw_data/untargeted/",
                                  "raw_fastqs/fastp_out_no_dedup/bwa_out/",
                                  "vs_vsp/mapped_and_low_complexity_removed")
unt_vs_vsp_bams <- list.files(unt_vs_vsp_bam_location, pattern = "*sorted.bam$", 
                              full.names = TRUE)
unt_vs_vsp_pileups <- lapply(unt_vs_vsp_bams, makePileup)
unt_vs_vsp_freqs <- lapply(unt_vs_vsp_pileups, freqs_from_pileup)
unt_vs_vsp_freqs <- bind_rows(unt_vs_vsp_freqs)
unt_vs_vsp_freqs$Enrichment <- "None"
unt_vs_vsp_h <- calculate_h(unt_vs_vsp_freqs)
unt_vs_vsp_pi <- calculate_pi_from_h(unt_vs_vsp_h)

unt_and_vsp_pi <- bind_rows(unt_vs_vsp_pi, vsp_pi) %>%
  parse_sample_ids() %>%
  parse_locations()

unt_and_vsp_pi <- addTax(unt_and_vsp_pi)

ggplot(unt_and_vsp_pi, aes(x = species, y = persitepi, colour = Enrichment)) + 
  geom_boxplot() + 
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))