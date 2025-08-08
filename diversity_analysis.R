BiocManager::install("Rsamtools")
library(Rsamtools)
bam_ <- Rsamtools::scanBam("input/link_to_raw_data/vsp_panels/raw_fastqs/fastp_out_no_dedup/bwa_out/mapped_and_low_complexity_removed/001-36397-FA-RP_S113_L007_bwa_vs_vsp_ref_mapped_and_filtered_sorted.bam")

pp_ <- PileupParam(max_depth = 999999, include_insertions = TRUE)

pileup_ <- Rsamtools::pileup("input/link_to_raw_data/vsp_panels/raw_fastqs/fastp_out_no_dedup/bwa_out/mapped_and_low_complexity_removed/012-32989-RA-RP_S124_L007_bwa_vs_vsp_ref_mapped_and_filtered_sorted.bam", 
                  index = "input/link_to_raw_data/vsp_panels/raw_fastqs/fastp_out_no_dedup/bwa_out/mapped_and_low_complexity_removed/012-32989-RA-RP_S124_L007_bwa_vs_vsp_ref_mapped_and_filtered_sorted.bam.bai", pileupParam = pp_)

freqs_from_pileup <- function(pileup) {
  pileup %>%
    as.data.frame() %>%
    group_by(seqnames, pos, nucleotide) %>%
    summarize(countfr = sum(count)) %>%
    mutate(freq = countfr / sum(countfr))
}

calculate_h <- function(freqdf) {
  freqdf %>%
    pivot_wider(names_from = nucleotide, values_from = c(freq, countfr), 
                values_fill = 0) %>%
    dplyr::rename(freq_Ins = `freq_+`, 
                  freq_Del = `freq_-`, 
                  countfr_Ins = `countfr_+`, 
                  countfr_Del = `countfr_-`) %>%
    mutate(total_nts = rowSums(across(c(freq_A, 
                                        freq_C,
                                        freq_G,
                                        freq_T,
                                        freq_Ins,
                                        freq_Del), 
                                      ~ .  > 0)),
           n = (countfr_Ins +
               countfr_Del +
               countfr_A +
               countfr_C +
               countfr_G +
               countfr_T),
           h = case_when(n > 1 ~ ((n  / (n - 1))*(1 - freq_A^2 -
                                                      freq_C^2 -
                                                      freq_G^2 -
                                                      freq_T^2 -
                                                      freq_Ins^2 -
                                                      freq_Del^2)),
                          .default = 0))
}

calculate_pi_from_h <- function(hdf) {
  hdf %>%
    group_by(seqnames) %>%
    summarise(num_covsites = n(),
              num_segsites = sum(h != 0),
              totalpi = sum(h),
              persitepi = totalpi/n())
}