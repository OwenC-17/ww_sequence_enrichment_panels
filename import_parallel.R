##For using multiple threads to import large kraken 2 output tables
##Input tables should be labeled for rRNA (filter_k2_out_by_rrna_status.R)

library(parallel)
library(foreach)
library(doParallel)
library(tidyverse)
library(zoo)
library(ggpubr)
library(gridExtra)
library(stringr)
library(microbenchmark)
library(vroom)

#Main dir that contains all targeted panel results
topdir <- "/projects/bios_microbe/cowen20/targeted_panels/"


#VSP Subdir
setwd(paste0(topdir, "vsp_panels/raw_fastqs/fastp_out_no_dedup/kraken2_out/",
             "k2_nt_20240530/rrna_labeled/"))

#Find all kraken 2 reports in current dir
vsp_file_list <- list.files(pattern = "*report.tsv")

#Define function to split content of SampleID column into more useful info
parse_sample_ids <- function(taxtable) {
  separate(taxtable,
           col = SampleID,
           into = c("QCSeqID", "LIMS_ID", "Treatment",
                    NA, NA, NA, NA, NA, 
                    "ribosomal", NA ), 
           sep="(-|_)")
}


#Define function to assign locations based on sample IDs
parse_locations <- function(taxtable) {
  mutate(taxtable, site = str_replace_all(
    LIMS_ID,
    c(
      "(36397|38156|41596|34970)" = "OBrien", #regex for any one of the given strings
      "(32976|32989|20040|22015)" = "OHare"
    )
    )
  )
}


#Start cluster for parallel computing
cl <- makeCluster(28)
registerDoParallel(cl = cl)

#Import all the VSP kraken 2 reports in parallel
#(generate_tax_table() and fill_tax_NAs() are from functions_for_tax_analysis.R)
imported_vsp_reports <- foreach(data = vsp_file_list, 
                          .packages = c("tidyverse", "zoo", "stringr")) %dopar%{
                            generate_tax_table(data) %>%
                              fill_tax_NAs() %>%
                              mutate(SampleID = data)
                            }

imported_vsp_reports <- parLapply(cl = cl, X = imported_vsp_reports, 
                                  fun = parse_sample_ids)

imported_vsp_reports <- parLapply(cl = cl, X = imported_vsp_reports, 
                                  fun = parse_locations)
stopCluster(cl)

################################################################################
#Now with RPIP

setwd(paste0(topdir, "rpip_panels/raw_fastqs/fastp_out_no_dedup/kraken2_out/k2_nt_20240530/rrna_labeled/"))

#Find all kraken 2 reports in current dir
rpip_file_list <- list.files(pattern = "*report.tsv")

#Start cluster for parallel computing
cl <- makeCluster(28)
registerDoParallel(cl = cl)

#Import all the RPIP kraken 2 reports in parallel
#(generate_tax_table() and fill_tax_NAs() are from functions_for_tax_analysis.R)
imported_rpip_reports <- foreach(data = rpip_file_list, 
                          .packages = c("tidyverse", "zoo", "stringr")) %dopar%{
                            generate_tax_table(data) %>%
                              fill_tax_NAs() %>%
                              mutate(SampleID = data)
                            }

imported_rpip_reports <- parLapply(cl = cl, X = imported_rpip_reports,
                                   fun = parse_sample_ids)
imported_rpip_reports <- parLapply(cl = cl, X = imported_rpip_reports,
                                   fun = parse_locations)

stopCluster(cl)

################################################################################
#And with untargeted

setwd(paste0(topdir, "untargeted/raw_fastqs/fastp_out_no_dedup/kraken2_out/k2_nt_20240530/rrna_labeled/"))

untargeted_file_list <- list.files(pattern = "*report.tsv")

#Start cluster for parallel computing
cl <- makeCluster(28)
registerDoParallel(cl = cl)

imported_untargeted_reports <- foreach(data = untargeted_file_list, 
                          .packages = c("tidyverse", "zoo", "stringr")) %dopar%{
                            generate_tax_table(data) %>%
                              fill_tax_NAs() %>%
                              mutate(SampleID = data)
                          }

imported_untargeted_reports <- parLapply(cl = cl, 
                                         X = imported_untargeted_reports,
                                         fun = parse_sample_ids)

imported_untargeted_reports <- parLapply(cl = cl,
                                         X = imported_untargeted_reports,
                                         fun = parse_locations)

################################################################################
#Combine the imported tables and write them to files:

vsp_big_df <- bind_rows(imported_vsp_reports) %>%
  mutate(SampleID = paste(.$LIMS_ID, .$Treatment, .$ribosomal, sep = "-"))

vroom_write(vsp_big_df,
            file = paste0(topdir, 
                      "vsp_panels/vsp_all_taxa_all_samples_rrna_separated.tsv"),
            col_names = TRUE,
            num_threads = 28)

rpip_big_df <- bind_rows(imported_rpip_reports) %>%
  mutate(SampleID = paste(.$LIMS_ID, .$Treatment, .$ribosomal, sep = "-"))

vroom_write(rpip_big_df, 
            file = paste0(topdir, 
                    "rpip_panels/rpip_all_taxa_all_samples_rrna_separated.tsv"),
            col_names = TRUE,
            num_threads = 28,)


untargeted_big_df <- bind_rows(imported_untargeted_reports) %>%
  mutate(SampleID = paste(.$LIMS_ID, .$Treatment, .$ribosomal, sep = "-"))

vroom_write(untargeted_big_df,
            file = paste0(topdir,
               "untargeted/untargeted_all_taxa_all_samples_rrna_separated.tsv"),
                col_names = TRUE,
                num_threads = 28)

################################################################################


ggplot(vsp_big_df, aes(x = as.factor(LIMS_ID), y = nodeOnly)) + geom_bar(position = 'stack', stat = 'identity') +
  ylab("Number of reads") +
  xlab("SampleID") + 
  theme_bw() +
  guides(col = "none") +
  theme(text = element_text(size = 15)) +
  guides(fill = guide_legend(ncol = 1)) +
  facet_grid(rows = vars(ribosomal), cols = vars(Treatment)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9))

ggplot(rpip_domains, aes(x = LIMS_ID, y = RA, fill = D)) + geom_bar(position = 'stack', stat = 'identity') +
  ylab("Relative abundance of reads") +
  xlab("SampleID") + 
  theme_bw() +
  labs(fill = "Domain") +
  guides(col = "none") +
  theme(text = element_text(size = 15)) +
  guides(fill = guide_legend(ncol = 1)) +
  facet_grid(rows = vars(ribosomal), cols = vars(Treatment)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9)) +
  ggtitle("RPIP")

ggplot(vsp_domains, aes(x = as.factor(LIMS_ID), y = RA, fill = D, color = D)) + geom_bar(position = 'stack', stat = 'identity') +
  ylab("Relative abundance of reads") +
  xlab("SampleID") + 
  theme_bw() +
  labs(fill = "Domain") +
  guides(col = "none") +
  theme(text = element_text(size = 15)) +
  guides(fill = guide_legend(ncol = 1)) +
  facet_grid(rows = vars(ribosomal), cols = vars(Treatment)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9)) +
  ggtitle("VSP")

ggplot(untargeted_domains, aes(x = LIMS_ID, y = RA, fill = D)) + geom_bar(position = 'stack', stat = 'identity') +
  ylab("Relative abundance of reads") +
  xlab("SampleID") + 
  theme_bw() +
  labs(fill = "Domain") +
  guides(col = "none") +
  theme(text = element_text(size = 15)) +
  guides(fill = guide_legend(ncol = 1)) +
  facet_grid(rows = vars(ribosomal), cols = vars(Treatment)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9)) +
  ggtitle("Untargeted")

ggplot(vsp_viral_phyla, aes(x = LIMS_ID, y = RA, fill = P)) + geom_bar(position = 'stack', stat = 'identity') +
  ylab("Relative abundance of reads") +
  xlab("SampleID") + 
  theme_bw() +
  labs(fill = "Phylum") +
  guides(col = "none") +
  theme(text = element_text(size = 15)) +
  guides(fill = guide_legend(ncol = 1)) +
  facet_grid(rows = vars(ribosomal), cols = vars(Treatment)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9))


rpip_mammalia <- filter(rpip_domains, C == "Mammalia") %>%
  filter(!startsWith(`F`, "unidentified"))

ggplot(rpip_mammalia, aes(x = LIMS_ID, y = RA, fill = `F`)) + geom_bar(position = 'stack', stat = 'identity') +
  ylab("Relative abundance of reads") +
  xlab("SampleID") + 
  theme_bw() +
  labs(fill = "Family") +
  guides(col = "none") +
  theme(text = element_text(size = 15)) +
  guides(fill = guide_legend(ncol = 3)) +
  facet_grid(rows = vars(ribosomal), cols = vars(Treatment)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9)) +
  ggtitle("Mammals?")

rpip_euk_phyla <- rpip_domains %>%
  filter(D == "Eukaryota") %>%
  filter(!startsWith(P, "uniden"))

test_read_vsp <- read_tsv("../../../../../../vsp_panels/vsp_all_taxa_all_samples_rrna_separated.tsv")
