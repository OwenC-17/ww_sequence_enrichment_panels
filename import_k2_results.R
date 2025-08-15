library(tidyverse)
library(zoo)

source("helper_functions.R")

#Import k2 reports with k2 confidence of 0.0 (default)
vsp_k2_reports_00conf <- import_kraken2_summary(vsp_dir_00conf, "report.tsv")
rpip_k2_reports_00conf <- import_kraken2_summary(rpip_dir_00conf, "report.tsv")
unt_k2_reports_00conf <- import_kraken2_summary(unt_dir_00conf, "report.tsv")

#Import deduped k2 reports with k2 confidence of 0.0
vsp_k2_reports_00conf_dedup <- import_kraken2_summary(vsp_dir_00conf_dedup,
                                                      "report.tsv")

#Import k2 reports with k2 confidence of 0.9
vsp_k2_reports_90conf <- import_kraken2_summary(vsp_dir_90conf, "report.tsv")
rpip_k2_reports_90conf <- import_kraken2_summary(rpip_dir_90conf, "report.tsv")
unt_k2_reports_90conf <- import_kraken2_summary(unt_dir_90conf, "report.tsv")


vsp_k2_reports_00conf <- vsp_k2_reports_00conf %>% 
  mutate(Enrichment = "VSP",
         Kraken2_confidence = "0.0", dedup = FALSE)
rpip_k2_reports_00conf <- rpip_k2_reports_00conf %>% 
  mutate(Enrichment = "RPIP",
         Kraken2_confidence = "0.0", dedup = FALSE)
unt_k2_reports_00conf <- unt_k2_reports_00conf %>% 
  mutate(Enrichment = "None",
         Kraken2_confidence = "0.0", dedup = FALSE)



vsp_k2_reports_90conf <- vsp_k2_reports_90conf %>% 
  mutate(Enrichment = "VSP",
         Kraken2_confidence = "0.9" dedup = FALSE)
rpip_k2_reports_90conf <- rpip_k2_reports_90conf %>% 
  mutate(Enrichment = "RPIP",
         Kraken2_confidence = "0.9", dedup = FALSE)
unt_k2_reports_90conf <- unt_k2_reports_90conf %>% 
  mutate(Enrichment = "None",
         Kraken2_confidence = "0.9", dedup = FALSE)


#Deduplicated versions
vsp_k2_reports_00conf_dedup <- vsp_k2_reports_00conf_dedup %>% 
  mutate(Enrichment = "VSP",
         Kraken2_confidence = "0.0", dedup = TRUE)
rpip_k2_reports_00conf_dedup <- rpip_k2_reports_00conf_dedup %>% 
  mutate(Enrichment = "RPIP",
         Kraken2_confidence = "0.0", dedup = TRUE)
unt_k2_reports_00conf_dedup <- unt_k2_reports_00conf_dedup %>% 
  mutate(Enrichment = "None",
         Kraken2_confidence = "0.0", dedup = TRUE)



vsp_k2_reports_90conf_dedup <- vsp_k2_reports_90conf_dedup %>% 
  mutate(Enrichment = "VSP",
         Kraken2_confidence = "0.9", dedup = TRUE)
rpip_k2_reports_90conf_dedup <- rpip_k2_reports_90conf_dedup %>% 
  mutate(Enrichment = "RPIP",
         Kraken2_confidence = "0.9", dedup = TRUE)
unt_k2_reports_90conf_dedup <- unt_k2_reports_90conf_dedup %>% 
  mutate(Enrichment = "None",
         Kraken2_confidence = "0.9", dedup = TRUE)



all_k2_reports_anyConf <- bind_rows(vsp_k2_reports_00conf,
                                    rpip_k2_reports_00conf,
                                    unt_k2_reports_00conf,
                                    vsp_k2_reports_90conf,
                                    rpip_k2_reports_90conf,
                                    unt_k2_reports_90conf)

rm(vsp_k2_reports_00conf, rpip_k2_reports_00conf, unt_k2_reports_00conf,
   vsp_k2_reports_90conf, rpip_k2_reports_90conf, unt_k2_reports_90conf)

dir.create("imported_k2_reports")
write_csv(all_k2_reports_anyConf, 
          "imported_k2_reports/all_k2_reports_anyConf.csv")

all_k2_reports_anyConf_dedup <- bins_rows(vsp_k2_reports_00conf_dedup,
                                          rpip_k2_reports_00conf_dedup,
                                          unt_k2_reports_00conf_dedup,
                                          vsp_k2_reports_90conf_dedup,
                                          rpip_k2_reports_90conf_dedup,
                                          unt_k2_reports_90conf_dedup)

rm(vsp_k2_reports_00conf_dedup, rpip_k2_reports_00conf_dedup, 
   unt_k2_reports_00conf_dedup, vsp_k2_reports_90conf_dedup, 
   rpip_k2_reports_90conf_dedup, unt_k2_reports_90conf_dedup)

write_csv(all_k2_reports_anyConf_dedup, 
          "imported_k2_reports/all_k2_reports_anyConf_dedup.csv")