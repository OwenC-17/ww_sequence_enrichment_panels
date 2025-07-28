vsp_nonrrna_loc <- "../../../cowen20/targeted_panels/vsp_panels/raw_fastqs/fastp_out_no_dedup/ribodetector_out/nonrrna/vsp_nonrrna_ids.tsv"
rpip_nonrrna_loc <- "../../../cowen20/targeted_panels/rpip_panels/raw_fastqs/fastp_out_no_dedup/ribodetector_out/nonrrna/rpip_nonrrna_ids.tsv"
untargeted_nonrrna_loc <- "../../../cowen20/targeted_panels/untargeted/raw_fastqs/fastp_out_no_dedup/ribodetector_out/nonrrna/untargeted_nonrrna_ids.tsv"

import_readlist <- function(list_file) {
  message("checking for needed libraries...")
  require(data.table) 
  require(dtplyr) 
  require(dplyr) 
  require(vroom)
  message("Done checking libraries!
          \nImporting and formatting table of read names...")
  readlist <- vroom(list_file, delim = "\t", col_names = c("readid")) %>%
    lazy_dt() %>%
    distinct() %>%
    as_tibble %>%
    lazy_dt()
  message("Done importing read names!
          \nCreating vector for read names...")
  namevec <- character(length = nrow(readlist))
  message("Naming vector...")
  names(namevec) <- readlist$parent$readid
  message("Populating vector...")
  namevec[1:length(namevec)] <- readlist$parent$readid
  message("Done creating vector of read names!
          \nWrapping up...")
  return(namevec)
}

label_k2_nonrrna <- function(k2_out_file, readlist) {
  message("checking for needed libraries...")
  require(data.table)
  require(dtplyr)
  require(dplyr) 
  require(vroom)
  require(tidyr)
  message("Done checking libraries!
          \nImporting and formatting kraken2 output table...")
  k2 <- vroom(k2_out_file, 
              col_names = c("status", "readid", "taxid", 
                            "length", "mapping")) %>%
    lazy_dt()
  
  message("Done importing kraken2 output table!
          \nLabelling non-rRNA reads...")
  k2 <- k2 %>%
    mutate(is_nonrrna = (readid == readlist[readid])) %>%
    mutate(is_nonrrna = replace_na(.$is_nonrrna, FALSE)) %>%
    as_tibble()
  
  message("Finished!")
  return(k2)
}



label_all_k2_outputs <- function(k2_out_dir, readlist_to_use, labeled_dest) {
  for (k2outfile in list.files(k2_out_dir, pattern = "*.out")) {
    return_filename <- paste(labeled_dest, 
                             gsub("\\.out", "_labeled.tsv", k2outfile),
                             sep = "")
    labeled_table <- label_k2_nonrrna(paste(k2_out_dir, k2outfile, sep = ""), 
                                      readlist = readlist_to_use)

    vroom_write(x = labeled_table, file = return_filename, delim = "\t", progress = TRUE)
  }
}


separate_all_k2_outputs <- function(k2_out_dir, readlist_to_use, labeled_dest) {
  for (k2outfile in list.files(k2_out_dir, pattern = "*.out")) {
    nonrrna_filename <- paste(labeled_dest, 
                             gsub("\\.out", "_nonrrna.out", k2outfile),
                             sep = "")
    rrnaOnly_filename <- paste(labeled_dest, 
                              gsub("\\.out", "_rrnaOnly.out", k2outfile),
                              sep = "")
    labeled_table <- label_k2_nonrrna(paste(k2_out_dir, k2outfile, sep = ""), 
                                      readlist = readlist_to_use) %>% lazy_dt()
    
    nonrrna_table <- labeled_table %>%
      filter(is_nonrrna == TRUE) %>%
      select(status:mapping) %>%
      as_tibble()
    vroom_write(x = nonrrna_table, file = nonrrna_filename, delim = "\t", col_names = FALSE, progress = TRUE)
    rm(nonrrna_table)
    
    rrnaOnly_table <- labeled_table %>%
      filter(is_nonrrna == FALSE) %>%
      select(status:mapping) %>%
      as_tibble()
    vroom_write(x = rrnaOnly_table, file = rrnaOnly_filename, delim = "\t", col_names = FALSE, progress = TRUE)
    rm(rrnaOnly_table)
  }
}

##############

#Untargeted
untargeted_nonr_readlist <- import_readlist(untargeted_nonrrna_loc)

untargeted_k2out_001 <- "../../../cowen20/targeted_panels/untargeted/raw_fastqs/fastp_out_no_dedup/kraken2_out/k2_nt_20240530/001-36397-FA-RP_S1_L006_k2out.out"
untargeted_k2out_001_90conf <- "../../../cowen20/targeted_panels/untargeted/raw_fastqs/fastp_out_no_dedup/kraken2_out/k2_nt_20240530/90conf/001-36397-FA-RP_S1_L006_90confk2out.out"


rm(labeled_untargeted_k2out_001)
labeled_untargeted_k2out_001 <- label_k2_nonrrna(untargeted_k2out_001, untargeted_nonr_readlist)
labeled_untargeted_k2out_001_90conf <- label_k2_nonrrna(untargeted_k2out_001_90conf, untargeted_nonr_readlist)


untargeted_k2out_dir <- "../../../cowen20/targeted_panels/untargeted/raw_fastqs/fastp_out_no_dedup/kraken2_out/k2_nt_20240530/"
untargeted_k2labeled_dir <- "../../../cowen20/targeted_panels/untargeted/raw_fastqs/fastp_out_no_dedup/kraken2_out/k2_nt_20240530/rrna_labeled/"

untargeted_k2out_dir_90conf <- "../../../cowen20/targeted_panels/untargeted/raw_fastqs/fastp_out_no_dedup/kraken2_out/k2_nt_20240530/90conf/"
untargeted_k2labeled_dir_90conf <- "../../../cowen20/targeted_panels/untargeted/raw_fastqs/fastp_out_no_dedup/kraken2_out/k2_nt_20240530/90conf/rrna_labeled/"


separate_all_k2_outputs(untargeted_k2out_dir, untargeted_nonr_readlist, untargeted_k2labeled_dir)
rm(untargeted_nonr_readlist)

separate_all_k2_outputs(untargeted_k2out_dir_90conf, untargeted_nonr_readlist, untargeted_k2labeled_dir_90conf)
rm(untargeted_nonr_readlist)

###########
#VSP
vsp_k2out_dir <- "../../../cowen20/targeted_panels/vsp_panels/raw_fastqs/fastp_out_no_dedup/kraken2_out/k2_nt_20240530/"
vsp_k2labeled_dir <- "../../../cowen20/targeted_panels/vsp_panels/raw_fastqs/fastp_out_no_dedup/kraken2_out/k2_nt_20240530/rrna_labeled/"

vsp_k2out_dir_90conf <- "../../../cowen20/targeted_panels/vsp_panels/raw_fastqs/fastp_out_no_dedup/kraken2_out/k2_nt_20240530/90conf/"
vsp_k2labeled_dir_90conf <- "../../../cowen20/targeted_panels/vsp_panels/raw_fastqs/fastp_out_no_dedup/kraken2_out/k2_nt_20240530/90conf/rrna_labeled/"


vsp_nonr_readlist <- import_readlist(vsp_nonrrna_loc)

separate_all_k2_outputs(vsp_k2out_dir, vsp_nonr_readlist, vsp_k2labeled_dir)
rm(vsp_nonr_readlist)

separate_all_k2_outputs(vsp_k2out_dir_90conf, vsp_nonr_readlist, vsp_k2labeled_dir_90conf)
rm(vsp_nonr_readlist)


rpip_k2out_dir <- "../../../cowen20/targeted_panels/rpip_panels/raw_fastqs/fastp_out_no_dedup/kraken2_out/k2_nt_20240530/"
rpip_k2labeled_dir <- "../../../cowen20/targeted_panels/rpip_panels/raw_fastqs/fastp_out_no_dedup/kraken2_out/k2_nt_20240530/rrna_labeled/"

rpip_k2out_dir_90conf <- "../../../cowen20/targeted_panels/rpip_panels/raw_fastqs/fastp_out_no_dedup/kraken2_out/k2_nt_20240530/90conf/"
rpip_k2labeled_dir_90conf <- "../../../cowen20/targeted_panels/rpip_panels/raw_fastqs/fastp_out_no_dedup/kraken2_out/k2_nt_20240530/90conf/rrna_labeled/"


rpip_nonr_readlist <- import_readlist(rpip_nonrrna_loc)
separate_all_k2_outputs(rpip_k2out_dir_90conf, rpip_nonr_readlist, rpip_k2labeled_dir_90conf)

