library(data.table)
library(dtplyr)
library(tidyverse)
library(vroom)
library(tidyr)

k2out001 <- vroom("../../csvtk_exp/header_added.tsv") %>% lazy_dt()

nonrrna_names <- vroom("../../targeted_panels/vsp_panels/raw_fastqs/fastp_out_no_dedup/ribodetector_out/rrnaOnly/vsp_rrnaOnly_ids.tsv", 
                       delim = "\t", 
                       col_names = c("readid")) %>%
  lazy_dt()

nonrrna_names <- nonrrna_names %>%
  distinct() %>%
  as_tibble %>%
  lazy_dt()

#nonrrna_names <- as.character(nonrrna_names$readid)
namevec <- character(length = nrow(nonrrna_names))
names(namevec) <- nonrrna_names$parent$readid
namevec[1:length(namevec)] <- nonrrna_names$parent$readid

###Testing
######################
library(microbenchmark)
miniout <- head(k2out001, n = 1000000) %>% as_tibble() %>% lazy_dt()

randids <- sample(1:nrow(nonrrna_names), size = 1000000, replace = FALSE)
mininonrrnanames <- nonrrna_names$parent[randids,]
mininamevec <- character(length = nrow(mininonrrnanames))
names(mininamevec) <- mininonrrnanames$readid

mininamevec[1:nrow(mininonrrnanames)] <- names(mininamevec)

#based on a vector
f1 <- function(a = miniout, b = mininamevec) {
  it <- a %>%
    mutate(isrrna = (readid %in% b)) %>%
    as_tibble()
  return(it)
}

#based on a column
f2 <- function(a = miniout, b = mininonrrnanames) {
  it <- a %>%
    mutate(isrrna = (readid %in% b$readid)) %>%
    as_tibble()
  return(it)
}

#based on indexed vector
f3 <- function(a = miniout, b = mininamevec) {
  it <- a %>%
    mutate(isrrna = (readid == b[readid])) %>%
    as_tibble()
  return(it)
}

test1 <- f1()

test2 <- f2()

test3 <- f3()

compare <- microbenchmark(f1, f2, f3, times = 100)
compare

classified_nonrrna <- f3(a = k2out001, b = namevec)

ggplot(test_single, aes(x = "SampleID", y = RA, fill = F)) + geom_bar(position = 'stack', stat = 'identity') +
  ylab("Relative abundance of reads") +
  xlab("sampleID") + 
  theme_bw() +
  labs(fill = "Domain") + 
  theme(text = element_text(size = 15)) +
  guides(fill = guide_legend(ncol = 1))

unique(merged_vsp_reports_full$LIMS_ID)

library(parallel)
detectCores()
