LIMS_IDa <- factor(group_data_90conf_no_rrna$LIMS_ID)
Fractiona <- factor(group_data_90conf_no_rrna$Fraction, levels = c("unfiltered", "retentate", "filtrate"))
Nanotrap_typea <- factor(group_data_90conf_no_rrna$Nanotrap_type, levels = c("none", "A", "A&B"))
Enrichmenta <- factor(group_data_90conf_no_rrna$Enrichment, levels = c("None", "RPIP", "VSP"))


treat_lista <- list()

for (x in levels(Fractiona)) {
  print(x)
  for (y in levels(Nanotrap_typea)) {
    print(y)
    for (z in levels(Enrichmenta)) {
      print(z)
      tx_name <- paste(x, y, z, sep = ".")
      print(tx_name)
      assign(tx_name, (Fractiona == x & Nanotrap_typea == y & Enrichmenta == z))
      print(get(tx_name))
      treat_lista <- c(treat_lista, tx_name)
    }
  }
}
print(treat_lista)
additionala <- mget(unlist(treat_lista))


treat_lista <- list()

for (x in levels(Fractiona)) {
#  print(x)
  for (y in levels(Nanotrap_typea)) {
#    print(y)
    for (z in levels(Enrichmenta)) {
#      print(z)
      tx_name <- paste(x, y, z, sep = ".")
      print(paste("txname =", tx_name))
      assign(tx_name, (Fractiona == x & Nanotrap_typea == y & Enrichmenta == z))
      print(paste("get(tx_name):", get(tx_name)))
      treat_lista <- c(treat_lista, setNames(list(get(tx_name)), tx_name))
    }
  }
}

identical(treat_lista, additionala)




cov90_example <- import_coverage_of_vsp_ref("001-36397-FA_filtered_paired_90cov.tsv")
cov90_example$Enrichment <- "VSP"
cov90_example <- cov90_example %>%
  group_by(Enrichment, LIMS_ID, Treatment, species, tname) %>%
  summarize(total_length = sum(endpos), numreads = sum(numreads), covbases = sum(covbases)) %>%
  mutate(coverage = (covbases / total_length) * 100)

cov00_example <- import_coverage_of_vsp_ref("001-36397-FA-RP_S113_L007_bwa_vs_vsp_ref_mapped_and_filtered_coverage.tsv")
cov00_example$Enrichment <- "VSP"
cov00_example <- cov00_example %>%
  group_by(Enrichment, LIMS_ID, Treatment, tname, species) %>%
  summarize(total_length = sum(endpos), numreads = sum(numreads), covbases = sum(covbases)) %>%
  mutate(coverage = (covbases / total_length) * 100)

cov00_example$readCov <- "Any"
cov90_example$readCov <- "90%"

cov0090comp <- bind_rows(cov00_example, cov90_example)

ggplot(cov0090comp, aes(y = species, x = numreads, colour = readCov)) + geom_point() +
  scale_x_log10()

coveragePath <- "001-36397-FA-RP_S113_L007_bwa_vs_vsp_ref_mapped_and_filtered_sorted_readcoverage_tagged_PR0.1.tsv"
taxPath = sqlPath

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

for (x in list.files()) {
  print(as.double(str_extract(x, "(?<=PR)[:digit:](\\.?[:digit:]){0,2}")))
}
