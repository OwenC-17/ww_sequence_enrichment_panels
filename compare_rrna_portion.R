library(tidyverse)
source("helper_functions.R")

all_k2_reports_anyConf_dedup <- read_csv(
  "imported_k2_reports/all_k2_reports_anyConf_dedup.csv"
  )


#Check how many entries there are for each:
table(all_k2_reports_anyConf_dedup$ribosomal)

#Confirm that RA adds to 1 in each category
all_k2_reports_anyConf_dedup %>%
  group_by(LIMS_ID, Treatment, Enrichment, ribosomal, dedup, Kraken2_confidence) %>%
  summarize(RA = sum(RA), nreads = sum(nodeOnly)) %>%
  View()


#Calculate the portion of all rRNA per sample:
rrna_portions <- all_k2_reports_anyConf_dedup %>% 
  #Add labels
  parse_sample_treatments() %>%
  #It doesn't matter which confidence level we use here, as long as only one is
  #included:
  filter(Kraken2_confidence == 0.0) %>%
  #Add IDs (unique for each sequence library)
  mutate(UniqueID = paste(LIMS_ID, Treatment, Enrichment, sep = "-")) %>%
  #Group by each variable that we want to keep:
  group_by(UniqueID, LIMS_ID, Treatment, Fraction, Nanotrap_type, Enrichment,
           ribosomal) %>%
  #Summarize rRNA and non-rRNA in each sample (without regard to taxa)
  summarize(nreads = sum(nodeOnly)) %>%
  #Make ribosomal and non-ribosomal cols:
  pivot_wider(names_from = ribosomal, values_from = nreads, id_cols = UniqueID:Enrichment) %>%
  mutate(totalreads = nonrrna + rrnaOnly) %>%
  mutate(rrna_portion = rrnaOnly / totalreads,
         nonrrna_portion = nonrrna / totalreads)


rrna_portion_plotdf <- rrna_portions %>%
  mutate(Fraction = str_replace_all(Fraction,
                                    c("filtrate" = "Fil",
                                      "retentate" = "Ret",
                                      "unfiltered" = "Unfil")),
         Nanotrap_type = str_replace_all(Nanotrap_type,
                                         c("^A$" = "NT-A",
                                           "A&B" = "NT-AB",
                                           "none" = "DirEx")),
         Enrichment = str_replace_all(Enrichment,
                                      c("None" = "Non-targeted"))
  )

ggplot(rrna_portion_plotdf,
       aes(x = Nanotrap_type, y = nonrrna_portion * 100, fill = Enrichment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3", "gold2")) +
  theme_bw() +
  ylab("Percent non-ribosomal") +
  xlab("Nanotrap_type")


ggplot(rrna_portion_plotdf,
       aes(x = Fraction, y = nonrrna_portion * 100, fill = Enrichment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3", "gold2")) +
  theme_bw() +
  ylab("Percent non-ribosomal") +
  xlab("Fraction")
