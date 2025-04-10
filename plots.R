library(tidyverse)

#Location where all the output tables are (in different subdirs):
topdir = "/projects/bios_microbe/cowen20/targeted_panels/"


#There are a lot of weird domains, like artificial sequences, plasmids etc. 
#This will help filter them.
desirable_domains <- c("Eukaryota", "Bacteria", "Archaea", "Viruses")

#Import data
################################################################################

#Read VSP df
vsp_big_df <- read_tsv(paste0(topdir, 
                      "vsp_panels/vsp_all_taxa_all_samples_rrna_separated.tsv"),
                      guess_max = Inf)
#Fix data types
vsp_big_df <- vsp_big_df %>%
  mutate(taxID = as.character(taxID),
         LIMS_ID = as.character(LIMS_ID))

########################################
#Read RPIP df
rpip_big_df <- read_tsv(paste0(topdir, 
                    "rpip_panels/rpip_all_taxa_all_samples_rrna_separated.tsv"),
                    guess_max = Inf)
#Fix data types
rpip_big_df <- rpip_big_df %>%
  mutate(taxID = as.character(taxID),
         LIMS_ID = as.character(LIMS_ID))

########################################
#Read untargeted df
untargeted_big_df <- read_tsv(paste0(topdir, 
               "untargeted/untargeted_all_taxa_all_samples_rrna_separated.tsv"),
               guess_max = Inf)
#Fix data types
untargeted_big_df <- untargeted_big_df %>%
  mutate(taxID = as.character(taxID),
         LIMS_ID = as.character(LIMS_ID))

################################################################################
#Filter to desired domains and group by domain
vsp_desired_domains <- vsp_big_df %>%
  group_by(Treatment, LIMS_ID, ribosomal, D) %>%
  filter(D %in% desirable_domains) %>%
  summarise(abundance = sum(RA))

#Barplot of desired domains
ggplot(vsp_desired_domains, aes(x = LIMS_ID, y = abundance, fill = D)) +
  geom_bar(position = 'stack', stat = 'identity') +
  ylab("Relative abundance of reads") +
  xlab("SampleID") +
  theme_bw() +
  labs(fill = "Domain") +
  guides(col = "none") +
  theme(text = element_text(size = 15)) +
  guides(fill = guide_legend(ncol = 1)) +
  facet_grid(rows = vars(ribosomal), cols = vars(Treatment)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9)) +
  ggtitle("VSP Domains")


#Filter to desired domains and group by domain
rpip_desired_domains <- rpip_big_df %>%
  group_by(Treatment, LIMS_ID, ribosomal, D) %>%
  filter(D %in% desirable_domains) %>%
  summarise(abundance = sum(RA))

#Barplot of desired domains
ggplot(rpip_desired_domains, aes(x = LIMS_ID, y = abundance, fill = D)) +
  geom_bar(position = 'stack', stat = 'identity') +
  ylab("Relative abundance of reads") +
  xlab("SampleID") +
  theme_bw() +
  labs(fill = "Domain") +
  guides(col = "none") +
  theme(text = element_text(size = 15)) +
  guides(fill = guide_legend(ncol = 1)) +
  facet_grid(rows = vars(ribosomal), cols = vars(Treatment)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9)) +
  ggtitle("RPIP Domains")


#Filter to desired domains and group by domain
untargeted_desired_domains <- untargeted_big_df %>%
  group_by(Treatment, LIMS_ID, ribosomal, D) %>%
  filter(D %in% desirable_domains) %>%
  summarise(abundance = sum(RA))

#Barplot of desired domains
ggplot(untargeted_desired_domains, aes(x = LIMS_ID, y = abundance, fill = D)) +
  geom_bar(position = 'stack', stat = 'identity') +
  ylab("Relative abundance of reads") +
  xlab("SampleID") +
  theme_bw() +
  labs(fill = "Domain") +
  guides(col = "none") +
  theme(text = element_text(size = 15)) +
  guides(fill = guide_legend(ncol = 1)) +
  facet_grid(rows = vars(ribosomal), cols = vars(Treatment)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9)) +
  ggtitle("Untargeted Domains")


vsp_viral <- vsp_big_df %>%
  filter(D == "Viruses")

vsp_viral_families <- vsp_viral %>%
  filter(!startsWith(`F`, "unidentified"))

vsp_viral_orders <- vsp_viral %>%
  filter(!startsWith(O, "unidentified"))

vsp_viral_classes <- vsp_viral %>%
  filter(!startsWith(C, "unidentified"))

vsp_viral_phyla <- vsp_viral %>%
  filter(!startsWith(P, "unidentified"))

vsp_coronaviridae <- vsp_viral %>%
  filter(`F` == "Coronaviridae")



#VSP barplot of classified domains
vsp_domains %>%
  group_by(Treatment, LIMS_ID, ribosomal, D) %>%
  summarize(abundance = sum(RA)) %>%
  ggplot(aes(x = LIMS_ID, y = abundance, fill = D)) + 
  geom_bar(position = 'stack', stat = 'identity') +
  ylab("Relative abundance of reads") +
  xlab("SampleID") +
  theme_bw() +
  labs(fill = "Domain") +
  guides(col = "none") +
  theme(text = element_text(size = 15)) +
  guides(fill = guide_legend(ncol = 1)) +
  facet_grid(rows = vars(ribosomal), cols = vars(Treatment)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9))

#VSP barplot, including weird domains
vsp_big_df %>%
  group_by(Treatment, LIMS_ID, ribosomal, D) %>%
  summarize(abundance = sum(RA)) %>%
  ggplot(aes(x = LIMS_ID, y = abundance, fill = D)) + 
  geom_bar(position = 'stack', stat = 'identity') +
  ylab("Relative abundance of reads") +
  xlab("SampleID") +
  theme_bw() +
  labs(fill = "Domain") +
  guides(col = "none") +
  theme(text = element_text(size = 15)) +
  guides(fill = guide_legend(ncol = 1)) +
  facet_grid(rows = vars(ribosomal), cols = vars(Treatment)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9))

vsp_weird_domains <- vsp_big_df %>%
  filter(!(D %in% desirable_domains))

vsp_weird_domains %>%
  group_by(Treatment, LIMS_ID, ribosomal, D) %>%
  summarize(abundance = sum(RA)) %>%
  ggplot(aes(x = LIMS_ID, y = abundance, fill = D)) + 
  geom_bar(position = 'stack', stat = 'identity') +
  ylab("Relative abundance of reads") +
  xlab("SampleID") +
  theme_bw() +
  labs(fill = "Domain") +
  guides(col = "none") +
  theme(text = element_text(size = 15)) +
  guides(fill = guide_legend(ncol = 1)) +
  facet_grid(rows = vars(ribosomal), cols = vars(Treatment)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9))


rpip_big_df <- bind_rows(imported_rpip_reports) %>%
  mutate(SampleID = paste(.$LIMS_ID, .$Treatment, .$ribosomal, sep = "-"))
rpip_domains <- rpip_big_df %>%
  filter(D %in% desirable_domains)


#rpip barplot of classified domains
ggplot(rpip_domains, aes(x = LIMS_ID, y = RA, fill = D)) + geom_bar(position = 'stack', stat = 'identity') +
  ylab("Relative abundance of reads") +
  xlab("SampleID") + 
  theme_bw() +
  labs(fill = "Domain") +
  guides(col = "none") +
  theme(text = element_text(size = 15)) +
  guides(fill = guide_legend(ncol = 1)) +
  facet_grid(rows = vars(ribosomal), cols = vars(Treatment)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9))


#rpip barplot of eukaryotes

rpip_euk <- rpip_big_df %>%
  filter(D == "Eukaryota")
ggplot(rpip_euk, aes(x = LIMS_ID, y = RA, fill = P)) + geom_bar(position = 'stack', stat = 'identity') +
  ylab("Relative abundance of reads") +
  xlab("SampleID") + 
  theme_bw() +
  labs(fill = "Phylum") +
  guides(col = "none") +
  theme(text = element_text(size = 15)) +
  guides(fill = guide_legend(ncol = 3)) +
  facet_grid(rows = vars(ribosomal), cols = vars(Treatment)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9))


rpip_fungi <- rpip_big_df %>%
  filter(D2 == "Fungi")

rpip_animals <- rpip_big_df %>%
  filter(D2 == "Metazoa")
ggplot(rpip_animals, aes(x = LIMS_ID, y = RA, fill = P)) + geom_bar(position = 'stack', stat = 'identity') +
  ylab("Relative abundance of reads") +
  xlab("SampleID") + 
  theme_bw() +
  labs(fill = "Phylum") +
  guides(col = "none") +
  theme(text = element_text(size = 15)) +
  guides(fill = guide_legend(ncol = 3)) +
  facet_grid(rows = vars(ribosomal), cols = vars(Treatment)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9))


rpip_arthropods <- rpip_big_df %>%
  filter(P == "Arthropoda")
rpip_insects <- rpip_big_df %>%
  filter(C == "Insecta") %>%
  arrange(desc(nodeOnly))
vsp_insects <- vsp_big_df %>%
  filter(C == "Insecta") %>%
  arrange(desc(nodeOnly))


ggplot(rpip_insects, aes(x = LIMS_ID, y = RA, fill = O)) + geom_bar(position = 'stack', stat = 'identity') +
  ylab("Relative abundance of reads") +
  xlab("SampleID") + 
  theme_bw() +
  labs(fill = "Order") +
  guides(col = "none") +
  theme(text = element_text(size = 15)) +
  guides(fill = guide_legend(ncol = 3)) +
  facet_grid(rows = vars(ribosomal), cols = vars(Treatment)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9))


vsp_corona <- vsp_big_df %>%
  filter(`F` == "Coronaviridae")
rpip_corona <- rpip_big_df %>%
  filter(`F` == "Coronaviridae")

vsp_mammalia <- vsp_big_df %>%
  filter(C == "Mammalia")
