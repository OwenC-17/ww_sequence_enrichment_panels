library(tidyverse)
library(edgeR)
library(vegan)
library(paletteer)

rpip_and_unt_counts_covstats_and_info <- read_rds(
  "input/modified/rpip_and_unt_counts_covstats_and_info.rds"
)

rpip_DGElist <- read_rds("input/modified/edgeR_bwamem_rpip_DGE.rds")

group_data <- read_csv("DS2_sample_metadata_w_rrna.csv") %>%
  arrange(by = UniqueID) %>%
  select(LIMS_ID, Treatment, site, Fraction, Nanotrap_type, Enrichment) %>%
  distinct() %>%
  mutate(UniqueID = paste(LIMS_ID, Treatment, Enrichment, sep=("-")),
         LIMS_ID = as.character(LIMS_ID))

#Get library sizes
rpip_libsizes = rpip_and_unt_counts_covstats_and_info %>%
  mutate(UniqueID = paste(LIMS_ID, Treatment, Enrichment, sep = "-")) %>%
  mutate(UniqueID = str_replace_all(UniqueID, "Non-targeted", "None")) %>%
  select(UniqueID, summary.after_dedup.total_reads) %>%
  distinct() %>%
  rename(libsize = summary.after_dedup.total_reads)

#Filter group data to RPIP and unenriched samples
rpip_group_data <- group_data %>%
  filter(Enrichment %in% c("RPIP", "None")) %>%
  left_join(rpip_libsizes)


#???
#edger_results_rpip <- read_rds("input/modified/rpip_bwamem_edger_results.rds")
#???
##################################
#####Read counts vs. Coverage#####
##################################
rpipunt_50plus_by_taxid <- rpip_and_unt_counts_covstats_and_info %>% 
  ungroup() %>%
  group_by(Enrichment, UniqueID, Treatment, LIMS_ID, Fraction, Nanotrap_type,
           across(taxid:domain), summary.after_dedup.total_bases, 
           summary.after_dedup.total_reads) %>%
  summarize(bases_covered = sum(bases_covered),
            nmapped = sum(nmapped),
            total_mapped_bases = sum(total_mapped_bases),
            comb_rlength = sum(rlength),
            mutation_count = sum(mutation_count)) %>%
  mutate(prcov = bases_covered / comb_rlength,
         RA_nmapped = nmapped / summary.after_dedup.total_reads,
         RA_bases = total_mapped_bases / summary.after_dedup.total_bases) %>%
  filter(bases_covered >= 50)


rpip_bacteria_df <- rpipunt_50plus_by_taxid %>%
  filter(domain == "Bacteria") %>%
  filter(nmapped > 0)


#Get the rpip_bacterias from the DGE object created in EdgeR_rgi_rpip_bacterialevel.R

#Transpose for use in Vegan functions, and select only RPIP samples
rpip_bacteria_community_matrix_noNA <- (rpip_DGElist$counts[
  rpip_DGElist$genes$Species %in% rpip_bacteria_df$species,
] %>% 
  t())[rpip_group_data$UniqueID,]

#Remove taxa with 0 total observed reads:
rpip_bacteria_community_matrix_noNA <- rpip_bacteria_community_matrix_noNA[
  , colSums(rpip_bacteria_community_matrix_noNA, na.rm = TRUE) > 0
]

#Make a version with NAs for some distance calculations:
rpip_bacteria_community_matrix <- rpip_bacteria_community_matrix_noNA
rpip_bacteria_community_matrix[rpip_bacteria_community_matrix == 0] <- NA_real_


#Calculate number of non-target sequences per sample:
Other <- rpip_group_data$libsize -
  rowSums(rpip_bacteria_community_matrix, na.rm = TRUE)
#Add column to community matrix
rpip_bacteria_community_matrix_complete <- cbind(
  rpip_bacteria_community_matrix, Other
)
rpip_bacteria_community_matrix_noNA_complete <- cbind(
  rpip_bacteria_community_matrix_noNA, Other
)

#Convert to relative abundance but don't include "Other" column
rpip_bacteria_relative <- 
  rpip_bacteria_community_matrix_noNA / rpip_group_data$libsize
#And again, this time including "Other" column
rpip_bacteria_relative_complete <- 
  rpip_bacteria_community_matrix_noNA_complete / rpip_group_data$libsize

#Wisconsin transformation
rpip_bacteria_wisconsin <- wisconsin(rpip_bacteria_community_matrix_noNA)

rpip_bacteria_wisconsin_complete <- wisconsin(
  rpip_bacteria_community_matrix_noNA_complete
)


#Bray-Curtis distance (delete if possible)
#rpip_bacteria_bray_relative <- vegdist(
#  rpip_bacteria_relative, method = "bray", na.rm = TRUE
#  )
#rpip_bacteria_bray_wisconsin <- vegdist(
#  rpip_bacteria_wisconsin, method = "bray", na.rm = TRUE)

#NMDS with 3 axes
#rpip_bacteria_bray_relative_nmds_3d <- metaMDS(
#  rpip_bacteria_relative,
#  distance = "bray",
#  k = 3,
#  try = 40,
#  trymax = 1000,
#  autotransform = FALSE
#)

#Get the NMDS coords:
#rpip_bacteria_bray_relative_nmds_3d_df <- rpip_bacteria_bray_relative_nmds_3d$points %>%
#  as.data.frame() %>%
#  mutate(sample_id = rownames(.)) %>%
#  left_join((rpip_group_data %>%
#               mutate(sample_id =  paste(
#                 LIMS_ID, Treatment, Enrichment, sep = "-"
#               ))), by = "sample_id")

#ggplot(rpip_bacteria_bray_relative_nmds_3d_df, aes(x = MDS2, y = MDS1, shape = Enrichment, color = Nanotrap_type)) + geom_point()


#NMDS with 2 axes
#rpip_bacteria_bray_relative_nmds_2d <- metaMDS(
#  rpip_bacteria_relative[rowSums(rpip_bacteria_relative) != 0, ],
#  distance = "bray",
#  k = 2,
#  try = 40,
#  trymax = 150,
#  autotransform = FALSE
#)

#Get NMDS coords:
#rpip_bacteria_bray_relative_nmds_2d_df <- rpip_bacteria_bray_relative_nmds_2d$points %>%
#  as.data.frame() %>%
#  mutate(sample_id = rownames(.)) %>%
#  left_join((rpip_group_data %>%
#               mutate(sample_id =  paste(
#                 LIMS_ID, Treatment, Enrichment, sep = "-"
#               ))), by = "sample_id")

#ggplot(rpip_bacteria_bray_relative_nmds_2d_df, aes(x = MDS1, y = MDS2, color = Fraction, shape = Enrichment)) + geom_point()

#Function to perform NMDS and add coords to a dataframe:
make_mds_df_rpip <- function(comm, group, k = 2, dist = "bray") {
  nmds <- metaMDS(
    comm[rowSums(comm, na.rm = TRUE) != 0,],
    distance = dist,
    k = k,
    try = 40,
    trymax = 400,
    autotransform = FALSE
  )
  
  nmds$points %>%
    as.data.frame() %>%
    mutate(UniqueID = rownames(.)) %>%
    left_join(group, by = "UniqueID")
}            

#Raw counts, not including "Other"
rpip_bacteria_bray_raw_nmds_2d_df_b <- make_mds_df_rpip(
  rpip_bacteria_community_matrix_noNA, rpip_group_data
)
ggplot(rpip_bacteria_bray_raw_nmds_2d_df_b, 
       aes(x = MDS1, y = MDS2, colour = site, shape = Enrichment)) +
  geom_point(size = 2) +
  geom_line(aes(group = paste(LIMS_ID, Treatment)), linetype = "dotted") +
  scale_colour_paletteer_d("ggsci::category10_d3") +
  theme_bw()


#Raw counts, including "Other"
rpip_bacteria_bray_raw_complete_nmds_2d_df_b <- make_mds_df_rpip(
  rpip_bacteria_community_matrix_noNA_complete, rpip_group_data
)
ggplot(rpip_bacteria_bray_raw_complete_nmds_2d_df_b, 
       aes(x = MDS1, y = MDS2, colour = site, shape = Enrichment)) +
  geom_point(size = 2) +
  geom_line(aes(group = paste(LIMS_ID, Treatment)), linetype = "dotted") +
  scale_colour_paletteer_d("ggsci::category10_d3") +
  theme_bw()


#Proportions, not including "Other"
rpip_bacteria_bray_relative_nmds_2d_df_b <- make_mds_df_rpip(rpip_bacteria_relative, rpip_group_data)
ggplot(rpip_bacteria_bray_relative_nmds_2d_df_b, 
       aes(x = MDS1, y = MDS2, colour = site, shape = Enrichment)) +
  geom_point(size = 2) +
  geom_line(aes(group = paste(LIMS_ID, Treatment)), linetype = "dotted") +
  scale_colour_paletteer_d("ggsci::category10_d3") +
  theme_bw()


#Wisconsin, not including "Other"
rpip_bacteria_bray_wisconsin_nmds_2d_df <- make_mds_df_rpip(
  rpip_bacteria_wisconsin, rpip_group_data
)
ggplot(rpip_bacteria_bray_wisconsin_nmds_2d_df, 
       aes(x = MDS1, y = MDS2, colour = site, shape = Enrichment)) +
  geom_point(size = 2) +
  geom_line(aes(group = paste(LIMS_ID, Treatment)), linetype = "dotted") +
  scale_colour_paletteer_d("ggsci::category10_d3") +
  theme_bw()


#Proportions, including "Other"
rpip_bacteria_bray_relative_complete_nmds_2d_df <- make_mds_df_rpip(
  rpip_bacteria_relative_complete, rpip_group_data
)
ggplot(rpip_bacteria_bray_relative_complete_nmds_2d_df, 
       aes(x = MDS1, y = MDS2, colour = site, shape = Enrichment)) +
  geom_point(size = 2) +
  geom_line(aes(group = paste(LIMS_ID, Treatment)), linetype = "dotted") +
  scale_colour_paletteer_d("ggsci::category10_d3") +
  theme_bw()

#Wisconsin, including "Other"
rpip_bacteria_bray_wisconsin_complete_nmds_2d_df <- make_mds_df_rpip(
  rpip_bacteria_wisconsin_complete, rpip_group_data
)
ggplot(rpip_bacteria_bray_wisconsin_complete_nmds_2d_df, 
       aes(x = MDS1, y = MDS2, colour = site, shape = Enrichment)) +
  geom_point(size = 2) +
  geom_line(aes(group = paste(LIMS_ID, Treatment)), linetype = "dotted") +
  scale_colour_paletteer_d("ggsci::category10_d3") +
  theme_bw()

#Jaccard, not including "Other" 
#(shouldn't need to include since it is present in all samples)
rpip_bacteria_jacc_nmds_2d_df <- make_mds_df_rpip(
  rpip_bacteria_community_matrix_noNA, rpip_group_data, dist = "jaccard"
)
ggplot(rpip_bacteria_jacc_nmds_2d_df,
       aes(x = MDS1, y = MDS2, colour = LIMS_ID, shape = Enrichment)) +
  geom_point(size = 2) +
  geom_line(aes(group = paste(LIMS_ID, Treatment)), linetype = "dotted") +
  scale_colour_paletteer_d("ggsci::category10_d3") +
  theme_bw()

#Jaccard-type Chao distance, not including "Other"
#(shouldn't need to include since it is present in all samples)
rpip_bacteria_chao_nmds_2d_df <- make_mds_df_rpip(
  rpip_bacteria_community_matrix_noNA, rpip_group_data, dist = "chao"
)
ggplot(rpip_bacteria_chao_nmds_2d_df,
       aes(x = MDS1, y = MDS2, colour = LIMS_ID, shape = Enrichment)) +
  geom_point(size = 2) +
  geom_line(aes(group = paste(LIMS_ID, Treatment)), linetype = "dotted") +
  scale_colour_paletteer_d("ggsci::category10_d3") +
  theme_bw()




####################
#####PERMANOVA######
####################

#Proportion, not including "Other"
adonis2(
  rpip_bacteria_relative[rowSums(rpip_bacteria_relative, na.rm = TRUE) != 0,] ~
    site + Fraction + Nanotrap_type + Enrichment, 
  data = rpip_group_data[rowSums(rpip_bacteria_relative, na.rm = TRUE) != 0,],
  by = "margin", na.rm = TRUE)


#Proportion, including "Other"
adonis2(rpip_bacteria_relative_complete ~ 
          site + Fraction + Nanotrap_type + Enrichment,
        data = rpip_group_data, by = "margin", na.rm = TRUE)

#Wisconsin, not including "Other"
set.seed(123)
adonis2(rpip_bacteria_wisconsin[
  rowSums(rpip_bacteria_wisconsin, na.rm = TRUE) != 0,] ~ 
    site + Fraction + Nanotrap_type + Enrichment,
  data = rpip_group_data[rowSums(rpip_bacteria_wisconsin, na.rm = TRUE) != 0,],
  by = "margin", na.rm = TRUE)

#Wisconsin, including "Other"
adonis2(rpip_bacteria_wisconsin_complete ~
          site + Fraction + Nanotrap_type + Enrichment,
        data = rpip_group_data, by = "margin", na.rm = TRUE)

#Jaccard
adonis2(rpip_bacteria_relative[rowSums(rpip_bacteria_relative, na.rm = TRUE) != 0,] ~ site + Fraction + Nanotrap_type + Enrichment, data = rpip_group_data[rowSums(rpip_bacteria_relative, na.rm = TRUE) != 0,], by = "margin", method = "jaccard", permutations = 9999)
adonis2(rpip_bacteria_relative_complete ~ site + Fraction + Nanotrap_type + Enrichment, data = rpip_group_data, by = "margin", method = "jaccard", permutations = 9999)


#Chao
adonis2(rpip_bacteria_community_matrix_noNA[rowSums(rpip_bacteria_community_matrix_noNA, na.rm = TRUE) != 0,] ~ site + Fraction + Nanotrap_type + Enrichment, data = rpip_group_data[rowSums(rpip_bacteria_community_matrix_noNA, na.rm = TRUE) != 0,], by = "margin", method = "chao")
adonis2(rpip_bacteria_community_matrix_noNA_complete ~ site + Fraction + Nanotrap_type + Enrichment, data = rpip_group_data, by = "margin", method = "chao")
