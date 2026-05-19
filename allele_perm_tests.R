#Perform ordinations and permutation tests on RGI output (ALLELE level)


library(tidyverse)
library(vegan)
library(edgeR)

##########################
#####Read in the data#####
##########################
rpip_and_unt_allele_and_info <- read_rds(
  "input/modified/rpip_and_unt_rgi_allele_protein_homolog.rds"
  )   %>%
  mutate(Fraction = str_replace_all(Fraction,
                                    c("filtrate" = "Fil",
                                      "retentate" = "Ret",
                                      "unfiltered" = "Unfil")),
         Nanotrap_type = str_replace_all(Nanotrap_type,
                                         c("^A$" = "NT-A",
                                           "A&B" = "NT-AB",
                                           "none" = "DirEx")),
         Enrichment = str_replace_all(Enrichment,
                                      c("None" = "Non-targeted")))

#Get group data for RPIP only
group_data <- read_csv("DS2_sample_metadata_w_rrna.csv") %>%
  arrange(by = UniqueID) %>%
  dplyr::select(LIMS_ID, Treatment, site, Fraction,
                Nanotrap_type, Enrichment) %>%
  distinct() %>%
  mutate(UniqueID = paste(LIMS_ID, Treatment, Enrichment, sep=("-")),
         LIMS_ID = as.character(LIMS_ID),
         Fraction = str_replace_all(Fraction,
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


#Get read counts
rpip_libsizes = rpip_and_unt_allele_and_info %>%
  select(UniqueID, summary.after_dedup.total_reads) %>%
  distinct() %>%
  rename("libsize" = "summary.after_dedup.total_reads")

#Remove VSP samples (irrelevant here)
rpip_allele_group_data <- group_data %>%
  filter(Enrichment %in% c("RPIP", "Non-targeted")) %>%
  left_join(rpip_libsizes)

#Get the alleles from the DGE object created in EdgeR_rgi_allelelevel.R:
rpip_allele_DGElist <- read_rds("input/modified/edgeR_rgi_allele_DGE.rds")


#######################################
#####Make community count matrices#####
#######################################

#Transpose for use in Vegan functions, and select only RPIP and untargeted
#samples
allele_community_matrix_noNA <- (rpip_allele_DGElist$counts %>%
                                   t())[rpip_allele_group_data$UniqueID, ]
allele_community_matrix_noNA <- allele_community_matrix_noNA[ ,
                  colSums(allele_community_matrix_noNA, na.rm = TRUE) > 0]
allele_community_matrix <- allele_community_matrix_noNA
allele_community_matrix[allele_community_matrix == 0] <- NA_real_


#Calculate number of non-ARG sequences per sample
Other <- rpip_allele_group_data$libsize - 
  rowSums(allele_community_matrix, na.rm = TRUE)
#Add column to community matrix
allele_community_matrix_complete <- cbind(
  allele_community_matrix, Other
  )
allele_community_matrix_noNA_complete <- cbind(
  allele_community_matrix_noNA, Other
  )

#Convert to relative abundance but don't include "Other" column
allele_relative <- allele_community_matrix_noNA / rpip_allele_group_data$libsize
#And again, this time including "Other" column
allele_relative_complete <- allele_community_matrix_noNA_complete /
  rpip_allele_group_data$libsize

#Wisconsin transformation, without "Other"
allele_wisconsin <- wisconsin(allele_community_matrix_noNA)

#Wisconsin transformation, with "Other"
allele_wisconsin_complete <- wisconsin(allele_community_matrix_noNA_complete)


#################################
#####Beta diversity analyses#####
#################################

#Function to perform NMDS and add coords to a dataframe:
make_mds_df <- function(comm, group, k = 2, dist = "bray") {
  nmds <- metaMDS(
    comm,
    distance = dist,
    k = k,
    try = 40,
    trymax = 400,
    autotransform = FALSE,
  )
  
  nmds$points %>%
    as.data.frame() %>%
    mutate(UniqueID = rownames(.)) %>%
    left_join(group, by = "UniqueID")
}            


#Raw counts, not including "Other"
allele_bray_raw_nmds_2d_df_b <- make_mds_df(allele_community_matrix_noNA,
                                            rpip_allele_group_data)
ggplot(allele_bray_raw_nmds_2d_df_b,
       aes(x = MDS1, y = MDS2, colour = Enrichment)) +
  geom_point()
#The odd grouping of the lower retentate sample (20040-RA-RPIP) can likely
#be attributed to low AMR read count and low diversity in that sample

#Proportions, not including "Other"
allele_bray_relative_nmds_2d_df <- make_mds_df(allele_relative,
                                                 rpip_allele_group_data)
ggplot(allele_bray_relative_nmds_2d_df,
       aes(x = MDS1, y = MDS2, colour = Enrichment)) +
  geom_point()


#Wisconsin, not including "Other"
allele_bray_wisconsin_nmds_2d_df <- make_mds_df(allele_wisconsin,
                                                rpip_allele_group_data)
ggplot(allele_bray_wisconsin_nmds_2d_df,
       aes(x = MDS1, y = MDS2, colour = Enrichment)) +
  geom_point()

#Jaccard, not including "Other" (shouldn't need to include since it is present in all samples)
allele_jacc_nmds_2d_df <- make_mds_df(allele_community_matrix_noNA,
                                      rpip_allele_group_data, dist = "jaccard")
ggplot(allele_jacc_nmds_2d_df,
       aes(x = MDS1, y = MDS2, colour = Enrichment)) +
  geom_point()

#Jaccard-type Chao distance, not including "Other" (shouldn't need to include since it is present in all samples)
allele_chao_nmds_2d_df <- make_mds_df(allele_community_matrix_noNA,
                                      rpip_allele_group_data, dist = "chao")
ggplot(allele_chao_nmds_2d_df,
       aes(x = MDS1, y = MDS2, colour = Enrichment)) +
  geom_point()


####################
#####PERMANOVA######
####################
#Proportion, not including "Other"
set.seed(123)
adonis2(allele_relative ~ Enrichment + Fraction + Nanotrap_type,
        data = rpip_allele_group_data, by = "margin", na.rm = TRUE)

#Wisconsin, not including "Other"
set.seed(123)
adonis2(allele_wisconsin ~ Enrichment + Fraction + Nanotrap_type,
        data = rpip_allele_group_data, by = "margin", na.rm = TRUE)


#Create presence/absence matrix for Jaccard distance:
allele_presabs_matrix <- allele_community_matrix_noNA
allele_presabs_matrix[allele_presabs_matrix > 0] <- 1
allele_presabs_matrix_complete <- allele_community_matrix_noNA_complete
allele_presabs_matrix_complete[allele_presabs_matrix_complete > 0] <- 1

#Jaccard
set.seed(123)
adonis2(allele_presabs_matrix ~ Enrichment + Fraction + Nanotrap_type,
        data = rpip_allele_group_data, by = "margin", method = "jaccard",
        permutations = 9999)

#Chao
#This is the only one with significant Nanotrap effect
set.seed(123)
adonis2(allele_community_matrix_noNA ~ Enrichment + Fraction + Nanotrap_type,
        data = rpip_allele_group_data, by = "margin", method = "chao")

##########################################
#####Restrict permutations by LIMS_ID#####
##########################################
library(permute)
rpip_perm_scheme <- how(within = Within(type = "free"),
                        plots = Plots(strata = rpip_allele_group_data$LIMS_ID,
                                      type = "none"),
                        nperm = 9999)

#Proportion, not including "Other"
set.seed(123)
adonis2(allele_relative ~ Fraction + Nanotrap_type + Enrichment,
        data = rpip_allele_group_data, by = "terms", na.rm = TRUE,
        permutations = rpip_perm_scheme)

#Wisconsin, not including "Other"
set.seed(123)
adonis2(allele_wisconsin ~ Fraction + Nanotrap_type + Enrichment,
        data = rpip_allele_group_data, by = "terms", na.rm = TRUE,
        permutations = rpip_perm_scheme)

#Jaccard, not including "Other"
set.seed(123)
adonis2(allele_presabs_matrix ~ Fraction + Nanotrap_type + Enrichment,
        data = rpip_allele_group_data, by = "terms", na.rm = TRUE,
        method = "jaccard",
        permutations = rpip_perm_scheme)

#Chao, not including "Other"
#Again the only distance measure with significant Nanotrap_type effect
set.seed(123)
adonis2(allele_community_matrix_noNA ~ Fraction + Nanotrap_type + Enrichment,
        data = rpip_allele_group_data, by = "terms", na.rm = TRUE,
        method = "chao",
        permutations = rpip_perm_scheme)

