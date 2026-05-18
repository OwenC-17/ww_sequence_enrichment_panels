library(tidyverse)
library(edgeR) #To be able to read DGE objects
library(vegan)

##########################
#####Read in the data#####
##########################

#Load RGI results (from rgi_ARG_import.R)
rpip_and_unt_gene_and_info <- read_rds(
  "input/modified/rpip_and_unt_rgi_gene_and_info.rds"
)  %>%
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

#Load the explanatory variables (group data) from Kraken analysis
#and modify for current data
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

#Remove VSP samples (irrelevant here)
rpip_gene_group_data <- group_data %>%
  filter(Enrichment %in% c("RPIP", "Non-targeted")) %>%
  left_join(rpip_libsizes)

#DGE object created in EdgeR_rgi_genelevel.R:
rpip_gene_DGElist <- read_rds("input/modified/edgeR_rgi_gene_DGE.rds")


#######################################
#####Make community count matrices#####
#######################################
#Get the genes from the DGE object created in EdgeR_rgi_genelevel.R
#Transpose for use in Vegan functions, and select only RPIP samples
gene_community_matrix_noNA <- (rpip_gene_DGElist$counts %>% 
                                 t())[rpip_gene_group_data$UniqueID, ]
gene_community_matrix_noNA <- gene_community_matrix_noNA[ ,
           colSums(gene_community_matrix_noNA, na.rm = TRUE) > 0]

#A version with NAs is useful for skipping pairs with one NA in some diversity
#calculations:
gene_community_matrix <- gene_community_matrix_noNA
gene_community_matrix[gene_community_matrix == 0] <- NA_real_


#Calculate number of non-ARG sequences per sample
Other <- rpip_only_group$libsize - rowSums(gene_community_matrix, na.rm = TRUE)
#Add column to community matrix
gene_community_matrix_complete <- cbind(gene_community_matrix, Other)
gene_community_matrix_noNA_complete <- cbind(gene_community_matrix_noNA, Other)

#Convert to relative abundance but don't include "Other" column
gene_relative <- gene_community_matrix_noNA / rpip_only_group$libsize

#And again, this time including "Other" column
gene_relative_complete <- gene_community_matrix_noNA_complete / 
  rpip_only_group$libsize

#Wisconsin transformation
gene_wisconsin <- wisconsin(gene_community_matrix_noNA)
gene_wisconsin_complete <- wisconsin(gene_community_matrix_noNA_complete)

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
    autotransform = FALSE
  )
  
  nmds$points %>%
    as.data.frame() %>%
    mutate(UniqueID = rownames(.)) %>%
    left_join(group, by = "UniqueID")
}            


#NMDS of raw ARG counts, not including "Other"
gene_bray_raw_nmds_2d_df_b <- make_mds_df(gene_community_matrix_noNA,
                                          rpip_gene_group_data)
ggplot(gene_bray_raw_nmds_2d_df_b,
       aes(x = MDS1, y = MDS2, colour = Enrichment)) +
  geom_point()


#Raw counts, including "Other"
#Including non-ARGs removes a lot of the variation, because ARGs make up 
#relatively few reads out of the whole sample.
gene_bray_raw_complete_nmds_2d_df_b <- make_mds_df(
  gene_community_matrix_noNA_complete, rpip_gene_group_data
  )
ggplot(gene_bray_raw_complete_nmds_2d_df_b,
       aes(x = MDS1, y = MDS2, colour = Fraction)) +
  geom_point()


#Proportions, not including "Other"
gene_bray_relative_nmds_2d_df <- make_mds_df(
  gene_relative, rpip_gene_group_data
  )
ggplot(gene_bray_relative_nmds_2d_df,
       aes(x = MDS1, y = MDS2, colour = Enrichment)) +
  geom_point()


#Wisconsin, not including "Other"
gene_bray_wisconsin_nmds_2d_df <- make_mds_df(
  gene_wisconsin, rpip_gene_group_data
  )
ggplot(gene_bray_wisconsin_nmds_2d_df,
       aes(x = MDS1, y = MDS2, colour = Enrichment)) +
  geom_point()

#Proportions, including "Other"
gene_bray_relative_complete_nmds_2d_df <- make_mds_df(
  gene_relative_complete, rpip_gene_group_data
  )
ggplot(gene_bray_relative_complete_nmds_2d_df,
       aes(x = MDS1, y = MDS2, colour = Enrichment)) +
  geom_point()

#Wisconsin, including "Other"
#This one is not so bad to look at
gene_bray_wisconsin_complete_nmds_2d_df <- make_mds_df(
  gene_wisconsin_complete, rpip_gene_group_data
  )
ggplot(gene_bray_wisconsin_complete_nmds_2d_df,
       aes(x = MDS1, y = MDS2, colour = Enrichment)) +
  geom_point()


#Jaccard, not including "Other" 
#(shouldn't need to include since "Other" is present in all samples)
#This one looks okay too
gene_jacc_nmds_2d_df <- make_mds_df(
  gene_community_matrix_noNA, rpip_gene_group_data, dist = "jaccard"
  )
ggplot(gene_jacc_nmds_2d_df,
       aes(x = MDS1, y = MDS2, colour = Enrichment)) +
  geom_point()

#Jaccard-type Chao distance, not including "Other"
#(shouldn't need to include since "Other" is present in all samples)
gene_chao_nmds_2d_df <- make_mds_df(
  gene_community_matrix_noNA, rpip_gene_group_data, dist = "chao"
  )
ggplot(gene_chao_nmds_2d_df,
       aes(x = MDS1, y = MDS2, colour = Enrichment)) +
  geom_point()


####################
#####PERMANOVA######
####################
#Proportion, not including "Other"
set.seed(123)
adonis2(gene_relative ~ Enrichment + Fraction + Nanotrap_type,
        data = rpip_gene_group_data, by = "margin", na.rm = TRUE)

#Wisconsin, not including "Other"
set.seed(123)
adonis2(gene_wisconsin ~ Enrichment + Fraction + Nanotrap_type,
        data = rpip_gene_group_data, by = "margin", na.rm = TRUE)


#Create presence/absence matrix for Jaccard distance:
gene_presabs_matrix <- gene_community_matrix_noNA
gene_presabs_matrix[gene_presabs_matrix > 0] <- 1
gene_presabs_matrix_complete <- gene_community_matrix_noNA_complete
gene_presabs_matrix_complete[gene_presabs_matrix_complete > 0] <- 1

#Jaccard
set.seed(123)
adonis2(gene_presabs_matrix ~ Enrichment + Fraction + Nanotrap_type,
        data = rpip_gene_group_data, by = "margin", method = "jaccard",
        permutations = 9999)

#Chao
#This is the only one with significant Nanotrap effect
set.seed(123)
adonis2(gene_community_matrix_noNA ~ Enrichment + Fraction + Nanotrap_type,
        data = rpip_gene_group_data, by = "margin", method = "chao")

##########################################
#####Restrict permutations by LIMS_ID#####
##########################################
library(permute)
rpip_perm_scheme <- how(within = Within(type = "free"),
                        plots = Plots(strata = rpip_gene_group_data$LIMS_ID,
                                      type = "none"),
                        nperm = 9999)

#Proportion, not including "Other"
set.seed(123)
adonis2(gene_relative ~ Fraction + Nanotrap_type + Enrichment,
        data = rpip_gene_group_data, by = "terms", na.rm = TRUE,
        permutations = rpip_perm_scheme)

#Wisconsin, not including "Other"
set.seed(123)
adonis2(gene_wisconsin ~ Fraction + Nanotrap_type + Enrichment,
        data = rpip_gene_group_data, by = "terms", na.rm = TRUE,
        permutations = rpip_perm_scheme)

#Jaccard, not including "Other"
set.seed(123)
adonis2(gene_presabs_matrix ~ Fraction + Nanotrap_type + Enrichment,
        data = rpip_gene_group_data, by = "terms", na.rm = TRUE,
        method = "jaccard",
        permutations = rpip_perm_scheme)

#Chao, not including "Other"
set.seed(123)
adonis2(gene_community_matrix_noNA ~ Fraction + Nanotrap_type + Enrichment,
        data = rpip_gene_group_data, by = "terms", na.rm = TRUE,
        method = "chao",
        permutations = rpip_perm_scheme)
