library(tidyverse)
library(vegan)

#Import data
rpip_and_unt_allele_and_info <- read_rds("input/modified/rpip_and_unt_rgi_allele_and_info.rds")
rpip_and_unt_gene_and_info <- read_rds("input/modified/rpip_and_unt_rgi_gene_and_info.rds")

#Get group data for RPIP only
rpip_only_group <- rpip_gene_group_data %>%
#  filter(Enrichment == "RPIP") %>%
  rename("libsize" = "summary.after_dedup.total_reads")

#Get the genes from the DGE object created in EdgeR_rgi_genelevel.R
#Transpose for use in Vegan functions, and select only RPIP samples
gene_community_matrix_noNA <- (rpip_gene_DGElist$counts %>% t())[rpip_only_group$UniqueID,]
gene_community_matrix_noNA <- gene_community_matrix_noNA[, colSums(gene_community_matrix_noNA, na.rm = TRUE) > 0]
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
gene_relative_complete <- gene_community_matrix_noNA_complete / rpip_only_group$libsize

#Wisconsin transformation
gene_wisconsin <- wisconsin(gene_community_matrix_noNA)

gene_wisconsin_complete <- wisconsin(gene_community_matrix_noNA_complete)


#Bray-Curtis distance
gene_bray_relative <- vegdist(gene_relative, method = "bray", na.rm = TRUE)
gene_bray_wisconsin <- vegdist(gene_wisconsin, method = "bray", na.rm = TRUE)

#NMDS with 3 axes
gene_bray_relative_nmds_3d <- metaMDS(
  gene_relative,
  distance = "bray",
  k = 3,
  try = 40,
  trymax = 150,
  autotransform = FALSE
)

#Get the NMDS coords:
gene_bray_relative_nmds_3d_df <- gene_bray_relative_nmds_3d$points %>%
  as.data.frame() %>%
  mutate(sample_id = rownames(.)) %>%
  left_join((group_data %>%
               mutate(sample_id =  paste(
                 LIMS_ID, Treatment, Enrichment, sep = "-"
               ))), by = "sample_id")

ggplot(gene_bray_relative_nmds_3d_df, aes(x = MDS2, y = MDS1, color = Fraction)) + geom_point()


#NMDS with 2 axes
gene_bray_relative_nmds_2d <- metaMDS(
  gene_relative,
  na.rm = TRUE,
  distance = "bray",
  k = 2,
  try = 40,
  trymax = 150,
  autotransform = FALSE
)

#Get NMDS coords:
gene_bray_relative_nmds_2d_df <- gene_bray_relative_nmds_2d$points %>%
  as.data.frame() %>%
  mutate(sample_id = rownames(.)) %>%
  left_join((group_data %>%
               mutate(sample_id =  paste(
                 LIMS_ID, Treatment, Enrichment, sep = "-"
               ))), by = "sample_id")

ggplot(gene_bray_relative_nmds_2d_df, aes(x = MDS1, y = MDS2, color = Fraction)) + geom_point() + geom_text(aes(label = sample_id))


#Function to make more NMDS's without typing so much:
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


#Raw counts, not including "Other"
gene_bray_raw_nmds_2d_df_b <- make_mds_df(gene_community_matrix_noNA, rpip_only_group)
ggplot(gene_bray_raw_nmds_2d_df_b, aes(x = MDS1, y = MDS2, colour = Fraction)) + geom_point()
#The odd grouping of the lower-left retentate sample (20040-RA-RPIP) can likely
#be attributed to low AMR read count and low diversity in that sample


#Raw counts, including "Other"
gene_bray_raw_complete_nmds_2d_df_b <- make_mds_df(gene_community_matrix_noNA_complete, rpip_only_group)
ggplot(gene_bray_raw_complete_nmds_2d_df_b, aes(x = MDS1, y = MDS2, colour = Fraction)) + geom_point()


#Proportions, not including "Other"
gene_bray_relative_nmds_2d_df_b <- make_mds_df(gene_relative, rpip_only_group)
ggplot(gene_bray_relative_nmds_2d_df_b, aes(x = MDS1, y = MDS2, colour = Fraction)) + geom_point()


#Wisconsin, not including "Other"
gene_bray_wisconsin_nmds_2d_df <- make_mds_df(gene_wisconsin, rpip_only_group)
ggplot(gene_bray_wisconsin_nmds_2d_df, aes(x = MDS1, y = MDS2, colour = Fraction)) + geom_point()

#Proportions, including "Other"
gene_bray_relative_complete_nmds_2d_df <- make_mds_df(gene_relative_complete, rpip_only_group)
ggplot(gene_bray_relative_nmds_2d_df, aes(x = MDS1, y = MDS2, colour = Fraction)) + geom_point()

#Wisconsin, including "Other"
gene_bray_wisconsin_complete_nmds_2d_df <- make_mds_df(gene_wisconsin_complete, rpip_only_group)
ggplot(gene_bray_wisconsin_nmds_2d_df, aes(x = MDS1, y = MDS2, colour = Fraction)) + geom_point()

#Jaccard, not including "Other" (shouldn't need to include since it is present in all samples)
gene_jacc_nmds_2d_df <- make_mds_df(gene_community_matrix_noNA, rpip_only_group, dist = "jaccard")
ggplot(gene_jacc_nmds_2d_df, aes(x = MDS1, y = MDS2, colour = Fraction, shape = site)) + geom_point(size = 5)

#Jaccard-type Chao distance, not including "Other" (shouldn't need to include since it is present in all samples)
gene_community_matrix_noNA <- gene_community_matrix
gene_community_matrix_noNA[is.na(gene_community_matrix_noNA)] <- 0
gene_chao_nmds_2d_df <- make_mds_df(gene_community_matrix_noNA, rpip_only_group, dist = "chao")
ggplot(gene_chao_nmds_2d_df, aes(x = MDS1, y = MDS2, colour = Fraction, shape = site)) + geom_point(size = 5)


####################
#####PERMANOVA######
####################

#Proportion, not including "Other"
adonis2(gene_relative ~ site + Fraction + Nanotrap_type, data = rpip_only_group, by = "margin", na.rm = TRUE)

#Proportion, including "Other"
adonis2(gene_relative_complete ~ site + Fraction + Nanotrap_type, data = rpip_only_group, by = "margin", na.rm = TRUE)

#Wisconsin, not including "Other"
set.seed(123)
adonis2(gene_wisconsin ~ site + Fraction + Nanotrap_type, data = rpip_only_group, by = "margin", na.rm = TRUE)

#Wisconsin, including "Other"
adonis2(gene_wisconsin_complete ~ site + Fraction + Nanotrap_type, data = rpip_only_group, by = "margin", na.rm = TRUE)


gene_presabs_matrix <- gene_community_matrix_noNA
gene_presabs_matrix[gene_presabs_matrix > 0] <- 1
gene_presabs_matrix_complete <- gene_community_matrix_noNA_complete
gene_presabs_matrix_complete[gene_presabs_matrix_complete > 0] <- 1

#Jaccard
adonis2(gene_presabs_matrix ~ site + Fraction + Nanotrap_type, data = rpip_only_group, by = "margin", method = "jaccard", permutations = 9999)
adonis2(gene_presabs_matrix_complete ~ site + Fraction + Nanotrap_type, data = rpip_only_group, by = "margin", method = "jaccard", permutations = 9999)


#Chao
adonis2(gene_community_matrix_noNA ~ site + Fraction + Nanotrap_type, data = rpip_only_group, by = "margin", method = "chao")


#####
#Alpha diversity
gene_relative_noNA <- gene_relative
gene_relative_noNA[is.na(gene_relative_noNA)] <- 0

gene_relative_noNA_complete <- gene_relative_complete
gene_relative_noNA_complete[is.na(gene_relative_noNA_complete)] <- 0


#Shannon by proportion, excluding "Other"
gene_shannon_RA <- diversity(gene_relative_noNA)
#Shannon prop including "Other"
gene_shannon_RA_complete <- diversity(gene_relative_noNA_complete)


gene_richness <- specnumber(gene_community_matrix_noNA)


rpip_only_group$richness <- specnumber(gene_community_matrix_noNA)

rpip_only_group <- rpip_only_group %>%
  mutate(relative_richness = richness / libsize) %>%
  mutate(shannon = diversity(gene_community_matrix_noNA, index = "shannon"),
         simpson = simpson.unb(gene_community_matrix_noNA),
         invsimpson = simpson.unb(gene_community_matrix_noNA, inverse = TRUE),
         fisher = fisher.alpha(gene_community_matrix_noNA),
         chao1 = estimateR(gene_community_matrix_noNA)["S.chao1",]) %>%
  mutate(Nanotrap_type = relevel(factor(Nanotrap_type, ordered = FALSE), ref = "none"),
         Fraction = relevel(factor(Fraction, ordered = FALSE), ref = "unfiltered"))

ggplot(rpip_only_group, aes(x = Treatment, y = richness, colour = site)) + geom_boxplot()
ggplot(rpip_only_group, aes(x = Treatment, y = shannon, colour = site)) + geom_boxplot()
ggplot(rpip_only_group, aes(x = Treatment, y = simpson, colour = site)) + geom_boxplot()
ggplot(rpip_only_group, aes(x = Treatment, y = invsimpson, colour = site)) + geom_boxplot()
ggplot(rpip_only_group, aes(x = Treatment, y = fisher, colour = site)) + geom_boxplot()
ggplot(rpip_only_group, aes(x = Treatment, y = chao1, colour = site)) + geom_boxplot()


library(spm2)

simpson_model <- glm(simpson ~ Nanotrap_type + Fraction + site, data = rpip_only_group)
summary(simpson_model)
richness_model <- glm(richness ~ Nanotrap_type + Fraction, data = rpip_only_group, family = "poisson")
summary(richness_model)
anova(richness_model)
shannon_model <- lm(shannon ~ Nanotrap_type + Fraction, data = rpip_only_group)
summary(shannon_model)


richness_model_cv <- glmcv(richness ~ Nanotrap_type + Fraction, rpip_only_group,rpip_only_group$richness, family = "poisson", validation = "LOO", cv.fold = 500)

library(randomForest)
set.seed(1234)
train <- sample(1:48, size = 36)
richness_model_rf <- randomForest(
  formula = richness ~ Nanotrap_type + Fraction + site,
  data = rpip_only_group[train,]
)

test <- rpip_only_group[-train, ]
test$pred <- predict(richness_model_rf, test)

ggplot(test, aes(x = pred, y = richness)) + geom_point()

summary(richness_model_rf)

ggplot(rpip_only_group, aes(x = log(simpson / (1 - simpson)))) + geom_density()


ggplot(rpip_only_group, aes(x = richness, colour = Treatment)) + geom_density(alpha = 0.5)

ggplot(rpip_only_group, aes(x = shannon)) + geom_density(alpha = 0.5)




##############GLS
library(nlme)
rpip_group <- rpip_only_group %>%
  mutate(LIMS_ID = factor(LIMS_ID), Treatment = factor(Treatment), site = factor(site), Fraction = factor(Fraction), Nanotrap_type = factor(Nanotrap_type), Enrichment = factor(Enrichment))


basic_richness_gls <- gls(relative_richness ~ Enrichment + Fraction + Nanotrap_type +
              Enrichment:Fraction +
              Enrichment:Nanotrap_type +
              Fraction:Nanotrap_type +
              Enrichment:Fraction:Nanotrap_type,
            weights = varIdent(form = ~ 1 | Enrichment * Fraction * Nanotrap_type),
            method = "ML", data = rpip_group)

anova(basic_richness_gls)
