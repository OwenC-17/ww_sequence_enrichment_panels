library(tidyverse)
library(vegan)

#Import data
rpip_and_unt_allele_and_info <- read_rds("input/modified/rpip_and_unt_rgi_allele_and_info.rds")

#Get group data for RPIP only
rpip_only_group <- rpip_gene_group_data %>%
  filter(Enrichment == "RPIP") %>%
  rename("libsize" = "summary.after_dedup.total_reads")

#Get the alleles from the DGE object created in EdgeR_rgi_allelelevel.R
#Transpose for use in Vegan functions, and select only RPIP samples
allele_community_matrix_noNA <- (rpip_allele_DGElist$counts %>% t())[rpip_only_group$UniqueID,]
allele_community_matrix_noNA <- allele_community_matrix_noNA[, colSums(allele_community_matrix_noNA, na.rm = TRUE) > 0]
allele_community_matrix <- allele_community_matrix_noNA
allele_community_matrix[allele_community_matrix == 0] <- NA_real_


#Calculate number of non-ARG sequences per sample
Other <- rpip_only_group$libsize - rowSums(allele_community_matrix, na.rm = TRUE)
#Add column to community matrix
allele_community_matrix_complete <- cbind(allele_community_matrix, Other)
allele_community_matrix_noNA_complete <- cbind(allele_community_matrix_noNA, Other)

#Convert to relative abundance but don't include "Other" column
allele_relative <- allele_community_matrix / rpip_only_group$libsize
#And again, this time including "Other" column
allele_relative_complete <- allele_community_matrix_complete / rpip_only_group$libsize

#Wisconsin transformation
allele_wisconsin <- wisconsin(allele_community_matrix, na.rm = TRUE)

allele_wisconsin_complete <- wisconsin(allele_community_matrix_complete, na.rm = TRUE)


#Bray-Curtis distance
allele_bray_relative <- vegdist(allele_relative, method = "bray", na.rm = TRUE)
allele_bray_wisconsin <- vegdist(allele_wisconsin, method = "bray", na.rm = TRUE)

#NMDS with 3 axes
allele_bray_relative_nmds_3d <- metaMDS(
  allele_relative,
  distance = "bray",
  k = 3,
  try = 40,
  trymax = 150,
  autotransform = FALSE
)

#Get the NMDS coords:
allele_bray_relative_nmds_3d_df <- allele_bray_relative_nmds_3d$points %>%
  as.data.frame() %>%
  mutate(sample_id = rownames(.)) %>%
  left_join((group_data %>%
               mutate(sample_id =  paste(
                 LIMS_ID, Treatment, Enrichment, sep = "-"
               ))), by = "sample_id")

ggplot(allele_bray_relative_nmds_3d_df, aes(x = MDS2, y = MDS1, color = Fraction)) + geom_point()


#NMDS with 2 axes
allele_bray_relative_nmds_2d <- metaMDS(
  allele_relative,
  distance = "bray",
  k = 2,
  try = 40,
  trymax = 150,
  autotransform = FALSE
)

#Get NMDS coords:
allele_bray_relative_nmds_2d_df <- allele_bray_relative_nmds_2d$points %>%
  as.data.frame() %>%
  mutate(sample_id = rownames(.)) %>%
  left_join((group_data %>%
               mutate(sample_id =  paste(
                 LIMS_ID, Treatment, Enrichment, sep = "-"
               ))), by = "sample_id")

ggplot(allele_bray_relative_nmds_2d_df, aes(x = MDS1, y = MDS2, color = Fraction)) + geom_point() + geom_text(aes(label = sample_id))


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
allele_bray_raw_nmds_2d_df_b <- make_mds_df(allele_community_matrix, rpip_only_group)
ggplot(allele_bray_raw_nmds_2d_df_b, aes(x = MDS1, y = MDS2, colour = Fraction)) + geom_point()
#The odd grouping of the lower-left retentate sample (20040-RA-RPIP) can likely
#be attributed to low AMR read count and low diversity in that sample


#Raw counts, including "Other"
allele_bray_raw_complete_nmds_2d_df_b <- make_mds_df(allele_community_matrix_complete, rpip_only_group)
ggplot(allele_bray_raw_complete_nmds_2d_df_b, aes(x = MDS1, y = MDS2, colour = Fraction)) + geom_point()


#Proportions, not including "Other"
allele_bray_relative_nmds_2d_df_b <- make_mds_df(allele_relative, rpip_only_group)
ggplot(allele_bray_relative_nmds_2d_df_b, aes(x = MDS1, y = MDS2, colour = Fraction)) + geom_point()


#Wisconsin, not including "Other"
allele_bray_wisconsin_nmds_2d_df <- make_mds_df(allele_wisconsin, rpip_only_group)
ggplot(allele_bray_wisconsin_nmds_2d_df, aes(x = MDS1, y = MDS2, colour = Fraction)) + geom_point()

#Proportions, including "Other"
allele_bray_relative_complete_nmds_2d_df <- make_mds_df(allele_relative_complete, rpip_only_group)
ggplot(allele_bray_relative_nmds_2d_df, aes(x = MDS1, y = MDS2, colour = Fraction)) + geom_point()

#Wisconsin, including "Other"
allele_bray_wisconsin_complete_nmds_2d_df <- make_mds_df(allele_wisconsin_complete, rpip_only_group)
ggplot(allele_bray_wisconsin_nmds_2d_df, aes(x = MDS1, y = MDS2, colour = Fraction)) + geom_point()

#Jaccard, not including "Other" (shouldn't need to include since it is present in all samples)
allele_jacc_nmds_2d_df <- make_mds_df(allele_community_matrix, rpip_only_group, dist = "jaccard")
ggplot(allele_jacc_nmds_2d_df, aes(x = MDS1, y = MDS2, colour = Fraction, shape = site)) + geom_point(size = 5)

#Jaccard-type Chao distance, not including "Other" (shouldn't need to include since it is present in all samples)
allele_community_matrix_noNA <- allele_community_matrix
allele_community_matrix_noNA[is.na(allele_community_matrix_noNA)] <- 0
allele_chao_nmds_2d_df <- make_mds_df(allele_community_matrix_noNA, rpip_only_group, dist = "chao")
ggplot(allele_chao_nmds_2d_df, aes(x = MDS1, y = MDS2, colour = Fraction, shape = site)) + geom_point(size = 5)


####################
#####PERMANOVA######
####################

#Proportion, not including "Other"
adonis2(allele_relative ~ site + Fraction + Nanotrap_type, data = rpip_only_group, by = "margin", na.rm = TRUE)

#Proportion, including "Other"
adonis2(allele_relative_complete ~ site + Fraction + Nanotrap_type, data = rpip_only_group, by = "margin", na.rm = TRUE)

#Wisconsin, not including "Other"
set.seed(123)
adonis2(allele_wisconsin ~ site + Fraction + Nanotrap_type, data = rpip_only_group, by = "margin", na.rm = TRUE)

#Wisconsin, including "Other"
adonis2(allele_wisconsin_complete ~ site + Fraction + Nanotrap_type, data = rpip_only_group, by = "margin", na.rm = TRUE)


allele_presabs_matrix <- allele_community_matrix_noNA
allele_presabs_matrix[allele_presabs_matrix > 0] <- 1
allele_presabs_matrix_complete <- allele_community_matrix_noNA_complete
allele_presabs_matrix_complete[allele_presabs_matrix_complete > 0] <- 1

#Jaccard
adonis2(allele_presabs_matrix ~ site + Fraction + Nanotrap_type, data = rpip_only_group, by = "margin", method = "jaccard", permutations = 9999)
adonis2(allele_presabs_matrix_complete ~ site + Fraction + Nanotrap_type, data = rpip_only_group, by = "margin", method = "jaccard", permutations = 9999)


#Chao
adonis2(allele_community_matrix_noNA ~ site + Fraction + Nanotrap_type, data = rpip_only_group, by = "margin", method = "chao")


#####
#Alpha diversity
allele_relative_noNA <- allele_relative
allele_relative_noNA[is.na(allele_relative_noNA)] <- 0

allele_relative_noNA_complete <- allele_relative_complete
allele_relative_noNA_complete[is.na(allele_relative_noNA_complete)] <- 0


#Shannon by proportion, excluding "Other"
allele_shannon_RA <- diversity(allele_relative_noNA)
#Shannon prop including "Other"
allele_shannon_RA_complete <- diversity(allele_relative_noNA_complete)


allele_richness <- specnumber(allele_community_matrix_noNA)


rpip_only_group$richness <- specnumber(allele_community_matrix_noNA)

rpip_only_group <- rpip_only_group %>%
  mutate(relative_richness = richness / libsize) %>%
  mutate(shannon = diversity(allele_community_matrix_noNA, index = "shannon"),
         simpson = simpson.unb(allele_community_matrix_noNA),
         invsimpson = simpson.unb(allele_community_matrix_noNA, inverse = TRUE),
         fisher = fisher.alpha(allele_community_matrix_noNA),
         chao1 = estimateR(allele_community_matrix_noNA)["S.chao1",]) %>%
  mutate(Nanotrap_type = relevel(factor(Nanotrap_type, ordered = FALSE), ref = "A"),
         Fraction = relevel(factor(Fraction, ordered = FALSE), ref = "unfiltered"))

ggplot(rpip_only_group, aes(x = Treatment, y = richness, colour = site)) + geom_boxplot()
ggplot(rpip_only_group, aes(x = Treatment, y = shannon, colour = site)) + geom_boxplot()
ggplot(rpip_only_group, aes(x = Treatment, y = simpson, colour = site)) + geom_boxplot()
ggplot(rpip_only_group, aes(x = Treatment, y = invsimpson, colour = site)) + geom_boxplot()
ggplot(rpip_only_group, aes(x = Treatment, y = fisher, colour = site)) + geom_boxplot()
ggplot(rpip_only_group, aes(x = Treatment, y = chao1, colour = site)) + geom_boxplot()


library(spm2)
library(lme4)
library(flexplot)

simpson_model <- lmer(simpson ~ 1 + Nanotrap_type + Fraction + (1 + Nanotrap_type + Fraction | LIMS_ID), data = rpip_only_group)
summary(simpson_model)
visualize(simpson_model, sample = 8)


invsimpson_model <- lmer(invsimpson ~ 1 + Nanotrap_type + Fraction + (1 + Nanotrap_type + Fraction | LIMS_ID), data = rpip_only_group)
summary(invsimpson_model)
visualize(invsimpson_model, sample = 8)

baseline_richness_me <- lmer(richness ~ 1 + (1 | LIMS_ID), data = rpip_only_group)
summary(baseline_richness_me)

richness_model_me_r <- lmer(richness ~ 1 + Nanotrap_type + Fraction + (1 + Nanotrap_type + Fraction | LIMS_ID), data = rpip_only_group)
summary(richness_model_me_r)
estimates(richness_model_me_r)
visualize(richness_model_me_r, plot = "model", formula = richness ~ Fraction + LIMS_ID, sample = 8)

richness_model_me_n <- lmer(richness ~ 1 + Nanotrap_type + Fraction + (1 | LIMS_ID), data = rpip_only_group)
summary(richness_model_me_n)
estimates(richness_model_me_n)
visualize(richness_model_me_n, plot = "model", formula = richness ~ Fraction + LIMS_ID, sample = 8)

visualize(richness_model_me_r, sample = 8)
visualize(richness_model_me_n, sample = 8)


compare.fits(richness ~ Nanotrap_type | LIMS_ID,
             data = rpip_only_group,
             model1 = richness_model_me_n,
             model2 = richness_model_me_r,
             re = T)

model.comparison(richness_model_me_r, richness_model_me_n)

richness_model_me 

richness_model_glme <- glmer(richness ~ 1 + Nanotrap_type + Fraction + (1 | LIMS_ID), data = rpip_only_group, 
                             family = "poisson")
estimates(richness_model_glme)
summary(richness_model_glme)
visualize(richness_model_glme, sample = 8)
visualize(richness_model_glme, sample = 8, plot = "residuals")

richness_model_binom_glme <- glmer(richness/libsize ~ 1 + Nanotrap_type + Fraction + (1 | LIMS_ID), data = rpip_only_group, 
                             family = "poisson")

richness_model_glme_r <- glmer(richness ~ 1 + Nanotrap_type + Fraction + (1 + Nanotrap_type + Fraction | LIMS_ID), data = rpip_only_group, 
                             family = "poisson")
estimates(richness_model_glme_r)
summary(richness_model_glme_r)
visualize(richness_model_glme_r, sample = 8)
visualize(richness_model_glme_r, sample = 8, plot = "residuals")


library(glmmTMB)
simpson_model_glme_r <- glmmTMB(simpson ~ Nanotrap_type + Fraction + (1 | LIMS_ID), data = rpip_only_group, family = beta_family())
estimates(simpson_model_glme_r)
summary(simpson_model_glme_r)
visualize(simpson_model_glme_r, sample = 8)
visualize(richness_model_glme_r, sample = 8, plot = "residuals")



compare.fits(richness ~ Fraction | Nanotrap_type + LIMS_ID,
             data = rpip_only_group,
             model1 = richness_model_glme,
             model2 = richness_model_glme_r,
             re = T, 
             clusters = 8)

model.comparison(richness_model_glme, richness_model_glme_r)

richness_model <- glm(richness ~ Nanotrap_type + Fraction + site, data = rpip_only_group, family = "poisson")
summary(richness_model)
anova(richness_model)


shannon_model <- lm(shannon ~ Nanotrap_type + Fraction + site, data = rpip_only_group)
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


richness_model_binom_glme <- glmer(richness/libsize ~ 1 + Nanotrap_type + Fraction + (1 | LIMS_ID), data = rpip_only_group, 
                                   family = "binomial")
richness_model_binom_glme_r <- glmer(richness/libsize ~ 1 + Nanotrap_type + Fraction + (1 + Nanotrap_type + Fraction | LIMS_ID), data = rpip_only_group, 
                               family = "binomial")
summary(richness_model_binom_glme_r)

