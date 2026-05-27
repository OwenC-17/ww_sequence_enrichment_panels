library(tidyverse)
library(edgeR)

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
rpip_bacteria_relative <- rpip_bacteria_community_matrix_noNA / rpip_group_data$libsize
#And again, this time including "Other" column
rpip_bacteria_relative_complete <- rpip_bacteria_community_matrix_noNA_complete / rpip_group_data$libsize

#Wisconsin transformation
rpip_bacteria_wisconsin <- wisconsin(rpip_bacteria_community_matrix_noNA, na.rm = TRUE)

rpip_bacteria_wisconsin_complete <- wisconsin(rpip_bacteria_community_matrix_noNA_complete, na.rm = TRUE)


#Bray-Curtis distance
rpip_bacteria_bray_relative <- vegdist(rpip_bacteria_relative, method = "bray", na.rm = TRUE)
rpip_bacteria_bray_wisconsin <- vegdist(rpip_bacteria_wisconsin, method = "bray", na.rm = TRUE)

#NMDS with 3 axes
rpip_bacteria_bray_relative_nmds_3d <- metaMDS(
  rpip_bacteria_relative,
  distance = "bray",
  k = 3,
  try = 40,
  trymax = 1000,
  autotransform = FALSE
)

#Get the NMDS coords:
rpip_bacteria_bray_relative_nmds_3d_df <- rpip_bacteria_bray_relative_nmds_3d$points %>%
  as.data.frame() %>%
  mutate(sample_id = rownames(.)) %>%
  left_join((rpip_group_data %>%
               mutate(sample_id =  paste(
                 LIMS_ID, Treatment, Enrichment, sep = "-"
               ))), by = "sample_id")

ggplot(rpip_bacteria_bray_relative_nmds_3d_df, aes(x = MDS2, y = MDS1, shape = Enrichment, color = Nanotrap_type)) + geom_point()


#NMDS with 2 axes
rpip_bacteria_bray_relative_nmds_2d <- metaMDS(
  rpip_bacteria_relative[rowSums(rpip_bacteria_relative) != 0, ],
  distance = "bray",
  k = 2,
  try = 40,
  trymax = 150,
  autotransform = FALSE
)

#Get NMDS coords:
rpip_bacteria_bray_relative_nmds_2d_df <- rpip_bacteria_bray_relative_nmds_2d$points %>%
  as.data.frame() %>%
  mutate(sample_id = rownames(.)) %>%
  left_join((rpip_group_data %>%
               mutate(sample_id =  paste(
                 LIMS_ID, Treatment, Enrichment, sep = "-"
               ))), by = "sample_id")

ggplot(rpip_bacteria_bray_relative_nmds_2d_df, aes(x = MDS1, y = MDS2, color = Fraction, shape = Enrichment)) + geom_point()


#Function to make more NMDS's without typing so much:
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
rpip_bacteria_bray_raw_nmds_2d_df_b <- make_mds_df_rpip(rpip_bacteria_community_matrix_noNA, rpip_group_data)
ggplot(rpip_bacteria_bray_raw_nmds_2d_df_b, aes(x = MDS1, y = MDS2, colour = Treatment, shape = Enrichment)) + geom_point()


#Raw counts, including "Other"
rpip_bacteria_bray_raw_complete_nmds_2d_df_b <- make_mds_df_rpip(rpip_bacteria_community_matrix_noNA_complete, rpip_group_data)
ggplot(rpip_bacteria_bray_raw_complete_nmds_2d_df_b, aes(x = MDS1, y = MDS2, colour = Fraction, shape = Enrichment)) + geom_point()


#Proportions, not including "Other"
rpip_bacteria_bray_relative_nmds_2d_df_b <- make_mds_df_rpip(rpip_bacteria_relative, rpip_group_data)
ggplot(rpip_bacteria_bray_relative_nmds_2d_df_b, aes(x = MDS1, y = MDS2, colour = Fraction, shape = Enrichment)) + geom_point(size = 2) + geom_line(aes(group = paste(LIMS_ID, Treatment)), linetype = "dotted") +
  scale_colour_paletteer_d("ggsci::default_igv") +
  theme_bw()


#Wisconsin, not including "Other"
rpip_bacteria_bray_wisconsin_nmds_2d_df <- make_mds_df_rpip(rpip_bacteria_wisconsin, rpip_group_data)
ggplot(rpip_bacteria_bray_wisconsin_nmds_2d_df, aes(x = MDS1, y = MDS2, colour = Fraction, shape = Enrichment)) + geom_point(size = 2) + geom_line(aes(group = paste(LIMS_ID, Treatment)), linetype = "dotted") +
  scale_colour_paletteer_d("ggsci::default_igv") +
  theme_bw()


#Proportions, including "Other"
rpip_bacteria_bray_relative_complete_nmds_2d_df <- make_mds_df_rpip(rpip_bacteria_relative_complete, rpip_group_data)
ggplot(rpip_bacteria_bray_relative_complete_nmds_2d_df, aes(x = MDS1, y = MDS2, colour = Fraction, shape = Enrichment)) + geom_point(size = 2) + geom_line(aes(group = paste(LIMS_ID, Treatment)), linetype = "dotted") +
  scale_colour_paletteer_d("ggsci::default_igv") +
  theme_bw()

#Wisconsin, including "Other"
rpip_bacteria_bray_wisconsin_complete_nmds_2d_df <- make_mds_df_rpip(rpip_bacteria_wisconsin_complete, rpip_group_data)
ggplot(rpip_bacteria_bray_wisconsin_complete_nmds_2d_df, aes(x = MDS1, y = MDS2, colour = Fraction, shape = Enrichment)) + geom_point(size = 2) + geom_line(aes(group = paste(LIMS_ID, Treatment)), linetype = "dotted") +
  scale_colour_paletteer_d("ggsci::default_igv") +
  theme_bw()

#Jaccard, not including "Other" (shouldn't need to include since it is present in all samples)
rpip_bacteria_jacc_nmds_2d_df <- make_mds_df_rpip(rpip_bacteria_community_matrix_noNA, rpip_group_data, dist = "jaccard")
ggplot(rpip_bacteria_jacc_nmds_2d_df, aes(x = MDS1, y = MDS2, colour = Fraction, shape = Enrichment)) + geom_point(size = 2) + geom_line(aes(group = paste(LIMS_ID, Treatment)), linetype = "dotted") +
  scale_colour_paletteer_d("ggsci::default_igv") +
  theme_bw()

#Jaccard-type Chao distance, not including "Other" (shouldn't need to include since it is present in all samples)
rpip_bacteria_chao_nmds_2d_df <- make_mds_df_rpip(rpip_bacteria_community_matrix_noNA, rpip_group_data, dist = "chao")
ggplot(rpip_bacteria_chao_nmds_2d_df, aes(x = MDS1, y = MDS2, colour = Fraction, shape = Enrichment)) + geom_point(size = 2) + geom_line(aes(group = paste(LIMS_ID, Treatment)), linetype = "dotted") +
  scale_colour_paletteer_d("ggsci::default_igv") +
  theme_bw()



####################
#####PERMANOVA######
####################

#Proportion, not including "Other"
adonis2(rpip_bacteria_relative[rowSums(rpip_bacteria_relative, na.rm = TRUE) != 0,] ~ site + Fraction + Nanotrap_type + Enrichment, data = rpip_group_data[rowSums(rpip_bacteria_relative, na.rm = TRUE) != 0,], by = "margin", na.rm = TRUE)


#Proportion, including "Other"
adonis2(rpip_bacteria_relative_complete ~ site + Fraction + Nanotrap_type + Enrichment, data = rpip_group_data, by = "margin", na.rm = TRUE)

#Wisconsin, not including "Other"
set.seed(123)
adonis2(rpip_bacteria_wisconsin[rowSums(rpip_bacteria_wisconsin, na.rm = TRUE) != 0,] ~ site + Fraction + Nanotrap_type + Enrichment, data = rpip_group_data[rowSums(rpip_bacteria_wisconsin, na.rm = TRUE) != 0,], by = "margin", na.rm = TRUE)

#Wisconsin, including "Other"
adonis2(rpip_bacteria_wisconsin_complete ~ site + Fraction + Nanotrap_type + Enrichment, data = rpip_group_data, by = "margin", na.rm = TRUE)


rpip_bacteria_presabs_matrix <- rpip_bacteria_community_matrix_noNA
rpip_bacteria_presabs_matrix[rpip_bacteria_presabs_matrix > 0] <- 1
rpip_bacteria_presabs_matrix_complete <- rpip_bacteria_community_matrix_noNA_complete
rpip_bacteria_presabs_matrix_complete[rpip_bacteria_presabs_matrix_complete > 0] <- 1

#Jaccard
adonis2(rpip_bacteria_presabs_matrix[rowSums(rpip_bacteria_presabs_matrix, na.rm = TRUE) != 0,] ~ site + Fraction + Nanotrap_type + Enrichment, data = rpip_group_data[rowSums(rpip_bacteria_presabs_matrix, na.rm = TRUE) != 0,], by = "margin", method = "jaccard", permutations = 9999)
adonis2(rpip_bacteria_presabs_matrix_complete ~ site + Fraction + Nanotrap_type + Enrichment, data = rpip_group_data, by = "margin", method = "jaccard", permutations = 9999)


#Chao
adonis2(rpip_bacteria_community_matrix_noNA[rowSums(rpip_bacteria_community_matrix_noNA, na.rm = TRUE) != 0,] ~ site + Fraction + Nanotrap_type + Enrichment, data = rpip_group_data[rowSums(rpip_bacteria_community_matrix_noNA, na.rm = TRUE) != 0,], by = "margin", method = "chao")
adonis2(rpip_bacteria_community_matrix_noNA_complete ~ site + Fraction + Nanotrap_type + Enrichment, data = rpip_group_data, by = "margin", method = "chao")


#####
#Alpha diversity
rpip_bacteria_relative_noNA <- rpip_bacteria_relative
rpip_bacteria_relative_noNA[is.na(rpip_bacteria_relative_noNA)] <- 0

rpip_bacteria_relative_noNA_complete <- rpip_bacteria_relative_complete
rpip_bacteria_relative_noNA_complete[is.na(rpip_bacteria_relative_noNA_complete)] <- 0


#Shannon by proportion, excluding "Other"
rpip_bacteria_shannon_RA <- diversity(rpip_bacteria_relative_noNA)
#Shannon prop including "Other"
rpip_bacteria_shannon_RA_complete <- diversity(rpip_bacteria_relative_noNA_complete)


rpip_bacteria_richness <- specnumber(rpip_bacteria_community_matrix_noNA)


rpip_group_data$richness <- specnumber(rpip_bacteria_community_matrix_noNA)
library(hillR)
rpip_group_data <- rpip_group_data %>%
  mutate(relative_richness = richness / libsize) %>%
  mutate(shannon = diversity(rpip_bacteria_community_matrix_noNA, index = "shannon"),
         simpson = simpson.unb(rpip_bacteria_community_matrix_noNA),
         simpson_complete = simpson.unb(rpip_bacteria_community_matrix_noNA_complete),
         invsimpson = simpson.unb(rpip_bacteria_community_matrix_noNA, inverse = TRUE),
         invsimpson_complete = simpson.unb(rpip_bacteria_community_matrix_noNA_complete, inverse = TRUE),
         gini_simpson = 1 - abs(simpson),
         gini_simpson_complete = 1 - simpson_complete,
         fisher = fisher.alpha(rpip_bacteria_community_matrix_noNA),
         chao1 = estimateR(rpip_bacteria_community_matrix_noNA)["S.chao1",],
         hill3 = hill_taxa(rpip_bacteria_community_matrix_noNA, q = 2)) %>%
  mutate(Nanotrap_type = relevel(factor(Nanotrap_type, ordered = FALSE), ref = "A"),
         Fraction = relevel(factor(Fraction, ordered = FALSE), ref = "unfiltered"),
         LIMS_ID = factor(LIMS_ID))

#rpip_group_data[rpip_group_data$simpson == -1,]$simpson <- 0
#rpip_group_data[rpip_group_data$gini_simpson == 1,]$gini_simpson <- 0
#rpip_group_data[rpip_group_data$gini_simpson_complete == 1,]$gini_simpson_complete <- 0


ggplot(rpip_group_data, aes(x = Treatment, y = richness, colour = site, fill = Enrichment)) + geom_boxplot()
ggplot(rpip_group_data, aes(x = Treatment, y = shannon, colour = site)) + geom_boxplot()
ggplot(rpip_group_data, aes(x = Treatment, y = simpson, colour = site)) + geom_boxplot()
ggplot(rpip_group_data, aes(x = Treatment, y = invsimpson, colour = site)) + geom_boxplot()
ggplot(rpip_group_data, aes(x = Treatment, y = fisher, colour = site)) + geom_boxplot()
ggplot(rpip_group_data, aes(x = Treatment, y = chao1, colour = site)) + geom_boxplot()


library(spm2)
library(lme4)
library(flexplot)


library(glmmTMB)
simpson_model_glme_r <- glmmTMB(simpson ~ Nanotrap_type + Fraction + (1 | LIMS_ID), data = rpip_group_data, family = beta_family())
estimates(simpson_model_glme_r)
summary(simpson_model_glme_r)
visualize(simpson_model_glme_r, sample = 8)
visualize(richness_model_glme_r, sample = 8, plot = "residuals")



compare.fits(richness ~ Fraction | Nanotrap_type + LIMS_ID,
             data = rpip_group_data,
             model1 = richness_model_glme,
             model2 = richness_model_glme_r,
             re = T, 
             clusters = 8)

model.comparison(richness_model_glme, richness_model_glme_r)

richness_model <- glm(richness ~ Nanotrap_type + Fraction + site, data = rpip_group_data, family = "poisson")
summary(richness_model)
anova(richness_model)


shannon_model <- lm(shannon ~ Nanotrap_type + Fraction + site, data = rpip_group_data)
summary(shannon_model)


richness_model_cv <- glmcv(richness ~ Nanotrap_type + Fraction, rpip_group_data,rpip_group_data$richness, family = "poisson", validation = "LOO", cv.fold = 500)

library(randomForest)
set.seed(1234)
train <- sample(1:48, size = 36)
richness_model_rf <- randomForest(
  formula = richness ~ Nanotrap_type + Fraction + site,
  data = rpip_group_data[train,]
)

test <- rpip_group_data[-train, ]
test$pred <- predict(richness_model_rf, test)

ggplot(test, aes(x = pred, y = richness)) + geom_point()

summary(richness_model_rf)

ggplot(rpip_group_data, aes(x = log(simpson / (1 - simpson)))) + geom_density()


ggplot(rpip_group_data, aes(x = richness, colour = Treatment)) + geom_density(alpha = 0.5)

ggplot(rpip_group_data, aes(x = shannon)) + geom_density(alpha = 0.5)



#Skillings-Mack
sm_rpip_group_data <- expand.grid(unique(rpip_group_data$Treatment), unique(rpip_group_data$LIMS_ID)) %>%
  rename(Treatment = Var1, LIMS_ID = Var2) %>%
  left_join(select(rpip_group_data, richness, Treatment, LIMS_ID))

Ski.Mack(sm_rpip_group_data$richness, groups = sm_rpip_group_data$Treatment, blocks = sm_rpip_group_data$LIMS_ID, simulate.p.value = TRUE)

library(PMCMRplus)
richness_smt <- skillingsMackTest(sm_rpip_group_data$richness, groups = sm_rpip_group_data$Treatment, blocks = sm_rpip_group_data$LIMS_ID, simulate.p.value = TRUE)
summary(richness_smt)

rpip_group_data$Enrichment <- factor(rpip_group_data$Enrichment)

library(ARTool)
complete_model_data <- expand.grid(unique(rpip_group_data$Fraction),
                                   unique(rpip_group_data$Nanotrap_type),
                                   unique(rpip_group_data$Enrichment),
                                   unique(rpip_group_data$LIMS_ID)) %>%
  rename(Fraction = Var1, 
         Nanotrap_type = Var2,
         Enrichment = Var3,
         LIMS_ID = Var4) %>%
  left_join(rpip_group_data)
armodel <- art(invsimpson ~ 1 + Nanotrap_type*Fraction*Enrichment + (1 + Enrichment | LIMS_ID), data = rpip_group_data)

summary(armodel)

anova(armodel)

art.con(armodel, "Nanotrap_type")



library(fitdistrplus)
library(glmmTMB)
library(DHARMa)
library(emmeans)

mean(rpip_group_data$gini_simpson)
var(rpip_group_data$gini_simpson) #Variance is < mean

mean(rpip_group_data$gini_simpson_complete)
var(rpip_group_data$gini_simpson_complete)


descdist(rpip_group_data$gini_simpson, boot = 1000)
descdist(rpip_group_data$fisher, boot = 1000)
descdist(rpip_group_data$chao1, boot = 1000)
descdist(rpip_group_data$simpson, boot = 1000)
descdist(rpip_group_data$shannon, boot = 1000)

descdist(rpip_group_data[rpip_group_data$Enrichment == "RPIP",]$invsimpson, boot = 1000)
descdist(rpip_group_data[rpip_group_data$Enrichment == "RPIP",]$gini_simpson, boot = 1000)
descdist(rpip_group_data[rpip_group_data$Enrichment == "RPIP",]$shannon, boot = 1000)
descdist(rpip_group_data[rpip_group_data$Enrichment == "RPIP",]$fisher, boot = 1000)
descdist(rpip_group_data[rpip_group_data$Enrichment == "RPIP",]$chao1, boot = 1000)
descdist(rpip_group_data[rpip_group_data$Enrichment == "RPIP",]$richness, boot = 1000)
descdist(rpip_group_data[rpip_group_data$Enrichment == "RPIP",]$relative_richness, boot = 1000)


shannon_gamma_model <- glmmTMB(shannon ~ Nanotrap_type*Enrichment + Fraction*Enrichment +
                                 Nanotrap_type * Fraction + (1 | LIMS_ID) + (1 | LIMS_ID:Nanotrap_type:Fraction),
                               family = Gamma,
                               data = rpip_group_data)
AIC(shannon_gamma_model)

sim_resid <- simulateResiduals(shannon_gamma_model, n = 1000)
plot(sim_resid)
testDispersion(sim_resid)

emmeans(shannon_gamma_model, pairwise ~ Nanotrap_type, adjust = "fdr")
emmeans(shannon_gamma_model, pairwise ~ Fraction, adjust = "fdr")
emmeans(shannon_gamma_model, pairwise ~ Nanotrap_type | Enrichment, adjust = "fdr")
emmeans(shannon_gamma_model, pairwise ~ Fraction | Enrichment, adjust = "fdr")
emmeans(shannon_gamma_model, pairwise ~ Fraction | Nanotrap_type | Enrichment, adjust = "fdr")

shannon_gaussian_model <- glmmTMB(shannon ~ Nanotrap_type*Enrichment + Fraction*Enrichment +
                                    Nanotrap_type * Fraction + (1 | LIMS_ID) + (1 | LIMS_ID:Nanotrap_type:Fraction),
                                  family = gaussian(link = "exp"),
                                  data = rpip_group_data)
AIC(shannon_gaussian_model)

sim_resid <- simulateResiduals(shannon_gaussian_model, n = 1000)
plot(sim_resid)
testDispersion(sim_resid)

emmeans(shannon_gaussian_model, pairwise ~ Nanotrap_type, adjust = "fdr")
emmeans(shannon_gaussian_model, pairwise ~ Fraction, adjust = "fdr")
emmeans(shannon_gaussian_model, pairwise ~ Nanotrap_type | Enrichment, adjust = "fdr")
emmeans(shannon_gaussian_model, pairwise ~ Fraction | Enrichment, adjust = "fdr")
emmeans(shannon_gaussian_model, pairwise ~ Fraction | Nanotrap_type | Enrichment, adjust = "fdr")


shannon_gamma_model <- glmmTMB(shannon ~ Nanotrap_type*Enrichment + Fraction*Enrichment +
                                  Nanotrap_type * Fraction + (1 | LIMS_ID) + (1 | LIMS_ID:Nanotrap_type:Fraction),
                                  family = Gamma(),
                                  data = rpip_group_data)
AIC(shannon_gamma_model)

sim_resid <- simulateResiduals(shannon_gamma_model, n = 1000)
plot(sim_resid)
testDispersion(sim_resid)

emmeans(shannon_lognormal_model, pairwise ~ Nanotrap_type, adjust = "fdr")
emmeans(shannon_lognormal_model, pairwise ~ Fraction, adjust = "fdr")
emmeans(shannon_lognormal_model, pairwise ~ Nanotrap_type | Enrichment, adjust = "fdr")
emmeans(shannon_lognormal_model, pairwise ~ Fraction | Enrichment, adjust = "fdr")
emmeans(shannon_lognormal_model, pairwise ~ Fraction | Nanotrap_type | Enrichment, adjust = "fdr")


richness_nbinom_model <- glmmTMB(richness ~ Nanotrap_type*Enrichment + Fraction*Enrichment +
                                   Nanotrap_type * Fraction + (1 | LIMS_ID) + (1 | LIMS_ID:Nanotrap_type:Fraction),
                                 family = nbinom1,
                                 data = rpip_group_data)
sim_resid <- simulateResiduals(richness_nbinom_model, n = 1000)
plot(sim_resid)
testDispersion(sim_resid)
testZeroInflation(sim_resid)

descdist(rpip_group_data$relative_richness, boot = 1000)
scaled_richness_beta_model <- glmmTMB(relative_richness ~ Nanotrap_type*Enrichment + Fraction*Enrichment +
                                     Nanotrap_type * Fraction + (1 | LIMS_ID) + (1 | LIMS_ID:Nanotrap_type:Fraction),
                                   family = beta_family(),
                                   data = rpip_group_data)
AIC(scaled_richness_beta_model)

sim_resid <- simulateResiduals(scaled_richness_beta_model, n = 1000)
plot(sim_resid)
testDispersion(sim_resid)
testZeroInflation(sim_resid)

emmeans(shannon_lognormal_model, pairwise ~ Nanotrap_type, adjust = "fdr")
emmeans(shannon_lognormal_model, pairwise ~ Fraction, adjust = "fdr")
emmeans(shannon_lognormal_model, pairwise ~ Nanotrap_type | Enrichment, adjust = "fdr")
emmeans(shannon_lognormal_model, pairwise ~ Fraction | Enrichment, adjust = "fdr")
emmeans(shannon_lognormal_model, pairwise ~ Fraction | Nanotrap_type | Enrichment, adjust = "fdr")



ginisimp_beta_model_zi <- glmmTMB(gini_simpson_complete ~ Nanotrap_type*Enrichment + Fraction*Enrichment +
                                    Nanotrap_type * Fraction + (1 | LIMS_ID) + (1 | LIMS_ID:Nanotrap_type:Fraction),
                                  family = beta_family(),
                                  data = rpip_group_data)
AIC(ginisimp_beta_model_zi)

sim_resid <- simulateResiduals(ginisimp_beta_model_zi, n = 1000)
plot(sim_resid)
testDispersion(sim_resid)
testZeroInflation(sim_resid)

emmeans(ginisimp_beta_model_zi, pairwise ~ Nanotrap_type, adjust = "fdr")
emmeans(ginisimp_beta_model_zi, pairwise ~ Fraction, adjust = "fdr")
emmeans(ginisimp_beta_model_zi, pairwise ~ Nanotrap_type | Enrichment, adjust = "fdr")
emmeans(ginisimp_beta_model_zi, pairwise ~ Fraction | Enrichment, adjust = "fdr")


chao_lognormal_model <- glmmTMB(shannon ~ Nanotrap_type*Enrichment + Fraction*Enrichment +
                                  Nanotrap_type * Fraction + (1 | LIMS_ID),
                                ziformula = ~Enrichment + LIMS_ID + Fraction,
                                family = lognormal,
                                data = rpip_group_data)
AIC(chao_lognormal_model)

sim_resid <- simulateResiduals(chao_lognormal_model, n = 1000)
plot(sim_resid)
testDispersion(sim_resid)

emmeans(chao_lognormal_model, pairwise ~ Nanotrap_type, adjust = "fdr")
emmeans(chao_lognormal_model, pairwise ~ Fraction, adjust = "fdr")
emmeans(chao_lognormal_model, pairwise ~ Nanotrap_type | Enrichment, adjust = "fdr")
emmeans(chao_lognormal_model, pairwise ~ Fraction | Enrichment, adjust = "fdr")
emmeans(chao_lognormal_model, pairwise ~ Fraction | Nanotrap_type | Enrichment, adjust = "fdr")

visualize(chao_lognormal_model)
