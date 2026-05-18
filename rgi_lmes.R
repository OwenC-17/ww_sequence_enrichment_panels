library(tidyverse)
library(vegan)

#Import data
rpip_and_unt_allele_and_info <- read_rds("input/modified/rpip_and_unt_rgi_allele_and_info.rds")
rpip_and_unt_gene_and_info <- read_rds("input/modified/rpip_and_unt_rgi_gene_and_info.rds")

#Get group data for RPIP only
rpip_gene_group_data <- rpip_gene_group_data %>%
  rename("libsize" = "summary.after_dedup.total_reads")

#Get the genes from the DGE object created in EdgeR_rgi_genelevel.R
#Transpose for use in Vegan functions, and select only RPIP samples
gene_community_matrix_noNA <- (rpip_gene_DGElist$counts %>% t())[rpip_gene_group_data$UniqueID,]
gene_community_matrix_noNA <- gene_community_matrix_noNA[, colSums(gene_community_matrix_noNA, na.rm = TRUE) > 0]
gene_community_matrix <- gene_community_matrix_noNA
gene_community_matrix[gene_community_matrix == 0] <- NA_real_


#Calculate number of non-ARG sequences per sample
Other <- rpip_gene_group_data$libsize - rowSums(gene_community_matrix, na.rm = TRUE)
#Add column to community matrix
gene_community_matrix_complete <- cbind(gene_community_matrix, Other)
gene_community_matrix_noNA_complete <- cbind(gene_community_matrix_noNA, Other)

#Convert to relative abundance but don't include "Other" column
gene_relative <- gene_community_matrix_noNA / rpip_gene_group_data$libsize
#And again, this time including "Other" column
gene_relative_complete <- gene_community_matrix_noNA_complete / rpip_gene_group_data$libsize

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
gene_bray_raw_nmds_2d_df_b <- make_mds_df(gene_community_matrix_noNA, rpip_gene_group_data)
ggplot(gene_bray_raw_nmds_2d_df_b, aes(x = MDS1, y = MDS2, colour = Fraction)) + geom_point()
#The odd grouping of the lower-left retentate sample (20040-RA-RPIP) can likely
#be attributed to low AMR read count and low diversity in that sample


#Raw counts, including "Other"
gene_bray_raw_complete_nmds_2d_df_b <- make_mds_df(gene_community_matrix_noNA_complete, rpip_gene_group_data)
ggplot(gene_bray_raw_complete_nmds_2d_df_b, aes(x = MDS1, y = MDS2, colour = Fraction)) + geom_point()


#Proportions, not including "Other"
gene_bray_relative_nmds_2d_df_b <- make_mds_df(gene_relative, rpip_gene_group_data)
ggplot(gene_bray_relative_nmds_2d_df_b, aes(x = MDS1, y = MDS2, colour = Fraction)) + geom_point()


#Wisconsin, not including "Other"
gene_bray_wisconsin_nmds_2d_df <- make_mds_df(gene_wisconsin, rpip_gene_group_data)
ggplot(gene_bray_wisconsin_nmds_2d_df, aes(x = MDS1, y = MDS2, colour = Fraction)) + geom_point()

#Proportions, including "Other"
gene_bray_relative_complete_nmds_2d_df <- make_mds_df(gene_relative_complete, rpip_gene_group_data)
ggplot(gene_bray_relative_nmds_2d_df, aes(x = MDS1, y = MDS2, colour = Fraction)) + geom_point()

#Wisconsin, including "Other"
gene_bray_wisconsin_complete_nmds_2d_df <- make_mds_df(gene_wisconsin_complete, rpip_gene_group_data)
ggplot(gene_bray_wisconsin_nmds_2d_df, aes(x = MDS1, y = MDS2, colour = Fraction)) + geom_point()

#Jaccard, not including "Other" (shouldn't need to include since it is present in all samples)
gene_jacc_nmds_2d_df <- make_mds_df(gene_community_matrix_noNA, rpip_gene_group_data, dist = "jaccard")
ggplot(gene_jacc_nmds_2d_df, aes(x = MDS1, y = MDS2, colour = Fraction, shape = site)) + geom_point(size = 5)

#Jaccard-type Chao distance, not including "Other" (shouldn't need to include since it is present in all samples)
gene_community_matrix_noNA <- gene_community_matrix
gene_community_matrix_noNA[is.na(gene_community_matrix_noNA)] <- 0
gene_chao_nmds_2d_df <- make_mds_df(gene_community_matrix_noNA, rpip_gene_group_data, dist = "chao")
ggplot(gene_chao_nmds_2d_df, aes(x = MDS1, y = MDS2, colour = Fraction, shape = site)) + geom_point(size = 5)


####################
#####PERMANOVA######
####################

#Proportion, not including "Other"
adonis2(gene_relative ~ site + Fraction + Nanotrap_type, data = rpip_gene_group_data, by = "margin", na.rm = TRUE)

#Proportion, including "Other"
adonis2(gene_relative_complete ~ site + Fraction + Nanotrap_type, data = rpip_gene_group_data, by = "margin", na.rm = TRUE)

#Wisconsin, not including "Other"
set.seed(123)
adonis2(gene_wisconsin ~ site + Fraction + Nanotrap_type, data = rpip_gene_group_data, by = "margin", na.rm = TRUE)

#Wisconsin, including "Other"
adonis2(gene_wisconsin_complete ~ site + Fraction + Nanotrap_type, data = rpip_gene_group_data, by = "margin", na.rm = TRUE)


gene_presabs_matrix <- gene_community_matrix_noNA
gene_presabs_matrix[gene_presabs_matrix > 0] <- 1
gene_presabs_matrix_complete <- gene_community_matrix_noNA_complete
gene_presabs_matrix_complete[gene_presabs_matrix_complete > 0] <- 1

#Jaccard
adonis2(gene_presabs_matrix ~ site + Fraction + Nanotrap_type, data = rpip_gene_group_data, by = "margin", method = "jaccard", permutations = 9999)
adonis2(gene_presabs_matrix_complete ~ site + Fraction + Nanotrap_type, data = rpip_gene_group_data, by = "margin", method = "jaccard", permutations = 9999)


#Chao
adonis2(gene_community_matrix_noNA ~ site + Enrichment + Fraction + Nanotrap_type, data = rpip_gene_group_data, by = "margin", method = "chao")


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


rpip_gene_group_data$richness <- specnumber(gene_community_matrix_noNA)

rpip_gene_group_data <- rpip_gene_group_data %>%
  mutate(relative_richness = richness / libsize) %>%
  mutate(shannon = diversity(gene_community_matrix_noNA, index = "shannon"),
         simpson = simpson.unb(gene_community_matrix_noNA),
         invsimpson = simpson.unb(gene_community_matrix_noNA, inverse = TRUE),
         fisher = fisher.alpha(gene_community_matrix_noNA),
         chao1 = estimateR(gene_community_matrix_noNA)["S.chao1",]) %>%
  mutate(Nanotrap_type = relevel(factor(Nanotrap_type, ordered = FALSE), ref = "none"),
         Fraction = relevel(factor(Fraction, ordered = FALSE), ref = "unfiltered"))

ggplot(rpip_gene_group_data, aes(x = Treatment, y = richness, colour = site)) + geom_boxplot()
ggplot(rpip_gene_group_data, aes(x = Treatment, y = shannon, colour = site)) + geom_boxplot()
ggplot(rpip_gene_group_data, aes(x = Treatment, y = simpson, colour = site)) + geom_boxplot()
ggplot(rpip_gene_group_data, aes(x = Treatment, y = invsimpson, colour = site)) + geom_boxplot()
ggplot(rpip_gene_group_data, aes(x = Treatment, y = fisher, colour = site)) + geom_boxplot()
ggplot(rpip_gene_group_data, aes(x = Treatment, y = chao1, colour = site)) + geom_boxplot()


library(spm2)

simpson_model <- glm(simpson ~ Nanotrap_type + Fraction + site, data = rpip_gene_group_data)
summary(simpson_model)
richness_model <- glm(richness ~ Nanotrap_type + Fraction, data = rpip_gene_group_data, family = "poisson")
summary(richness_model)
anova(richness_model)
shannon_model <- lm(shannon ~ Nanotrap_type + Fraction, data = rpip_gene_group_data)
summary(shannon_model)


richness_model_cv <- glmcv(richness ~ Nanotrap_type + Fraction, rpip_gene_group_data,rpip_gene_group_data$richness, family = "poisson", validation = "LOO", cv.fold = 500)

library(randomForest)
set.seed(1234)
train <- sample(1:48, size = 36)
richness_model_rf <- randomForest(
  formula = richness ~ Nanotrap_type + Fraction + site,
  data = rpip_gene_group_data[train,]
)

test <- rpip_gene_group_data[-train, ]
test$pred <- predict(richness_model_rf, test)

ggplot(test, aes(x = pred, y = richness)) + geom_point()

summary(richness_model_rf)

ggplot(rpip_gene_group_data, aes(x = log(simpson / (1 - simpson)))) + geom_density()


ggplot(rpip_gene_group_data, aes(x = richness, colour = Treatment)) + geom_density(alpha = 0.5)

ggplot(rpip_gene_group_data, aes(x = shannon)) + geom_density(alpha = 0.5)




##############GLS
library(nlme)
var1 <- varIdent(form = ~ 1 | Enrichment * Nanotrap_type * Fraction)
var2 <- varIdent(form = ~ 1 | Enrichment * Nanotrap_type)



basic_richness_lme <- lme(relative_richness ~ 1 + Fraction*Nanotrap_type*Enrichment, random = ~ 1 | LIMS_ID, data = rpip_gene_group_data)
basic_richness_lme_weighted <- lme(relative_richness ~ 1 + Fraction*Nanotrap_type*Enrichment, random = ~ 1 | LIMS_ID, data = rpip_gene_group_data, weights = var1)
basic_richness_lme_weighted2 <- lme(relative_richness ~ 1 + Fraction*Nanotrap_type*Enrichment, random = ~ 1 | LIMS_ID, data = rpip_gene_group_data, weights = var2)


summary(basic_richness_lme)
plot(basic_richness_lme)
anova(basic_richness_lme, basic_richness_lme_weighted, basic_richness_lme_weighted2)
#Full weighted version is better

rpip_gene_group_data$LIMS_ID <- factor(rpip_gene_group_data$LIMS_ID)
rpip_gene_group_data$site <- factor(rpip_gene_group_data$site)
rpip_gene_group_data$Fraction <- factor(rpip_gene_group_data$Fraction)
rpip_gene_group_data$Nanotrap_type <- factor(rpip_gene_group_data$Nanotrap_type)
rpip_gene_group_data$Enrichment <- factor(rpip_gene_group_data$Enrichment)



#No random effect:
basic_richness_gls <- gls(relative_richness ~ 1 + Fraction*Nanotrap_type*Enrichment, data = rpip_gene_group_data)
basic_richness_gls_weighted <- gls(relative_richness ~ 1 + Fraction*Nanotrap_type*Enrichment, data = rpip_gene_group_data, weights = var1)


summary(basic_richness_gls)
plot(basic_richness_gls)
anova(basic_richness_gls)

anova(basic_richness_lme_weighted, basic_richness_gls_weighted)
#The lme is better; the random effect is very significant (p < 0.0001)

#Drop 3-way interaction
#First the long-form model:
R.Base <- lme(relative_richness ~ 1 + Fraction + Nanotrap_type + Enrichment +
                Fraction:Nanotrap_type +
                Fraction:Enrichment +
                Nanotrap_type:Enrichment +
                Fraction:Nanotrap_type:Enrichment,
              random = ~ 1 | LIMS_ID,
              data = rpip_gene_group_data,
              method = "ML",
              weights = var1)

#No 3-way interaction:
R.Drop1 <- lme(relative_richness ~ 1 + Fraction + Nanotrap_type + Enrichment +
                Fraction:Nanotrap_type +
                Fraction:Enrichment +
                Nanotrap_type:Enrichment,
              random = ~ 1 | LIMS_ID,
              data = rpip_gene_group_data,
              method = "ML",
              weights = var1)

anova(R.Base, R.Drop1)
#p = 0.76 so 3-way interaction is not significant; R.Drop1 is new base. 

#Now try dropping the other interactions.
R1.Base <- R.Drop1
R1.Drop1 <- update(R1.Base, .~. -Nanotrap_type:Enrichment)
anova(R1.Base, R1.Drop1)
#p = 0.0057; signif.

R1.Drop2 <- update(R1.Base, .~. -Fraction:Enrichment)
anova(R1.Base, R1.Drop2)
#p = .0035; signif.

R1.Drop3 <- update(R1.Base, .~. -Fraction:Nanotrap_type, control = lmeControl(msMaxIter = 1000, msMaxEval = 1000))
anova(R1.Base, R1.Drop3)
#p = 0.16; not signif. R1.Drop3 is the new best model.

#Stage 2:
R2.Base <- R1.Drop3
anova(R2.Base)
#Now all terms are significant.
E <- resid(R2.Base, type = "normalized")
plot(x = fitted(R2.Base), y = E) #Looks pretty good

hist(E) #Looks great!

anova(R2.Base)
summary(R2.Base)
plot(x = rpip_gene_group_data$Nanotrap_type, y = E)
plot(x = rpip_gene_group_data$LIMS_ID, y = E)


rpip_gene_group_data_plot <- rpip_gene_group_data %>%
  mutate(Fraction = str_replace_all(Fraction, c("filtrate" = "Fil", "retentate" = "Ret", "unfiltered" = "Unfil")),
         Nanotrap_type = str_replace_all(Nanotrap_type, c("^A$" = "NT-A", "A&B" = "NT-AB", "none" = "DirEx")),
         Enrichment = str_replace_all(Enrichment, "None", "Non-targeted"))

ggplot(rpip_gene_group_data_plot, aes(x = Nanotrap_type, y = relative_richness*10^6, fill = Enrichment)) + geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3")) +
  theme_bw() +
  ylab("Unique ARGs per million reads") +
  xlab("Concentration method")


ggplot(rpip_gene_group_data_plot, aes(x = Fraction, y = relative_richness*10^6, fill = Enrichment)) + geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3")) +
  theme_bw() +
  ylab("Unique ARGs per million reads") +
  xlab("Fraction")


SMod <- lme(shannon ~ 1 + Fraction + Nanotrap_type + Enrichment +
                         Fraction:Enrichment +
                         Nanotrap_type:Enrichment,
                       random = ~ 1 | LIMS_ID,
                       data = rpip_gene_group_data,
                       method = "ML",
                       weights = var1,
            control = lmeControl(msMaxIter = 1000, msMaxEval = 1000)
            )
plot(SMod)


