#Fit linear mixed effects models to diversity measures of reads mapped to RPIP
#bacterial targets

library(tidyverse)
library(vegan)
library(edgeR)
library(nlme)

rpip_and_unt_counts_covstats_and_info <- read_rds(
  "input/modified/rpip_and_unt_counts_covstats_and_info.rds"
)

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
                                         c("^A$" = "NTM-A",
                                           "A&B" = "NTM-AB",
                                           "none" = "DirEx")),
         Enrichment = str_replace_all(Enrichment,
                                      c("None" = "Non-targeted"))
  )

#Get library sizes
rpip_libsizes = rpip_and_unt_counts_covstats_and_info %>%
  mutate(UniqueID = paste(LIMS_ID, Treatment, Enrichment, sep = "-")) %>%
  mutate(UniqueID = str_replace_all(UniqueID, "Non-targeted", "None")) %>%
  select(UniqueID, summary.after_dedup.total_reads) %>%
  distinct() %>%
  rename(libsize = summary.after_dedup.total_reads)


#Filter group data to RPIP and unenriched samples
rpip_bact_group_data <- group_data %>%
  filter(Enrichment != "VSP") %>%
  left_join(rpip_libsizes)

#DGE object created in EdgeR_rgi_bactlevel.R:
rpip_bact_DGElist <- read_rds("input/modified/edgeR_bwamem_rpip_DGE.rds")

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


#Transpose for use in Vegan functions, and select only RPIP samples
rpip_bacteria_community_matrix_noNA <- (rpip_bact_DGElist$counts[
  rpip_bact_DGElist$genes$Species %in% rpip_bacteria_df$species,
] %>% 
  t())[rpip_bact_group_data$UniqueID,]

#Remove taxa with 0 total observed reads:
rpip_bacteria_community_matrix_noNA <- rpip_bacteria_community_matrix_noNA[
  , colSums(rpip_bacteria_community_matrix_noNA, na.rm = TRUE) > 0
]

#Make a version with NAs for some distance calculations:
rpip_bacteria_community_matrix <- rpip_bacteria_community_matrix_noNA
rpip_bacteria_community_matrix[rpip_bacteria_community_matrix == 0] <- NA_real_


#Calculate number of non-target sequences per sample:
Other <- rpip_bact_group_data$libsize -
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
  rpip_bacteria_community_matrix_noNA / rpip_bact_group_data$libsize
#And again, this time including "Other" column
rpip_bacteria_relative_complete <- 
  rpip_bacteria_community_matrix_noNA_complete / rpip_bact_group_data$libsize

#Wisconsin transformation
rpip_bacteria_wisconsin <- wisconsin(rpip_bacteria_community_matrix_noNA)

rpip_bacteria_wisconsin_complete <- wisconsin(
  rpip_bacteria_community_matrix_noNA_complete
)


#####
#Alpha diversity

rpip_bact_group_data$richness <- specnumber(rpip_bacteria_community_matrix_noNA)
rpip_bact_group_data <- rpip_bact_group_data %>%
  mutate(relative_richness = richness / libsize) %>%
  mutate(shannon = diversity(rpip_bacteria_community_matrix_noNA, index = "shannon"),
         simpson = simpson.unb(rpip_bacteria_community_matrix_noNA),
         simpson_complete = simpson.unb(rpip_bacteria_community_matrix_noNA_complete),
         invsimpson = simpson.unb(rpip_bacteria_community_matrix_noNA, inverse = TRUE),
         invsimpson_complete = simpson.unb(rpip_bacteria_community_matrix_noNA_complete, inverse = TRUE),
         gini_simpson = 1 - abs(simpson),
         gini_simpson_complete = 1 - simpson_complete,
         fisher = fisher.alpha(rpip_bacteria_community_matrix_noNA),
         chao1 = estimateR(rpip_bacteria_community_matrix_noNA)["S.chao1",]
  ) %>%
  mutate(Nanotrap_type = relevel(factor(Nanotrap_type, ordered = FALSE), ref = "NTM-A"),
         Fraction = relevel(factor(Fraction, ordered = FALSE), ref = "Unfil"),
         LIMS_ID = factor(LIMS_ID))



rpip_bact_group_data$LIMS_ID <- factor(rpip_bact_group_data$LIMS_ID)
rpip_bact_group_data$site <- factor(rpip_bact_group_data$site)
rpip_bact_group_data$Fraction <- factor(rpip_bact_group_data$Fraction)
rpip_bact_group_data$Nanotrap_type <- factor(rpip_bact_group_data$Nanotrap_type)
rpip_bact_group_data$Enrichment <- factor(rpip_bact_group_data$Enrichment)

####################
#####Fit Models#####
####################

#Potential variance structures:
var1 <- varIdent(form = ~ 1 | Enrichment * Nanotrap_type * Fraction)
var2 <- varIdent(form = ~ 1 | Enrichment * Nanotrap_type)
var3 <- varIdent(form = ~ 1 | Enrichment * Fraction)
var4 <- varIdent(form = ~ 1 | Enrichment)
var5 <- varIdent(form = ~ 1 | Nanotrap_type * Fraction)
var6 <- varIdent(form = ~ 1 | Nanotrap_type)
var7 <- varIdent(form = ~ 1 | Fraction)


#Try basic LME baselines:
basic_richness_lme <- lme(
  relative_richness ~ 1 + Fraction*Nanotrap_type*Enrichment,
  random = ~ 1 | LIMS_ID,
  data = rpip_bact_group_data)
richness_lme_weighted1 <- lme(
  relative_richness ~ 1 + Fraction*Nanotrap_type*Enrichment,
  random = ~ 1 | LIMS_ID,
  data = rpip_bact_group_data,
  weights = var1)
richness_lme_weighted2 <- lme(
  relative_richness ~ 1 + Fraction*Nanotrap_type*Enrichment,
  random = ~ 1 | LIMS_ID,
  data = rpip_bact_group_data,
  weights = var2)
richness_lme_weighted3 <- lme(
  relative_richness ~ 1 + Fraction*Nanotrap_type*Enrichment,
  random = ~ 1 | LIMS_ID,
  data = rpip_bact_group_data,
  weights = var3)
richness_lme_weighted4 <- lme(
  relative_richness ~ 1 + Fraction*Nanotrap_type*Enrichment,
  random = ~ 1 | LIMS_ID,
  data = rpip_bact_group_data,
  weights = var4)
richness_lme_weighted5 <- lme(
  relative_richness ~ 1 + Fraction*Nanotrap_type*Enrichment,
  random = ~ 1 | LIMS_ID,
  data = rpip_bact_group_data,
  weights = var5)
richness_lme_weighted6 <- lme(
  relative_richness ~ 1 + Fraction*Nanotrap_type*Enrichment,
  random = ~ 1 | LIMS_ID,
  data = rpip_bact_group_data,
  weights = var6)
richness_lme_weighted7 <- lme(
  relative_richness ~ 1 + Fraction*Nanotrap_type*Enrichment,
  random = ~ 1 | LIMS_ID,
  data = rpip_bact_group_data,
  weights = var7)



#Compare the three models: full weighted version is better
anova(basic_richness_lme, 
      richness_lme_weighted1)

#1 is better than 2
anova(richness_lme_weighted1,
      richness_lme_weighted2)

#1 is better than 3
anova(richness_lme_weighted1,
      richness_lme_weighted3)

#1 is better than 4
anova(richness_lme_weighted1,
      richness_lme_weighted4)

#1 is better than 5
anova(richness_lme_weighted1,
      richness_lme_weighted5)

#1 is better than 6
anova(richness_lme_weighted1,
      richness_lme_weighted6)


#1 is better than 7
anova(richness_lme_weighted1,
      richness_lme_weighted7)


#Try without random effect to be sure it's necessary:
basic_richness_gls <- gls(
  relative_richness ~ 1 + Fraction * Nanotrap_type * Enrichment, 
  data = rpip_bact_group_data)
richness_gls_weighted1 <- gls(
  relative_richness ~ 1 + Fraction*Nanotrap_type*Enrichment,
  data = rpip_bact_group_data,
  weights = var1)

anova(basic_richness_gls, richness_gls_weighted1)

######################################
#####Select variables for the LME#####
######################################

###Drop 3-way interaction###
#First the long-form model:
R.Base <- lme(relative_richness ~ 1 + Fraction + Nanotrap_type + Enrichment +
                Fraction:Nanotrap_type +
                Fraction:Enrichment +
                Nanotrap_type:Enrichment +
                Fraction:Nanotrap_type:Enrichment,
              random = ~ 1 | LIMS_ID,
              data = rpip_bact_group_data,
              method = "ML",
              weights = var1,
              control = lmeControl(msMaxIter = 1000, msMaxEval = 1000)
              )

#No 3-way interaction:
R.Drop1 <- lme(relative_richness ~ 1 + Fraction + Nanotrap_type + Enrichment +
                 Fraction:Nanotrap_type +
                 Fraction:Enrichment +
                 Nanotrap_type:Enrichment,
               random = ~ 1 | LIMS_ID,
               data = rpip_bact_group_data,
               method = "ML",
               weights = var1)

anova(R.Base, R.Drop1)
#p = 0.71 so 3-way interaction is not significant; R.Drop1 is new base. 

###Now try dropping the other interactions###
R1.Base <- R.Drop1
R1.Drop1 <- update(R1.Base, .~. -Nanotrap_type:Enrichment)
anova(R1.Base, R1.Drop1)
#p = 1e-04; signif.

R1.Drop2 <- update(R1.Base, .~. -Fraction:Enrichment)
anova(R1.Base, R1.Drop2)
#p < 0.0001; signif.

R1.Drop3 <- update(R1.Base, .~. -Fraction:Nanotrap_type,
                   control = lmeControl(msMaxIter = 1000, msMaxEval = 1000))
anova(R1.Base, R1.Drop3)
#p = 0.46; not signif. R1.Drop3 is the new best model.


R2.Base <- R1.Drop3

#Check if any insignificant terms remain:
anova(R2.Base)
#Now all terms are significant.

#Update the model method:
R2.REML <- update(R2.Base, method = "REML")

#Plot residuals for final model:
E <- resid(R2.REML, type = "normalized")
plot(x = fitted(R2.REML), y = E) #Looks pretty good

hist(E) #Looks good!

#Plot residuals by categorical variables:
plot(x = rpip_bact_group_data$Nanotrap_type, y = E)
plot(x = rpip_bact_group_data$Fraction, y = E)
plot(x = rpip_bact_group_data$Enrichment, y = E)
plot(x = rpip_bact_group_data$LIMS_ID, y = E)


#Get model info:
anova(R2.REML)
summary(R2.REML)



###############
#####Plots#####
###############
rpip_bact_group_data_plot <- rpip_bact_group_data %>%
  mutate(Fraction = factor(Fraction, levels = c("Fil", "Ret", "Unfil")))

ggplot(rpip_bact_group_data,
       aes(x = Nanotrap_type, y = relative_richness*10^6, fill = Enrichment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3")) +
  theme_bw() +
  ylab("Unique bacterial targets per million reads") +
  xlab("Concentration method")


ggplot(rpip_bact_group_data,
       aes(x = Fraction, y = relative_richness*10^6, fill = Enrichment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3")) +
  theme_bw() +
  ylab("Unique bacterial targets per million reads") +
  xlab("Fraction") +
  scale_y_log10()

######################
#####Other models#####
######################

#Fraction and nanotrap type are significant for Shannon
SMod <- lme(shannon ~ 1 + Fraction + Nanotrap_type + Enrichment +
              Fraction:Enrichment +
              Nanotrap_type:Enrichment,
            random = ~ 1 | LIMS_ID,
            data = rpip_bact_group_data,
            method = "REML",
            weights = var1,
            control = lmeControl(msMaxIter = 1000, msMaxEval = 1000)
)
plot(SMod)
anova(SMod)

ggplot(rpip_bact_group_data_plot,
       aes(x = Fraction, y = shannon, fill = Enrichment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3")) +
  theme_bw() +
  ylab("Shannon") +
  xlab("Fraction")

ggplot(rpip_bact_group_data_plot,
       aes(x = Nanotrap_type, y = shannon, fill = Enrichment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3")) +
  theme_bw() +
  ylab("Shannon") +
  xlab("Nanotrap_type")


#Nanotrap_type:Enrichment and Fraction significant for Simpson
SiMod <- lme(simpson ~ 1 + Fraction + Nanotrap_type + Enrichment +
               Fraction:Enrichment +
               Nanotrap_type:Enrichment,
             random = ~ 1 | LIMS_ID,
             data = rpip_bact_group_data,
             method = "REML",
             weights = var1,
             control = lmeControl(msMaxIter = 1000, msMaxEval = 1000)
)

plot(SiMod)
anova(SiMod)

ggplot(rpip_bact_group_data_plot,
       aes(x = Fraction, y = simpson, fill = Enrichment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3")) +
  theme_bw() +
  ylab("Simpson") +
  xlab("Fraction")

ggplot(rpip_bact_group_data_plot,
       aes(x = Nanotrap_type, y = simpson, fill = Enrichment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3")) +
  theme_bw() +
  ylab("Simpson") +
  xlab("Nanotrap_type")


#All three main effects but no interactions are significant for Fisher
FiMod <- lme(fisher ~ 1 + Fraction + Nanotrap_type + Enrichment +
               Fraction:Enrichment +
               Nanotrap_type:Enrichment,
             random = ~ 1 | LIMS_ID,
             data = rpip_bact_group_data,
             method = "REML",
             weights = var1,
             control = lmeControl(msMaxIter = 1000, msMaxEval = 1000)
)

plot(FiMod)
anova(FiMod)

ggplot(rpip_bact_group_data_plot,
       aes(x = Fraction, y = fisher, fill = Enrichment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3")) +
  theme_bw() +
  ylab("Fisher") +
  xlab("Fraction")

ggplot(rpip_bact_group_data_plot,
       aes(x = Nanotrap_type, y = fisher, fill = Enrichment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3")) +
  theme_bw() +
  ylab("Fisher") +
  xlab("Nanotrap_type")


#Fraction and Nanotrap type significant for Chao1
ChMod <- lme(chao1 ~ 1 + Fraction + Nanotrap_type + Enrichment +
               Fraction:Enrichment +
               Nanotrap_type:Enrichment,
             random = ~ 1 | LIMS_ID,
             data = rpip_bact_group_data,
             method = "REML",
             weights = var1,
             control = lmeControl(msMaxIter = 1000, msMaxEval = 1000)
)

plot(ChMod)
anova(ChMod)

ggplot(rpip_bact_group_data_plot,
       aes(x = Fraction, y = chao1, fill = Enrichment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3")) +
  theme_bw() +
  ylab("Chao1") +
  xlab("Fraction")

ggplot(rpip_bact_group_data_plot,
       aes(x = Nanotrap_type, y = chao1, fill = Enrichment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3")) +
  theme_bw() +
  ylab("Fisher") +
  xlab("Nanotrap_type")
