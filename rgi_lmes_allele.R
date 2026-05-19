#Fit linear mixed effects models to RGI diversity measures

library(tidyverse)
library(vegan)
library(edgeR)
library(nlme)


#Import data
rpip_and_unt_allele_and_info <- read_rds(
  "input/modified/rpip_and_unt_rgi_allele_and_info.rds"
)

rpip_and_unt_allele_protein_homolog <- rpip_and_unt_allele_and_info %>%
  filter(`Reference Model Type` == "protein homolog model")

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

rpip_libsizes = rpip_and_unt_allele_protein_homolog %>%
  select(UniqueID, summary.after_dedup.total_reads) %>%
  distinct()

rpip_allele_group_data <- group_data %>%
  filter(Enrichment != "VSP") %>%
  left_join(rpip_libsizes) %>%
  rename("libsize" = "summary.after_dedup.total_reads")


#DGE object created in EdgeR_rgi_allelelevel.R:
rpip_allele_DGElist <- read_rds("input/modified/edgeR_rgi_allele_DGE.rds")

#Get the alleles from the DGE object created in EdgeR_rgi_allelelevel.R
#Transpose for use in Vegan functions, and select only RPIP samples
allele_community_matrix_noNA <- (rpip_allele_DGElist$counts %>%
                                 t())[rpip_allele_group_data$UniqueID, ]
allele_community_matrix_noNA <- allele_community_matrix_noNA[ ,
                colSums(allele_community_matrix_noNA, na.rm = TRUE) > 0]
allele_community_matrix <- allele_community_matrix_noNA
allele_community_matrix[allele_community_matrix == 0] <- NA_real_


#Calculate number of non-ARG sequences per sample
Other <- rpip_allele_group_data$libsize - rowSums(allele_community_matrix,
                                                na.rm = TRUE)
#Add column to community matrix
allele_community_matrix_complete <- cbind(allele_community_matrix, Other)
allele_community_matrix_noNA_complete <- cbind(allele_community_matrix_noNA,
                                               Other)

#Convert to relative abundance but don't include "Other" column
allele_relative <- allele_community_matrix_noNA / rpip_allele_group_data$libsize
#And again, this time including "Other" column
allele_relative_complete <- allele_community_matrix_noNA_complete /
  rpip_allele_group_data$libsize

#Wisconsin transformation
allele_wisconsin <- wisconsin(allele_community_matrix_noNA)

allele_wisconsin_complete <- wisconsin(allele_community_matrix_noNA_complete)

###########################################
#####Calculate alpha diversity metrics#####
###########################################

rpip_allele_group_data <- rpip_allele_group_data %>%
  mutate(richness = specnumber(allele_community_matrix_noNA),
         relative_richness = richness / libsize) %>%
  mutate(shannon = diversity(allele_community_matrix_noNA, index = "shannon"),
         simpson = simpson.unb(allele_community_matrix_noNA),
         invsimpson = simpson.unb(allele_community_matrix_noNA, inverse = TRUE),
         gini_simpson = 1 - simpson,
         fisher = fisher.alpha(allele_community_matrix_noNA),
         chao1 = estimateR(allele_community_matrix_noNA)["S.chao1",]) %>%
  mutate(Nanotrap_type = relevel(factor(Nanotrap_type, ordered = FALSE), ref = "DirEx"),
         Fraction = relevel(factor(Fraction, ordered = FALSE), ref = "Unfil"))


rpip_allele_group_data$LIMS_ID <- factor(rpip_allele_group_data$LIMS_ID)
rpip_allele_group_data$site <- factor(rpip_allele_group_data$site)
rpip_allele_group_data$Fraction <- factor(rpip_allele_group_data$Fraction)
rpip_allele_group_data$Nanotrap_type <- factor(
  rpip_allele_group_data$Nanotrap_type)
rpip_allele_group_data$Enrichment <- factor(rpip_allele_group_data$Enrichment)


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
  data = rpip_allele_group_data)
richness_lme_weighted1 <- lme(
  relative_richness ~ 1 + Fraction*Nanotrap_type*Enrichment,
  random = ~ 1 | LIMS_ID,
  data = rpip_allele_group_data,
  weights = var1)
richness_lme_weighted2 <- lme(
  relative_richness ~ 1 + Fraction*Nanotrap_type*Enrichment,
  random = ~ 1 | LIMS_ID,
  data = rpip_allele_group_data,
  weights = var2)
richness_lme_weighted3 <- lme(
  relative_richness ~ 1 + Fraction*Nanotrap_type*Enrichment,
  random = ~ 1 | LIMS_ID,
  data = rpip_allele_group_data,
  weights = var3)
richness_lme_weighted4 <- lme(
  relative_richness ~ 1 + Fraction*Nanotrap_type*Enrichment,
  random = ~ 1 | LIMS_ID,
  data = rpip_allele_group_data,
  weights = var4)
richness_lme_weighted5 <- lme(
  relative_richness ~ 1 + Fraction*Nanotrap_type*Enrichment,
  random = ~ 1 | LIMS_ID,
  data = rpip_allele_group_data,
  weights = var5)
richness_lme_weighted6 <- lme(
  relative_richness ~ 1 + Fraction*Nanotrap_type*Enrichment,
  random = ~ 1 | LIMS_ID,
  data = rpip_allele_group_data,
  weights = var6)
richness_lme_weighted7 <- lme(
  relative_richness ~ 1 + Fraction*Nanotrap_type*Enrichment,
  random = ~ 1 | LIMS_ID,
  data = rpip_allele_group_data,
  weights = var7)


#View residuals 
#Hard to infer anything from these tbh
plot(basic_richness_lme)
plot(richness_lme_weighted1)
plot(richness_lme_weighted2)
plot(richness_lme_weighted3)
plot(richness_lme_weighted4)
plot(richness_lme_weighted5)
plot(richness_lme_weighted6)
plot(richness_lme_weighted7)


#Compare the three models: weighted 2 (Enrichment*Nanotrap_type)
#is better based on AIC
anova(basic_richness_lme, 
      richness_lme_weighted1,
      richness_lme_weighted2)

#2 is better than 3
anova(richness_lme_weighted2,
      richness_lme_weighted3)

#4 (variance structured by Enrichment only) is better than 2
anova(richness_lme_weighted2,
      richness_lme_weighted4)

#4 is better than 5
anova(richness_lme_weighted4,
      richness_lme_weighted5)

#4 is better than 6
anova(richness_lme_weighted4,
      richness_lme_weighted6)


#4 is better than 7
anova(richness_lme_weighted4,
      richness_lme_weighted7)


#Try without random effect to be sure it's necessary:
basic_richness_gls <- gls(
  relative_richness ~ 1 + Fraction * Nanotrap_type * Enrichment, 
  data = rpip_allele_group_data)
richness_gls_weighted4 <- gls(
  relative_richness ~ 1 + Fraction*Nanotrap_type*Enrichment,
  data = rpip_allele_group_data,
  weights = var4)


anova(richness_lme_weighted1, richness_gls_weighted4)
#The LME model is better; the random effect is very significant (p < 0.0001)


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
              data = rpip_allele_group_data,
              method = "ML",
              weights = var4)

#No 3-way interaction:
R.Drop1 <- lme(relative_richness ~ 1 + Fraction + Nanotrap_type + Enrichment +
                 Fraction:Nanotrap_type +
                 Fraction:Enrichment +
                 Nanotrap_type:Enrichment,
               random = ~ 1 | LIMS_ID,
               data = rpip_allele_group_data,
               method = "ML",
               weights = var4)

anova(R.Base, R.Drop1)
#p = 0.804 so 3-way interaction is not significant; R.Drop1 is new base. 

###Now try dropping the other interactions###
R1.Base <- R.Drop1
R1.Drop1 <- update(R1.Base, .~. -Nanotrap_type:Enrichment)
anova(R1.Base, R1.Drop1)
#p = 0.61; not signif.

R1.Drop2 <- update(R1.Base, .~. -Fraction:Enrichment)
anova(R1.Base, R1.Drop2)
#p = .001; signif.

R1.Drop3 <- update(R1.Base, .~. -Fraction:Nanotrap_type,
                   control = lmeControl(msMaxIter = 1000, msMaxEval = 1000))
anova(R1.Base, R1.Drop3)
#p = 0.11; not signif. but Nanotrap_type:Enrichment has higher p value;
#R1.Drop1 is the new best model.


R2.Base <- R1.Drop1

#Check if any insignificant terms remain:
anova(R2.Base)
#Fraction:Nanotrap_type and Nanotrap_type main effect are not signif.

###Try dropping interactions first:
R2.Drop1 <- update(R2.Base, .~. -Fraction:Nanotrap_type)
anova(R2.Base, R2.Drop1)
#p = 0.1168; not signif. 

R2.Drop2 <- update(R2.Base, .~. -Fraction:Enrichment)
anova(R2.Base, R2.Drop2)
#p = 0.0011; signif.
#R2.Drop1 is new best

R3.Base <- R2.Drop1

#Check remaining terms:
anova(R3.Base)
#Nanotrap_type is still not signif (p = 0.09)

###Try dropping Nanotrap_type main effect:
R3.Drop1 <- update(R3.Base, .~. -Nanotrap_type)
anova(R3.Base, R3.Drop1)
#p = 0.07; not signif
#R3.Drop1 is new best

R4.Base = R3.Drop1

#Check remaining terms:
anova(R4.Base)
#All have p < 0.05

#Update the model method:
R4.REML <- update(R4.Base, method = "REML")

#Plot residuals for final model:
E <- resid(R4.REML, type = "normalized")
plot(x = fitted(R4.REML), y = E) #Looks pretty good

hist(E) #Looks okay

#Plot residuals by categorical variables:
plot(x = rpip_allele_group_data$Nanotrap_type, y = E)
plot(x = rpip_allele_group_data$Fraction, y = E)
plot(x = rpip_allele_group_data$Enrichment, y = E)
plot(x = rpip_allele_group_data$LIMS_ID, y = E)


#Get model info:
anova(R4.REML)
summary(R4.REML)

rpip_allele_group_data_plot <- rpip_allele_group_data %>%
  mutate(Fraction = factor(Fraction, levels = c("Fil", "Ret", "Unfil")))

ggplot(rpip_allele_group_data,
       aes(x = Nanotrap_type, y = relative_richness*10^6, fill = Enrichment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3")) +
  theme_bw() +
  ylab("Unique ARGs per million reads") +
  xlab("Concentration method")


ggplot(rpip_allele_group_data_plot,
       aes(x = Fraction, y = relative_richness*10^6, fill = Enrichment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3")) +
  theme_bw() +
  ylab("Unique ARGs per million reads") +
  xlab("Fraction")


#Nanotrap_type ans Nanotrap_type:Enrichment may be insignificant for Shannon
SMod <- lme(shannon ~ 1 + Fraction + Nanotrap_type + Enrichment +
              Fraction:Enrichment +
              Nanotrap_type:Enrichment,
            random = ~ 1 | LIMS_ID,
            data = rpip_allele_group_data,
            method = "REML",
            weights = var4,
            control = lmeControl(msMaxIter = 1000, msMaxEval = 1000)
)
plot(SMod)
anova(SMod)

ggplot(rpip_allele_group_data_plot,
       aes(x = Fraction, y = shannon, fill = Enrichment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3")) +
  theme_bw() +
  ylab("Shannon") +
  xlab("Fraction")

ggplot(rpip_allele_group_data_plot,
       aes(x = Nanotrap_type, y = shannon, fill = Enrichment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3")) +
  theme_bw() +
  ylab("Shannon") +
  xlab("Nanotrap_type")


#Nanotrap_type:Enrichment insignificant for Simpson
SiMod <- lme(simpson ~ 1 + Fraction + Nanotrap_type + Enrichment +
               Fraction:Enrichment +
               Nanotrap_type:Enrichment,
             random = ~ 1 | LIMS_ID,
             data = rpip_allele_group_data,
             method = "REML",
             weights = var4,
             control = lmeControl(msMaxIter = 1000, msMaxEval = 1000)
)

plot(SiMod)
anova(SiMod)

ggplot(rpip_allele_group_data_plot,
       aes(x = Fraction, y = simpson, fill = Enrichment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3")) +
  theme_bw() +
  ylab("Simpson") +
  xlab("Fraction")

ggplot(rpip_allele_group_data_plot,
       aes(x = Nanotrap_type, y = simpson, fill = Enrichment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3")) +
  theme_bw() +
  ylab("Simpson") +
  xlab("Nanotrap_type")


#Only Fraction is significant for Fisher
FiMod <- lme(fisher ~ 1 + Fraction + Nanotrap_type + Enrichment +
               Fraction:Enrichment +
               Nanotrap_type:Enrichment,
             random = ~ 1 | LIMS_ID,
             data = rpip_allele_group_data,
             method = "REML",
             weights = var4,
             control = lmeControl(msMaxIter = 1000, msMaxEval = 1000)
)

plot(FiMod)
anova(FiMod)

ggplot(rpip_allele_group_data_plot,
       aes(x = Fraction, y = fisher, fill = Enrichment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3")) +
  theme_bw() +
  ylab("Fisher") +
  xlab("Fraction")

ggplot(rpip_allele_group_data_plot,
       aes(x = Nanotrap_type, y = fisher, fill = Enrichment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3")) +
  theme_bw() +
  ylab("Fisher") +
  xlab("Nanotrap_type")


#Nanotrap_type:Enrichment is insignificant for Chao1
ChMod <- lme(chao1 ~ 1 + Fraction + Nanotrap_type + Enrichment +
               Fraction:Enrichment +
               Nanotrap_type:Enrichment,
             random = ~ 1 | LIMS_ID,
             data = rpip_allele_group_data,
             method = "REML",
             weights = var4,
             control = lmeControl(msMaxIter = 1000, msMaxEval = 1000)
)

plot(ChMod)
anova(ChMod)

ggplot(rpip_allele_group_data_plot,
       aes(x = Fraction, y = chao1, fill = Enrichment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3")) +
  theme_bw() +
  ylab("Chao1") +
  xlab("Fraction")

ggplot(rpip_allele_group_data_plot,
       aes(x = Nanotrap_type, y = chao1, fill = Enrichment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3")) +
  theme_bw() +
  ylab("Fisher") +
  xlab("Nanotrap_type")
