#Fit linear mixed effects models to RGI diversity measures

library(tidyverse)
library(vegan)
library(edgeR)
library(nlme)


#Import data
rpip_and_unt_allele_and_info <- read_rds(
  "input/modified/rpip_and_unt_rgi_allele_and_info.rds"
  )
rpip_and_unt_gene_and_info <- read_rds(
  "input/modified/rpip_and_unt_rgi_gene_and_info.rds"
  )

rpip_and_unt_gene_protein_homolog <- rpip_and_unt_gene_and_info %>%
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

rpip_libsizes = rpip_and_unt_gene_protein_homolog %>%
  select(UniqueID, summary.after_dedup.total_reads) %>%
  distinct()

rpip_gene_group_data <- group_data %>%
  filter(Enrichment != "VSP") %>%
  left_join(rpip_libsizes) %>%
  rename("libsize" = "summary.after_dedup.total_reads")


#DGE object created in EdgeR_rgi_genelevel.R:
rpip_gene_DGElist <- read_rds("input/modified/edgeR_rgi_gene_DGE.rds")

#Get the genes from the DGE object created in EdgeR_rgi_genelevel.R
#Transpose for use in Vegan functions, and select only RPIP samples
gene_community_matrix_noNA <- (rpip_gene_DGElist$counts %>%
                                 t())[rpip_gene_group_data$UniqueID, ]
gene_community_matrix_noNA <- gene_community_matrix_noNA[ ,
                colSums(gene_community_matrix_noNA, na.rm = TRUE) > 0]
gene_community_matrix <- gene_community_matrix_noNA
gene_community_matrix[gene_community_matrix == 0] <- NA_real_


#Calculate number of non-ARG sequences per sample
Other <- rpip_gene_group_data$libsize - rowSums(gene_community_matrix,
                                                na.rm = TRUE)
#Add column to community matrix
gene_community_matrix_complete <- cbind(gene_community_matrix, Other)
gene_community_matrix_noNA_complete <- cbind(gene_community_matrix_noNA, Other)

#Convert to relative abundance but don't include "Other" column
gene_relative <- gene_community_matrix_noNA / rpip_gene_group_data$libsize
#And again, this time including "Other" column
gene_relative_complete <- gene_community_matrix_noNA_complete /
  rpip_gene_group_data$libsize

#Wisconsin transformation
gene_wisconsin <- wisconsin(gene_community_matrix_noNA)

gene_wisconsin_complete <- wisconsin(gene_community_matrix_noNA_complete)

###########################################
#####Calculate alpha diversity metrics#####
###########################################
gene_relative_noNA <- gene_relative
gene_relative_noNA[is.na(gene_relative_noNA)] <- 0

rpip_gene_group_data$richness <- specnumber(gene_community_matrix_noNA)

rpip_gene_group_data <- rpip_gene_group_data %>%
  mutate(richness = specnumber(gene_community_matrix_noNA),
         relative_richness = richness / libsize) %>%
  mutate(shannon = diversity(gene_community_matrix_noNA, index = "shannon"),
         simpson = simpson.unb(gene_community_matrix_noNA),
         invsimpson = simpson.unb(gene_community_matrix_noNA, inverse = TRUE),
         gini_simpson = 1 - simpson,
         fisher = fisher.alpha(gene_community_matrix_noNA),
         chao1 = estimateR(gene_community_matrix_noNA)["S.chao1",]) %>%
  mutate(Nanotrap_type = relevel(factor(Nanotrap_type, ordered = FALSE), ref = "DirEx"),
         Fraction = relevel(factor(Fraction, ordered = FALSE), ref = "Unfil"))


rpip_gene_group_data$LIMS_ID <- factor(rpip_gene_group_data$LIMS_ID)
rpip_gene_group_data$site <- factor(rpip_gene_group_data$site)
rpip_gene_group_data$Fraction <- factor(rpip_gene_group_data$Fraction)
rpip_gene_group_data$Nanotrap_type <- factor(rpip_gene_group_data$Nanotrap_type)
rpip_gene_group_data$Enrichment <- factor(rpip_gene_group_data$Enrichment)


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
  data = rpip_gene_group_data)
richness_lme_weighted1 <- lme(
  relative_richness ~ 1 + Fraction*Nanotrap_type*Enrichment,
  random = ~ 1 | LIMS_ID,
  data = rpip_gene_group_data,
  weights = var1)
richness_lme_weighted2 <- lme(
  relative_richness ~ 1 + Fraction*Nanotrap_type*Enrichment,
  random = ~ 1 | LIMS_ID,
  data = rpip_gene_group_data,
  weights = var2)
richness_lme_weighted3 <- lme(
  relative_richness ~ 1 + Fraction*Nanotrap_type*Enrichment,
  random = ~ 1 | LIMS_ID,
  data = rpip_gene_group_data,
  weights = var3)
richness_lme_weighted4 <- lme(
  relative_richness ~ 1 + Fraction*Nanotrap_type*Enrichment,
  random = ~ 1 | LIMS_ID,
  data = rpip_gene_group_data,
  weights = var4)
richness_lme_weighted5 <- lme(
  relative_richness ~ 1 + Fraction*Nanotrap_type*Enrichment,
  random = ~ 1 | LIMS_ID,
  data = rpip_gene_group_data,
  weights = var5)
richness_lme_weighted6 <- lme(
  relative_richness ~ 1 + Fraction*Nanotrap_type*Enrichment,
  random = ~ 1 | LIMS_ID,
  data = rpip_gene_group_data,
  weights = var6)
richness_lme_weighted7 <- lme(
  relative_richness ~ 1 + Fraction*Nanotrap_type*Enrichment,
  random = ~ 1 | LIMS_ID,
  data = rpip_gene_group_data,
  weights = var7)



#Basic version has clearly heterogeneous variance.
#Weighted versions do not; but they look similar to each other. 
plot(basic_richness_lme)
plot(richness_lme_weighted1)
plot(richness_lme_weighted2)
plot(richness_lme_weighted3)
plot(richness_lme_weighted4)
plot(richness_lme_weighted5)
plot(richness_lme_weighted6)
plot(richness_lme_weighted7)


#Compare the three models: full weighted version is better
anova(basic_richness_lme, 
      richness_lme_weighted1,
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
  data = rpip_gene_group_data)
richness_gls_weighted1 <- gls(
  relative_richness ~ 1 + Fraction*Nanotrap_type*Enrichment,
  data = rpip_gene_group_data,
  weights = var1)


anova(richness_lme_weighted1, richness_gls_weighted1)
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

###Now try dropping the other interactions###
R1.Base <- R.Drop1
R1.Drop1 <- update(R1.Base, .~. -Nanotrap_type:Enrichment)
anova(R1.Base, R1.Drop1)
#p = 0.0057; signif.

R1.Drop2 <- update(R1.Base, .~. -Fraction:Enrichment)
anova(R1.Base, R1.Drop2)
#p = .0035; signif.

R1.Drop3 <- update(R1.Base, .~. -Fraction:Nanotrap_type,
                   control = lmeControl(msMaxIter = 1000, msMaxEval = 1000))
anova(R1.Base, R1.Drop3)
#p = 0.16; not signif. R1.Drop3 is the new best model.


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
plot(x = rpip_gene_group_data$Nanotrap_type, y = E)
plot(x = rpip_gene_group_data$Fraction, y = E)
plot(x = rpip_gene_group_data$Enrichment, y = E)
plot(x = rpip_gene_group_data$LIMS_ID, y = E)


#Get model info:
anova(R2.REML)
summary(R2.REML)

rpip_gene_group_data_plot <- rpip_gene_group_data %>%
  mutate(Fraction = factor(Fraction, levels = c("Fil", "Ret", "Unfil")))

ggplot(rpip_gene_group_data, aes(x = Nanotrap_type, y = relative_richness*10^6, fill = Enrichment)) + geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3")) +
  theme_bw() +
  ylab("Unique ARGs per million reads") +
  xlab("Concentration method")


ggplot(rpip_gene_group_data_plot, aes(x = Fraction, y = relative_richness*10^6, fill = Enrichment)) + geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3")) +
  theme_bw() +
  ylab("Unique ARGs per million reads") +
  xlab("Fraction")


#Nanotrap_type ans Nanotrap_type:Enrichment may be insignificant for Shannon
SMod <- lme(shannon ~ 1 + Fraction + Nanotrap_type + Enrichment +
                         Fraction:Enrichment +
                         Nanotrap_type:Enrichment,
                       random = ~ 1 | LIMS_ID,
                       data = rpip_gene_group_data,
                       method = "REML",
                       weights = var1,
            control = lmeControl(msMaxIter = 1000, msMaxEval = 1000)
            )
plot(SMod)
anova(SMod)

ggplot(rpip_gene_group_data_plot, aes(x = Fraction, y = shannon, fill = Enrichment)) + geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3")) +
  theme_bw() +
  ylab("Shannon") +
  xlab("Fraction")

ggplot(rpip_gene_group_data_plot, aes(x = Nanotrap_type, y = shannon, fill = Enrichment)) + geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3")) +
  theme_bw() +
  ylab("Shannon") +
  xlab("Nanotrap_type")


#Nanotrap_type:Enrichment insignificant for Simpson
SiMod <- lme(simpson ~ 1 + Fraction + Nanotrap_type + Enrichment +
               Fraction:Enrichment +
               Nanotrap_type:Enrichment,
             random = ~ 1 | LIMS_ID,
             data = rpip_gene_group_data,
             method = "REML",
             weights = var1,
             control = lmeControl(msMaxIter = 1000, msMaxEval = 1000)
)

plot(SiMod)
anova(SiMod)

ggplot(rpip_gene_group_data_plot, aes(x = Fraction, y = simpson, fill = Enrichment)) + geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3")) +
  theme_bw() +
  ylab("Simpson") +
  xlab("Fraction")

ggplot(rpip_gene_group_data_plot, aes(x = Nanotrap_type, y = simpson, fill = Enrichment)) + geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3")) +
  theme_bw() +
  ylab("Simpson") +
  xlab("Nanotrap_type")


#Only Fraction is significant for Fisher
FiMod <- lme(fisher ~ 1 + Fraction + Nanotrap_type + Enrichment +
               Fraction:Enrichment +
               Nanotrap_type:Enrichment,
             random = ~ 1 | LIMS_ID,
             data = rpip_gene_group_data,
             method = "REML",
             weights = var1,
             control = lmeControl(msMaxIter = 1000, msMaxEval = 1000)
)

plot(FiMod)
anova(FiMod)

ggplot(rpip_gene_group_data_plot, aes(x = Fraction, y = fisher, fill = Enrichment)) + geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3")) +
  theme_bw() +
  ylab("Fisher") +
  xlab("Fraction")

ggplot(rpip_gene_group_data_plot, aes(x = Nanotrap_type, y = fisher, fill = Enrichment)) + geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3")) +
  theme_bw() +
  ylab("Fisher") +
  xlab("Nanotrap_type")


#Nanotrap_type:Enrichment is insignificant for Chao1
ChMod <- lme(chao1 ~ 1 + Fraction + Nanotrap_type + Enrichment +
               Fraction:Enrichment +
               Nanotrap_type:Enrichment,
             random = ~ 1 | LIMS_ID,
             data = rpip_gene_group_data,
             method = "REML",
             weights = var1,
             control = lmeControl(msMaxIter = 1000, msMaxEval = 1000)
)

plot(ChMod)
anova(ChMod)

ggplot(rpip_gene_group_data_plot, aes(x = Fraction, y = chao1, fill = Enrichment)) + geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3")) +
  theme_bw() +
  ylab("Chao1") +
  xlab("Fraction")

ggplot(rpip_gene_group_data_plot, aes(x = Nanotrap_type, y = chao1, fill = Enrichment)) + geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3")) +
  theme_bw() +
  ylab("Fisher") +
  xlab("Nanotrap_type")


