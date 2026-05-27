library(tidyverse)
library(nlme)

source("helper_functions.R")

all_k2_reports_anyConf_dedup <- read_csv(
  "imported_k2_reports/all_k2_reports_anyConf_dedup.csv"
  )


#Check how many entries there are for each:
table(all_k2_reports_anyConf_dedup$ribosomal)

#Confirm that RA adds to 1 in each category
all_k2_reports_anyConf_dedup %>%
  group_by(LIMS_ID, Treatment, Enrichment, ribosomal, dedup,
           Kraken2_confidence) %>%
  summarize(RA = sum(RA), nreads = sum(nodeOnly)) %>%
  View()


#Calculate the portion of all rRNA per sample:
rrna_portions <- all_k2_reports_anyConf_dedup %>% 
  #Add labels
  parse_sample_treatments() %>%
  #It doesn't matter which confidence level we use here, as long as only one is
  #included:
  filter(Kraken2_confidence == 0.0) %>%
  #Add IDs (unique for each sequence library)
  mutate(UniqueID = paste(LIMS_ID, Treatment, Enrichment, sep = "-")) %>%
  #Group by each variable that we want to keep:
  group_by(UniqueID, LIMS_ID, Treatment, Fraction, Nanotrap_type, Enrichment,
           ribosomal) %>%
  #Summarize rRNA and non-rRNA in each sample (without regard to taxa)
  summarize(nreads = sum(nodeOnly)) %>%
  #Make ribosomal and non-ribosomal cols:
  pivot_wider(names_from = ribosomal,
              values_from = nreads,
              id_cols = UniqueID:Enrichment) %>%
  mutate(totalreads = nonrrna + rrnaOnly) %>%
  mutate(rrna_portion = rrnaOnly / totalreads,
         nonrrna_portion = nonrrna / totalreads)


rrna_portion_plotdf <- rrna_portions %>%
  mutate(Fraction = str_replace_all(Fraction,
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


rrna_modeldata <- rrna_portion_plotdf %>%
  ungroup() %>%
  mutate(Fraction = relevel(factor(Fraction, ordered = FALSE), ref = "Unfil"),
         Nanotrap_type = relevel(factor(Nanotrap_type, ordered = FALSE), ref = "NTM-A"),
         Enrichment = relevel(factor(Enrichment, ordered = FALSE), ref = "Non-targeted"),
         LIMS_ID = as.factor(LIMS_ID)) %>%
  mutate(scaled_weights = (totalreads - min(totalreads)) / (max(totalreads) - min(totalreads)))

ggplot(rrna_modeldata, aes(x = rrna_portion, colour = Treatment)) + geom_density()


###############
####Models#####
###############

#Potential variance structures:
var1 <- varIdent(form = ~ 1 | Enrichment * Nanotrap_type * Fraction)
var2 <- varIdent(form = ~ 1 | Enrichment * Nanotrap_type)
var3 <- varIdent(form = ~ 1 | Enrichment * Fraction)
var4 <- varIdent(form = ~ 1 | Enrichment)
var5 <- varIdent(form = ~ 1 | Nanotrap_type * Fraction)
var6 <- varIdent(form = ~ 1 | Nanotrap_type)
var7 <- varIdent(form = ~ 1 | Fraction)


#Try basic LME baselines:
basic_rrna_lme <- lme(
  nonrrna_portion ~ 1 + Fraction*Nanotrap_type*Enrichment,
  random = ~ 1 | LIMS_ID,
  data = rrna_modeldata)
rrna_lme_weighted1 <- lme(
  nonrrna_portion ~ 1 + Fraction*Nanotrap_type*Enrichment,
  random = ~ 1 | LIMS_ID,
  data = rrna_modeldata,
  weights = var1,
  control = lmeControl(msMaxIter = 1000, msMaxEval = 1000))
rrna_lme_weighted2 <- lme(
  nonrrna_portion ~ 1 + Fraction*Nanotrap_type*Enrichment,
  random = ~ 1 | LIMS_ID,
  data = rrna_modeldata,
  weights = var2)
rrna_lme_weighted3 <- lme(
  nonrrna_portion ~ 1 + Fraction*Nanotrap_type*Enrichment,
  random = ~ 1 | LIMS_ID,
  data = rrna_modeldata,
  weights = var3)
rrna_lme_weighted4 <- lme(
  nonrrna_portion ~ 1 + Fraction*Nanotrap_type*Enrichment,
  random = ~ 1 | LIMS_ID,
  data = rrna_modeldata,
  weights = var4)
rrna_lme_weighted5 <- lme(
  nonrrna_portion ~ 1 + Fraction*Nanotrap_type*Enrichment,
  random = ~ 1 | LIMS_ID,
  data = rrna_modeldata,
  weights = var5)
rrna_lme_weighted6 <- lme(
  nonrrna_portion ~ 1 + Fraction*Nanotrap_type*Enrichment,
  random = ~ 1 | LIMS_ID,
  data = rrna_modeldata,
  weights = var6)
rrna_lme_weighted7 <- lme(
  nonrrna_portion ~ 1 + Fraction*Nanotrap_type*Enrichment,
  random = ~ 1 | LIMS_ID,
  data = rrna_modeldata,
  weights = var7)



#Compare the models: basic version is better
anova(basic_rrna_lme, 
      rrna_lme_weighted1)

#Weighted 2 is better than basic
anova(basic_rrna_lme,
      rrna_lme_weighted2)

#3 is better than 2 (very slightly)
anova(rrna_lme_weighted2,
      rrna_lme_weighted3)

#3 is better than 4
anova(rrna_lme_weighted3,
      rrna_lme_weighted4)

#5 is better than 3
anova(rrna_lme_weighted3,
      rrna_lme_weighted5)

#5 is better than 6
anova(rrna_lme_weighted5,
      rrna_lme_weighted6)


#5 is better than 7
anova(rrna_lme_weighted5,
      rrna_lme_weighted7)


#Try without random effect to be sure it's necessary:
basic_rrna_gls <- gls(
  nonrrna_portion ~ 1 + Fraction * Nanotrap_type * Enrichment, 
  data = rrna_modeldata)
rrna_gls_weighted5 <- gls(
  nonrrna_portion ~ 1 + Fraction*Nanotrap_type*Enrichment,
  data = rrna_modeldata,
  weights = var5)

anova(basic_rrna_gls, rrna_gls_weighted1)

######################################
#####Select variables for the LME#####
######################################

###Drop 3-way interaction###
#First the long-form model:
R.Base <- lme(nonrrna_portion ~ 1 + Fraction + Nanotrap_type + Enrichment +
                Fraction:Nanotrap_type +
                Fraction:Enrichment +
                Nanotrap_type:Enrichment +
                Fraction:Nanotrap_type:Enrichment,
              random = ~ 1 | LIMS_ID,
              data = rrna_modeldata,
              method = "ML",
              weights = var5,
              control = lmeControl(msMaxIter = 1000, msMaxEval = 1000)
)

#No 3-way interaction:
R.Drop1 <- lme(nonrrna_portion ~ 1 + Fraction + Nanotrap_type + Enrichment +
                 Fraction:Nanotrap_type +
                 Fraction:Enrichment +
                 Nanotrap_type:Enrichment,
               random = ~ 1 | LIMS_ID,
               data = rrna_modeldata,
               method = "ML",
               weights = var5)

anova(R.Base, R.Drop1)
#p = 0.52 so 3-way interaction is not significant; R.Drop1 is new base. 

###Now try dropping the other interactions###
R1.Base <- R.Drop1
R1.Drop1 <- update(R1.Base, .~. -Nanotrap_type:Enrichment)
anova(R1.Base, R1.Drop1)
#p = 0.3979; insignif.

R1.Drop2 <- update(R1.Base, .~. -Fraction:Enrichment)
anova(R1.Base, R1.Drop2)
#p < 0.1236; insignif.

R1.Drop3 <- update(R1.Base, .~. -Fraction:Nanotrap_type,
                   control = lmeControl(msMaxIter = 1000, msMaxEval = 1000))
anova(R1.Base, R1.Drop3)
#p = 0.0001; signif.
#R1.Drop1 has least significant difference; it is the new best model


R2.Base <- R1.Drop1

#Check if any insignificant terms remain (Fraction:Enrichment has p = 0.15):
anova(R2.Base)

R2.Drop1 <- update(R2.Base, .~. -Fraction:Enrichment)
anova(R2.Base, R2.Drop1)
#p = 0.14; insignif.

R2.Drop2 <- update(R2.Base, .~. -Fraction:Nanotrap_type)
anova(R2.Base, R2.Drop2)
#p < 0.0001; signif.
#R2.Drop1 is new best

R3.Base <- R2.Drop1
anova(R3.Base)
#All remaining terms are significant.


#Update the model method:
R3.REML <- update(R3.Base, method = "REML")

#Plot residuals for final model:
E <- resid(R3.REML, type = "normalized")
plot(x = fitted(R3.REML), y = E) #Looks pretty good

hist(E) #Looks good!

#Plot residuals by categorical variables:
plot(x = rrna_modeldata$Nanotrap_type, y = E)
plot(x = rrna_modeldata$Fraction, y = E)
plot(x = rrna_modeldata$Enrichment, y = E)
plot(x = rrna_modeldata$LIMS_ID, y = E)


#Get model info:
anova(R3.REML)
summary(R3.REML)

###############
#####Plots#####
###############

ggplot(rrna_portion_plotdf,
       aes(x = Nanotrap_type, y = nonrrna_portion * 100, fill = Enrichment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3", "gold2")) +
  theme_bw() +
  ylab("Percent non-ribosomal") +
  xlab("Nanotrap_type")


ggplot(rrna_portion_plotdf,
       aes(x = Fraction, y = nonrrna_portion * 100, fill = Enrichment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3", "gold2")) +
  theme_bw() +
  ylab("Percent non-ribosomal") +
  xlab("Fraction")