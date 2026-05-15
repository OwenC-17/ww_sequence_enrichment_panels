library(tidyverse)
library(nlme)
library(emmeans)
library(cowplot)

#Data from read_fastp_reports.R
all_dedup_reports <- read_rds("input/modified/all_fastp_reports_dupdedup.rds")

#Create an ID that contains treatments before enrichment:
all_dedup_reports$sample_id <- paste(
  all_dedup_reports$LIMS_ID,
  all_dedup_reports$Fraction,
  all_dedup_reports$Nanotrap_type,
  sep = "-"
)

#Change labels for consistency with figures:
all_dedup_reports <- all_dedup_reports %>%
  mutate(Fraction = str_replace_all(Fraction,
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

#Factorize explanatory variables:
all_dedup_reports$LIMS_ID <- factor(all_dedup_reports$LIMS_ID, ordered = FALSE)
all_dedup_reports$Nanotrap_type <- relevel(
  factor(all_dedup_reports$Nanotrap_type, ordered = FALSE),
  ref = "DirEx")
all_dedup_reports$Fraction <- relevel(
  factor(all_dedup_reports$Fraction, ordered = FALSE),
  ref = "Unfil")
all_dedup_reports$Enrichment <- relevel(
  factor(all_dedup_reports$Enrichment, ordered = FALSE),
  ref = "Non-targeted")


#Change "rate" variable because it's not very intuitive:
all_dedup_reports <- all_dedup_reports %>%
  rename("duprate" = "rate")


#Look at the distribution of duprate:
hist(all_dedup_reports$duprate)
#It seems "close enough" to normal overall. 

#Try by category:
ggplot(all_dedup_reports, aes(x = duprate, colour = Enrichment)) + geom_density()
#They aren't quite normal, but normal fit may be fine.

##########################
#####Fit linear model#####
##########################
library(nlme)
var1 <- varIdent(form = ~ 1 | Enrichment * Nanotrap_type * Fraction)
var2 <- varIdent(form = ~ 1 | Enrichment * Nanotrap_type)
var3 <- varIdent(form = ~ 1 | Enrichment * Fraction)
var4 <- varIdent(form = ~ 1 | Enrichment)
var5 <- varIdent(form = ~ 1 | Nanotrap_type * Fraction)
var6 <- varIdent(form = ~ 1 | Nanotrap_type)
var7 <- varIdent(form = ~ 1 | Fraction)




basic_duprate_lme <- lme(duprate ~ 1 + Fraction*Nanotrap_type*Enrichment, random = ~ 1 | LIMS_ID, data = all_dedup_reports)
basic_duprate_lme_weighted1 <- lme(duprate ~ 1 + Fraction*Nanotrap_type*Enrichment, random = ~ 1 | LIMS_ID, data = all_dedup_reports, weights = var1, control = lmeControl(msMaxIter = 1000, msMaxEval = 1000))
basic_duprate_lme_weighted2 <- lme(duprate ~ 1 + Fraction*Nanotrap_type*Enrichment, random = ~ 1 | LIMS_ID, data = all_dedup_reports, weights = var2)
basic_duprate_lme_weighted3 <- lme(duprate ~ 1 + Fraction*Nanotrap_type*Enrichment, random = ~ 1 | LIMS_ID, data = all_dedup_reports, weights = var3)
basic_duprate_lme_weighted4 <- lme(duprate ~ 1 + Fraction*Nanotrap_type*Enrichment, random = ~ 1 | LIMS_ID, data = all_dedup_reports, weights = var4)
basic_duprate_lme_weighted5 <- lme(duprate ~ 1 + Fraction*Nanotrap_type*Enrichment, random = ~ 1 | LIMS_ID, data = all_dedup_reports, weights = var5)
basic_duprate_lme_weighted6 <- lme(duprate ~ 1 + Fraction*Nanotrap_type*Enrichment, random = ~ 1 | LIMS_ID, data = all_dedup_reports, weights = var6)
basic_duprate_lme_weighted7 <- lme(duprate ~ 1 + Fraction*Nanotrap_type*Enrichment, random = ~ 1 | LIMS_ID, data = all_dedup_reports, weights = var7)


plot(basic_duprate_lme) #Residuals seem decent, maybe a bit heterogeneous
plot(basic_duprate_lme_weighted1) #Definitely more visually homogeneous than unweighted
plot(basic_duprate_lme_weighted2) #Seems in between previous two
plot(basic_duprate_lme_weighted3) #Also in between
plot(basic_duprate_lme_weighted4) #Seems pretty heterogeneous
plot(basic_duprate_lme_weighted5) #Similar
plot(basic_duprate_lme_weighted6) #Similar to unweighted
plot(basic_duprate_lme_weighted7) #Also similar to unweighted

anova(basic_duprate_lme,
      basic_duprate_lme_weighted1
      ) #Unweighted better

anova(basic_duprate_lme,
      basic_duprate_lme_weighted2
) #Weighted 2 better

anova(basic_duprate_lme_weighted2,
      basic_duprate_lme_weighted3
) #Weighted 3 better

anova(basic_duprate_lme_weighted3,
      basic_duprate_lme_weighted4
) #Weighted 3 better

anova(basic_duprate_lme_weighted3,
      basic_duprate_lme_weighted5
) #Weighted 3 better

anova(basic_duprate_lme_weighted3,
      basic_duprate_lme_weighted6
) #Weighted 3 better

anova(basic_duprate_lme_weighted3,
      basic_duprate_lme_weighted7
) #Weighted 7 better but not significantly (p = 0.134)


#No random effect:
basic_duprate_gls <- gls(duprate ~ 1 + Fraction*Nanotrap_type*Enrichment, data = all_dedup_reports)
basic_duprate_gls_weighted1 <- gls(duprate ~ 1 + Fraction*Nanotrap_type*Enrichment, data = all_dedup_reports, weights = var1)
basic_duprate_gls_weighted2 <- gls(duprate ~ 1 + Fraction*Nanotrap_type*Enrichment, data = all_dedup_reports, weights = var2)
basic_duprate_gls_weighted3 <- gls(duprate ~ 1 + Fraction*Nanotrap_type*Enrichment, data = all_dedup_reports, weights = var3)
basic_duprate_gls_weighted4 <- gls(duprate ~ 1 + Fraction*Nanotrap_type*Enrichment, data = all_dedup_reports, weights = var4)
basic_duprate_gls_weighted5 <- gls(duprate ~ 1 + Fraction*Nanotrap_type*Enrichment, data = all_dedup_reports, weights = var5)
basic_duprate_gls_weighted6 <- gls(duprate ~ 1 + Fraction*Nanotrap_type*Enrichment, data = all_dedup_reports, weights = var6)
basic_duprate_gls_weighted7 <- gls(duprate ~ 1 + Fraction*Nanotrap_type*Enrichment, data = all_dedup_reports, weights = var7)



anova(basic_duprate_gls, 
      basic_duprate_gls_weighted1
      ) #Unweighted better 

anova(basic_duprate_gls, 
      basic_duprate_gls_weighted2
) #Unweighted better

anova(basic_duprate_gls, 
      basic_duprate_gls_weighted3
) #Unweighted better

anova(basic_duprate_gls, 
      basic_duprate_gls_weighted4
) #Weighted 4 better

anova(basic_duprate_gls_weighted4, 
      basic_duprate_gls_weighted5
) #Weighted 5 better

anova(basic_duprate_gls_weighted5, 
      basic_duprate_gls_weighted6
) #Weighted 5 better


anova(basic_duprate_gls_weighted5, 
      basic_duprate_gls_weighted7
) #Weighted 7 better


plot(basic_duprate_gls) #Now there does seem to be uneven variance in the residuals
plot(x = basic_duprate_gls$fitted, y = basic_duprate_gls$residuals, col = all_dedup_reports$Enrichment)
plot(x = basic_duprate_lme$fitted, y = basic_duprate_lme$residuals, col = all_dedup_reports$Enrichment)
#Visually it looks like variance depends on enrichment

anova(basic_duprate_gls_weighted7, basic_duprate_lme_weighted7)
#The lme is better; the random effect is very significant (p < 0.0001)

#Drop 3-way interaction
#First the long-form model:
Dpr.Base <- lme(duprate ~ 1 + Fraction + Nanotrap_type + Enrichment +
                Fraction:Nanotrap_type +
                Fraction:Enrichment +
                Nanotrap_type:Enrichment +
                Fraction:Nanotrap_type:Enrichment,
              random = ~ 1 | LIMS_ID,
              data = all_dedup_reports,
              method = "ML",
              weights = var7)

#No 3-way interaction:
Dpr.Drop1 <- lme(duprate ~ 1 + Fraction + Nanotrap_type + Enrichment +
                 Fraction:Nanotrap_type +
                 Fraction:Enrichment +
                 Nanotrap_type:Enrichment,
               random = ~ 1 | LIMS_ID,
               data = all_dedup_reports,
               method = "ML",
               weights = var7)

anova(Dpr.Base, Dpr.Drop1)
#p = 0.0937 so 3-way interaction is not significant; Dpr.Drop1 is new base. 

#Now try dropping the other interactions.
Dpr1.Base <- Dpr.Drop1
anova(Dpr1.Base) #Find least significant interaction (Nanotrap_type:Enrichment)

Dpr1.Drop1 <- update(Dpr1.Base, .~. -Nanotrap_type:Enrichment)
anova(Dpr1.Base, Dpr1.Drop1)
#p = 0.7316; not signif.

Dpr1.Drop2 <- update(Dpr1.Base, .~. -Fraction:Enrichment)
anova(Dpr1.Base, Dpr1.Drop2)
#p = .001; marginally signif.

Dpr1.Drop3 <- update(Dpr1.Base, .~. -Fraction:Nanotrap_type)
anova(Dpr1.Base, Dpr1.Drop3)
#p = 0.0116; marginally signif. Dpr1.Drop1 is the new best model.

#Stage 2:
Dpr2.Base <- Dpr1.Drop1
anova(Dpr2.Base)
#Now all terms are significant.
E <- resid(Dpr2.Base, type = "normalized")
plot(x = fitted(Dpr2.Base), y = E) #Looks pretty good

hist(E) #Looks okay

anova(Dpr2.Base)
summary(Dpr2.Base)
plot(x = all_dedup_reports$Nanotrap_type, y = E)
plot(x = all_dedup_reports$LIMS_ID, y = E)

contrast(emmeans(Dpr2.Base, ~ Enrichment), method = "pairwise")
contrast(emmeans(Dpr2.Base, ~ Fraction), method = "pairwise")
contrast(emmeans(Dpr2.Base, ~ Nanotrap_type), method = "pairwise")
contrast(emmeans(Dpr2.Base, ~ Fraction | Enrichment), method = "pairwise")
contrast(emmeans(Dpr2.Base, ~ Nanotrap_type | Enrichment), method = "pairwise")
contrast(emmeans(Dpr2.Base, ~ Nanotrap_type | Enrichment), method = "pairwise")
contrast(emmeans(Dpr2.Base, ~ Nanotrap_type | Fraction), method = "pairwise")
contrast(emmeans(Dpr2.Base, ~ Fraction | Nanotrap_type), method = "pairwise")

contrast(emmeans(Dpr2.Base, ~ Fraction:Nanotrap_type | Enrichment), method = "pairwise")

################################
#####Plot duplication rates#####
################################
nt_dpr_boxplot <- ggplot(all_dedup_reports, aes(x = Nanotrap_type, y = duprate, fill = Enrichment)) + geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3", "goldenrod")) +
  theme_bw() +
  ylab("Read duplication rate") +
  xlab("Concentration method") +
  theme(legend.title = element_blank())

fraction_dpr_boxplot <- ggplot(all_dedup_reports, aes(x = Fraction, y = duprate, fill = Enrichment)) + geom_boxplot() +
  scale_fill_manual(values = c("deeppink3", "royalblue3", "goldenrod")) +
  theme_bw() +
  ylab("Read duplication rate") +
  xlab("Fraction") +
  theme(legend.position = "none")

cowplot::plot_grid(fraction_dpr_boxplot, nt_dpr_boxplot,
                   nrow = 2,
                   align = "hv")