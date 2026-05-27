

#####
#Alpha diversity

rpip_group_data$richness <- specnumber(rpip_bacteria_community_matrix_noNA)
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
         chao1 = estimateR(rpip_bacteria_community_matrix_noNA)["S.chao1",]
         ) %>%
  mutate(Nanotrap_type = relevel(factor(Nanotrap_type, ordered = FALSE), ref = "A"),
         Fraction = relevel(factor(Fraction, ordered = FALSE), ref = "unfiltered"),
         LIMS_ID = factor(LIMS_ID))




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
