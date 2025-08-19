library(tidyverse)
library(ComplexHeatmap)
library(pheatmap)

get_contrast_table <- function(edger_glm) {
  bind_cols(edger_glm$genes, edger_glm$table)
}

get_contrast_table(VSP_filtrate_NanoAB_no_rrna_90conf_lf_removed_dedup) %>% View()

#All of these contrasts are for 0.9 confidence, deduplicated, 
#rRNA-removed samples:

###VSP, Direct Extraction
VSP_DE_Unf <- glmQLFTest(
  edger_family_fit_lfRemoved_90conf_noRrna_dedup, 
  contrast = (is.VSP & is.DirectExt & is.Unfiltered) - 
    (is.Untargeted & is.DirectExt & is.Unfiltered))

VSP_DE_Fil <- glmQLFTest(
  edger_family_fit_lfRemoved_90conf_noRrna_dedup, 
  contrast = (is.VSP & is.DirectExt & is.Filtrate) - 
    (is.Untargeted & is.DirectExt & is.Filtrate))

VSP_DE_Ret <- glmQLFTest(
  edger_family_fit_lfRemoved_90conf_noRrna_dedup, 
  contrast = (is.VSP & is.DirectExt & is.Retentate) - 
    (is.Untargeted & is.DirectExt & is.Retentate))


###VSP, Nanotrap A
VSP_NTA_Unf <- glmQLFTest(
  edger_family_fit_lfRemoved_90conf_noRrna_dedup, 
  contrast = (is.VSP & is.NanoA & is.Unfiltered) - 
    (is.Untargeted & is.NanoA & is.Unfiltered))

VSP_NTA_Fil <- glmQLFTest(
  edger_family_fit_lfRemoved_90conf_noRrna_dedup, 
  contrast = (is.VSP & is.NanoA & is.Filtrate) - 
    (is.Untargeted & is.NanoA & is.Filtrate))

VSP_NTA_Ret <- glmQLFTest(
  edger_family_fit_lfRemoved_90conf_noRrna_dedup, 
  contrast = (is.VSP & is.NanoA & is.Retentate) - 
    (is.Untargeted & is.NanoA & is.Retentate))

###VSP, Nanotrap A and B
VSP_NTAB_Unf <- glmQLFTest(
  edger_family_fit_lfRemoved_90conf_noRrna_dedup, 
  contrast = (is.VSP & is.NanoB & is.Unfiltered) - 
    (is.Untargeted & is.NanoB & is.Unfiltered))

VSP_NTAB_Fil <- glmQLFTest(
  edger_family_fit_lfRemoved_90conf_noRrna_dedup, 
  contrast = (is.VSP & is.NanoB & is.Filtrate) - 
    (is.Untargeted & is.NanoB & is.Filtrate))

VSP_NTAB_Ret <- glmQLFTest(
  edger_family_fit_lfRemoved_90conf_noRrna_dedup, 
  contrast = (is.VSP & is.NanoB & is.Retentate) - 
    (is.Untargeted & is.NanoB & is.Retentate))


#######RPIP
###RPIP, Direct Extraction
RPIP_DE_Unf <- glmQLFTest(
  edger_family_fit_lfRemoved_90conf_noRrna_dedup, 
  contrast = (is.RPIP & is.DirectExt & is.Unfiltered) - 
    (is.Untargeted & is.DirectExt & is.Unfiltered))

RPIP_DE_Fil <- glmQLFTest(
  edger_family_fit_lfRemoved_90conf_noRrna_dedup, 
  contrast = (is.RPIP & is.DirectExt & is.Filtrate) - 
    (is.Untargeted & is.DirectExt & is.Filtrate))

RPIP_DE_Ret <- glmQLFTest(
  edger_family_fit_lfRemoved_90conf_noRrna_dedup, 
  contrast = (is.RPIP & is.DirectExt & is.Retentate) - 
    (is.Untargeted & is.DirectExt & is.Retentate))


###RPIP, Nanotrap A
RPIP_NTA_Unf <- glmQLFTest(
  edger_family_fit_lfRemoved_90conf_noRrna_dedup, 
  contrast = (is.RPIP & is.NanoA & is.Unfiltered) - 
    (is.Untargeted & is.NanoA & is.Unfiltered))

RPIP_NTA_Fil <- glmQLFTest(
  edger_family_fit_lfRemoved_90conf_noRrna_dedup, 
  contrast = (is.RPIP & is.NanoA & is.Filtrate) - 
    (is.Untargeted & is.NanoA & is.Filtrate))

RPIP_NTA_Ret <- glmQLFTest(
  edger_family_fit_lfRemoved_90conf_noRrna_dedup, 
  contrast = (is.RPIP & is.NanoA & is.Retentate) - 
    (is.Untargeted & is.NanoA & is.Retentate))

###RPIP, Nanotrap A and B
RPIP_NTAB_Unf <- glmQLFTest(
  edger_family_fit_lfRemoved_90conf_noRrna_dedup, 
  contrast = (is.RPIP & is.NanoB & is.Unfiltered) - 
    (is.Untargeted & is.NanoB & is.Unfiltered))

RPIP_NTAB_Fil <- glmQLFTest(
  edger_family_fit_lfRemoved_90conf_noRrna_dedup, 
  contrast = (is.RPIP & is.NanoB & is.Filtrate) - 
    (is.Untargeted & is.NanoB & is.Filtrate))

RPIP_NTAB_Ret <- glmQLFTest(
  edger_family_fit_lfRemoved_90conf_noRrna_dedup, 
  contrast = (is.RPIP & is.NanoB & is.Retentate) - 
    (is.Untargeted & is.NanoB & is.Retentate))

VSP_results <- list(
  VSP_DE_Unf = VSP_DE_Unf,
  VSP_DE_Fil = VSP_DE_Fil,
  VSP_DE_Ret = VSP_DE_Ret,
  VSP_NTA_Unf = VSP_NTA_Unf,
  VSP_NTA_Fil = VSP_NTA_Fil,
  VSP_NTA_Ret = VSP_NTA_Ret,
  VSP_NTAB_Unf = VSP_NTAB_Unf,
  VSP_NTAB_Fil = VSP_NTAB_Fil,
  VSP_NTAB_Ret = VSP_NTAB_Ret
)
 
VSP_results <- lapply(VSP_results, get_contrast_table) %>%
  imap(~ mutate(.x, ID = .y)) %>%
  bind_rows()

VSP_results <- VSP_results %>%
  separate(ID, into = c("Enrichment", "Nanotrap_type", "Fraction"),
           sep = "_")  %>% 
  mutate(FDR = p.adjust(PValue, method = "BH"))


signif <- VSP_results %>%
  filter(FDR <= 0.05)
most_important <- signif %>%
  group_by(Family) %>%
  summarize(mean_logfc = mean(logFC)) %>%
  arrange(desc(mean_logfc))
most_important <- most_important[str_detect(most_important$Family, pattern = "^*viridae$"),]




VSP_results_important <- VSP_results %>%
  filter(Family %in% most_important$Family)


ggplot(VSP_results_important, aes(x = Nanotrap_type, 
                                  y = Fraction, 
                                  fill = logFC)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkred") + 
  facet_wrap(~ Family, ncol = 3) +
  theme_minimal()

VSP_results_important <- VSP_results_important %>% 
  mutate(FDR = p.adjust(PValue, method = "BH"))

VSP_results_important$stars <- cut(
  VSP_results_important$FDR,
  breaks = c(-Inf, 0.00001, 0.0001, 0.001, 0.05, Inf),
  labels = c("****", "***", "**", "*", "-")
)



ggplot(VSP_results_important, aes(x = Nanotrap_type, 
                                  y = Fraction, 
                                  fill = logFC)) +
  geom_tile(color = "black", linewidth = .5) +
  scale_fill_gradient(low = "white", high = "darkred") + 
  facet_wrap(~ Family, ncol = 3) +
  theme_minimal() +
  theme(strip.text = element_text(size = 12), 
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 10))+
  geom_text(aes(label = stars), size = 5, vjust = 0.5, fontface = "bold") +
  xlab("Concentration method")



###########################33
############################3
#############################
RPIP_results <- list(
  RPIP_DE_Unf = RPIP_DE_Unf,
  RPIP_DE_Fil = RPIP_DE_Fil,
  RPIP_DE_Ret = RPIP_DE_Ret,
  RPIP_NTA_Unf = RPIP_NTA_Unf,
  RPIP_NTA_Fil = RPIP_NTA_Fil,
  RPIP_NTA_Ret = RPIP_NTA_Ret,
  RPIP_NTAB_Unf = RPIP_NTAB_Unf,
  RPIP_NTAB_Fil = RPIP_NTAB_Fil,
  RPIP_NTAB_Ret = RPIP_NTAB_Ret
)

RPIP_results <- lapply(RPIP_results, get_contrast_table) %>%
  imap(~ mutate(.x, ID = .y)) %>%
  bind_rows()

RPIP_results <- RPIP_results %>%
  separate(ID, into = c("Enrichment", "Nanotrap_type", "Fraction"),
           sep = "_")  %>% 
  mutate(FDR = p.adjust(PValue, method = "BH"))


RPIP_signif <- RPIP_results %>%
  filter(FDR <= 0.05)
most_important <- signif %>%
  group_by(Family) %>%
  summarize(mean_logfc = mean(logFC)) %>%
  arrange(desc(mean_logfc))
most_important <- most_important[str_detect(most_important$Family, pattern = "^*viridae$"),]




RPIP_results_important <- RPIP_results %>%
  filter(Family %in% RPIP_signif$Family)


ggplot(RPIP_results_important, aes(x = Nanotrap_type, 
                                  y = Fraction, 
                                  fill = logFC)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkred") + 
  facet_wrap(~ Family, ncol = 3) +
  theme_minimal()

RPIP_results_important <- RPIP_results_important %>% 
  mutate(FDR = p.adjust(PValue, method = "BH"))

RPIP_results_important$stars <- cut(
  RPIP_results_important$FDR,
  breaks = c(-Inf, 0.00001, 0.0001, 0.001, 0.05, Inf),
  labels = c("****", "***", "**", "*", "-")
)



ggplot(RPIP_results_important, aes(x = Nanotrap_type, 
                                  y = Fraction, 
                                  fill = logFC)) +
  geom_tile(color = "black", linewidth = .5) +
  scale_fill_gradient(low = "white", high = "darkred") + 
  facet_wrap(~ Family, ncol = 3) +
  theme_minimal() +
  theme(strip.text = element_text(size = 12), 
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 10))+
  geom_text(aes(label = stars), size = 5, vjust = 0.5, fontface = "bold") +
  xlab("Concentration method")
