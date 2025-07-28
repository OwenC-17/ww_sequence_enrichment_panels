#Untargeted vs VSP:
if(exists("read_coverage_vs_ref_coverage_unt_vs_vsp")){
  rm(read_coverage_vs_ref_coverage_unt_vs_vsp)
}

for (cvt in list.files(path = unt_vs_vsp_dir, pattern = "*.tsv$")) {
  x = import_readportion_coverage_of_vsp_ref(paste0(unt_vs_vsp_dir, cvt))
  if(exists("read_coverage_vs_ref_coverage_unt_vs_vsp")) {
    read_coverage_vs_ref_coverage_unt_vs_vsp <- bind_rows(read_coverage_vs_ref_coverage_unt_vs_vsp, x)
  } else {
    read_coverage_vs_ref_coverage_unt_vs_vsp <- x
  }
}

read_coverage_vs_ref_coverage_unt_vs_vsp$Enrichment <- "None"

read_coverage_vs_ref_coverage_both_vs_vsp <- bind_rows(
  read_coverage_vs_ref_coverage_vsp_vs_vsp,
  read_coverage_vs_ref_coverage_unt_vs_vsp
)


group_data_90_conf_with_rrna_for_maps <- group_data_90conf

libsizes <- all_paired_fastp_summaries %>%
  select(LIMS_ID, Treatment, Enrichment, after_filtering.total_reads) %>%
  mutate(LIMS_ID = as.integer(LIMS_ID))

group_data_90_conf_with_rrna_for_maps <- group_data_90_conf_with_rrna_for_maps %>%
  left_join(libsizes, by = c("LIMS_ID", "Treatment", "Enrichment"))

group_data_90_conf_with_rrna_for_maps <- group_data_90_conf_with_rrna_for_maps %>%
  rename(lib.size = after_filtering.total_reads)

group_data_90_conf_with_rrna_for_maps$UniqueID <- str_remove(group_data_90_conf_with_rrna_for_maps$UniqueID,
                                                             "\\-0\\.9$")

view(group_data_90_conf_with_rrna_for_maps)

rpip_v_unt_mapped_group_data <- group_data_90_conf_with_rrna_for_maps %>%
  filter(Enrichment != "VSP")

vsp_v_unt_mapped_group_data <- group_data_90_conf_with_rrna_for_maps %>%
  filter(Enrichment != "RPIP")


colnames(all_taxa_90rc_both_vs_vsp)
edger_counts_all_taxa_90rc_both_vs_vsp <- all_taxa_90rc_both_vs_vsp %>%
  unite("UniqueID", LIMS_ID, Treatment, Enrichment, sep = "-", remove = FALSE) %>%
  select(UniqueID, rname, numreads) %>%
#filter(edger_counts_all_taxa_90rc_both_vs_vsp, str_starts(UniqueID, "20040-FD-VSP"))
  pivot_wider(names_from = UniqueID, values_from = numreads, id_cols = rname,
              values_fill = 0)

#colSums(edger_counts_all_taxa_90rc_both_vs_vsp[,2:ncol(edger_counts_all_taxa_90rc_both_vs_vsp)])

sum(!vsp_v_unt_mapped_group_data$UniqueID %in% colnames(edger_counts_all_taxa_90rc_both_vs_vsp))

edger_counts_all_taxa_90rc_both_vs_vsp <- edger_counts_all_taxa_90rc_both_vs_vsp[,c("rname", vsp_v_unt_mapped_group_data$UniqueID)]

sum(vsp_v_unt_mapped_group_data$UniqueID != colnames(edger_counts_all_taxa_90rc_both_vs_vsp[2:ncol(edger_counts_all_taxa_90rc_both_vs_vsp)]))

dge_counts_all_taxa_90rc_both_vs_vsp <- DGEList(counts = edger_counts_all_taxa_90rc_both_vs_vsp, lib.size = vsp_v_unt_mapped_group_data$lib.size)

generate_levels_for_readmappings <- function(group_df) {
  LIMS_ID <- factor(group_df$LIMS_ID)
  Fraction <- factor(group_df$Fraction, levels = c("unfiltered", "retentate", "filtrate"))
  Nanotrap_type <- factor(group_df$Nanotrap_type, levels = c("none", "A", "A&B"))
  Enrichment <- factor(group_df$Enrichment)
  Enrichment <- relevel(Enrichment, ref = "None")
  treat_list <- list()
  for (x in levels(Fraction)) {
    for (y in levels(Nanotrap_type)) {
      for (z in levels(Enrichment)) {
        tx_name <- paste(x, y, z, sep = ".")
        assign(tx_name, (Fraction == x & Nanotrap_type == y & Enrichment == z))
        treat_list <- c(treat_list, tx_name)
      }
    }
  }
  additional_groups <- mget(unlist(treat_list))
  design <- model.matrix(~LIMS_ID)
  for (group in 1:length(additional_groups)) {
    design <- cbind(design, unlist(additional_groups[group]))
  }
  colnames(design)[(length(levels(LIMS_ID)) + 1):ncol(design)] <- names(additional_groups)
  design <- design[, colnames(design) != "unfiltered.none.None"]
  return(design)
}

design_mapped_reads <- generate_levels_for_readmappings(vsp_v_unt_mapped_group_data)

expressed_mapped <- filterByExpr(dge_counts_all_taxa_90rc_both_vs_vsp, design = design_mapped_reads)
mapped_expressed <- dge_counts_all_taxa_90rc_both_vs_vsp[expressed_mapped, , keep.lib.sizes = TRUE]

mapped_expressed$genes$taxid <- accessionToTaxa(mapped_expressed$genes$rname, sqlFile = sqlPath)

genes <- mapped_expressed$genes

unique_accessions <- unique(genes$rname)
  
taxonomyDf <- data.frame(rname = unique_accessions)
taxonomyDf$taxid <- accessionToTaxa(taxonomyDf$rname, sqlPath)
taxonomyDf$taxid <- as.character(taxonomyDf$taxid)
  
uniqueTaxId <- unique(taxonomyDf$taxid)
  
rawTaxa <- getRawTaxonomy(uniqueTaxId, sqlPath)
rawTaxaDf <- bind_rows(rawTaxa, .id = "taxid")
rawTaxaDf$taxid <- str_remove_all(rawTaxaDf$taxid, "\\s")
  
accessions2taxids <- taxonomyDf %>%
  filter(!is.na(taxid))
  
rawTaxaDf <- rawTaxaDf %>%
  full_join(accessions2taxids, by = "taxid")
  
rawTaxaDf <- rawTaxaDf %>%
  unite("tname",c("species", "no rank", "no rank.1", 
                    "no rank.2", "clade", 
                    "clade.1"), sep = " / ", na.rm = TRUE, 
          remove = FALSE) 
  
genes <- genes %>%
  left_join(rawTaxaDf, by = c("rname"))
  
  

mapped_expressed$genes <- genes

disp_mapped <- estimateDisp(y = mapped_expressed, design = design_mapped_reads)
disp_mapped$common.dispersion
fit_mapped <- glmQLFit(disp_mapped, design_mapped_reads)
plotBCV(disp_mapped)

is.VSP <- str_detect(colnames(design_mapped_reads), "VSP")
is.Untargeted <- str_detect(colnames(design_mapped_reads), "None")
is.Filtrate <- str_detect(colnames(design_mapped_reads), "filtrate")
is.Retentate <- str_detect(colnames(design_mapped_reads), "retentate")
is.Unfiltered <- str_detect(colnames(design_mapped_reads), "unfiltered")
is.NanoA <- str_detect(colnames(design_mapped_reads), "A\\.")
is.NanoB <- str_detect(colnames(design_mapped_reads), "A&B")
is.DirectExt <- str_detect(colnames(design_mapped_reads), "none")


#####VSP - Untargeted
#Main effect:
VSP_all_treatments_mapped <- glmQLFTest(fit_mapped, contrast = is.VSP - is.Untargeted)
top25(VSP_all_treatments_mapped)$table %>% View()

#No NT, no filtration:
VSP_vs_Unt_DirExt_Unfiltered_mapped <- glmQLFTest(fit_mapped, contrast = (is.VSP & is.Unfiltered & is.DirectExt) - (is.Untargeted & is.Unfiltered & is.DirectExt))
top25(VSP_vs_Unt_DirExt_Unfiltered_mapped)$table %>% View()

#NTA, no filtration:
VSP_vs_Unt_NanoA_Unfiltered_mapped <- glmQLFTest(fit_mapped, contrast = (is.VSP & is.Unfiltered & is.NanoA) - (is.Untargeted & is.Unfiltered & is.NanoA))
top25(VSP_vs_Unt_NanoA_Unfiltered_mapped)$table %>% View()

#NTB, no filtration:
VSP_vs_Unt_NanoBA_Unfiltered_mapped <- glmQLFTest(fit_mapped, contrast = (is.VSP & is.NanoB & is.Unfiltered) - (is.Untargeted & is.NanoB & is.Unfiltered))
top25(VSP_vs_Unt_NanoBA_Unfiltered_mapped)$table %>% View()

#No NT, filtrate:
VSP_vs_Unt_DirExt_Filtrate_mapped <- glmQLFTest(fit_mapped, contrast = (is.VSP & is.Filtrate & is.DirectExt) - (is.Untargeted & is.Filtrate & is.DirectExt))
top25(VSP_vs_Unt_DirExt_Filtrate_mapped)$table %>% View()

#NTA, filtrate:
VSP_vs_Unt_NanoA_Filtrate_mapped <- glmQLFTest(fit_mapped, contrast = (is.VSP & is.NanoA & is.Filtrate) - (is.Untargeted & is.NanoA & is.Filtrate))
top25(VSP_vs_Unt_NanoA_Filtrate_mapped)$table %>% View()

#NTB, filtrate:
VSP_vs_Unt_NanoBA_Filtrate_mapped <- glmQLFTest(fit_mapped, contrast = (is.VSP & is.NanoB & is.Filtrate) - (is.Untargeted & is.NanoB & is.Filtrate))
top25(VSP_vs_Unt_NanoBA_Filtrate_mapped)$table %>% View()

#No NT, retentate:
VSP_vs_Unt_DirExt_Retentate_mapped <- glmQLFTest(fit_mapped, contrast = (is.VSP & is.Retentate & is.DirectExt) - (is.Untargeted & is.Retentate & is.DirectExt))
top25(VSP_vs_Unt_DirExt_Retentate_mapped)$table %>% View()

#NTA, retentate:
VSP_vs_Unt_NanoA_Retentate_mapped <- glmQLFTest(fit_mapped, contrast = (is.VSP & is.NanoA & is.Retentate) - (is.Untargeted & is.NanoA & is.Retentate))
top25(VSP_vs_Unt_NanoA_Retentate_mapped)$table %>% View()

#NTB, retentate:
VSP_vs_Unt_NanoBA_Retentate_mapped <- glmQLFTest(fit_mapped, contrast = (is.VSP & is.NanoB & is.Retentate) - (is.Untargeted & is.NanoB & is.Retentate))
top25(VSP_vs_Unt_NanoBA_Retentate_mapped)$table %>% View()


########################
##########Nanotrap contrasts
########################
#NanoA main effect:
NanoA_vs_DirExt_main_effect_mapped <- glmQLFTest(fit_mapped, contrast = is.NanoA - is.DirectExt)
top25(NanoA_vs_DirExt_main_effect_mapped)

#NanoA, unfiltered, untargeted:
NanoA_vs_DirExt_Unfiltered_Untargeted_mapped <- glmQLFTest(fit_mapped, contrast = (is.NanoA & is.Unfiltered & is.Untargeted) - (is.DirectExt & is.Unfiltered & is.Untargeted))
top25(NanoA_vs_DirExt_Unfiltered_Untargeted_mapped)

#NanoA, unfiltered, VSP:
NanoA_vs_DirExt_Unfiltered_VSP_mapped <- glmQLFTest(fit_mapped, contrast = (is.NanoA & is.Unfiltered & is.VSP) - (is.DirectExt & is.Unfiltered & is.VSP))
top25(NanoA_vs_DirExt_Unfiltered_VSP_mapped)

#NanoA, unfiltered, RPIP:
NanoA_vs_DirExt_Unfiltered_RPIP_mapped <- glmQLFTest(fit_mapped, contrast = (is.NanoA & is.Unfiltered & is.RPIP) - (is.DirectExt & is.Unfiltered & is.RPIP))
top25(NanoA_vs_DirExt_Unfiltered_RPIP_mapped)

#NanoA, filtrate, untargeted:
NanoA_vs_DirExt_Filtrate_Untargeted_mapped <- glmQLFTest(fit_mapped, contrast = (is.NanoA & is.Filtrate & is.Untargeted) - (is.DirectExt & is.Filtrate & is.Untargeted))
top25(NanoA_vs_DirExt_Filtrate_Untargeted_mapped)

#NanoA, filtrate, VSP:
NanoA_vs_DirExt_Filtrate_VSP_mapped <- glmQLFTest(fit_mapped, contrast = (is.NanoA & is.Filtrate & is.VSP) - (is.DirectExt & is.Filtrate & is.VSP))
top25(NanoA_vs_DirExt_Filtrate_VSP_mapped)

#NanoA, filtrate, RPIP:
NanoA_vs_DirExt_Filtrate_RPIP_mapped <- glmQLFTest(fit_mapped, contrast = (is.NanoA & is.Filtrate & is.RPIP) - (is.DirectExt & is.Filtrate & is.RPIP))
top25(NanoA_vs_DirExt_Filtrate_RPIP_mapped)

#NanoA, Retentate, untargeted:
NanoA_vs_DirExt_Retentate_Untargeted_mapped <- glmQLFTest(fit_mapped, contrast = (is.NanoA & is.Retentate & is.Untargeted) - (is.DirectExt & is.Retentate & is.Untargeted))
top25(NanoA_vs_DirExt_Retentate_Untargeted_mapped)

#NanoA, Retentate, VSP:
NanoA_vs_DirExt_Retentate_VSP_mapped <- glmQLFTest(fit_mapped, contrast = (is.NanoA & is.Retentate & is.VSP) - (is.DirectExt & is.Retentate & is.VSP))
top25(NanoA_vs_DirExt_Retentate_VSP_mapped)

#NanoA, Retentate, RPIP:
NanoA_vs_DirExt_Retentate_RPIP_mapped <- glmQLFTest(fit_mapped, contrast = (is.NanoA & is.Retentate & is.RPIP) - (is.DirectExt & is.Retentate & is.RPIP))
top25(NanoA_vs_DirExt_Retentate_RPIP_mapped)

#NanoB main effect:
NanoB_vs_DirExt_main_effect_mapped <- glmQLFTest(fit_mapped, contrast = is.NanoB - is.DirectExt)
top25(NanoB_vs_DirExt_main_effect_mapped)

#NanoB, unfiltered, untargeted:
NanoB_vs_DirExt_Unfiltered_Untargeted_mapped <- glmQLFTest(fit_mapped, contrast = (is.NanoB & is.Unfiltered & is.Untargeted) - (is.DirectExt & is.Unfiltered & is.Untargeted))
top25(NanoB_vs_DirExt_Unfiltered_Untargeted_mapped)

#NanoB, unfiltered, VSP:
NanoB_vs_DirExt_Unfiltered_VSP_mapped <- glmQLFTest(fit_mapped, contrast = (is.NanoB & is.Unfiltered & is.VSP) - (is.DirectExt & is.Unfiltered & is.VSP))
top25(NanoB_vs_DirExt_Unfiltered_VSP_mapped)

#NanoB, unfiltered, RPIP:
NanoB_vs_DirExt_Unfiltered_RPIP_mapped <- glmQLFTest(fit_mapped, contrast = (is.NanoB & is.Unfiltered & is.RPIP) - (is.DirectExt & is.Unfiltered & is.RPIP))
top25(NanoB_vs_DirExt_Unfiltered_RPIP_mapped)

#NanoB, filtrate, untargeted:
NanoB_vs_DirExt_Filtrate_Untargeted_mapped <- glmQLFTest(fit_mapped, contrast = (is.NanoB & is.Filtrate & is.Untargeted) - (is.DirectExt & is.Filtrate & is.Untargeted))
top25(NanoB_vs_DirExt_Filtrate_Untargeted_mapped)

#NanoB, filtrate, VSP:
NanoB_vs_DirExt_Filtrate_VSP_mapped <- glmQLFTest(fit_mapped, contrast = (is.NanoB & is.Filtrate & is.VSP) - (is.DirectExt & is.Filtrate & is.VSP))
top25(NanoB_vs_DirExt_Filtrate_VSP_mapped)

#NanoB, filtrate, RPIP:
NanoB_vs_DirExt_Filtrate_RPIP_mapped <- glmQLFTest(fit_mapped, contrast = (is.NanoB & is.Filtrate & is.RPIP) - (is.DirectExt & is.Filtrate & is.RPIP))
top25(NanoB_vs_DirExt_Filtrate_RPIP_mapped)

#NanoB, Retentate, untargeted:
NanoB_vs_DirExt_Retentate_Untargeted_mapped <- glmQLFTest(fit_mapped, contrast = (is.NanoB & is.Retentate & is.Untargeted) - (is.DirectExt & is.Retentate & is.Untargeted))
top25(NanoB_vs_DirExt_Retentate_Untargeted_mapped)

#NanoB, Retentate, VSP:
NanoB_vs_DirExt_Retentate_VSP_mapped <- glmQLFTest(fit_mapped, contrast = (is.NanoB & is.Retentate & is.VSP) - (is.DirectExt & is.Retentate & is.VSP))
top25(NanoB_vs_DirExt_Retentate_VSP_mapped)

#NanoB, Retentate, RPIP:
NanoB_vs_DirExt_Retentate_RPIP_mapped <- glmQLFTest(fit_mapped, contrast = (is.NanoB & is.Retentate & is.RPIP) - (is.DirectExt & is.Retentate & is.RPIP))
top25(NanoB_vs_DirExt_Retentate_RPIP_mapped)

##################
#Compare Nanotrap Types
##################

#NanoB vs NanoA main effect:
NanoB_vs_NanoA_main_effect_mapped <- glmQLFTest(fit_mapped, contrast = is.NanoB - is.NanoA)
top25(NanoB_vs_NanoA_main_effect_mapped)

#NanoB vs NanoA, unfiltered, untargeted:
NanoB_vs_NanoA_Unfiltered_Untargeted_mapped <- glmQLFTest(fit_mapped, contrast = (is.NanoB & is.Unfiltered & is.Untargeted) - (is.NanoA & is.Unfiltered & is.Untargeted))
top25(NanoB_vs_NanoA_Unfiltered_Untargeted_mapped)

#NanoB vs NanoA, unfiltered, VSP:
NanoB_vs_NanoA_Unfiltered_VSP_mapped <- glmQLFTest(fit_mapped, contrast = (is.NanoB & is.Unfiltered & is.VSP) - (is.NanoA & is.Unfiltered & is.VSP))
top25(NanoB_vs_NanoA_Unfiltered_VSP_mapped)

#NanoB vs NanoA, unfiltered, RPIP:
NanoB_vs_NanoA_Unfiltered_RPIP_mapped <- glmQLFTest(fit_mapped, contrast = (is.NanoB & is.Unfiltered & is.RPIP) - (is.NanoA & is.Unfiltered & is.RPIP))
top25(NanoB_vs_NanoA_Unfiltered_RPIP_mapped)

#NanoB vs NanoA, filtrate, untargeted:
NanoB_vs_NanoA_Filtrate_Untargeted_mapped <- glmQLFTest(fit_mapped, contrast = (is.NanoB & is.Filtrate & is.Untargeted) - (is.NanoA & is.Filtrate & is.Untargeted))
top25(NanoB_vs_NanoA_Filtrate_Untargeted_mapped)

#NanoB vs NanoA, filtrate, VSP:
NanoB_vs_NanoA_Filtrate_VSP_mapped <- glmQLFTest(fit_mapped, contrast = (is.NanoB & is.Filtrate & is.VSP) - (is.NanoA & is.Filtrate & is.VSP))
top25(NanoB_vs_NanoA_Filtrate_VSP_mapped)

#NanoB vs NanoA, filtrate, RPIP:
NanoB_vs_NanoA_Filtrate_RPIP_mapped <- glmQLFTest(fit_mapped, contrast = (is.NanoB & is.Filtrate & is.RPIP) - (is.NanoA & is.Filtrate & is.RPIP))
top25(NanoB_vs_NanoA_Filtrate_RPIP_mapped)

#NanoB vs NanoA, Retentate, untargeted:
NanoB_vs_NanoA_Retentate_Untargeted_mapped <- glmQLFTest(fit_mapped, contrast = (is.NanoB & is.Retentate & is.Untargeted) - (is.NanoA & is.Retentate & is.Untargeted))
top25(NanoB_vs_NanoA_Retentate_Untargeted_mapped)

#NanoB vs NanoA, Retentate, VSP:
NanoB_vs_NanoA_Retentate_VSP_mapped <- glmQLFTest(fit_mapped, contrast = (is.NanoB & is.Retentate & is.VSP) - (is.NanoA & is.Retentate & is.VSP))
top25(NanoB_vs_NanoA_Retentate_VSP_mapped)

#NanoB vs NanoA, Retentate, RPIP:
NanoB_vs_NanoA_Retentate_RPIP_mapped <- glmQLFTest(fit_mapped, contrast = (is.NanoB & is.Retentate & is.RPIP) - (is.NanoA & is.Retentate & is.RPIP))
top25(NanoB_vs_NanoA_Retentate_RPIP_mapped)

##################
###Effect of Filtration
##################
#Filtrate vs. unfiltered main effect:
Filtrate_vs_Unfiltered_main_effect_mapped <- glmQLFTest(fit_mapped, contrast = is.Filtrate - is.Unfiltered)
top25(Filtrate_vs_Unfiltered_main_effect_mapped)

#Filtrate vs. unfiltered, no Nanotrap, Untargeted
Filtrate_vs_Unfiltered_DirExt_Untargeted <- glmQLFTest(fit_mapped, contrast = (is.Filtrate & is.DirectExt & is.Untargeted) - (is.Unfiltered & is.DirectExt & is.Untargeted))
top25(Filtrate_vs_Unfiltered_DirExt_Untargeted)

#Filtrate vs. unfiltered, no Nanotrap, VSP
Filtrate_vs_Unfiltered_DirExt_VSP <- glmQLFTest(fit_mapped, contrast = (is.Filtrate & is.DirectExt & is.VSP) - (is.Unfiltered & is.DirectExt & is.VSP))
top25(Filtrate_vs_Unfiltered_DirExt_VSP)

#Filtrate vs. unfiltered, no Nanotrap, RPIP
Filtrate_vs_Unfiltered_DirExt_RPIP <- glmQLFTest(fit_mapped, contrast = (is.Filtrate & is.DirectExt & is.RPIP) - (is.Unfiltered & is.DirectExt & is.RPIP))
top25(Filtrate_vs_Unfiltered_DirExt_RPIP)

#Filtrate vs. unfiltered, NanoA, Untargeted
Filtrate_vs_Unfiltered_NanoA_Untargeted <- glmQLFTest(fit_mapped, contrast = (is.Filtrate & is.NanoA & is.Untargeted) - (is.Unfiltered & is.NanoA & is.Untargeted))
top25(Filtrate_vs_Unfiltered_NanoA_Untargeted)

#Filtrate vs. unfiltered, NanoA, VSP
Filtrate_vs_Unfiltered_NanoA_VSP <- glmQLFTest(fit_mapped, contrast = (is.Filtrate & is.NanoA & is.VSP) - (is.Unfiltered & is.NanoA & is.VSP))
top25(Filtrate_vs_Unfiltered_NanoA_VSP)

#Filtrate vs. unfiltered, NanoA, RPIP
Filtrate_vs_Unfiltered_NanoA_RPIP <- glmQLFTest(fit_mapped, contrast = (is.Filtrate & is.NanoA & is.RPIP) - (is.Unfiltered & is.NanoA & is.RPIP))
top25(Filtrate_vs_Unfiltered_NanoA_RPIP)

#Filtrate vs. unfiltered, NanoB, Untargeted
Filtrate_vs_Unfiltered_NanoB_Untargeted <- glmQLFTest(fit_mapped, contrast = (is.Filtrate & is.NanoB & is.Untargeted) - (is.Unfiltered & is.NanoB & is.Untargeted))
top25(Filtrate_vs_Unfiltered_NanoB_Untargeted)

#Filtrate vs. unfiltered, NanoB, VSP
Filtrate_vs_Unfiltered_NanoB_VSP <- glmQLFTest(fit_mapped, contrast = (is.Filtrate & is.NanoB & is.VSP) - (is.Unfiltered & is.NanoB & is.VSP))
top25(Filtrate_vs_Unfiltered_NanoB_VSP)

#Filtrate vs. unfiltered, NanoB, RPIP
Filtrate_vs_Unfiltered_NanoB_RPIP <- glmQLFTest(fit_mapped, contrast = (is.Filtrate & is.NanoB & is.RPIP) - (is.Unfiltered & is.NanoB & is.RPIP))
top25(Filtrate_vs_Unfiltered_NanoB_RPIP)

##################
#Retentate vs. unfiltered main effect:
Retentate_vs_Unfiltered_main_effect_mapped <- glmQLFTest(fit_mapped, contrast = is.Retentate - is.Unfiltered)
top25(Retentate_vs_Unfiltered_main_effect_mapped)

#Retentate vs. unfiltered, no Nanotrap, Untargeted
Retentate_vs_Unfiltered_DirExt_Untargeted <- glmQLFTest(fit_mapped, contrast = (is.Retentate & is.DirectExt & is.Untargeted) - (is.Unfiltered & is.DirectExt & is.Untargeted))
top25(Retentate_vs_Unfiltered_DirExt_Untargeted)

#Retentate vs. unfiltered, no Nanotrap, VSP
Retentate_vs_Unfiltered_DirExt_VSP <- glmQLFTest(fit_mapped, contrast = (is.Retentate & is.DirectExt & is.VSP) - (is.Unfiltered & is.DirectExt & is.VSP))
top25(Retentate_vs_Unfiltered_DirExt_VSP)

#Retentate vs. unfiltered, no Nanotrap, RPIP
Retentate_vs_Unfiltered_DirExt_RPIP <- glmQLFTest(fit_mapped, contrast = (is.Retentate & is.DirectExt & is.RPIP) - (is.Unfiltered & is.DirectExt & is.RPIP))
top25(Retentate_vs_Unfiltered_DirExt_RPIP)

#Retentate vs. unfiltered, NanoA, Untargeted
Retentate_vs_Unfiltered_NanoA_Untargeted <- glmQLFTest(fit_mapped, contrast = (is.Retentate & is.NanoA & is.Untargeted) - (is.Unfiltered & is.NanoA & is.Untargeted))
top25(Retentate_vs_Unfiltered_NanoA_Untargeted)

#Retentate vs. unfiltered, NanoA, VSP
Retentate_vs_Unfiltered_NanoA_VSP <- glmQLFTest(fit_mapped, contrast = (is.Retentate & is.NanoA & is.VSP) - (is.Unfiltered & is.NanoA & is.VSP))
top25(Retentate_vs_Unfiltered_NanoA_VSP)

#Retentate vs. unfiltered, NanoA, RPIP
Retentate_vs_Unfiltered_NanoA_RPIP <- glmQLFTest(fit_mapped, contrast = (is.Retentate & is.NanoA & is.RPIP) - (is.Unfiltered & is.NanoA & is.RPIP))
top25(Retentate_vs_Unfiltered_NanoA_RPIP)

#Retentate vs. unfiltered, NanoB, Untargeted
Retentate_vs_Unfiltered_NanoB_Untargeted <- glmQLFTest(fit_mapped, contrast = (is.Retentate & is.NanoB & is.Untargeted) - (is.Unfiltered & is.NanoB & is.Untargeted))
top25(Retentate_vs_Unfiltered_NanoB_Untargeted)

#Retentate vs. unfiltered, NanoB, VSP
Retentate_vs_Unfiltered_NanoB_VSP <- glmQLFTest(fit_mapped, contrast = (is.Retentate & is.NanoB & is.VSP) - (is.Unfiltered & is.NanoB & is.VSP))
top25(Retentate_vs_Unfiltered_NanoB_VSP)

#Retentate vs. unfiltered, NanoB, RPIP
Retentate_vs_Unfiltered_NanoB_RPIP <- glmQLFTest(fit_mapped, contrast = (is.Retentate & is.NanoB & is.RPIP) - (is.Unfiltered & is.NanoB & is.RPIP))
top25(Retentate_vs_Unfiltered_NanoB_RPIP)



