library(tidyverse)
library(cowplot)

#Data from read_fastp_reports.R
all_fastp_summaries_unmerged <- read_rds("input/modified/all_fastp_reports_dupdedup.rds")
summary(all_fastp_summaries_unmerged)

#############################################################
#####Show how many reads and bases were removed by FASTP#####
#############################################################
#Portion of reads removed during QC filtering and deduplication combined:
portion_reads_removed_box <- ggplot(all_fastp_summaries_unmerged, 
                                    aes(x = Enrichment, 
                                        y = portion_reads_removed_qcANDdedup)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  ylab("Portion reads removed")

portion_reads_removed_box

#Portion bases removed during QC filtering and deduplication combined:
portion_bases_removed_box <- ggplot(all_fastp_summaries_unmerged, 
                                    aes(x = Enrichment, 
                                        y = portion_bases_removed_qcANDdedup)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  ylab("Portion bases removed")

portion_bases_removed_box

plot_grid(portion_bases_removed_box,
          portion_reads_removed_box,
          align = "v",
          ncol = 1)


##########################################################
#####Generate summary statistics table for supplement#####
##########################################################
QC_summary_statistics <- all_fastp_summaries_unmerged %>%
  group_by(Enrichment) %>%
  summarise(mean_prefilter_readcount = mean(summary.before_filtering.total_reads),
            sd_prefilter_readcount = sd(summary.before_filtering.total_reads),
            min_prefilter_readcount = min(summary.before_filtering.total_reads),
            max_prefilter_readcount = max(summary.before_filtering.total_reads),
            mean_prefilter_basecount = mean(summary.before_filtering.total_bases),
            sd_prefilter_basecount = sd(summary.before_filtering.total_bases),
            min_prefilter_basecount = min(summary.before_filtering.total_bases),
            max_prefilter_basecount = max(summary.before_filtering.total_bases),
            mean_postfilter_readcount = mean(summary.after_filtering.total_reads),
            sd_postfilter_readcount = sd(summary.after_filtering.total_reads),
            min_postfilter_readcount = min(summary.after_filtering.total_reads),
            max_postfilter_readcount = max(summary.after_filtering.total_reads),
            mean_postfilter_basecount = mean(summary.after_filtering.total_bases),
            sd_postfilter_basecount = sd(summary.after_filtering.total_bases),
            min_postfilter_basecount = min(summary.after_filtering.total_bases),
            max_postfilter_basecount = max(summary.after_filtering.total_bases)
            )

View(QC_summary_statistics)

###############################
#####Read length summaries#####
###############################
QC_length_ins_summary <- all_fastp_summaries_unmerged %>%
  group_by(Enrichment) %>%
  summarize(mean_mean_pre_read1length = mean(summary.before_filtering.read1_mean_length),
            mean_mean_post_read1length = mean(summary.after_filtering.read1_mean_length),
            mean_mean_pre_read2length = mean(summary.before_filtering.read2_mean_length),
            mean_mean_post_read2length = mean(summary.after_filtering.read2_mean_length),
            mean_peak_post_ins_size = mean(insert_size_peak))


#########################################
#####Miscellaneous boxplot summaries#####
#########################################

#Portion of bases removed by QC before deduplication:
ggplot(all_fastp_summaries_unmerged, aes(x = Enrichment, 
                                         y = portion_bases_removed)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle("Bases removed by QC filters")


#Portion of reads removed due to to many N bases:
ggplot(all_fastp_summaries_unmerged, aes(x = Enrichment, 
                                         y = filtering_result.too_many_N_reads /
                                           summary.before_filtering.total_reads)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle("Reads removed: too many N")

#Reads removed due to low complexity:
ggplot(all_fastp_summaries_unmerged, 
       aes(x = Enrichment, 
       y = filtering_result.low_complexity_reads /
         summary.before_filtering.total_reads)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle("Reads removed: low complexity")

#Reads corrected by fastp:
ggplot(all_fastp_summaries_unmerged,
       aes(x = Enrichment, 
       y = filtering_result.corrected_reads /
         summary.before_filtering.total_reads)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle("Portion of reads corrected")


#Reads removed by length filter:
ggplot(all_fastp_summaries_unmerged,
       aes(x = Enrichment, 
       y = filtering_result.too_short_reads /
         summary.before_filtering.total_reads)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle("Reads removed: too short")

#Insert size peaks:
ggplot(all_fastp_summaries_unmerged,
       aes(x = Enrichment, 
           y = insert_size_peak)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle("Insert size peak")
