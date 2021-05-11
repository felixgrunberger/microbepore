# >> Mapped read analysis 2 << #

# load libraries ----
library(here)
source(here("Rscripts/load_libraries.R"))

# functions & defs ----

## sample_names ====
bc_to_sample <- data.table(seq_run = c("RNA001_Ecoli",
                                       rep("PCB109_PCR15_Ecoli",4),
                                       rep("PCB109_PCR12_Ecoli",2),
                                       rep("DCS109_Ecoli",2),
                                       "RNA002_Ecoli_run1",
                                       "RNA002_Ecoli_run2"),
                           barcode = c(rep("no_barcode",1),
                                       paste0(rep("barcode0",4),1:6),
                                       paste0(rep("barcode0",2),1:2),
                                       rep("no_barcode",2)),
                           sample  = c("RNA001_Ecoli_TEX_replicate1",
                                       "PCB109_PCR15_Ecoli_NOTEX_replicate4",
                                       "PCB109_PCR15_Ecoli_NOTEX_replicate5",
                                       "PCB109_PCR15_Ecoli_TEX_replicate4",
                                       "PCB109_PCR15_Ecoli_TEX_replicate5",
                                       "PCB109_PCR12_Ecoli_NOTEX_replicate4",
                                       "PCB109_PCR12_Ecoli_TEX_replicate4",
                                       "DCS109_Ecoli_NOTEX_replicate2",
                                       "DCS109_Ecoli_NOTEX_replicate3",
                                       "RNA002_Ecoli_NOTEX_replicate2",
                                       "RNA002_Ecoli_NOTEX_replicate3"))

## calculate number of mapped reads per gene and join with gff file ====
mutate_bam_files <- function(inputBAMtable){
  inputBAMtable %>%
    group_by(gene, sample) %>%
    dplyr::mutate(reads = n(),
                  bases = sum(aligned_reads)) %>%
    dplyr::select(gene, reads, bases, sample) %>%
    distinct(gene, sample,reads, bases,.keep_all = T) %>%
    arrange(desc(reads)) %>%
    left_join(gff_table, by = "gene") %>%
    ungroup()
}


## calculate total counts per feature ====
calc_bam_counts <- function(inputBAMtable, mode){
  inputBAMtable %>%
    group_by(sample) %>%
    mutate(total_reads = sum(reads),
           total_bases = sum(bases)) %>%
    dplyr::mutate(type = ifelse(type == "rRNA" & mode == "rRNA", locus_name, as.character(type)),
                  type = ifelse(type == "tRNA", "ncRNA", as.character(type))) %>%
    group_by(type, sample) %>%
    summarise(number_of_reads = sum(reads),
              number_of_bases = sum(bases),
              number_of_reads_p = number_of_reads/total_reads*100,
              number_of_bases_p = number_of_bases/total_bases*100) %>%
    distinct(type, number_of_reads, number_of_bases,sample, .keep_all = T) %>%
    ungroup() %>%
    dplyr::mutate(type = ifelse(is.na(type) == T, "unknown", as.character(type))) 
}

raw_reads_plotting <- function(mydf, myx, myy, myfill, mypalette){
  ggplot(data = mydf, aes(x = {{myx}}, y = {{myy}}, fill = {{myfill}})) +
    theme_Publication_white() +
    ylab("") +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(linetype = "dashed", color = "black")) +
    scale_fill_manual(values = mypalette)
}

# load & tidy data ----

## load saved mapped data table ====
dir <- here()
total_frame <- vroom(paste0(dir, "/data/mapped_data_no_trimming.tsv"), num_threads = 8)

total_frame_plot <- total_frame %>%
  dplyr::mutate(mode = substr(sample, 1,3),
                type = ifelse(type == "CDS", "mRNA",
                              ifelse(type == "rRNA", "rRNA","other_ncRNA"))) %>%
  dplyr::filter(!is.na(type))
  

## prepare quantification summary ====
# > mutate bam frames
total_frame_q <- mutate_bam_files(total_frame)

# > calc bam counts 
total_frame_c <- calc_bam_counts(total_frame_q, mode = "other") %>%
  dplyr::mutate(type = factor(type, levels=c("unknown","ncRNA", "rRNA","CDS"))) %>%
  mutate(mode = substr(sample, 1,3))

## combine with raw read stats ====
summary_frame_sample <- vroom(paste0(dir, "/data/summary_data_overview.tsv"),num_threads = 8)
total_frame_all <- left_join(total_frame_plot, summary_frame_sample %>%
                               dplyr::select(read_id, sequence_length_template, mean_qscore_template), 
                             by = c("minion_read_name" = "read_id") )

# calculate stats ----
summary_stats <- total_frame_c %>%
  dplyr::filter(type == "CDS") %>% 
  mutate(number_of_bases_p = round(number_of_bases_p, digits = 2),
         sequencing_depth = round(number_of_bases/length(ecoli_fasta$chr), digits = 2)) %>%
  arrange(factor(sample, levels = bc_to_sample$sample[c(1,10,11,8,9,2,4,3,5,6,7)]))

# PLOTS ----

## reorder levels ====
total_frame_c$sample <- factor(total_frame_c$sample,
                              levels = rev(bc_to_sample$sample[c(1,10,11,8,9,2,4,3,5,6,7)]))
total_frame_c$mode <- factor(total_frame_c$mode,
                            levels = c("RNA", "DCS", "PCB"))

total_frame_all$mode <- factor(total_frame_all$mode,
                               levels = c("RNA", "DCS", "PCB"))

total_frame_plot$sample <- factor(total_frame_plot$sample,
                               levels = rev(bc_to_sample$sample[c(1,10,11,8,9,2,4,3,5,6,7)]))

total_frame_plot$type <- factor(total_frame_plot$type,
                                levels = c("mRNA", "rRNA", "other_ncRNA"))


## color palette ====
cbf1_high <- c("#EFEAFF", "#ABC2DA","#F6B2FB","#648FFF")

## plotting ==== 

### Proportion of mapped reads to features - Supplementary Fig. 7A ####
raw_reads_plotting(total_frame_c, number_of_reads_p, sample, type, cbf1_high) +
  geom_bar(stat = "identity", color = "black") +
  scale_x_continuous(expand = c(0,0), limits = c(0,100)) +
  xlab("Proportion of mapped reads to features (%)") 

### Proportion of mapped reads to features - Supplementary Fig. 7A ####
raw_reads_plotting(total_frame_c, number_of_bases_p, sample, type, cbf1_high) +
  geom_bar(stat = "identity", color = "black") +
  scale_x_continuous(expand = c(0,0), limits = c(0,100)) +
  xlab("Proportion of mapped bases to features (%)") 

### Aligned read length distribution - Supplementary Fig. 8A ####
raw_reads_plotting(total_frame_plot, 
                   aligned_reads, sample, mode, cbf1[c(2,5,3)]) +
  facet_grid(cols = vars(type)) +
  geom_density_ridges(aes(height =..ndensity..), scale = 0.9, color = "black") +
  scale_x_continuous(limits = c(0,3500), expand = c(0,0)) +
  scale_y_discrete(expand = c(0.01,0.1)) +
  xlab("Aligned bases (nt)") 

### Read identity distribution - Supplementary Fig. 8B ####
raw_reads_plotting(total_frame_plot, 
                   identity, sample, mode, cbf1[c(2,5,3)]) +
  facet_grid(cols = vars(type)) +
  geom_density_ridges(aes(height =..ndensity..), scale = 0.9, color = "black") +
  scale_x_continuous(limits = c(50, 100), expand = c(0,0)) +
  scale_y_discrete(expand = c(0.01,0.1)) +
  xlab("Aligned read identity (%)") 

### Read length vs aligned bases - Supplementary Fig. 9A ####
raw_reads_plotting(total_frame_all %>% dplyr::filter(type == "mRNA") %>% group_by(mode) %>% sample_n(5000), 
                   sequence_length_template, aligned_reads, mode, cbf1[c(2,5,3)]) +
  geom_abline(linetype = "dashed", slope = 1) +
  geom_abline(linetype = "dashed", slope = .5) +
  facet_grid(rows = vars(mode)) +
  geom_point(aes(color = mode), alpha = 0.1, size = 0.5) +
  stat_density2d(aes(alpha=..level.., fill = mode),color = NA,
                 bins=10, geom="polygon") +
  geom_density2d(color = "black", contour_var = "ndensity", bins = 10, aes(alpha = ..level..)) +
  xlab("Read length (bases)") +
  ylab("Aligned length (nt)") +
  guides(alpha = F) +
  scale_x_continuous(limits = c(0,3500), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,3500), expand = c(0,0)) +
  scale_color_manual(values = cbf1[c(3,2,5)]) +
  theme_Publication_white()

### Read qscore vs identity - Supplementary Fig. 9B ####
raw_reads_plotting(total_frame_all %>% dplyr::filter(type == "mRNA") %>% group_by(mode) %>% sample_n(5000), 
                   mean_qscore_template, identity, mode, cbf1[c(2,5,3)]) +
  geom_abline(linetype = "dashed", slope = 1) +
  geom_abline(linetype = "dashed", slope = .5) +
  facet_grid(rows = vars(mode)) +
  geom_point(aes(color = mode), alpha = 0.1, size = 0.5) +
  stat_density2d(aes(alpha=..level.., fill = mode),color = NA,
                 bins=10, geom="polygon") +
  geom_density2d(color = "black", contour_var = "ndensity", bins = 10, aes(alpha = ..level..)) +
  xlab("Read length (bases)") +
  ylab("Aligned length (nt)") +
  guides(alpha = F) +
  scale_x_continuous(limits = c(0,20), expand = c(0,0)) +
  scale_y_continuous(limits = c(50,100), expand = c(0,0)) +
  scale_color_manual(values = cbf1[c(3,2,5)]) +
  theme_Publication_white()

### Qscore per group - Supplementary Fig. 9C ####
ggplot(data = total_frame_all %>% dplyr::filter(type == "mRNA"), 
       aes(x = mode, 
           y = mean_qscore_template, 
           factor = sequence_length_template <= 1.75*aligned_reads,
           alpha  = sequence_length_template <= 1.75*aligned_reads)) +
  geom_split_violin(trim = T, scale = "width", aes(fill = mode)) +
  coord_flip() +
  theme_Publication_white() +
  scale_y_continuous(limits = c(0,20), expand = c(0,0)) +
  scale_fill_manual(values = cbf1[c(5,2,3)]) 

### Identity per group - Supplementary Fig. 9D ####
ggplot(data = total_frame_all %>% dplyr::filter(type == "mRNA"), 
       aes(x = mode, 
           y = identity, 
           factor = sequence_length_template <= 1.75*aligned_reads,
           alpha  = sequence_length_template <= 1.75*aligned_reads)) +
  geom_split_violin(trim = T, scale = "width", aes(fill = mode)) +
  coord_flip() +
  theme_Publication_white() +
  scale_y_continuous(limits = c(50,100), expand = c(0,0)) +
  scale_fill_manual(values = cbf1[c(5,2,3)]) 

### Summary group n - Supplementary Fig. 9E ####
total_frame_all %>% 
  dplyr::filter(type == "mRNA") %>%
  mutate(strandswitch = ifelse(sequence_length_template <= 1.75*aligned_reads, "correct", "not_correct")) %>%
  group_by(mode, strandswitch) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = strandswitch, values_from = n) %>%
  mutate(c_p = correct/(correct+not_correct),
         n_p = not_correct/(correct+not_correct)) %>%
  dplyr::select(-correct, -not_correct) %>%
  pivot_longer(c_p:n_p, values_to = "perc") %>%
  ggplot(aes(y = mode, x = perc, fill = mode, alpha = name == "c_p")) +
    geom_col(color = "black") +
    scale_fill_manual(values = cbf1[c(5,2,3)]) +
    theme_Publication_white() +
    scale_x_continuous(limits = c(0,1), expand = c(0,0)) +
    ylab("") +
    xlab("Proportion of 1D/2D reads (%)") +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linetype = "dashed", color = "black")) 


