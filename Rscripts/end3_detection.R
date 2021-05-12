# >> 3´ end analysis << #

# load libraries ----
library(here)
source(here("Rscripts/load_libraries.R"))

# functions & defs ----
mod_termseq_peaks <- function(input, strand_s, cov_min = 3, merge_w = 20){
  suppressMessages(vroom(input, col_names = F, num_threads = 8)) %>%
    dplyr::rename(chr = X1, start_peak = X2, end_peak = X3, 
                  prominence = X5, strand_peak = X6, width = X10,
                  start_cov = X12, end_cov = X13, cov = X14, width_cov = X15) %>%
    dplyr::select(-X4, -X7, -X8, -X11) %>%
    group_by(start_peak, end_peak) %>%
    dplyr::filter(cov == max(cov)) %>%
    dplyr::mutate(decision_v = ifelse(strand_s == "+", 
                                      max(end_cov), min(end_cov))) %>%
    dplyr::filter(end_cov == decision_v) %>%
    ungroup() %>%
    arrange(end_cov) %>%
    mutate(index = lag(end_cov, default = 1) + as.integer(merge_w),
           index1 = cumsum(ifelse(index >= start_cov, 0, 1))+1) %>%
    dplyr::group_by(index1) %>%
    dplyr::filter(cov == max(cov),
                  cov >= cov_min) 
}


merge_peaks <- function(input_TTS, input_genome_gff){
  
  # prepare TTS data
  TTS_p <- input_TTS %>%
    ungroup() %>%
    dplyr::rename(seqname = chr, start = start_cov, end = end_cov, strand = strand_peak) %>%
    dplyr::select(seqname, start, end, strand, cov) %>%
    makeGRangesFromDataFrame()
  
  # prepare annotation file
  gff_p <- input_genome_gff %>%
    dplyr::filter(strand_feature == levels(as.factor(input_TTS$strand_peak))) %>%
    mutate(start_feature = ifelse(strand_feature == "-",start_feature - 300, start_feature),
           end_feature = ifelse(strand_feature == "+",end_feature + 300, end_feature)) %>%
    dplyr::rename(seqnames = seqid, start = start_feature, end = end_feature, names = id_name) %>%
    distinct(start, end, .keep_all = T) %>%
    makeGRangesFromDataFrame()
  
  # find overlapping peaks
  ol <- findOverlapsOfPeaks(TTS_p, gff_p)
  
  # select overlapping features
  overlaps <- ol$overlappingPeaks[["TTS_p///gff_p"]]
}

annotate_peaks <- function(merged_peaks, input_TTS, input_genome_gff = ecoli_gff, input_genome_fasta = ecoli_fasta){
  
  # set colnames
  colnames(merged_peaks) <- c(colnames(merged_peaks)[1:6],
                              paste0(colnames(merged_peaks)[7:12], "_2"),
                              colnames(merged_peaks)[13:14]) 
  
  final_peaks <- merged_peaks %>%
    as_tibble() %>%
    dplyr::rename(TTS = end) %>%
    dplyr::select(seqnames, TTS,start_2, end_2,strand) %>%
    dplyr::mutate(end_2 = ifelse(strand == "+", end_2 - 300, end_2),
                  start_2 = ifelse(strand == "-", start_2 + 300, start_2)) %>%
    dplyr::rename(start = start_2,end = end_2) %>%
    left_join(input_genome_gff, by = c("start" = "start_feature", "end" = "end_feature")) %>%
    rowwise() %>%
    dplyr::mutate(tts_sequence = ifelse(strand == "+", as.character(input_genome_fasta$`NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome`[(TTS - 20):(TTS + 10 )]),
                                        as.character(reverseComplement(input_genome_fasta$`NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome`[(TTS - 10):(TTS + 20 )]))),
                  UTR3 = ifelse(strand == "+", TTS - end, start - TTS)) %>%
    left_join(input_TTS %>% 
                ungroup() %>% 
                dplyr::rename(TTS = end_cov) %>% 
                dplyr::select(TTS, cov), by = "TTS") %>%
    ungroup() %>%
    group_by(id_name) %>%
    mutate(TTS_type = ifelse(cov == max(cov) & strand == "+" & TTS >= end - 5, "primary", 
                             ifelse(cov != max(cov) & strand == "+" & TTS >= end, "secondary", 
                                    ifelse(cov == max(cov) & strand == "-" & TTS <= start + 5, "primary", 
                                           ifelse(cov != max(cov) & strand == "-" & TTS <= start, "secondary", 
                                                  ifelse(TTS < end & TTS > start, "internal", "rest")))))) %>%
    group_by(TTS, TTS_type) %>%
    dplyr::slice(which.min(UTR3)) %>%
    ungroup()
  return(final_peaks)
}

peaks_pipe <- function(file_dir, names, map_method, p_f, m_f, gff = ecoli_gff){
  
  # > set frames
  tts_frame   <- data.table()
  partial_frame <- data.table()
  
  # > loop through files
  for (i in seq_along(names)){
    
    f <- paste0(dir, "/data/tts_data/", file_dir, "/")
    
    # modify termseq-peaks results
    tic("mod peaks")
    peaks_p     <- mod_termseq_peaks(paste0(f,p_f[i]), strand_s = "+")
    peaks_m     <- mod_termseq_peaks(paste0(f,m_f[i]), strand_s = "-")
    toc()
    
    # merge peaks with genome
    tic("merge peaks")
    peaks_tts_p     <- merge_peaks(peaks_p,gff)
    peaks_tts_m     <- merge_peaks(peaks_m,gff)
    toc()
    
    # add gene information
    tic("add info")
    peaks_tts_p_anno <- annotate_peaks(peaks_tts_p, peaks_p) 
    peaks_tts_m_anno <- annotate_peaks(peaks_tts_m, peaks_m) 
    toc()
    
    # bind plus and minus
    tic("bind frames")
    partial_frame <- rbind(peaks_tts_p_anno, peaks_tts_m_anno) %>%
      mutate(dataset = names[i]) %>%
      mutate(method = map_method)
    toc()
    
    # export
    tts_frame <- rbind(tts_frame, partial_frame) %>%
      as_tibble() 
    
  }
  return(tts_frame)
}

# load & tidy data ----

## genome data ====
dir <- here()
ecoli_gff   <- read_in_gff(paste0(dir, "/data/genome_data/NC_000913.3.gff3"))
ecoli_fasta <- readDNAStringSet(paste0(dir, "/data/genome_data/NC_000913.3.fasta"))

## Peak tables ====
### read in (narrowpeaks files if without replicates) of full-length > polyA.trimmed > polySSPtrimmed > soft/hard-clip removed ####
files <- list.files(paste0(dir, "/data/tts_data/tts_data_pychopper_auto_cutadapt_SSP_clipped/"), recursive = T, pattern = ".narrowPeak.counts")
plus_files     <- files[which(1:length(files) %% 2 == 0)] 
minus_files    <- files[which(1:length(files) %% 2 == 1)] 
dataset_names  <- str_split_fixed(str_split_fixed(str_split_fixed(plus_files, "\\/", n = 3)[,3],"_fu",2)[,1], "\\.",2)[,1]

tts_ssp <- peaks_pipe("tts_data_pychopper_auto_cutadapt_SSP_clipped",dataset_names,"clipped",plus_files, minus_files)
tts_data_pychopper_auto_cutadapt_SSP_clipped <- tts_ssp
fwrite(tts_data_pychopper_auto_cutadapt_SSP_clipped, paste0(dir, "/data/tts_data/tts_data_pychopper_auto_cutadapt_SSP_clipped.tsv"), col.names = T, sep = "\t")

## read in (narrowpeaks files if without replicates) of raw reads
files <- list.files(paste0(dir, "/data/tts_data/tts_data_notrimming/"), recursive = T, pattern = ".narrowPeak.counts")
plus_files     <- files[which(1:length(files) %% 2 == 0)] 
minus_files    <- files[which(1:length(files) %% 2 == 1)] 
dataset_names  <- str_split_fixed(str_split_fixed(str_split_fixed(plus_files, "\\/", n = 3)[,3],"_fu",2)[,1], "\\.",2)[,1]

### calc TTS ##############################
tts_raw <- peaks_pipe("tts_data_notrimming",dataset_names,"raw",plus_files, minus_files)
tts_data_notrimming <- tts_raw
fwrite(tts_data_notrimming, paste0(dir, "/data/tts_data/tts_data_notrimming.tsv"), col.names = T, sep = "\t")

# comparison to other datasets -------------------------------------------------------------------------------
## SMRT-cap seq ======================================================================================
smrt_tts <- read_xlsx(here("data/comparison_data/TTS/smrt_cap-41467_2018_5997_MOESM4_ESM.xlsx"),sheet = "Rich_TTS", skip = 1) %>%
  dplyr::mutate(TSS = ifelse(strand == "forward",start,end),
                TTS = ifelse(strand == "forward", end, start),
                gene = fully_covered_genes) %>%
  separate_rows(gene, sep = "\\|") %>%
  left_join(ecoli_gff, by = "gene") %>%
  mutate(UTR3 = ifelse(strand == "forward", TTS - end_feature, start_feature - TTS)) %>%
  group_by(TTS) %>%
  dplyr::slice(which.max(number_of_reads_TTS)) %>%
  ungroup() %>%
  group_by(TTS) %>%
  dplyr::slice(which.min(UTR3)) %>%
  dplyr::select(TTS, strand, gene, UTR3) 

## Elife-Termseq-results ======================================================================================
elife_tts <- read_xlsx(here("data/comparison_data/TTS/elife-62438-supp1-v1.xlsx"),sheet = "LB 0.4", skip = 1) %>%
  dplyr::rename(TTS = `3´ end position`) %>%
  separate_rows(classification, sep = ",") %>%
  dplyr::filter(classification != "") %>%
  dplyr::filter(!is.na(details), classification %in% c("primary", " primary")) %>%
  dplyr::mutate(ecogene = substring(str_split_fixed(str_split_fixed(details, "; gene_name", 2)[,1], "gene_id",2)[,2],3,9)) %>%
  left_join(ecoli_gff, by = "ecogene") %>%
  distinct(TTS, gene, .keep_all = T) %>%
  dplyr::filter(!is.na(gene)) %>%
  dplyr::mutate(UTR3 = ifelse(strand_feature == "+", TTS - (end_feature), start_feature - TTS)) %>%
  dplyr::select(TTS, strand, gene, UTR3)


levels(as.factor(tts_data_notrimming$dataset))
levels(as.factor(tts_total_comparison$method))
## merge

# STATS ---------------------------------------------------------------
## calc ===============================================================
summary_tts_all <-  rbind(tts_data_pychopper_auto_cutadapt_SSP_clipped,
                          tts_data_notrimming) %>%
  group_by(method, id_name, dataset) %>%
  dplyr::mutate(total_n = n(),
                group_sec = ifelse(total_n >1, "has_secondary", "single")) %>%
  ungroup() %>%
  dplyr::filter(type == "CDS") %>%
  dplyr::filter(dataset %in% bc_to_sample$sample[c(1,5:12,15:16)]) %>%
  group_by(dataset, TTS_type, method) %>%
  summarise(n = n()) %>%
  mutate(mode = substr(dataset, 8,10)) %>%
  dplyr::filter(map_method == "clipped") %>%
  remove_missing(vars = "TSS_type")

## reorder ===============================================================
summary_tts_all$dataset <- factor(summary_tts_all$dataset,
                                  levels = rev((bc_to_sample$sample[c(1,15,16,11,12,5,7,6,8:10)])))

summary_tts_all$TTS_type <- factor(summary_tts_all$TTS_type,
                                   levels = rev(c("primary", "secondary", "internal", "rest")))

# The palette with highlight:
cbf1_high <- c("#EFEAFF", "#ABC2DA","#F6B2FB","#648FFF")

### plot #############
pdf(here("figures/R_plots/210416_summary_TTS.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = summary_tts_all,
       aes(x = dataset, y = n, fill = TTS_type)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,4000)) +
  geom_bar(stat = "identity", color = "black") +
  facet_grid(cols = vars(method)) +
  theme_Publication_white() +
  coord_flip() +
  scale_fill_manual(values = cbf1_high)
dev.off()

# CORRELATION ANALYSIS ----------------------------------------------------------------------------
tts_total_comparison <- rbind(tts_data_pychopper_auto_cutadapt_SSP_clipped,
                              tts_data_notrimming) %>%
  dplyr::filter(TTS_type %in% c("primary", "secondary")) %>%
  group_by(method, id_name, dataset) %>%
  dplyr::mutate(total_n = n(),
                group_sec = ifelse(total_n >1, "has_secondary", "single")) %>%
  dplyr::filter(cov >= 5) %>%
  ungroup() %>%
  dplyr::filter(TTS_type == "primary", type == "CDS") %>%
  dplyr::select(gene, dataset, method, UTR3) %>%
  dplyr::filter(dataset %in% bc_to_sample$sample[c(1,5:12,15:16)]) %>%
  rbind(elife_tts %>% mutate(dataset = "Elife", 
                             method = "Elife") %>%
          dplyr::select(gene, dataset, method, UTR3)) %>%
  rbind(dar_tts %>% mutate(dataset = "illumina", 
                           method = "illumina") %>%
          dplyr::select(gene, dataset, method, UTR3)) %>%
  rbind(smrt_tts %>% ungroup () %>% mutate(dataset = "Smrt_CAP", 
                                           method = "Smrt_CAP") %>%
          dplyr::select(gene, dataset, method, UTR3)) %>%
  dplyr::filter(UTR3 >= 0 & !is.na(UTR3) & is.finite(UTR3), UTR3 <= 300) %>%
  mutate(method2 = paste0(dataset, "_", method))


tts_frame_wide <- tts_total_comparison %>%
  dplyr::select(-method, -dataset) %>%
  pivot_wider(names_from = method2, values_from = UTR3, values_fn = {max}) %>%
  left_join(ecoli_gff) %>%
  dplyr::filter(!is.na(gene)) 

# >  reorder
tts_total_comparison$dataset <- factor(tts_total_comparison$dataset,
                                       levels = (bc_to_sample$sample[c(1,15,16,11,12,5,7,6,8:10)]))

## REPLICATE COMPARISON =================================================================
### cDNA vs PCR-cDNA
levels(as.factor(tts_total_comparison$dataset))

# > make frame
tts_total_comparison_cDNA <- tts_total_comparison %>%
  dplyr::select(dataset,UTR3, gene) %>%
  dplyr::filter(dataset %in% c("201210_PCB109_Ecoli_NOTEX_replicate1", "210317_DCS109_Ecoli_NOTEX_replicate2")) %>%
  pivot_wider(names_from = dataset, values_from = UTR3, values_fn = {max}) %>%
  dplyr::rename(PCR_cDNA = 2, cDNA = 3) %>%
  remove_missing() 

# > how many pairwise comparisons
nrow(tts_total_comparison_cDNA)

# > calc_density
tts_total_comparison_cDNA$density <- get_density(tts_total_comparison_cDNA$PCR_cDNA,tts_total_comparison_cDNA$cDNA)

# > plot PCR-cDNA vs cDNA
pdf(here("figures/R_plots/210415_tts_PCR-cDNA_cDNA.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = tts_total_comparison_cDNA, 
       aes(x = PCR_cDNA, 
           y = cDNA,
           fill = density)) +
  geom_abline(linetype = "dashed", slope = 1) +
  geom_point(alpha = 0.8, size = 5, shape = 21, color = "black") +
  scale_fill_gradientn(colours = brewer.pal(name = "Blues", n = 9)) +
  theme_Publication_white() +
  stat_cor(method = "pearson", label.x = 1) +
  coord_equal() 
dev.off()  



# DISTANCE BETWEEN FOUND SITES -------------------------------------------------------
tts_total_comparison$dataset <- factor(tts_total_comparison$dataset,
                                       levels = (bc_to_sample$sample[c(1,15,16,11,12,5,7,6,8:10)]))

pdf(here("figures/R_plots/210415_UTR3_distances_dRNA_Elife.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = tts_total_comparison, aes(x = (UTR3.x - UTR3.y.y), fill = mode)) +
  facet_grid(cols = vars(method), rows = vars(dataset)) +
  geom_histogram(binwidth = 1, aes(y=..density..), color = "black") +
  scale_x_continuous(limits = c(-30,30), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
  scale_fill_manual(values = cbf1[c(2,5,3)]) +
  theme_Publication_white()
dev.off()

levels(as.factor(tts_total_comparison$dataset))
pdf(here("figures/R_plots/210415_UTR3_distances_dRNA_Elife.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = subset(tts_total_comparison,dataset %in% "201210_PCB109_Ecoli_NOTEX_replicate1"),
       aes(x = UTR3.x, y = UTR3.x.x)) +
  facet_grid(cols = vars(method), rows = vars(dataset)) +
  geom_abline(linetype = "dashed", slope = 1) +
  geom_point(alpha = 0.8, size = 5, shape = 21, color = "black") +
  scale_fill_gradientn(colours = brewer.pal(name = "Blues", n = 9)) +
  scale_x_continuous(limits = c(0,300)) +
  scale_y_continuous(limits = c(0,300)) +
  theme_Publication_white() +
  stat_cor(method = "pearson", label.x = 1) +
  coord_equal() 
dev.off()


## STATS =======================================================
### > WHAT´s the number/percentage of sites found at the 0 site
tts_total_comparison_ONT <- rbind(tts_data_pychopper_auto_cutadapt_SSP_clipped,
                                  tts_data_notrimming) %>%
  dplyr::filter(TTS_type %in% c("primary", "secondary")) %>%
  group_by(method, id_name, dataset) %>%
  dplyr::mutate(total_n = n(),
                group_sec = ifelse(total_n >1, "has_secondary", "single")) %>%
  ungroup() %>%
  dplyr::filter(TTS_type == "primary", type == "CDS") %>%
  dplyr::select(gene, dataset, method, UTR3) %>%
  dplyr::filter(dataset %in% bc_to_sample$sample[c(1,5:12,15:16)])

stats_tts_distance_only <- tts_total_comparison_ONT %>%
  left_join(term_tts, by = "gene") %>%
  mutate(distance = UTR3.x - UTR3.y) %>%
  remove_missing(vars = "distance") %>%
  remove_missing(vars = "distance") %>%
  mutate(mode = substr(dataset, 8,10)) %>%
  group_by(dataset, method) %>%
  mutate(distance_p = sum(distance == 0)/n()*100) %>%
  distinct(dataset,method, distance_p, mode) 

stats_tts_distance_only$dataset <- factor(stats_tts_distance_only$dataset,
                                          levels = rev(bc_to_sample$sample[c(1,15,16,11,12,5,7,6,8:10)]))

stats_tts_distance_only$method <- factor(stats_tts_distance_only$method,
                                         levels = (c("raw", "clipped")))


pdf(here("figures/R_plots/210416_UTR3_raw_clipped_ONT_TERM.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = stats_tts_distance_only, aes(x = distance_p, y = dataset, fill = mode)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge2(width = 2), alpha = 1) +
  theme_Publication_white() +
  #facet_grid(cols = vars(method)) +
  scale_x_continuous(limits = c(0,50), expand = c(0,0)) +
  scale_fill_manual(values = cbf1[c(2,5,3)]) 
dev.off()




# DISTANCE BETWEEN FOUND SITES -------------------------------------------------------
stats_tts_distance$dataset <- factor(stats_tts_distance$dataset,
                                     levels = (bc_to_sample$sample[c(1,15,16,11,12,5,7,6,8:10)]))
stats_tts_distance$method <- factor(stats_tts_distance$method,
                                    levels = (c("raw", "clipped")))

stats_tts_distance %>%
  dplyr::filter(distance > -17 & distance < +17) %>%
  group_by(dataset) %>%
  distinct(gene) %>%
  summarise(n = n())

pdf(here("figures/R_plots/210415_UTR3_distances_ONT_SMRT.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = stats_tts_distance, aes(x = distance, fill = mode)) +
  facet_grid(cols = vars(method), rows = vars(dataset)) +
  geom_histogram(binwidth = 1, aes(y=..density..), color = "black") +
  scale_x_continuous(limits = c(-17,17), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.75)) +
  scale_fill_manual(values = cbf1[c(2,5,3)]) +
  theme_Publication_white()+
  theme(panel.grid.major.y = element_blank()) +
  geom_vline(xintercept = -12, linetype = "dashed")
dev.off()


## Correlation matrix ================================================================
### PEARSON PAIRWISE #################
tts_frame_wide <- rbind(tts_data_pychopper_auto_cutadapt_SSP_clipped) %>%
  dplyr::filter(TTS_type %in% c("primary", "secondary")) %>%
  group_by(method, id_name, dataset) %>%
  dplyr::mutate(total_n = n(),
                group_sec = ifelse(total_n >1, "has_secondary", "single")) %>%
  dplyr::filter(cov >= 5) %>%
  ungroup() %>%
  dplyr::filter(TTS_type == "primary", type == "CDS") %>%
  dplyr::select(gene, dataset, method, UTR3) %>%
  dplyr::filter(dataset %in% bc_to_sample$sample[c(1,5:12,15:16)]) %>%
  rbind(elife_tts %>% mutate(dataset = "Elife", 
                             method = "Elife") %>%
          dplyr::select(gene, dataset, method, UTR3)) %>%
  rbind(dar_tts %>% mutate(dataset = "illumina", 
                           method = "illumina") %>%
          dplyr::select(gene, dataset, method, UTR3)) %>%
  rbind(smrt_tts %>% ungroup () %>% mutate(dataset = "Smrt_CAP", 
                                           method = "Smrt_CAP") %>%
          dplyr::select(gene, dataset, method, UTR3)) %>%
  dplyr::filter(UTR3 >= 0 & !is.na(UTR3) & is.finite(UTR3), UTR3 <= 300) %>%
  mutate(method2 = paste0(dataset, "_", method)) %>%
  dplyr::select(-method, -dataset) %>%
  pivot_wider(names_from = method2, values_from = UTR3, values_fn = {max}) %>%
  left_join(ecoli_gff) %>%
  dplyr::filter(!is.na(gene)) 

colnames(tts_frame_wide)
res    <- cor(tts_frame_wide[c(2,11,12,9,10,3,5,4,6,7,8,14,15)],  method = "pearson", use = "pairwise.complete.obs")
#res    <- cor(tts_frame_wide[c(2:15)],  method = "pearson", use = "pairwise.complete.obs")
res_gg <- melt(get_upper_tri(res))

pdf(here("figures/R_plots/210416_TTS_CORR_MAT.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = res_gg, aes(Var2, Var1, fill = value, color = value, size = value)) +
  geom_tile(color = "grey", size = 0.3, fill = "white") +
  geom_tile(color = "black", size = 0.3, aes(width = value, height = value)) +
  theme_void() +
  scale_fill_gradientn(colours = brewer.pal(name = "YlGnBu", n = 9),
                       limit = c(0.5,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  geom_text(aes(label=round(value, digits = 2)), color = "white", size = 4) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1)) +
  coord_fixed() +
  scale_y_discrete(limits=rev)
dev.off()

### NUMBERS PAIRWISE ################
colnames(tss_frame_wide)
res_counts <- pairwiseCount(tts_frame_wide[c(2,11,12,9,10,3,5,4,6,7,8,14,15)], diagonal = F)
res_counts_gg <- melt(get_lower_tri(res_counts))

pdf(here("figures/R_plots/210416_TTS_CORR_MAT_NUMBERS.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = res_counts_gg, aes(as.factor(Var2), as.factor(Var1), fill = value, color = value, size = value)) +
  geom_tile(color = "grey", size = 0.3, fill = "white") +
  geom_point(color = "black", shape = 21) +
  #geom_tile(color = "black", size = 0.3, aes(width = value/500, height = value/500)) +
  theme_void() +
  scale_fill_gradientn(colours = brewer.pal(name = "Blues", n = 9),
                       space = "Lab", 
                       name="Overlap\nCorrelation") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1)) +
  coord_fixed() +
  scale_y_discrete(limits=rev)
dev.off()


# plotting .............................................................
## logos ========================================================
tts_total_comparison_ONT <- tts_data_pychopper_auto_cutadapt_SSP_clipped %>%
  group_by(method, id_name, dataset) %>%
  dplyr::mutate(total_n = n(),
                group_sec = ifelse(total_n >1, "has_secondary", "single")) %>%
  ungroup() %>%
  dplyr::filter(TTS_type == "primary", type == "CDS") %>%
  #dplyr::select(gene, dataset, method, UTR5) %>%
  dplyr::filter(dataset %in% bc_to_sample$sample[c(1,5:12,15:16)]) %>%
  mutate(mode = substr(dataset, 8,10)) %>%
  dplyr::select(mode, dataset, UTR3) %>%
  ungroup() %>%
  rbind(dar_tts %>% 
          ungroup() %>%
          dplyr::mutate(mode = "ILL", dataset = "ILL") %>%
          dplyr::select(mode, dataset, UTR3)) %>%
  rbind(smrt_tts %>% 
          ungroup() %>%
          dplyr::mutate(mode = "SMRT-Cap", dataset = "SMRT-Cap") %>%
          dplyr::select(mode, dataset, UTR3)) %>%
  dplyr::filter(UTR3 >= 0 & !is.na(UTR3) & is.finite(UTR3), UTR3 <= 300) 

### reorder #############
tts_total_comparison_ONT$dataset <- factor(tts_total_comparison_ONT$dataset,
                                           levels = c("SMRT-Cap","ILL", rev((bc_to_sample$sample[c(1,15,16,11,12,5,7,6,8:10)]))))


## 3´UTR ========================================================
pdf(here("figures/R_plots/210416_UTR3_all.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = tts_total_comparison_ONT, aes(y = dataset, fill = mode)) +
  geom_density_ridges(stat = "binline",binwidth = 4,
                      aes(x = UTR3, height =..ndensity..), 
                      scale = 0.9, alpha = 1) +
  theme_Publication_white() +
  #scale_x_continuous(limits = c(0,300), expand = c(0,0)) +
  scale_fill_manual(values = cbf1[c(2,1,5,3,4)])
dev.off()

## SEQUENCE LOGO ===============================================
acgt_color_scale = make_col_scheme(chars=c('A', 'T', 'C', 'G'), 
                                   cols=c("#8591B3","#6CB8D4","#A1CEC2","#AD9D86"))

# TSS-40 :TSS+10
seqs <- tts_data_pychopper_auto_cutadapt_SSP_clipped %>%
  dplyr::filter(dataset %in% "201210_PCB109_Ecoli_NOTEX_replicate1", TTS_type %in% "primary", type == "CDS", cov >= 3) %>%
  distinct(gene, .keep_all = T) %>%
  dplyr::filter(UTR3 >= 0 & !is.na(UTR3) & is.finite(UTR3), UTR3 <= 300)  %>%
  rowwise() %>%
  dplyr::mutate(tts_sequence2 = ifelse(strand == "+", as.character(ecoli_fasta$`NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome`[(TTS-20):(TTS+10)]),
                                       as.character(reverseComplement(ecoli_fasta$`NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome`[(TTS-10):(TTS+20)])))) %>%
  dplyr::select(tts_sequence2) 

pdf(here("figures/R_plots/210416_geom_logo_TTS.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot() + 
  geom_logo(seqs, font = "helvetica_bold", col_scheme = acgt_color_scale, seq_type = "dna") + 
  theme_logo() +
  theme_Publication_white() +
  theme(panel.grid.major = element_line(colour = NA),
        axis.ticks.x = element_line(colour = NA), 
        axis.text.x = element_text(size = 0)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0,0.4), expand = c(0,0))
dev.off()





## 5´UTR comparison ========================================================

utr3_table <- rbind(read_xlsx(here("data/comparison_data/TTS/elife-62438-supp1-v1.xlsx"),sheet = "LB 0.4", skip = 1) %>%
                      dplyr::rename(TTS = `3´ end position`) %>%
                      separate_rows(classification, sep = ",") %>%
                      dplyr::filter(classification != "") %>%
                      dplyr::filter(!is.na(details), classification %in% c("primary", " primary")) %>%
                      dplyr::mutate(ecogene = substring(str_split_fixed(str_split_fixed(details, "; gene_name", 2)[,1], "gene_id",2)[,2],3,9)) %>%
                      left_join(ecoli_gff, by = "ecogene") %>%
                      distinct(TTS, gene, .keep_all = T) %>%
                      dplyr::filter(!is.na(gene)) %>%
                      dplyr::mutate(UTR3 = ifelse(strand_feature == "+", TTS - (end_feature), start_feature - TTS)) %>%
                      rowwise() %>%
                      dplyr::mutate(tts_sequence = ifelse(strand == "+", as.character(ecoli_fasta$`NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome`[(TTS - 20):(TTS + 10 )]),
                                                          as.character(reverseComplement(ecoli_fasta$`NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome`[(TTS - 10):(TTS + 20 )])))) %>%
                      dplyr::select(UTR3, tts_sequence) %>% mutate(method = "Elife"),
                    read_xlsx(here("data/comparison_data/TTS/smrt_cap-41467_2018_5997_MOESM4_ESM.xlsx"),sheet = "Rich_TTS", skip = 1) %>%
                      dplyr::mutate(TSS = ifelse(strand == "forward", start, end), 
                                    TTS = ifelse(strand == "forward", end, start),
                                    gene = fully_covered_genes) %>%
                      dplyr::select(TTS, strand, gene) %>%
                      left_join(ecoli_gff, by = "gene") %>%
                      dplyr::mutate(UTR3 = ifelse(strand == "forward", TTS - (end_feature), start_feature - TTS)) %>%
                      rowwise() %>%
                      dplyr::mutate(tts_sequence = ifelse(strand == "+", as.character(ecoli_fasta$`NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome`[(TTS - 20):(TTS + 10 )]),
                                                          as.character(reverseComplement(ecoli_fasta$`NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome`[(TTS - 10):(TTS + 20 )])))) %>%
                      dplyr::select(UTR3, tts_sequence) %>% mutate(method = "SMRT-CAP"),
                    tts_data_pychopper_auto_cutadapt_clipped %>%
                      dplyr::filter(TTS_type %in% c("primary"),type == "CDS",
                                    dataset == "201210_PCB109_Ecoli_NOTEX_replicate1",
                                    cov > 3) %>%
                      distinct(id_name, .keep_all = T) %>%
                      dplyr::select(UTR3,tts_sequence) %>%
                      mutate(method = "clipped")) %>%
  #tts_data_pychopper_auto %>%
  #  dplyr::filter(TTS_type %in% c("primary"),type == "CDS",
  #                dataset == "201210_PCB109_Ecoli_NOTEX_replicate1",
  #                cov > 3) %>%
  #  distinct(id_name, .keep_all = T) %>%
  #  dplyr::select(UTR3,tts_sequence) %>%
  #  mutate(method = "pychopper"),
  #tts_data_notrimming %>%
  #  dplyr::filter(TTS_type %in% c("primary"),type == "CDS",
  #                dataset == "201210_PCB109_Ecoli_NOTEX_replicate1",
  #                cov > 3) %>%
#  distinct(id_name, .keep_all = T) %>%
#  dplyr::select(UTR3,tts_sequence) %>%
#  mutate(method = "raw")) %>%
group_by(method) %>%
  dplyr::filter(!is.na(UTR3)) %>%
  mutate(dens = approxfun(density(UTR3))(UTR3))


## geom seq logos ===================================================
color_scale = make_col_scheme(chars=c('A', 'T', 'C', 'G'), 
                              cols=ibm_colors[c(2,4,1,5)])
utr3_table %>% group_by(method) %>% summarise(n = n())
levels(as.factor(utr3_table$method))[2]
motif_set <-   c(as.character(utr3_table$tts_sequence[utr3_table$method == levels(as.factor(utr3_table$method))[2]]))

pdf(here("figures/R_plots/210217_geom_logo_TTS_Termseq.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot() + 
  geom_logo(motif_set, font = "helvetica_bold", col_scheme = color_scale, seq_type = "dna") + 
  theme_logo() +
  theme_Publication_white() +
  theme(panel.grid.major = element_line(colour = NA),
        axis.ticks.x = element_line(colour = NA), 
        axis.text.x = element_text(size = 0)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0,0.5), expand = c(0,0))
dev.off()


## UTR3 comparison ==================================================
pdf(here("figures/R_plots/210217_UTR3_comparison.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = utr3_table %>% dplyr::filter(method %in% c("SMRT-CAP", "Elife", "clipped")), 
       aes(x = method, y = UTR3, fill = dens)) +
  geom_half_violin(draw_quantiles = c(0.25,0.5,0.75), side = "r", scale = "width") +
  geom_point(position = position_jitter(width = 0.1, height = 1), size = 2, alpha = 0.5, shape = 21) +
  coord_flip() +
  theme_Publication_white() +
  scale_y_continuous(limits = c(0,300)) +
  scale_fill_viridis_c(option = "cividis")
dev.off()


pdf(here("figures/R_plots/210222_UTR3_comparison.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = utr3_table %>% dplyr::filter(method %in% c("SMRT-CAP", "Elife", "clipped")), 
       aes(x = UTR3, fill = dens, y = method)) +
  stat_density_ridges(geom = "density_ridges_gradient",bandwidth = 5, na.rm = T,from = 0,to = 300,
                      scale = 0.75, aes(height =..ndensity..),
                      jittered_points = TRUE,quantile_lines = T,quantiles = c(0.5),
                      position = position_points_jitter(yoffset = -0.025,width = 0, height = 0),
                      point_shape = '|', point_size = 3, point_alpha = 1)  +
  theme_Publication_white() +
  scale_x_continuous(limits = c(0,300), expand = c(0,0))
dev.off()



