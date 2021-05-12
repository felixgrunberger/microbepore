# >> 5´ end analysis << #

# load libraries ----
library(here)
source(here("Rscripts/load_libraries.R"))

# functions & defs ----

mod_tss_peaks <- function(input, strand_s, cov_min = 3, merge_w = 20){
  
  suppressMessages(vroom(input, col_names = F, num_threads = 8)) %>%
    dplyr::rename(chr = X1, start_peak = X2, end_peak = X3, 
                  prominence = X5, strand_peak = X6, width = X10,
                  start_cov = X12, end_cov = X13, cov = X14, width_cov = X15) %>%
    dplyr::select(-X4, -X7, -X8, -X11) %>%
    group_by(start_peak, end_peak) %>%
    dplyr::filter(cov == max(cov)) %>%
    dplyr::mutate(decision_v = ifelse(strand_s == "+", 
                                      min(end_cov), max(end_cov))) %>%
    dplyr::filter(end_cov == decision_v) %>%
    ungroup() %>%
    arrange(end_cov) %>%
    mutate(index = lag(end_cov, default = 1) + as.integer(merge_w),
           index1 = cumsum(ifelse(index >= start_cov, 0, 1))+1) %>%
    dplyr::group_by(index1) %>%
    dplyr::filter(cov == max(cov),
                  cov >= cov_min) 
}


merge_tss_peaks <- function(input_TTS, input_genome_gff){
  
  # prepare TTS data
  TTS_p <- input_TTS %>%
    ungroup() %>%
    dplyr::rename(seqname = chr, start = start_cov, end = end_cov, strand = strand_peak) %>%
    dplyr::select(seqname, start, end, strand, cov) %>%
    makeGRangesFromDataFrame()
  
  # prepare annotation file
  gff_p <- input_genome_gff %>%
    dplyr::filter(strand_feature == levels(as.factor(input_TTS$strand_peak))) %>%
    mutate(start_feature = ifelse(strand_feature == "+",start_feature - 300, start_feature),
           end_feature = ifelse(strand_feature == "-",end_feature + 300, end_feature)) %>%
    dplyr::rename(seqnames = seqid, start = start_feature, end = end_feature, names = id_name) %>%
    distinct(start, end, .keep_all = T) %>%
    makeGRangesFromDataFrame()
  
  # find overlapping peaks
  ol <- findOverlapsOfPeaks(TTS_p, gff_p)
  
  # select overlapping features
  overlaps <- ol$overlappingPeaks[["TTS_p///gff_p"]]
}

annotate_tss_peaks <- function(merged_peaks, input_TTS, input_genome_gff = ecoli_gff, input_genome_fasta = ecoli_fasta, type_select){
  
  # set colnames
  colnames(merged_peaks) <- c(colnames(merged_peaks)[1:6],
                              paste0(colnames(merged_peaks)[7:12], "_2"),
                              colnames(merged_peaks)[13:14]) 
  
  final_peaks <- merged_peaks %>%
    as_tibble() %>%
    dplyr::rename(TSS = end) %>%
    dplyr::select(seqnames, TSS,start_2, end_2,strand) %>%
    dplyr::mutate(end_2 = ifelse(strand == "-", end_2 - 300, end_2),
                  start_2 = ifelse(strand == "+", start_2 + 300, start_2)) %>%
    dplyr::rename(start = start_2,end = end_2) %>%
    left_join(input_genome_gff, by = c("start" = "start_feature", "end" = "end_feature")) %>%
    rowwise() %>%
    dplyr::mutate(tss_sequence = ifelse(strand == "+", as.character(input_genome_fasta$`NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome`[(TSS - 40):(TSS + 10 )]),
                                        as.character(reverseComplement(input_genome_fasta$`NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome`[(TSS - 10):(TSS + 40 )]))),
                  UTR5 = ifelse(strand == "+", start - TSS, TSS - end)) %>%
    left_join(input_TTS %>% 
                ungroup() %>% 
                dplyr::rename(TSS = end_cov) %>% 
                dplyr::select(TSS, cov), by = "TSS") %>%
    ungroup() %>%
    mutate(type_set = type_select) %>%
    group_by(id_name) %>%
    mutate(TSS_type = ifelse(cov == max(cov) & strand == "+" & TSS <= start & TSS >= start-300 & type_set == "cDNA", "primary", 
                             ifelse(cov != max(cov) & strand == "+" & TSS <= start & TSS >= start-300 & type_set == "cDNA", "secondary", 
                                    ifelse(cov == max(cov) & strand == "-" & TSS >= end & TSS <= end+300 & type_set == "cDNA", "primary", 
                                           ifelse(cov != max(cov) & strand == "-" & TSS >= end & TSS <= end+300 & type_set == "cDNA", "secondary", 
                                                  ifelse(cov == max(cov) & strand == "+" & TSS <= start+12 & TSS >= start-300+12 & type_set == "RNA", "primary", 
                                                         ifelse(cov != max(cov) & strand == "+" & TSS <= start+12 & TSS >= start-300+12 & type_set == "RNA", "secondary", 
                                                                ifelse(cov == max(cov) & strand == "-" & TSS >= end-12 & TSS <= end+300-12 & type_set == "RNA", "primary", 
                                                                       ifelse(cov != max(cov) & strand == "-" & TSS >= end-12 & TSS <= end+300-12 & type_set == "RNA", "secondary", 
                                                                              ifelse(TSS > start & TSS < end & type_set == "cDNA", "internal",
                                                                                     ifelse(TSS > start+12 & TSS < end-12 & type_set == "RNA", "internal",
                                                                                            "rest"))))))))))) %>%
    
    group_by(TSS, TSS_type) %>%
    dplyr::slice(which.min(UTR5)) %>%
    ungroup()
  
  
  return(final_peaks)
}


tss_peaks_pipe <- function(file_dir, names, map_method, p_f, m_f, gff = ecoli_gff){
  
  # > set frames
  tss_frame   <- data.table()
  partial_frame <- data.table()
  
  # > loop through files
  for (i in seq_along(names)){
    
    f <- paste0(dir, "/data/tss_data/", file_dir, "/")
    
    type_f <- ifelse(str_detect(names[i], "RNA") == T, "RNA", "cDNA")
    
    # modify termseq-peaks results
    tic("mod peaks")
    peaks_p     <- mod_tss_peaks(paste0(f,p_f[i]), strand_s = "+")
    peaks_m     <- mod_tss_peaks(paste0(f,m_f[i]), strand_s = "-")
    toc()
    
    # merge peaks with genome
    tic("merge peaks")
    peaks_tts_p     <- merge_tss_peaks(peaks_p,gff)
    peaks_tts_m     <- merge_tss_peaks(peaks_m,gff)
    toc()
    
    # add gene information
    tic("add info")
    peaks_tts_p_anno <- annotate_tss_peaks(peaks_tts_p, peaks_p, type_select = type_f) 
    peaks_tts_m_anno <- annotate_tss_peaks(peaks_tts_m, peaks_m, type_select = type_f) 
    toc()
    
    # bind plus and minus
    tic("bind frames")
    partial_frame <- rbind(peaks_tts_p_anno, peaks_tts_m_anno) %>%
      mutate(dataset = names[i]) %>%
      mutate(method = map_method)
    toc()
    
    # export
    tss_frame <- rbind(tss_frame, partial_frame) %>%
      as_tibble() 
    
  }
  return(tss_frame)
}


# load & tidy data ----

## genome data ====
dir <- here()
ecoli_gff   <- read_in_gff(paste0(dir, "/data/genome_data/NC_000913.3.gff3"))
ecoli_fasta <- readDNAStringSet(paste0(dir, "/data/genome_data/NC_000913.3.fasta"))

## Peak tables ====

### read in pychopper auto trimmed data > cutadapt polyA > cutadapt SSP > clipped removed ####
files <- list.files(paste0(dir,"/data/tss_data/tss_data_pychopper_auto_cutadapt_SSP_clipped/"), recursive = T, pattern = ".narrowPeak.counts")
tss_data_trimmed <- tss_peaks_pipe("tss_data_pychopper_auto_cutadapt_SSP_clipped",
                                   str_split_fixed(str_split_fixed(str_split_fixed(files[which(1:length(files) %% 2 == 0)], "\\/", n = 3)[,3],".plus",2)[,1],"_fu",2)[,1],
                                   "trimmed",
                                   files[which(1:length(files) %% 2 == 0)],
                                   files[which(1:length(files) %% 2 == 1)])


fwrite(tss_data_trimmed, 
       paste0(dir, "/data/tss_data/tss_data_pychopper_auto_cutadapt_SSP_clipped.tsv"), col.names = T, sep = "\t")


### read in untrimmed - raw mapped data ####
files <- list.files(paste0(dir,"/data/tss_data/tss_data_notrimming/"), recursive = T, pattern = ".narrowPeak.counts")
tss_data_untrimmed <- tss_peaks_pipe("tss_data_notrimming",
                                     str_split_fixed(str_split_fixed(str_split_fixed(files[which(1:length(files) %% 2 == 0)], "\\/", n = 3)[,3],".plus",2)[,1],"_fu",2)[,1],
                                     "untrimmed",
                                     files[which(1:length(files) %% 2 == 0)],
                                     files[which(1:length(files) %% 2 == 1)])

fwrite(tss_data_untrimmed, paste0(dir, "/data/tss_data/tss_data_notrimming.tsv"), col.names = T, sep = "\t")

## TSS from other studies ====
### TEX-TSS-results ####
dir <- "/Volumes/EX_SSD/"
tex_tss <- read_xlsx(paste0(dir, "data/comparison_data/TSS/2014_zjb999093409sd1.xlsx"),sheet = "TSS Map MasterTable", skip = 2) %>%
  dplyr::filter(Condition == "LB_0.4", Primary == 1) %>%
  dplyr::rename(TSS = Pos, 
                gene = Locus_tag,
                strand = Strand) %>%
  dplyr::select(TSS, strand, gene, UTRlength, enrichmentFactor, `Sequence -50 nt upstream + TSS (51nt)`) %>%
  left_join(ecoli_gff, by = "gene") %>%
  dplyr::select(TSS, gene, start_feature, end_feature, strand, UTRlength, enrichmentFactor, `Sequence -50 nt upstream + TSS (51nt)`) %>%
  mutate(UTR5 = as.numeric(UTRlength)) %>%
  group_by(gene, strand) %>%
  dplyr::filter(UTR5 >= 0, UTR5 <= 300) %>%
  dplyr::slice(which.max(enrichmentFactor)) %>% 
  ungroup() %>%
  dplyr::select(-enrichmentFactor)


### SMRT-Cap-results ####
smrt_tss <- read_xlsx(paste0(dir,"data/comparison_data/TTS/smrt_cap-41467_2018_5997_MOESM4_ESM.xlsx"),sheet = "Rich_TTS", skip = 1) %>%
  dplyr::mutate(TSS = ifelse(strand == "forward", start, end), 
                TTS = ifelse(strand == "forward", end, start),
                gene = fully_covered_genes) %>%
  
  separate_rows(gene, sep = "\\|") %>%
  dplyr::select(TSS, strand, gene,number_of_reads_TSS) %>%
  left_join(ecoli_gff, by = "gene") %>%
  dplyr::select(TSS, gene, start_feature, end_feature, strand,number_of_reads_TSS) %>%
  dplyr::mutate(UTR5 = ifelse(strand == "forward", (start_feature - TSS), TSS - end_feature)) %>%
  group_by(gene, strand) %>%
  dplyr::filter(UTR5 >= 0, UTR5 <= 300) %>%
  dplyr::slice(which.max(number_of_reads_TSS)) %>% 
  ungroup() %>%
  dplyr::select(-number_of_reads_TSS)


## Comparison ====
summary_tss_ONT <-  rbindlist(list(tss_data_trimmed,
                          tss_data_untrimmed)) %>%
  dplyr::filter(type == "CDS") %>%
  group_by(sample, TSS_type, method) %>%
  summarise(n = n()) %>%
  mutate(mode = substr(sample, 1,3)) %>%
  remove_missing(vars = "TSS_type")

# PLOTS ----

## reorder levels ====
summary_tss_ONT$sample <- factor(summary_tss_ONT$sample,
                               levels = rev(bc_to_sample$sample[c(1,10,11,8,9,2,4,3,5,6,7)]))
summary_tss_ONT$mode <- factor(summary_tss_ONT$mode,
                               levels = c("RNA", "DCS", "PCB"))
summary_tss_ONT$TSS_type <- factor(summary_tss_ONT$TSS_type,
                                   levels = rev(c("primary", "secondary", "internal", "rest")))
summary_tss_ONT$method <- factor(summary_tss_ONT$method,
                                   levels = c("untrimmed", "trimmed"))

## plotting ==== 

### Number of 5´ends in category - Supplementary Fig. 14A #### 
raw_reads_plotting(summary_tss_ONT, n, sample, TSS_type, cbf1_high) +
  geom_bar(stat = "identity", color = "black", position = position_stack()) +
  scale_x_continuous(expand = c(0,0), limits = c(0,6000)) +
  facet_grid(cols = vars(method)) +
  xlab("Number of 5´ends in category") 



## STATS =======================================================
### > WHAT´s the number/percentage of sites found at the 0 site
tss_total_comparison_ONT <- rbind(tss_data_pychopper_auto_cutadapt_SSP_clipped,
                                  tss_data_notrimming) %>%
  dplyr::filter(TSS_type %in% c("primary", "secondary")) %>%
  group_by(method, id_name, dataset) %>%
  dplyr::mutate(total_n = n(),
                group_sec = ifelse(total_n >1, "has_secondary", "single")) %>%
  ungroup() %>%
  dplyr::filter(TSS_type == "primary", type == "CDS") %>%
  dplyr::select(gene, dataset, method, UTR5) %>%
  dplyr::filter(dataset %in% bc_to_sample$sample[c(1,5:12,15:16)])

stats_tss_distance <- tss_total_comparison_ONT %>%
  left_join(tex_tss, by = "gene") %>%
  mutate(distance = UTR5.x - UTR5.y) %>%
  remove_missing(vars = "distance") %>%
  mutate(mode = substr(dataset, 8,10)) %>%
  group_by(dataset, method) %>%
  mutate(distance_p = ifelse(mode != "RNA", sum(distance == 0)/n()*100, sum(distance == -12)/n()*100)) %>%
  distinct(dataset,method, distance_p, mode) 

ggplot(data = stats_tss_distance, aes(x = distance_p, y = dataset, fill = mode, factor = method)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge2(width = 2), alpha = 1) +
  facet_grid(cols = vars(method)) +
  theme_Publication_white() +
  scale_fill_manual(values = cbf1[c(2,5,3)]) 



## REPLICATE comparison ====================================================================
tss_total_comparison <- rbind(tss_data_pychopper_auto_cutadapt_SSP_clipped,
                              tss_data_notrimming) %>%
  dplyr::filter(cov >= 5) %>%
  group_by(method, id_name, dataset) %>%
  dplyr::mutate(total_n = n(),
                group_sec = ifelse(total_n >1, "has_secondary", "single")) %>%
  ungroup() %>%
  dplyr::filter(TSS_type == "primary", type == "CDS", method == "raw") %>%
  dplyr::select(gene, dataset, method, UTR5) %>%
  rbind(tex_tss %>% mutate(dataset = "illumina", 
                           method = "illumina") %>%
          dplyr::select(gene, dataset, method, UTR5)) %>%
  rbind(smrt_tss %>% ungroup () %>% mutate(dataset = "Smrt_CAP", 
                                           method = "Smrt_CAP") %>%
          dplyr::select(gene, dataset, method, UTR5)) %>%
  dplyr::filter(UTR5 >= 0 & !is.na(UTR5) & is.finite(UTR5), UTR5 <= 300) %>%
  mutate(method2 = paste0(dataset, "_", method))


## STATS =======================================================
### > WHAT´s the number/percentage of sites found at the 0 site
tss_total_comparison_ONT <- rbind(tss_data_pychopper_auto_cutadapt_SSP_clipped %>% mutate(method = "pychopper_auto_cutadapt_clipped"),
                                  tss_data_notrimming %>% mutate(method = "notrimming")) %>%
  group_by(method, id_name, dataset) %>%
  dplyr::mutate(total_n = n(),
                group_sec = ifelse(total_n >1, "has_secondary", "single")) %>%
  ungroup() %>%
  dplyr::filter(TSS_type == "primary", type == "CDS") %>%
  dplyr::select(gene, dataset, method, UTR5) %>%
  dplyr::filter(dataset %in% bc_to_sample$sample[c(1,5:12,15:16)])

stats_tss_distance <- tss_total_comparison_ONT %>%
  dplyr::filter(method == "raw") %>%
  left_join(tex_tss, by = "gene") %>%
  mutate(distance = UTR5.x - UTR5.y) %>%
  remove_missing(vars = "distance") %>%
  mutate(mode = substr(dataset, 8,10)) 

# DISTANCE BETWEEN FOUND SITES -------------------------------------------------------
stats_tss_distance$dataset <- factor(stats_tss_distance$dataset,
                                     levels = (bc_to_sample$sample[c(1,15,16,11,12,5,7,6,8:10)]))
stats_tss_distance$method <- factor(stats_tss_distance$method,
                                    levels = (c("raw", "clipped")))

stats_tss_distance %>%
  dplyr::filter(distance > -17 & distance < +17) %>%
  group_by(dataset) %>%
  distinct(gene) %>%
  summarise(n = n())

## plot distance histogram =====================================
pdf(here("figures/R_plots/210415_UTR5_distances_ONT_SMRT.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = stats_tss_distance, aes(x = distance, fill = mode)) +
  facet_grid(cols = vars(method), rows = vars(dataset)) +
  geom_histogram(binwidth = 1, aes(y=..density..), color = "black") +
  scale_x_continuous(limits = c(-17,17), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.75)) +
  scale_fill_manual(values = cbf1[c(2,5,3)]) +
  theme_Publication_white()+
  theme(panel.grid.major.y = element_blank()) +
  geom_vline(xintercept = -12, linetype = "dashed")
dev.off()

## calc fractions =====================================
stats_tss_distance_only <- tss_total_comparison_ONT %>%
  left_join(tex_tss, by = "gene") %>%
  mutate(distance = UTR5.x - UTR5.y) %>%
  remove_missing(vars = "distance") %>%
  mutate(mode = substr(dataset, 8,10)) %>%
  group_by(dataset, method) %>%
  mutate(distance_p = ifelse(mode != "RNA", sum(distance == 0)/n()*100, sum(distance == -12)/n()*100)) %>%
  distinct(dataset,method, distance_p, mode) 

stats_tss_distance_only$dataset <- factor(stats_tss_distance_only$dataset,
                                          levels = rev(bc_to_sample$sample[c(1,15,16,11,12,5,7,6,8:10)]))

stats_tss_distance_only$method <- factor(stats_tss_distance_only$method,
                                         levels = c("raw", "clipped"))


pdf(here("figures/R_plots/210415_UTR5_raw_clipped_ONT_TEX.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = stats_tss_distance_only, aes(x = distance_p, y = dataset, fill = mode, factor = method)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge2(width = 2), alpha = 1) +
  #facet_grid(cols = vars(method)) +
  theme_Publication_white() +
  scale_x_continuous(limits = c(0,50), expand = c(0,0)) +
  scale_fill_manual(values = cbf1[c(2,5,3)]) 
dev.off()

tss_frame_wide <- tss_total_comparison %>%
  dplyr::select(-method, -dataset) %>%
  pivot_wider(names_from = method2, values_from = UTR5, values_fn = {max}) %>%
  left_join(ecoli_gff) %>%
  dplyr::filter(!is.na(gene)) 

# >  reorder
tss_total_comparison$dataset <- factor(tss_total_comparison$dataset,
                                       levels = (bc_to_sample$sample[c(1,15,16,11,12,5,7,6,8:10)]))

## REPLICATE COMPARISON =================================================================
levels(as.factor(tss_total_comparison$dataset))
### PCR-cDNA vs SMRT-CAP
# > detect secondary start sites that are primary in SMRT CAP
tss_total_comparison_SMRT_sec <- tss_data_notrimming %>%
  dplyr::filter(dataset %in% "201210_PCB109_Ecoli_NOTEX_replicate1", 
                method == "raw", 
                TSS_type %in% c("primary", "secondary")) %>%
  dplyr::filter(cov >= 5) %>%
  group_by(id_name) %>%
  dplyr::mutate(total_n = n(),
                group_sec = ifelse(total_n >1, "has_secondary", "single")) %>%
  ungroup() %>%
  dplyr::filter(group_sec == "has_secondary", TSS_type == "secondary", type == "CDS") %>%
  dplyr::select(gene, dataset, method, UTR5, group_sec) %>%
  left_join(smrt_tss, by = "gene") %>%
  dplyr::filter(abs(UTR5.x-UTR5.y) < 5) %>%
  dplyr::select(gene, UTR5.x)

# > label those
tss_total_comparison_SMRT <- tss_data_notrimming %>%
  dplyr::filter(dataset %in% "201210_PCB109_Ecoli_NOTEX_replicate1", 
                method == "raw", 
                TSS_type %in% c("primary", "secondary")) %>%
  dplyr::filter(cov >= 5) %>%
  group_by(id_name) %>%
  dplyr::mutate(total_n = n(),
                group_sec = ifelse(total_n >1, "has_secondary", "single")) %>%
  ungroup() %>%
  dplyr::filter(TSS_type == "primary", type == "CDS") %>%
  dplyr::mutate(group = ifelse(group_sec == "single", "single",
                               ifelse(group_sec == "has_secondary" & gene %in% tss_total_comparison_SMRT_sec$gene, "SMRT_first", "rest"))) %>%
  left_join(tss_total_comparison_SMRT_sec, by = "gene") %>%
  dplyr::mutate(UTR5 = ifelse(group == "SMRT_first", UTR5.x, UTR5)) %>%
  dplyr::select(-UTR5.x) %>%
  dplyr::select(gene, dataset, method, UTR5, group) %>%
  left_join(smrt_tss, by = "gene") %>%
  remove_missing(name = "UTR5.y")

tss_total_comparison_SMRT %>%
  group_by(group) %>%
  summarise(n = n())

# > how many pairwise comparisons
nrow(tss_total_comparison_SMRT)

# > calc_density
tss_total_comparison_SMRT$density <- get_density(tss_total_comparison_SMRT$PCR_cDNA,tss_total_comparison_SMRT$SMRT)

# > plot PCR-cDNA vs SMRT
pdf(here("figures/R_plots/210415_TSS_PCR-cDNA_SMRT.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = tss_total_comparison_SMRT, 
       aes(x = UTR5.x, 
           y = UTR5.y,
           fill = group)) +
  geom_abline(linetype = "dashed", slope = 1) +
  geom_point(alpha = 0.8, size = 5, shape = 21, color = "black") +
  #scale_fill_gradientn(colours = brewer.pal(name = "Blues", n = 9)) +
  theme_Publication_white() +
  stat_cor(method = "pearson", label.x = 1) +
  scale_fill_manual(values = cbf1_high[c(3,2,4)]) +
  coord_equal() 
dev.off()  


### cDNA vs SMRT
# > make frame
tss_total_comparison_SMRT1 <- tss_total_comparison %>%
  dplyr::filter(method %in% c("raw", "Smrt_CAP")) %>%
  dplyr::select(dataset,UTR5, gene) %>%
  dplyr::filter(dataset %in% c("201210_PCB109_Ecoli_NOTEX_replicate1", "Smrt_CAP")) %>%
  pivot_wider(names_from = dataset, values_from = UTR5, values_fn = {max}) %>%
  dplyr::rename(PCR_cDNA = 2, SMRT = 3) %>%
  remove_missing() 

# > how many pairwise comparisons
nrow(tss_total_comparison_SMRT1)

# > calc_density
tss_total_comparison_SMRT1$density <- get_density(tss_total_comparison_SMRT1$PCR_cDNA,tss_total_comparison_SMRT1$SMRT)

pdf(here("figures/R_plots/210415_TSS_PCR-cDNA_SMRT2.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = tss_total_comparison_SMRT1, 
       aes(x = PCR_cDNA, 
           y = SMRT,
           fill = density)) +
  geom_abline(linetype = "dashed", slope = 1) +
  geom_point(alpha = 0.8, size = 5, shape = 21, color = "black") +
  scale_fill_gradientn(colours = brewer.pal(name = "Blues", n = 9)) +
  theme_Publication_white() +
  stat_cor(method = "pearson", label.x = 1) +
  coord_equal() 
dev.off()  

### cDNA vs PCR-cDNA
# > make frame
tss_total_comparison_cDNA <- tss_total_comparison %>%
  dplyr::filter(method == "raw") %>%
  dplyr::select(dataset,UTR5, gene) %>%
  dplyr::filter(dataset %in% c("201210_PCB109_Ecoli_NOTEX_replicate1", "210317_DCS109_Ecoli_NOTEX_replicate2")) %>%
  pivot_wider(names_from = dataset, values_from = UTR5, values_fn = {max}) %>%
  dplyr::rename(PCR_cDNA = 2, cDNA = 3) %>%
  remove_missing() 

# > how many pairwise comparisons
nrow(tss_total_comparison_cDNA)

# > calc_density
tss_total_comparison_cDNA$density <- get_density(tss_total_comparison_cDNA$PCR_cDNA,tss_total_comparison_cDNA$cDNA)

# > plot PCR-cDNA vs cDNA
pdf(here("figures/R_plots/210415_TSS_PCR-cDNA_cDNA.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = tss_total_comparison_cDNA, 
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

# > plot PCR-cDNA vs cDNA
pdf(here("figures/R_plots/210415_TSS_PCR-cDNA_cDNA.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = tss_total_comparison_cDNA, 
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

### RNA vs cDNA
# > make frame
tss_total_comparison_RNA <- tss_total_comparison %>%
  group_by(dataset,gene) %>%
  #dplyr::slice(which.min(UTR5)) %>%
  dplyr::select(dataset,UTR5, gene) %>%
  dplyr::filter(dataset %in% c("190123_RNA001_Ecoli_TEX_replicate1", "201210_PCB109_Ecoli_NOTEX_replicate1")) %>%
  pivot_wider(names_from = dataset, values_from = UTR5, values_fn = {max}) %>%
  dplyr::rename(RNA = 2, PCR_cDNA = 3) %>%
  remove_missing() 

# > how many pairwise comparisons
nrow(tss_total_comparison_RNA)

# > calc_density
tss_total_comparison_RNA$density <- get_density(tss_total_comparison_RNA$RNA,tss_total_comparison_RNA$PCR_cDNA)

# > plot PCR-cDNA vs cDNA
pdf(here("figures/R_plots/210415_TSS_PCR-cDNA_RNA.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = tss_total_comparison_RNA, 
       aes(x = PCR_cDNA, 
           y = RNA,
           fill = density)) +
  geom_abline(linetype = "dashed", slope = 1) +
  geom_point(alpha = 0.8, size = 5, shape = 21, color = "black") +
  scale_fill_gradientn(colours = brewer.pal(name = "Blues", n = 9)) +
  theme_Publication_white() +
  stat_cor(method = "pearson", label.x = 1) +
  coord_equal() 
dev.off()  

## Correlation matrix: 
### PEARSON PAIRWISE #################
colnames(tss_frame_wide)
res    <- cor(tss_frame_wide[c(2,16,17,12,13,6,8,7,9,10,11)],  method = "pearson", use = "pairwise.complete.obs")
res_gg <- melt(get_upper_tri(res))

pdf(here("figures/R_plots/210415_TSS_CORR_MAT.pdf"), 
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
res_counts <- pairwiseCount(tss_frame_wide[c(2,16,17,12,13,6,8,7,9,10,11)], diagonal = F)
res_counts_gg <- melt(get_lower_tri(res_counts))

pdf(here("figures/R_plots/210415_TSS_CORR_MAT_NUMBERS.pdf"), 
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
tss_total_comparison_ONT <- rbind(tss_data_pychopper_auto_cutadapt_SSP_clipped %>% mutate(method = "pychopper_auto_cutadapt_clipped"),
                                  tss_data_notrimming %>% mutate(method = "notrimming")) %>%
  group_by(method, id_name, dataset) %>%
  dplyr::mutate(total_n = n(),
                group_sec = ifelse(total_n >1, "has_secondary", "single")) %>%
  ungroup() %>%
  dplyr::filter(TSS_type == "primary", type == "CDS") %>%
  #dplyr::select(gene, dataset, method, UTR5) %>%
  dplyr::filter(dataset %in% bc_to_sample$sample[c(1,5:12,15:16)]) %>%
  mutate(mode = substr(dataset, 8,10)) %>%
  mutate(UTR5 = ifelse(mode == "RNA", UTR5+12, UTR5)) %>%
  dplyr::select(mode, dataset, UTR5) %>%
  rbind(tex_tss %>% 
          dplyr::mutate(mode = "ILL", dataset = "ILL") %>%
          dplyr::select(mode, dataset, UTR5)) %>%
  rbind(smrt_tss %>% 
          dplyr::mutate(mode = "SMRT-Cap", dataset = "SMRT-Cap") %>%
          dplyr::select(mode, dataset, UTR5))

### reorder #############
tss_total_comparison_ONT$dataset <- factor(tss_total_comparison_ONT$dataset,
                                           levels = c("SMRT-Cap","ILL", rev((bc_to_sample$sample[c(1,15,16,11,12,5,7,6,8:10)]))))


## 5´UTR ========================================================
pdf(here("figures/R_plots/210416_UTR5_all.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = tss_total_comparison_ONT, aes(y = dataset, fill = mode)) +
  geom_density_ridges(stat = "binline",binwidth = 4,
                      aes(x = UTR5, height =..ndensity..), 
                      scale = 0.9, alpha = 1) +
  theme_Publication_white() +
  #scale_x_continuous(limits = c(0,300), expand = c(0,0)) +
  scale_fill_manual(values = cbf1[c(2,1,5,3,4)])
dev.off()

## SEQUENCE LOGO ===============================================
acgt_color_scale = make_col_scheme(chars=c('A', 'T', 'C', 'G'), 
                                   cols=c("#8591B3","#6CB8D4","#A1CEC2","#AD9D86"))

# TSS-40 :TSS+10
TSS=50
start=10
end=60

seqs <- tss_data_notrimming %>%
  dplyr::filter(dataset %in% "201210_PCB109_Ecoli_NOTEX_replicate1", TSS_type %in% "primary", type == "CDS", cov >= 3) %>%
  distinct(gene, .keep_all = T) %>%
  rowwise() %>%
  dplyr::mutate(tss_sequence2 = ifelse(strand == "+", as.character(ecoli_fasta$`NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome`[(TSS-40):(TSS)]),
                                       as.character(reverseComplement(ecoli_fasta$`NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome`[(TSS):(TSS+40)])))) %>%
  dplyr::select(tss_sequence2) 

pdf(here("figures/R_plots/210416_geom_logo_TSS.pdf"), 
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




## more important to see effect of sequencing depth on TSS detection ==========================================
tss_data_pychopper_auto_cutadapt_clipped %>%
  dplyr::filter(TSS_type %in% c("primary"),type == "CDS",
                dataset == "201210_PCB109_Ecoli_NOTEX_replicate1") %>%
  mutate(map_method = "clipped")


g <- tss_data_pychopper_auto_cutadapt_clipped %>% 
  dplyr::filter(TSS_type %in% c("primary","secondary")) %>%
  group_by(dataset, id_name) %>%
  dplyr::mutate(total_n = n(),
                group_sec = ifelse(total_n >1, "has_secondary", "single")) %>%
  #dplyr::mutate(UTR5_new = ifelse(group_sec == "single" & TSS_type == "primary", UTR5, 
  #                                ifelse(group_sec != "single", UTR5[TSS_type == "secondary"], NA))) %>%
  ungroup() %>%
  dplyr::filter(TSS_type == "primary") %>%
  rowwise() %>%
  mutate(polyG_detected = str_count(string = substring(tss_sequence,31,41),pattern = "GGG+") >= 1) %>%
  dplyr::filter(polyG_detected == F) %>%
  dplyr::select(UTR5, dataset, gene) %>%
  pivot_wider(names_from = dataset, values_from = UTR5, values_fn = {sum}) %>%
  left_join(subset(tss_data_pychopper_auto_cutadapt_clipped, dataset %in% "201208_PCB109_Ecoli_NOTEX_replicate1") %>%
              group_by(id_name) %>%
              dplyr::mutate(total_n = n(),
                            group_sec = ifelse(total_n >1, "has_secondary", "single")), by = "gene") %>%
  distinct(gene, .keep_all = T)
g
colnames(g)
ggplot(data = g, aes(x = `201208_PCB109_Ecoli_NOTEX_replicate1` - `201208_PCB109_Ecoli_NOTEX_replicate2`, fill = cov > 10)) +
  geom_histogram() +
  scale_x_continuous(limits = c(-30,30))

cg <- g %>% dplyr::filter(!is.na(`201208_PCB109_Ecoli_NOTEX_replicate1`),
                          !is.na(`201208_PCB109_Ecoli_NOTEX_replicate2`)) %>%
  dplyr::filter(cov > 3) %>%
  group_by(group_sec) %>%
  summarise(n = n())
cg
sum()
pdf(here("figures/R_plots/210217_UTR5_cDNA-PCR_replicates.pdf"), 
    width = 3.2, height = 3.2, paper = "special", onefile=FALSE)
ggplot(data = cg, aes(x = `201208_PCB109_Ecoli_NOTEX_replicate1`, 
                      y  = `201208_PCB109_Ecoli_NOTEX_replicate2`, fill =group_sec)) +
  geom_point(size = 3, alpha = 0.7, shape = 21) +
  geom_abline(linetype = "dashed") +
  scale_x_continuous(limits = c(0,300)) +
  scale_y_continuous(limits = c(0,300)) +
  stat_cor(method = "pearson", label.x = 1) +
  xlab("peakcaller") +
  ylab("compare set") +
  theme_Publication_white() +
  scale_fill_manual(values = ibm_colors[c(4,1)]) +
  coord_equal()
dev.off()

# DISTANCE TO ILLUMINA DATA ---------------------------------
tss_total_to_illumina <- rbind(tss_data_pychopper_auto_cutadapt_clipped %>% mutate(method = "pychopper_auto_cutadapt_clipped"),
                               tss_data_pychopper_auto_cutadapt_SSP_clipped %>% mutate(method = "pychopper_auto_cutadapt_SSP_clipped"),
                               tss_data_pychopper_auto %>% mutate(method = "pychopper_auto"),
                               tss_data_notrimming %>% mutate(method = "notrimming")) %>%
  group_by(method, id_name, dataset) %>%
  dplyr::mutate(total_n = n(),
                group_sec = ifelse(total_n >1, "has_secondary", "single")) %>%
  ungroup() %>%
  dplyr::filter(TSS_type == "primary", type == "CDS") %>%
  mutate(polyG_detected = str_count(string = substring(tss_sequence,31,41),pattern = "GGG+") >= 1) %>%
  dplyr::select(gene, dataset, method, UTR5) 

tss_total_to_illumina2 <- left_join(tss_total_to_illumina, 
                                    smrt_tss %>% dplyr::rename(UTR5_Ill = UTR5) %>% dplyr::select(gene,UTR5_Ill), by = "gene")   %>%
  mutate(distance = UTR5 - UTR5_Ill,
         mode = substr(dataset, 8,10)) %>%
  dplyr::filter(!is.na(dataset)) %>%
  dplyr::filter(dataset %in% bc_to_sample$sample[c(1,5:12,15:16)]) 

levels(as.factor(tss_total_to_illumina2$method))
tss_total_to_illumina2$dataset <- factor(tss_total_to_illumina2$dataset,
                                         levels = (bc_to_sample$sample[c(1,15,16,11,12,5,7,6,8:10)]))

ggplot(data = tss_total_to_illumina2, aes(x = distance, fill= mode)) +
  facet_grid(rows = vars(dataset), cols = vars(method)) +
  scale_x_continuous(limits = c(-30,30), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, .65)) +
  geom_histogram(binwidth = 1, aes(y=..density..), color = "black") +
  scale_fill_manual(values = cbf1[c(2,5,3)]) +
  geom_vline(xintercept = -11, linetype = "dashed") +
  theme_Publication_white() 





