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
    dplyr::mutate(tts_sequence = ifelse(strand == "+", as.character(input_genome_fasta$chr[(TTS - 20):(TTS + 10 )]),
                                        as.character(reverseComplement(input_genome_fasta$chr[(TTS - 10):(TTS + 20 )]))),
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

tts_peaks_pipe <- function(file_dir, names, map_method, p_f, m_f, gff = ecoli_gff){
  
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

plot_3end_distance <- function(trimtype = c("untrimmed", "trimmed"), compset = c("diff", "SMRT"), output = c("plot", "stats")){
  
  inputdf <- if(trimtype == "untrimmed"){tts_data_untrimmed}else{tts_data_trimmed}
  compdf  <- if(compset == "diff"){dar_tts}else{smrt_tts}
  
  # > calc
  outputdf <- inputdf %>%
    dplyr::filter(TTS_type == "primary", type == "CDS") %>%
    dplyr::select(gene, sample, method, UTR3) %>%
    left_join(compdf, by = "gene") %>%
    mutate(distance = UTR3.x - UTR3.y) %>%
    remove_missing(vars = "distance") %>%
    mutate(mode = substr(sample, 1,3)) 
  
  # > reorder levels
  outputdf$sample <- factor(outputdf$sample,
                            levels = (bc_to_sample$sample[c(1,10,11,8,9,2,4,3,5,6,7)]))
  
  # > plot
  if(output == "plot"){
    raw_reads_plotting(outputdf, distance, sample, mode, cbf1[c(2,5,3)]) +
      facet_grid(rows = vars(sample)) +
      geom_histogram(binwidth = 1, aes(y=..density..), color = "black") +
      scale_x_continuous(limits = c(-17,17), expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0), limits = c(0,0.75)) 
  }else{
    outputdf %>%
      dplyr::filter(distance > -17 & distance < +17) %>%
      group_by(sample) %>%
      distinct(gene) %>%
      summarise(n = n())
  }
}


# load & tidy data ----

## Peak tables ====
### read in pychopper auto trimmed data > cutadapt polyA > cutadapt SSP > clipped removed ####
files <- list.files(paste0(dir,"/data/tts_data/tts_data_pychopper_auto_cutadapt_SSP_clipped/"), recursive = T, pattern = ".narrowPeak.counts")
tts_data_trimmed <- tts_peaks_pipe("tts_data_pychopper_auto_cutadapt_SSP_clipped",
                                   str_split_fixed(str_split_fixed(str_split_fixed(files[which(1:length(files) %% 2 == 0)], "\\/", n = 3)[,3],".plus",2)[,1],"_fu",2)[,1],
                                   "trimmed",
                                   files[which(1:length(files) %% 2 == 0)],
                                   files[which(1:length(files) %% 2 == 1)])

fwrite(tts_data_trimmed, paste0(dir, "/tables/tts_tables/tts_data_trimmed.tsv"), col.names = T, sep = "\t")

### read in untrimmed - raw mapped data ####
files <- list.files(paste0(dir,"/data/tts_data/tts_data_notrimming/"), recursive = T, pattern = ".narrowPeak.counts")
tts_data_untrimmed <- tts_peaks_pipe("tts_data_notrimming",
                                     str_split_fixed(str_split_fixed(str_split_fixed(files[which(1:length(files) %% 2 == 0)], "\\/", n = 3)[,3],".plus",2)[,1],"_fu",2)[,1],
                                     "untrimmed",
                                     files[which(1:length(files) %% 2 == 0)],
                                     files[which(1:length(files) %% 2 == 1)])

fwrite(tts_data_untrimmed, paste0(dir, "/tables/tts_tables/tts_data_untrimmed.tsv"), col.names = T, sep = "\t")

## TTS from other studies ====
### Term-seq (Dar et al) ####
dar_tts <- read_xlsx(paste0(dir,"data/comparison_data/TTS/dar_gky274_supplemental_files.xlsx"),sheet = "Table S1", skip = 11) %>%
  dplyr::rename(TTS = `primary 3' end position`,
                strand = `Gene strand`) %>%
  mutate(gene = paste0("b",str_split_fixed(`Locus tag`, "_",2)[,2]),
         short_gene = `Gene name`) %>%
  left_join(ecoli_gff, by = "short_gene") %>%
  distinct(TTS, short_gene, .keep_all = T) %>%
  dplyr::mutate(UTR3 = ifelse(strand == "+", TTS - `gene to`, `gene fr` - TTS)) %>%
  dplyr::rename(gene = gene.y) %>%
  dplyr::select(TTS, strand, gene, UTR3)

### SMRT-Cap-results ####
smrt_tts <- read_xlsx(paste0(dir,"data/comparison_data/TTS/smrt_cap-41467_2018_5997_MOESM4_ESM.xlsx"),sheet = "Rich_TTS", skip = 1) %>%
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

## Elife-Termseq-results ====
elife_tts <- read_xlsx(paste0(dir,"data/comparison_data/TTS/elife-62438-supp1-v1.xlsx"),sheet = "LB 0.4", skip = 1) %>%
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

## Comparison ====
summary_tts_ONT <-  rbindlist(list(tts_data_trimmed,
                                   tts_data_untrimmed)) %>%
  dplyr::filter(type == "CDS") %>%
  group_by(sample, TTS_type, method) %>%
  summarise(n = n()) %>%
  mutate(mode = substr(sample, 1,3)) %>%
  remove_missing(vars = "TTS_type")

## Distance between diff RNA-seq & ONT ====
end3_comp <- rbindlist(list(tts_data_trimmed,
                            tts_data_untrimmed)) %>%
  dplyr::filter(TTS_type == "primary", type == "CDS") %>%
  dplyr::select(gene, sample, method, UTR3) %>%
  left_join(dar_tts, by = "gene") %>%
  mutate(distance = UTR3.x - UTR3.y) %>%
  remove_missing(vars = "distance") %>%
  mutate(mode = substr(sample, 1,3)) %>%
  group_by(sample, method) %>%
  mutate(distance_p = sum(distance == 0)/n()*100) %>%
  distinct(sample,method, distance_p, mode) 

## REPLICATE comparison ====
### for point plots ####
tts_total_comparison <- tts_data_trimmed %>%
  dplyr::filter(cov >= 5) %>%
  dplyr::filter(TTS_type == "primary", type == "CDS") %>%
  dplyr::select(gene, sample, method, UTR3) %>%
  rbind(dar_tts %>% mutate(sample = "illumina", 
                           method = "illumina") %>%
          dplyr::select(gene, sample, method, UTR3)) %>%
  rbind(smrt_tts %>% ungroup () %>% mutate(sample = "Smrt_CAP", 
                                           method = "Smrt_CAP") %>%
          dplyr::select(gene, sample, method, UTR3)) %>%
  dplyr::filter(UTR3 >= 0 & !is.na(UTR3) & is.finite(UTR3), UTR3 <= 300) %>%
  dplyr::select(-method) %>%
  pivot_wider(names_from = sample, values_from = UTR3, values_fn = {max}) %>%
  left_join(ecoli_gff) %>%
  dplyr::filter(!is.na(gene)) 

### correlation matrix ####
#### pairwise complete Pearson correlation ####
res             <- cor(tts_total_comparison[c(2,11,12,9,10,3,5,4,6,7,8,13,14)],  
                       method = "pearson", use = "pairwise.complete.obs")
res_gg          <- reshape2::melt(get_upper_tri(res))

#### pairwise complete Pearson observations ####
res_counts      <- pairwiseCount(tts_total_comparison[c(2,11,12,9,10,3,5,4,6,7,8,13,14)], diagonal = F)
res_counts_gg   <- reshape2::melt(get_lower_tri(res_counts))

### for utr5 comparison ####
utr3 <- tts_total_comparison %>%
  dplyr::filter(type == "CDS") %>%
  pivot_longer(cols = 2:14,names_to = "sample", values_to = "UTR3") %>% 
  mutate(mode = str_sub(sample, 1, 3)) %>%
  distinct(UTR3, gene, sample, .keep_all = T)

utr3_logo <- tts_data_trimmed %>%
  dplyr::filter(sample %in% "PCB109_PCR12_Ecoli_NOTEX_replicate4", 
                TTS_type %in% "primary", 
                type == "CDS", cov >= 3) %>%
  distinct(gene, .keep_all = T) %>%
  rowwise() %>%
  dplyr::mutate(tts_sequence = ifelse(strand == "+", as.character(ecoli_fasta$chr[(TTS-20):(TTS+10)]),
                                      as.character(reverseComplement(ecoli_fasta$chr[(TTS-10):(TTS+20)])))) %>%
  dplyr::select(tts_sequence) 

# PLOTS ----

## reorder levels ====
summary_tts_ONT$sample <- factor(summary_tts_ONT$sample,
                                 levels = rev(bc_to_sample$sample[c(1,10,11,8,9,2,4,3,5,6,7)]))
summary_tts_ONT$mode <- factor(summary_tts_ONT$mode,
                               levels = c("RNA", "DCS", "PCB"))
summary_tts_ONT$TTS_type <- factor(summary_tts_ONT$TTS_type,
                                   levels = rev(c("primary", "secondary", "internal", "rest")))
summary_tts_ONT$method <- factor(summary_tts_ONT$method,
                                 levels = c("untrimmed", "trimmed"))

end3_comp$sample <- factor(end3_comp$sample,
                           levels = rev(bc_to_sample$sample[c(1,10,11,8,9,2,4,3,5,6,7)]))
end3_comp$method <- factor(end3_comp$method,
                           levels = c("untrimmed", "trimmed"))

utr3$sample <- factor(utr3$sample,
                      levels = c("Smrt_CAP", "illumina",rev(bc_to_sample$sample[c(1,10,11,8,9,2,4,3,5,6,7)])))


## plotting ==== 

### Number of 3´ends in category - Supplementary Fig. 18A #### 
raw_reads_plotting(summary_tts_ONT, n, sample, TTS_type, cbf1_high) +
  geom_bar(stat = "identity", color = "black", position = position_stack()) +
  scale_x_continuous(expand = c(0,0), limits = c(0,4000)) +
  facet_grid(cols = vars(method)) +
  xlab("Number of 3´ends in category") 

### Proportion of ONT 3´ends with 0 distance to Tern-seq 3´ends - Supplementary Fig. 18B #### 
raw_reads_plotting(end3_comp, distance_p, sample, mode, cbf1[c(2,5,3)]) +
  geom_bar(aes(group = method, alpha = method),
           stat = "identity", color = "black", 
           position = position_dodge2(width = 2)) +
  xlab("Proportion of ONT 3´ends with 0 distance to Term-seq 3´ends (%)") +
  scale_x_continuous(limits = c(0,50), expand = c(0,0))

### Correlation matrix - Supplementary Fig. 19 #### 
#### Part 1 ####
corr_matrix_plot(res_gg, Var2, Var1, value) +
  geom_tile(color = "black", size = 0.3, aes(width = value, height = value)) +
  geom_text(aes(label=round(value, digits = 2)), color = "white", size = 4) 

#### Part 2 ####  
corr_matrix_plot(res_counts_gg, Var2, Var1, value) +
  geom_point(shape = 21, color = "black") +
  scale_fill_gradientn(colours = brewer.pal(name = "Blues", n = 9),
                       limit = c(0,600), space = "Lab", 
                       name="Pearson\nCorrelation")

### 3´UTR length distribution - Supplementary Fig. 20A #### 
ggplot(data = utr3, aes(y = sample, fill = mode)) +
  geom_density_ridges(stat = "binline",binwidth = 4,
                      aes(x = UTR3, height =..ndensity..), 
                      scale = 0.9, alpha = 1) +
  theme_Publication_white() +
  scale_fill_manual(values = cbf1[c(2,1,5,3,4)])

### Promoter logo - Supplementary Fig. 20B #### 
ggplot() + 
  geom_logo(utr3_logo, font = "helvetica_bold", col_scheme = acgt_color_scale, seq_type = "dna") + 
  theme_logo() +
  theme_Publication_white() +
  theme(panel.grid.major = element_line(colour = NA),
        axis.ticks.x = element_line(colour = NA), 
        axis.text.x = element_text(size = 0)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0,0.4), expand = c(0,0))

### 3´end histograms accuracy (including n) - Fig. 3C #### 
plot_3end_distance(trimtype = "trimmed", compset = "diff", output = "plot")  
plot_3end_distance(trimtype = "trimmed", compset = "diff", output = "stats")  
plot_3end_distance(trimtype = "trimmed", compset = "SMRT", output = "plot")  
plot_3end_distance(trimtype = "trimmed", compset = "SMRT", output = "stats")  
