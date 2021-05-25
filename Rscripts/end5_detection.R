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
    dplyr::mutate(tss_sequence = ifelse(strand == "+", as.character(input_genome_fasta$chr[(TSS - 40):(TSS + 10 )]),
                                        as.character(reverseComplement(input_genome_fasta$chr[(TSS - 10):(TSS + 40 )]))),
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


plot_5end_distance <- function(trimtype = c("untrimmed", "trimmed"), compset = c("diff", "SMRT"), output = c("plot", "stats")){
  
  inputdf <- if(trimtype == "untrimmed"){tss_data_untrimmed}else{tss_data_trimmed}
  compdf  <- if(compset == "diff"){tex_tss}else{smrt_tss}
  
  # > calc
  outputdf <- inputdf %>%
    dplyr::filter(TSS_type == "primary", type == "CDS") %>%
    dplyr::select(gene, sample, method, UTR5) %>%
    left_join(compdf, by = "gene") %>%
    mutate(distance = UTR5.x - UTR5.y) %>%
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

point_cor_ends <-  function(mydf, myx, myy, myfill){
  df <- {{mydf}} %>%
    dplyr::filter(!is.na({{myx}}),
                  !is.na({{myy}})) %>%
    mutate(density = get_density(({{myx}}), ({{myy}}))) 
  
  ggplot(data = df, aes(x = {{myx}}, y = {{myy}}, fill = {{myfill}})) +
    geom_abline(linetype = "dashed", slope = 1) +
    geom_point(alpha = 1, shape = 21, color = "black", size = 4) +
    scale_fill_gradientn(colours = brewer.pal(name = "Blues", n = 9)) +
    theme_Publication_white() +
    scale_x_continuous(limits = c(0, 300)) +
    scale_y_continuous(limits = c(0, 300)) +
    stat_cor(method = "pearson", label.x = 1) +
    coord_equal() 
}

detect_secondary <- function(inputdf,sample_s){
  inputdf %>%
    dplyr::filter(sample %in% sample_s, TSS_type %in% c("primary", "secondary")) %>%
    dplyr::filter(cov >= 5) %>%
    group_by(id_name) %>%
    dplyr::mutate(total_n = n(),
                  group_sec = ifelse(total_n >1, "has_secondary", "single")) %>%
    ungroup()
}


# load & tidy data ----

dir <- here()

## Peak tables ====

### read in pychopper auto trimmed data > cutadapt polyA > cutadapt SSP > clipped removed ####
files <- list.files(paste0(dir,"/data/tss_data/tss_data_pychopper_auto_cutadapt_SSP_clipped/"), recursive = T, pattern = ".narrowPeak.counts")
tss_data_trimmed <- tss_peaks_pipe("tss_data_pychopper_auto_cutadapt_SSP_clipped",
                                   str_split_fixed(str_split_fixed(str_split_fixed(files[which(1:length(files) %% 2 == 0)], "\\/", n = 3)[,3],".plus",2)[,1],"_fu",2)[,1],
                                   "trimmed",
                                   files[which(1:length(files) %% 2 == 0)],
                                   files[which(1:length(files) %% 2 == 1)])

fwrite(tss_data_trimmed, paste0(dir, "/tables/tss_tables/tss_data_trimmed.tsv"), col.names = T, sep = "\t")

### read in untrimmed - raw mapped data ####
files <- list.files(paste0(dir,"/data/tss_data/tss_data_notrimming/"), recursive = T, pattern = ".narrowPeak.counts")
tss_data_untrimmed <- tss_peaks_pipe("tss_data_notrimming",
                                     str_split_fixed(str_split_fixed(str_split_fixed(files[which(1:length(files) %% 2 == 0)], "\\/", n = 3)[,3],".plus",2)[,1],"_fu",2)[,1],
                                     "untrimmed",
                                     files[which(1:length(files) %% 2 == 0)],
                                     files[which(1:length(files) %% 2 == 1)])

fwrite(tss_data_untrimmed, paste0(dir, "/tables/tss_tables/tss_data_untrimmed.tsv"), col.names = T, sep = "\t")

### write to Supplementary Table 4 ####
tss_data_untrimmed %>%
  dplyr::rename(end5 = TSS, end5_type = TSS_type) %>%
  dplyr::select(gene, end5, end5_type, sample) %>%
  dplyr::filter(!is.na(gene)) %>%
  arrange(gene) %>%
  write_xlsx(path = here("tables/Supplementary_Table4.xlsx"))



## TSS from other studies ====
### TEX-TSS-results ####
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

## Distance between diff RNA-seq & ONT ====
end5_comp <- rbindlist(list(tss_data_trimmed,
                            tss_data_untrimmed)) %>%
  dplyr::filter(TSS_type == "primary", type == "CDS") %>%
  dplyr::select(gene, sample, method, UTR5) %>%
  left_join(tex_tss, by = "gene") %>%
  mutate(distance = UTR5.x - UTR5.y) %>%
  remove_missing(vars = "distance") %>%
  mutate(mode = substr(sample, 1,3)) %>%
  group_by(sample, method) %>%
  mutate(distance_p = ifelse(mode != "RNA", sum(distance == 0)/n()*100, sum(distance == -12)/n()*100)) %>%
  distinct(sample,method, distance_p, mode) 

## REPLICATE comparison ====
### for point plots ####
tss_total_comparison <- tss_data_untrimmed %>%
  dplyr::filter(cov >= 5) %>%
  dplyr::filter(TSS_type == "primary", type == "CDS") %>%
  dplyr::select(gene, sample, method, UTR5) %>%
  rbind(tex_tss %>% mutate(sample = "illumina", 
                           method = "illumina") %>%
          dplyr::select(gene, sample, method, UTR5)) %>%
  rbind(smrt_tss %>% ungroup () %>% mutate(sample = "Smrt_CAP", 
                                           method = "Smrt_CAP") %>%
          dplyr::select(gene, sample, method, UTR5)) %>%
  dplyr::filter(UTR5 >= 0 & !is.na(UTR5) & is.finite(UTR5), UTR5 <= 300) %>%
  dplyr::select(-method) %>%
  pivot_wider(names_from = sample, values_from = UTR5, values_fn = {max}) %>%
  left_join(ecoli_gff) %>%
  dplyr::filter(!is.na(gene)) 

### for utr5 comparison ####
utr5 <- tss_total_comparison %>%
  dplyr::filter(type == "CDS") %>%
  pivot_longer(cols = 2:14,names_to = "sample", values_to = "UTR5") %>% 
  mutate(mode = str_sub(sample, 1, 3)) %>%
  mutate(UTR5 = ifelse(mode == "RNA", UTR5+12, UTR5)) %>%
  distinct(UTR5, gene, sample, .keep_all = T)

utr5_logo <- tss_data_untrimmed %>%
  dplyr::filter(sample %in% "PCB109_PCR12_Ecoli_NOTEX_replicate4", 
                TSS_type %in% "primary", 
                type == "CDS", cov >= 3) %>%
  distinct(gene, .keep_all = T) %>%
  rowwise() %>%
  dplyr::mutate(tss_sequence = ifelse(strand == "+", as.character(ecoli_fasta$chr[(TSS-40):(TSS)]),
                                      as.character(reverseComplement(ecoli_fasta$chr[(TSS):(TSS+40)])))) %>%
  dplyr::select(tss_sequence) 

### correlation matrix ####
#### pairwise complete Pearson correlation ####
res             <- cor(tss_total_comparison[c(2,11,12,9,10,3,5,4,6,7,8)],  
                       method = "pearson", use = "pairwise.complete.obs")
res_gg          <- reshape2::melt(get_upper_tri(res))

#### pairwise complete Pearson observations ####
res_counts      <- pairwiseCount(tss_total_comparison[c(2,11,12,9,10,3,5,4,6,7,8)], diagonal = F)
res_counts_gg   <- reshape2::melt(get_lower_tri(res_counts))

### compare secondary & primary ONT to SMRT ####
c_SMRT_sec <- detect_secondary(tss_data_untrimmed,
                               "PCB109_PCR12_Ecoli_NOTEX_replicate4") %>% 
  dplyr::filter(group_sec == "has_secondary", TSS_type == "secondary", type == "CDS") %>%
  dplyr::select(gene, sample, method, UTR5, group_sec) %>%
  left_join(smrt_tss, by = "gene") %>%
  dplyr::filter(abs(UTR5.x-UTR5.y) < 5) %>%
  dplyr::select(gene, UTR5.x)

c_SMRT <- detect_secondary(tss_data_untrimmed,
                           "PCB109_PCR12_Ecoli_NOTEX_replicate4") %>% 
  dplyr::filter(TSS_type == "primary", type == "CDS") %>%
  dplyr::mutate(group = ifelse(group_sec == "single", "single",
                               ifelse(group_sec == "has_secondary" & gene %in% c_SMRT_sec$gene, "SMRT_first", "rest"))) %>%
  left_join(c_SMRT_sec, by = "gene") %>%
  dplyr::mutate(UTR5 = ifelse(group == "SMRT_first", UTR5.x, UTR5)) %>%
  dplyr::select(-UTR5.x) %>%
  dplyr::select(gene, sample, method, UTR5, group) %>%
  left_join(smrt_tss, by = "gene") %>%
  remove_missing(name = "UTR5.y")

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

end5_comp$sample <- factor(end5_comp$sample,
                                 levels = rev(bc_to_sample$sample[c(1,10,11,8,9,2,4,3,5,6,7)]))
end5_comp$method <- factor(end5_comp$method,
                                 levels = c("untrimmed", "trimmed"))

utr5$sample <- factor(utr5$sample,
                           levels = c("Smrt_CAP", "illumina",rev(bc_to_sample$sample[c(1,10,11,8,9,2,4,3,5,6,7)])))

## plotting ==== 

### Number of 5´ends in category - Supplementary Fig. 14A #### 
raw_reads_plotting(summary_tss_ONT, n, sample, TSS_type, cbf1_high) +
  geom_bar(stat = "identity", color = "black", position = position_stack()) +
  scale_x_continuous(expand = c(0,0), limits = c(0,6000)) +
  facet_grid(cols = vars(method)) +
  xlab("Number of 5´ends in category") 

### Proportion of ONT 5´ends with 0 distance to diff RNA-seq 5´ends - Supplementary Fig. 14B #### 
raw_reads_plotting(end5_comp, distance_p, sample, mode, cbf1[c(2,5,3)]) +
  geom_bar(aes(group = method, alpha = method),
           stat = "identity", color = "black", 
           position = position_dodge2(width = 2)) +
  xlab("Proportion of ONT 5´ends with 0 distance to diff RNA-seq 5´ends (%)") +
  scale_x_continuous(limits = c(0,50), expand = c(0,0))

### 5´end histograms accuracy (including n) - Fig. 3B #### 
plot_5end_distance(trimtype = "untrimmed", compset = "diff", output = "plot")  
plot_5end_distance(trimtype = "untrimmed", compset = "diff", output = "stats")  
plot_5end_distance(trimtype = "untrimmed", compset = "SMRT", output = "plot")  
plot_5end_distance(trimtype = "untrimmed", compset = "SMRT", output = "stats")  

### Correlation 5´end detection reproducibility 1 - Supplementary Fig. 15A #### 
point_cor_ends(tss_total_comparison, PCB109_PCR12_Ecoli_NOTEX_replicate4, DCS109_Ecoli_NOTEX_replicate2, density) 

### Correlation 5´end detection reproducibility 2 - Supplementary Fig. 15B #### 
point_cor_ends(tss_total_comparison, PCB109_PCR12_Ecoli_NOTEX_replicate4, RNA001_Ecoli_TEX_replicate1, density) 

### Correlation matrix - Supplementary Fig. 15C #### 
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

### Correlation 5´end detection ONT vs SMRT - Supplementary Fig. 16A #### 
point_cor_ends(tss_total_comparison, PCB109_PCR12_Ecoli_NOTEX_replicate4, Smrt_CAP, density) 

### Correlation 5´end detection ONT vs SMRT secondary sites - Supplementary Fig. 16B #### 
ggplot(data = c_SMRT, 
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

### 5´UTR length distribution - Supplementary Fig. 17A #### 
ggplot(data = utr5, aes(y = sample, fill = mode)) +
  geom_density_ridges(stat = "binline",binwidth = 4,
                      aes(x = UTR5, height =..ndensity..), 
                      scale = 0.9, alpha = 1) +
  theme_Publication_white() +
  scale_fill_manual(values = cbf1[c(2,1,5,3,4)])

### Promoter logo - Supplementary Fig. 17B #### 
ggplot() + 
  geom_logo(utr5_logo, font = "helvetica_bold", col_scheme = acgt_color_scale, seq_type = "dna") + 
  theme_logo() +
  theme_Publication_white() +
  theme(panel.grid.major = element_line(colour = NA),
        axis.ticks.x = element_line(colour = NA), 
        axis.text.x = element_text(size = 0)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0,0.4), expand = c(0,0))
