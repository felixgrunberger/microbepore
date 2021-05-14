# >> gene body coverage << #

# load libraries ----
library(here)
source(here("Rscripts/load_libraries.R"))

# functions & defs ----
filter_bed_files <- function(input_file){
  
  strand_f <- str_split_fixed(str_split_fixed(string = input_file, pattern = "\\.coverage", n = 2)[,1], "\\.", 2)[,2]
  
  suppressMessages(vroom(input_file,col_names = c("seqid","TSS", "TTS", "gene","rel_pos", "counts"), num_threads = 8, progress = F)) %>%
    as_tibble() %>%
    mutate(strand = strand_f) %>%
    group_by(gene) %>%
    mutate(perc_pos = ifelse(strand == "plus", round(scales::rescale(rel_pos, to=c(0,100)), digits = 0),
                             round(scales::rescale(rel_pos, to=c(100,0)), digits= 0))) %>%
    dplyr::select(-strand) 
}

merge_bed_files <- function(input_plus, input_minus){
  
  # Suppress summarise info
  options(dplyr.summarise.inform = FALSE)
  
  rbind(input_plus, input_minus) %>%
    group_by(gene) %>%
    mutate(perc_coverage = scales::rescale(counts, to = c(0,100))) %>%
    group_by(perc_pos, gene) %>%
    summarise(C_mean = mean(perc_coverage),
              C_max = max(perc_coverage),
              C_min = min(perc_coverage)) %>%
    group_by(perc_pos) %>%
    summarise(C_mean_sum = mean(C_mean),
              C_m_sum = min(C_mean)) %>%
    mutate(median_C = median(C_mean_sum),
           IQR = IQR(C_mean_sum),
           QCoV = IQR/median_C)
}

modify_coverage_files <- function(folder, minimum_wanted_seqs = 10, output = c("normal", "genesizes")){
  
  coverage_files   <- list.files(folder, recursive = T, pattern = ".coverage$")
  plus_c_files     <- coverage_files[which(1:length(coverage_files) %% 2 == 0)] 
  minus_c_files    <- coverage_files[which(1:length(coverage_files) %% 2 == 1)] 
  dataset_names    <- str_split_fixed(str_split_fixed(str_split_fixed(plus_c_files, "\\/", n = 3)[,3],"_fu",2)[,1], "\\.", 2)[,1]
  
  coverage_frame   <- data.table()
  for(i in seq_along(plus_c_files)){
    
    # files
    f <- folder
    
    tic("normalise coverage files")
    print(paste0("file number ", i, " of ", length(plus_c_files)))
    
    if(output == "normal"){
      # coverage to normalised coverage
      p_t <- filter_bed_files(paste0(f,plus_c_files[i])) %>%
        group_by(gene) %>% dplyr::filter(min(counts) >= minimum_wanted_seqs) %>% ungroup()
      m_t <- filter_bed_files(paste0(f,minus_c_files[i])) %>%
        group_by(gene) %>% dplyr::filter(min(counts) >= minimum_wanted_seqs) %>% ungroup()
      a_t <- merge_bed_files(p_t, m_t) %>%
        mutate(dataset = dataset_names[i])
      coverage_frame <- rbind(coverage_frame, a_t)
    }else{
      # coverage to normalised coverage
      p_t <- filter_bed_files_size(paste0(f,plus_c_files[i])) %>%
        group_by(gene) %>% dplyr::filter(min(counts) >= minimum_wanted_seqs) %>% ungroup()
      m_t <- filter_bed_files_size(paste0(f,minus_c_files[i])) %>%
        group_by(gene) %>% dplyr::filter(min(counts) >= minimum_wanted_seqs) %>% ungroup()
      a_t <- merge_bed_files_sizes(p_t, m_t) %>%
        mutate(dataset = dataset_names[i])
      coverage_frame <- rbind(coverage_frame, a_t)
    }
    
    toc()
  }
  return(coverage_frame)
}

keep_highest_site <- function(inputdf, selected_end, merge_w = 20, cov_min = 3){
  
  inputdf %>% 
    distinct({{selected_end}}, .keep_all = T) %>%
    arrange({{selected_end}}) %>%
    mutate(index = lag({{selected_end}}, default = 1) + as.integer(merge_w),
           index1 = cumsum(ifelse(index >= {{selected_end}}, 0, 1))+1) %>%
    dplyr::group_by(index1) %>%
    dplyr::filter(cov == max(cov),
                  cov >= cov_min) %>%
    ungroup() %>%
    group_by(gene) %>%
    dplyr::slice(which.max(cov)) %>%
    ungroup() %>%
    dplyr::select(id_name, {{selected_end}}, strand)
  
}
  
filter_bed_files_size <- function(input_file, min_wanted = 10){

  strand_f <- str_split_fixed(str_split_fixed(string = input_file, pattern = "\\.coverage", n = 2)[,1], "\\.", 2)[,2]
  
  suppressMessages(vroom(input_file,col_names = c("seqid","TSS", "TTS", "gene","rel_pos", "counts"), num_threads = 8, progress = F)) %>%
    as_tibble() %>%
    mutate(strand = strand_f) %>%
    group_by(gene) %>%
    mutate(perc_pos = ifelse(strand == "plus", round(scales::rescale(rel_pos, to=c(0,100)), digits = 0),
                             round(scales::rescale(rel_pos, to=c(100,0)), digits= 0))) %>%
    dplyr::select(-strand) %>%
    group_by(gene) %>% 
    dplyr::filter(min(counts) >= min_wanted) %>% 
    ungroup() %>%
    left_join(ecoli_gff %>% dplyr::select(id_name, width) %>% dplyr::rename(gene = id_name), by = c("gene")) %>%
    mutate(read_group = ifelse(width <= 500, "sub500",
                               ifelse(width > 500 & width <= 1000, "sub1000",
                                      ifelse(width > 1000 & width <= 1500, "sub1500",
                                             ifelse(width > 1500 & width <= 2000, "sub1500",
                                                    ifelse(width > 2000, "big2000",NA)))))) 
  
}

merge_bed_files_sizes <- function(input_plus, input_minus){
  
  # Suppress summarise info
  options(dplyr.summarise.inform = FALSE)
  
  f <- rbind(input_plus, input_minus) %>%
    group_by(gene) %>%
    mutate(perc_coverage = scales::rescale(counts, to = c(0,100))) %>%
    group_by(perc_pos, gene,read_group) %>%
    summarise(C_mean = mean(perc_coverage),
              C_max = max(perc_coverage),
              C_min = min(perc_coverage)) %>%
    group_by(perc_pos,read_group) %>%
    summarise(C_mean_sum = mean(C_mean),
              C_m_sum = min(C_mean)) %>%
    ungroup() %>%
    group_by(read_group) %>%
    mutate(median_C = median(C_mean_sum),
           IQR = IQR(C_mean_sum),
           QCoV = IQR/median_C)
  
  gene_per_group <- rbind(input_plus, input_minus) %>%
    distinct(gene, read_group) %>%
    group_by(read_group) %>%
    summarise(n = n())
  
  f2 <- left_join(f, gene_per_group, by = "read_group")
  
  return(f2)
}



# load & tidy data ----

## prepare 5´-3´end tables for coverage calculations ====

### primary 5´end ####
dir <- here()

tss <- vroom(file = paste0(dir,"/tables/tss_tables/tss_data_untrimmed.tsv"), num_threads = 8, progress = F) %>%
  mutate(mode = str_sub(sample, 1,3),
         TSS = ifelse(mode == "RNA" & strand == "+", TSS - 12, 
                      ifelse(mode == "RNA" & strand == "-", TSS + 12,TSS))) %>%
  dplyr::filter(TSS_type == "primary", type == "CDS") %>%
  keep_highest_site(inputdf = .,selected_end = TSS)

### primary 3´end ####
tts <- vroom(file = paste0(dir,"/tables/tts_tables/tts_data_trimmed.tsv"), num_threads = 8, progress = F) %>%
  mutate(mode = str_sub(sample, 1,3)) %>%
  dplyr::filter(TTS_type == "primary", type == "CDS") %>%
  keep_highest_site(inputdf = .,selected_end = TTS)

### find gens with annotated primary 5´and 3´end ####
w <- left_join(tss, tts, by = c("id_name", "strand")) %>%
  dplyr::filter(!is.na(TSS), !is.na(TTS)) %>% 
  mutate(seqnames = ecoli_gff$seqid[1]) %>%
  dplyr::select(seqnames, TSS, TTS, id_name, strand)

### write to bed-like file which can be used with bedtools coverage
fwrite(w %>% dplyr::filter(strand == "+") %>% dplyr::select(-strand), 
       paste0(dir,"/tables/transcript_tables/transcripts.plus.bedgraph"), sep = "\t", col.names = F, quote = F)
fwrite(w %>% dplyr::filter(strand == "-") %>% dplyr::select(-strand), 
       paste0(dir,"/tables/transcript_tables/transcripts.minus.bedgraph"), sep = "\t", col.names = F, quote = F)

## read in files from bedtools coverage ====
dir <- "/Volumes/EX_SSD/"
### full-length > polyA-trimmed > polyA & SSP adapter trimmed > clipping removed > stranded ####
cov_trimmed   <- modify_coverage_files(folder = paste0(dir, "/data/coverage_data/coverage_data_pychopper_auto_cutadapt_SSP_clipped_stranded/"),
                                       output = "normal")

### notrimming > stranded ####
cov_untrimmed <- modify_coverage_files(folder = paste0(dir,"/data/coverage_data/coverage_data_notrimming_stranded/"),
                                       output = "normal")

### merge datasets ####
cov_sets <- rbind(cov_untrimmed %>% mutate(method = "untrimmed"),
                  cov_trimmed %>% mutate(method = "trimmed")) %>%
  dplyr::left_join(old_new, by = c("dataset" = "old_name")) %>%
  mutate(sample = new_name) %>%
  dplyr::select(-new_name, -dataset) %>%
  dplyr::filter(!is.na(sample))

## calculate QCoV sizes
cov_trimmed_sizes   <- modify_coverage_files(folder = paste0(dir, "/data/coverage_data/coverage_data_pychopper_auto_cutadapt_SSP_clipped_stranded/"),
                                             output = "genesizes")

###

# PLOTS ----

## reorder levels ====
cov_sets$sample <- factor(cov_sets$sample,
                          levels = (bc_to_sample$sample[c(1,10,11,8,9,2,4,3,5,6,7)]))

cov_sets$method <- factor(cov_sets$method,
                          levels = rev(c("untrimmed", "trimmed")))

cov_trimmed_sizes$sample <- factor(cov_trimmed_sizes$sample,
                                   levels = rev(bc_to_sample$sample[c(1,10,11,8,9,2,4,3,5,6,7)]))
cov_trimmed_sizes$read_group <- factor(cov_trimmed_sizes$read_group,
                                       levels = (c("sub500", 
                                                   "sub1000",
                                                   "sub1500", 
                                                   "big2000")))
## plotting ==== 

### Gene body coverage - Fig. 4A #### 
ggplot(data = cov_sets ,
       aes(x = perc_pos, y = C_mean_sum)) +
  geom_line(size = 1.2,aes(linetype = method), color = "black") +
  facet_grid(cols = vars(sample)) +
  scale_x_continuous(limits = c(0,100), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,100), expand = c(0,0)) +
  geom_ribbon(aes(fill = method, ymin = 0, ymax = C_mean_sum), alpha = 0.5, color = NA) +
  scale_fill_manual(values = c("#AD9D86","#A1CEC2")) +
  theme_Publication_white() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "black")) +
  ylab("Mean gene body coverage (%)") +
  xlab("Relative 5´to 3´gene body position (%)")

### Gene body coverage - Fig. 4C #### 
cov_sets %>%
  group_by(sample,method) %>% 
  summarise(QCoV = max(QCoV)) %>%
  mutate(mode = str_sub(sample, 1,3),
         sample = factor(sample, levels = rev(bc_to_sample$sample[c(1,10,11,8,9,2,4,3,5,6,7)]))) %>%
  ggplot() +
    geom_bar(aes(y = sample, x = QCoV, group = method, fill = method),position = position_dodge(),
             stat = "identity", color = "black") +
    theme_Publication_white() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(color = "black", linetype = "dashed")) +
    scale_x_continuous(limits = c(0.0,0.3, expand = c(0,0))) +
    scale_fill_manual(values = c("#AD9D86","#A1CEC2")) +
    ylab("") +
    xlab("QCoV")

### Cov5 - Fig. 4D #### 
cov_sets %>%
  dplyr::filter(perc_pos <= 10, method == "trimmed") %>%
  dplyr::group_by(sample) %>%
  summarise(prime5 = mean(C_mean_sum, na.rm = T)/median_C*100) %>%
  distinct(sample, .keep_all = T) %>%
  mutate(mode = str_sub(sample, 1,3),
         sample = factor(sample, levels = rev(bc_to_sample$sample[c(1,10,11,8,9,2,4,3,5,6,7)]))) %>%
  ggplot() +
    geom_bar(aes(y = sample, x = prime5 - 100, fill = mode), 
             stat = "identity", color = "black") +
    scale_x_continuous(limits = c(-35,35), expand = c(0,0)) +
    theme_Publication_white() +
    scale_fill_manual(values = cbf1[c(2,5,3)]) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(color = "black", linetype = "dashed")) +
    ylab("") +
    xlab("CoV5") 

### Cov3 - Fig. 4E #### 
cov_sets %>%
  dplyr::filter(perc_pos <= 90, method == "trimmed") %>%
  dplyr::group_by(sample) %>%
  summarise(prime3 = mean(C_mean_sum, na.rm = T)/median_C*100) %>%
  distinct(sample, .keep_all = T) %>%
  mutate(mode = str_sub(sample, 1,3),
         sample = factor(sample, levels = rev(bc_to_sample$sample[c(1,10,11,8,9,2,4,3,5,6,7)]))) %>%
  ggplot() +
    geom_bar(aes(y = sample, x = prime3 - 100, fill = mode), 
             stat = "identity", color = "black") +
    scale_x_continuous(limits = c(-35,35), expand = c(0,0)) +
    theme_Publication_white() +
    scale_fill_manual(values = cbf1[c(2,5,3)]) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(color = "black", linetype = "dashed")) +
    ylab("") +
    xlab("CoV3")  

### QCoV gene size - Fig. 4F #### 
ggplot(data = cov_trimmed_sizes %>% mutate(mode = str_sub(sample,1,3)) %>% distinct(sample, read_group, QCoV,n, .keep_all =T),
       aes(x = QCoV, y = sample, fill = mode, factor = read_group)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black", 
           aes(alpha = read_group),size = 0.5) +
  scale_alpha_manual(values = rev(c(0,0.33,0.66,1))) +
  scale_x_continuous(limits = c(0.0,2.5),expand = c(0,0)) +
  theme_Publication_white() +
  scale_fill_manual(values = cbf1[c(2,5,3)]) 

