# >> transcriptional unit analysis << #

# load libraries ----
library(here)
source(here("Rscripts/load_libraries.R"))

# functions & defs ----
perform_analysis <- function(dataset_choice, 
                             tss_input = paste0(dir,"/data/tss_data/tss_data_notrimming.tsv"), 
                             tts_input = paste0(dir,"/data/tts_data/tts_data_pychopper_auto_cutadapt_SSP_clipped.tsv"), 
                             read_table = paste0(dir, "/data/mapped_data_pychopper_auto_cutadapt_clipped.tsv"), 
                             merge_w = 20,
                             cov_min = 3){
  
  # > get TSS
  tss_p <- vroom(file = tss_input, num_threads = 8, progress = F) %>%
    mutate(mode = str_sub(sample, 1,3),
           TSS = ifelse(mode == "RNA" & strand == "+", TSS - 12, 
                        ifelse(mode == "RNA" & strand == "-", TSS + 12,TSS))) %>%
    dplyr::filter(TSS_type == "primary", type == "CDS") %>%
    dplyr::filter(sample %in% dataset_choice) %>%
    distinct(TSS, .keep_all = T) %>%
    arrange(TSS) %>%
    mutate(index = lag(TSS, default = 1) + as.integer(merge_w),
           index1 = cumsum(ifelse(index >= TSS, 0, 1))+1) %>%
    dplyr::group_by(index1) %>%
    dplyr::filter(cov == max(cov),
                  cov >= cov_min) %>%
    ungroup() %>%
    group_by(gene) %>%
    dplyr::slice(which.max(cov)) %>%
    ungroup() %>%
    dplyr::select(id_name, TSS, TSS_type)
  
  # > get TTS
  tts_p <- vroom(file = tts_input, num_threads = 8, progress = F) %>%
      dplyr::left_join(old_new, by = c("dataset" = "old_name")) %>%
      mutate(sample = new_name) %>%
      dplyr::select(-new_name, -dataset) %>%
      dplyr::filter(!is.na(sample)) %>%
    mutate(mode = str_sub(sample, 1,3)) %>%
    dplyr::filter(TTS_type == "primary", type == "CDS") %>%
    dplyr::filter(sample %in% dataset_choice) %>%
    distinct(TTS, .keep_all = T) %>%
    arrange(TTS) %>%
    mutate(index = lag(TTS, default = 1) + as.integer(merge_w),
           index1 = cumsum(ifelse(index >= TTS, 0, 1))+1) %>%
    dplyr::group_by(index1) %>%
    dplyr::filter(cov == max(cov),
                  cov >= cov_min) %>%
    ungroup() %>%
    group_by(gene) %>%
    dplyr::slice(which.max(cov)) %>%
    ungroup() %>%
    dplyr::select(id_name, TTS, TTS_type)
  
  # > mod read table
  working_table <- vroom(read_table, num_threads = 8) %>%
    dplyr::left_join(old_new, by = c("method" = "old_name")) %>%
    dplyr::rename(sample = new_name) %>%
    dplyr::select(-method) %>%
    dplyr::filter(!is.na(sample)) %>%
    dplyr::filter(sample %in% dataset_choice) %>%
    distinct(minion_read_name,.keep_all = T) %>%
    left_join(tss_p, by = "id_name") %>%
    left_join(tts_p, by = "id_name") 
  
  f <- working_table %>%
    dplyr::filter(soft_l < 10 & soft_r < 10 & soft_r < 10 & hard_r < 10 & width.x < (aligned_reads * 1.5)) %>%
    distinct(start, end, gene, .keep_all = T) %>%
    mutate(sample = dataset_choice)
  return(f)
}

find_reads_in_operons <- function(dataset_choice, gff_table = ecoli_gff_cds){
  
  map_table_f <- perform_analysis(dataset_choice) %>%
    dplyr::filter(strand == strand_feature) %>%
    arrange(desc(end)) %>%
    distinct(minion_read_name, .keep_all = T) %>%
    rownames_to_column("id") %>%
    mutate(id = as.numeric(id))
  
  seen_at_least <-  nrow(map_table_f)/21923.9
  
  # > reads 
  x1 <- map_table_f %>% dplyr::select(seqnames,start,end, minion_read_name) %>% distinct(minion_read_name, .keep_all = T) 
  gr1 <- makeGRangesFromDataFrame(x1)
  
  # > CDS & ncRNA
  x2 <- ecoli_gff_cds %>% dplyr::rename(seqnames = 2, start = 5, end = 6) %>% dplyr::select(seqnames,start,end, gene)
  gr2 <- makeGRangesFromDataFrame(x2)  
  
  # > within range --> overlapping by at least 100 bases
  type3 <- findOverlaps(query = gr2, subject = gr1, type = "any", minoverlap = 100)
  type3_df <- data.frame(x1[subjectHits(type3),],x2[queryHits(type3),])
  
  # > each read is at the beginning a single operon
  read_name_to_operon <- data.table(minion_read_name = x1$minion_read_name, operon_name = 1:nrow(x1))
  
  # > only a percentage of those reads covers a read fully, calc number of total reads, and how often each gene is covered by a read
  type3_df_operon <- type3_df %>%
    mutate(total_reads = length(unique(minion_read_name))) %>%
    distinct(minion_read_name, gene, .keep_all = T) %>%
    group_by(gene) %>%
    mutate(covered_by = n()) %>%
    left_join(read_name_to_operon, by = "minion_read_name")
  
  # > operons are named by collapsing all the genes in the list
  type3_df_operon_names <- type3_df_operon %>%
    distinct(minion_read_name, operon_name, .keep_all = T) %>%
    group_by(operon_name) %>%
    mutate(genes_in_operon = paste(gene, collapse = ",")) 
  
  # > calculate how many times an operon is "seen" by a single read, keep in mind that every read name is a single operon name at this point
  type3_df_operon_names_D <- type3_df_operon_names %>%
    distinct(genes_in_operon, operon_name, .keep_all = F) %>%
    group_by(genes_in_operon) %>%
    mutate(operon_seen = n()) %>%
    dplyr::filter(operon_seen >= seen_at_least)
  
  # > calc number of reads per gene
  cov_frame <- type3_df_operon %>%
    dplyr::filter(operon_name %in% type3_df_operon_names_D$operon_name) %>%
    distinct(minion_read_name, operon_name, .keep_all = T) %>%
    distinct(gene, covered_by, .keep_all = F)
  
  # > start with reads as rows, genes containing and how often seen, seperate rows by "," but contain operon
  type3_df_operon_names_D_collapse <- type3_df_operon_names_D %>% 
    distinct(genes_in_operon, .keep_all = T) %>%
    dplyr::mutate(size_operon = 1 + str_count(genes_in_operon, ",")) %>%
    dplyr::select(genes_in_operon, size_operon, operon_seen) %>%
    rownames_to_column("operon_name") %>%
    separate_rows(genes_in_operon,sep = ",") %>%
    left_join(x2, by = c("genes_in_operon" = "gene")) %>%
    left_join(cov_frame, by = c("genes_in_operon" = "gene")) %>%
    arrange(genes_in_operon) %>%
    group_by(operon_name) %>%
    mutate(genes_in_operon_all = paste(genes_in_operon, collapse = ",")) %>%
    mutate(perc_seen = operon_seen/covered_by*100) 
  
  return(type3_df_operon_names_D_collapse)
}

write_operon <- function(inputDF){
  inputDF %>%
    ungroup() %>%
    dplyr::rename(gene = genes_in_operon) %>%
    dplyr::select(seqnames,genes_in_operon_all, size_operon, operon_seen, covered_by, perc_seen, gene) %>%
    write.table(x = ., file = paste0(dir, "/tables/operon_tables/", sample, ".operons.tsv"), 
                sep = "\t", quote = F,col.names = T, row.names = F)
}

get_context <- function(input, dataset_name){
  input %>%
    dplyr::rename(operon = genes_in_operon_all) %>%
    distinct(operon) %>%
    arrange(operon) %>%
    mutate(c = str_split(operon, '\\|'), 
           c = map_chr(c, ~toString(sort(.x))),
           c = str_replace_all(c, ",|;", "|"),
           c = str_replace_all(c, " ", "")) %>%
    dplyr::select(c) %>%
    separate_rows(c,sep = "\\|") %>%
    group_by(c) %>%
    summarise(n = n()) %>%
    #group_by(n) %>%
    #summarise(total = n()) %>%
    mutate(n = ifelse(n >=6, 6, n)) %>%
    group_by(n) %>%
    summarise(total = n()) %>%
    ungroup() %>%
    mutate(perc = total/sum(total)*100) %>%
    mutate(sample = dataset_name)
}

modify_TU_output <- function(input_set){
  input_set %>%
    distinct(genes_in_operon_all, .keep_all =T) %>%
    arrange(desc(size_operon)) %>%
    
    dplyr:::select(size_operon,genes_in_operon_all, operon_seen) %>%
    
    mutate(operon = str_replace_all(genes_in_operon_all, ",", "\\|")) %>%
    arrange(operon) %>%
    arrange(desc(size_operon))
}

# load & tidy data ----

## calc transcriptional units for example data sets ====
sets <- c("RNA001_Ecoli_TEX_replicate1",
          "PCB109_PCR12_Ecoli_TEX_replicate4",
          "DCS109_Ecoli_NOTEX_replicate2")

RNA_operon <- find_reads_in_operons(dataset_choice = sets[1]) %>% 
  mutate(sample = sets[1])
PCB_operon <- find_reads_in_operons(dataset_choice = sets[2]) %>% 
  mutate(sample = sets[2])
DCS_operon <- find_reads_in_operons(dataset_choice = sets[3]) %>% 
  mutate(sample = sets[3])

## write to output tables ====
write_operon(RNA_operon)
write_operon(PCB_operon)
write_operon(DCS_operon)

## for Upset comparison ====

### prepare data ####
pcb_tu <- modify_TU_output(PCB_operon)
dcs_tu <- modify_TU_output(DCS_operon)
rna_tu <- modify_TU_output(RNA_operon)

### calc intersections ####
inter_DCS_PCB    <- intersect(dcs_tu$operon, pcb_tu$operon)
inter_DCS_RNA    <- intersect(dcs_tu$operon, rna_tu$operon)
inter_RNA_PCB    <- intersect(rna_tu$operon, pcb_tu$operon)
inter_all        <- intersect(rna_tu$operon, pcb_tu$operon) %>%
  intersect(dcs_tu$operon)

### Upset table ####
upset_table <- c(`DCS` = length(dcs_tu$operon) - length(inter_DCS_RNA) + length(inter_all)  - length(inter_DCS_PCB),
                 `PCB`= length(pcb_tu$operon) - length(inter_RNA_PCB) + length(inter_all)  - length(inter_DCS_PCB),
                 `RNA` = length(rna_tu$operon) - length(inter_RNA_PCB) + length(inter_all)  - length(inter_DCS_RNA),
                 `DCS&PCB`= length(inter_DCS_PCB) - length(inter_all),
                 `DCS&RNA` = length(inter_DCS_RNA) - length(inter_all),
                 `ONT&RNA` = length(inter_RNA_PCB) - length(inter_all),
                 `DCS&PCB&RNA` = length(inter_all))


## calc transcriptional unit context ====
ONT_context_all <- rbindlist(list(get_context(PCB_operon, "PCB"),
                                  get_context(DCS_operon, "DCS"),
                                  get_context(RNA_operon, "RNA")))

# PLOTS ----

## reorder levels ====
ONT_context_all$sample <- factor(ONT_context_all$sample, 
                                 levels=c("PCB", "DCS", "RNA"))

cov_sets$method <- factor(cov_sets$method,
                          levels = rev(c("untrimmed", "trimmed")))

## plotting ==== 

### Operon comparison ONT methods - Fig. 5B #### 
upset(data = fromExpression(upset_table), 
      sets = c("DCS", "PCB", "RNA"),
      keep.order = TRUE, order.by = c("degree"),
      sets.x.label = "Number of TUs", mainbar.y.label = "Intersections", 
      line.size = 3, point.size = 11, 
      sets.bar.color = "gray60",
      matrix.color = "gray60")

### Transcriptional contexts of a gene per method - Fig. 5C #### 
ggplot(data = ONT_context_all, aes(y = sample, x = perc, fill = as.factor(n))) +
  geom_bar(stat = "identity", color = "black") +
  theme_Publication_white() +
  scale_fill_manual(values = cbf1_high2) +
  coord_cartesian(xlim = c(0,100), expand = c(0,0)) +
  ylab("") +
  xlab("Distribution of transcriptional contexts") +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "black", linetype = "dashed")) 
