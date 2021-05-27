# >> transcriptional unit analysis << #

# load libraries ----
library(here)
source(here("Rscripts/load_libraries.R"))

# functions & defs ----
p <- function(v) {
  Reduce(f=paste, x = v)
}

perform_analysis <- function(dataset_choice, 
                             read_table = paste0(dir, "/data/mapped_data_pychopper_auto_cutadapt_clipped.tsv")){
  
  # > mod read table
  f <- vroom(read_table, num_threads = 8) %>%
    distinct(minion_read_name,.keep_all = T) %>%
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
  s <- inputDF %>%
    ungroup() %>%
    dplyr::rename(gene = genes_in_operon) %>%
    dplyr::select(seqnames,genes_in_operon_all, size_operon, operon_seen, covered_by, perc_seen, gene, start, end, sample, operon_name) 
    write.table(x = s, file = paste0(dir, "/tables/operon_tables/", s$sample[1], ".operons.tsv"), 
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

get_intersect <- function(inputdf, dataset){
  inputdf %>%
    as_tibble() %>%
    mutate(size_operon = str_count(operon, "\\|") + 1) %>%
    group_by(size_operon) %>%
    summarise(how_many = n()) %>%
    mutate(intersection = dataset) %>%
    ungroup()
}


# load & tidy data ----

## calc transcriptional units for example data sets ====
sets <- c("RNA001_Ecoli_TEX_replicate1",
          "PCB109_PCR12_Ecoli_NOTEX_replicate4",
          "DCS109_Ecoli_NOTEX_replicate2")

RNA_operon <- find_reads_in_operons(dataset_choice = sets[1]) %>% 
  mutate(sample = sets[1])
PCB_operon <- find_reads_in_operons(dataset_choice = sets[2]) %>% 
  mutate(sample = sets[2])
DCS_operon <- find_reads_in_operons(dataset_choice = sets[3]) %>% 
  mutate(sample = sets[3])

### write to Supplementary Table 6 ####
rbindlist(list(RNA_operon,
               PCB_operon,
               DCS_operon)) %>%
  ungroup() %>%
  dplyr::rename(gene = genes_in_operon) %>%
  dplyr::select(seqnames,gene, genes_in_operon_all, size_operon, operon_seen, covered_by, perc_seen, sample) %>%
  write_xlsx(path = here("tables/Supplementary_Table6.xlsx"))


  

## write to output tables ====
write_operon(RNA_operon)
write_operon(PCB_operon)
write_operon(DCS_operon)

## for Upset comparison ONT ====

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

## for Upset comparison ALL ====

### load other data ####

#### regulon db ####
OperonSet <- vroom(paste0(dir ,"data/comparison_data/Operon/OperonSet.txt"), comment = "#", col_names = F) %>%
  separate_rows(X6, sep = ",") %>%
  left_join(ecoli_gff, by = c("X6" = "short_gene")) %>%
  dplyr::filter(!is.na(seqid)) %>% 
  dplyr::rename(operon_name = 1, strand = 4, evidence = 8) %>%
  dplyr::select(operon_name, strand, evidence, gene) %>%
  group_by(operon_name) %>% 
  arrange(gene) %>%
  mutate(operon = str_replace_all(p(gene), " ", "|")) %>%
  mutate(operon = gsub('*^\\|', '', operon)) %>%
  mutate(operon = gsub('*^\\|', '', operon)) %>%
  mutate(operon = gsub('*^\\|', '', operon)) %>%
  dplyr::select(operon, evidence) %>%
  arrange(operon) %>%
  distinct(operon, .keep_all = T) 

#### SMRT CAP operons ####
smrt_operons <- read_xlsx(paste0(dir ,"data/comparison_data/Operon/smrt_cap-41467_2018_5997_MOESM5_ESM.xlsx"), skip = 1) %>% 
  dplyr::rename(operon = 1) %>%
  select(operon, type) %>%
  distinct(operon, .keep_all = T) %>%
  mutate(operon = str_split(operon, '\\|'), 
         operon = map_chr(operon, ~toString(sort(.x))),
         operon = str_replace_all(operon, ",|;", "|"),
         operon = str_replace_all(operon, " ", ""))

### calc intersection ####
inter_DB_NANO    <- intersect(OperonSet$operon, pcb_tu$operon)
inter_DB_SMRT    <- intersect(OperonSet$operon, smrt_operons$operon)
inter_SMRT_NANO  <- intersect(smrt_operons$operon, pcb_tu$operon)
inter_all        <- intersect(smrt_operons$operon, pcb_tu$operon) %>%
  intersect(OperonSet$operon)

### Upset table ####
upset_table_all <- c(`db` = length(smrt_operons$operon) - length(inter_DB_SMRT) + length(inter_all)  - length(inter_DB_NANO),
                 `ONT`= length(pcb_tu$operon) - length(inter_SMRT_NANO) + length(inter_all)  - length(inter_DB_NANO),
                 `smrt` = length(smrt_operons$operon) - length(inter_SMRT_NANO) + length(inter_all)  - length(inter_DB_SMRT),
                 `db&ONT`= length(inter_DB_NANO) - length(inter_all),
                 `db&smrt` = length(inter_DB_SMRT) - length(inter_all),
                 `ONT&smrt` = length(inter_SMRT_NANO) - length(inter_all),
                 `db&ONT&smrt` = length(inter_all))

### Operon size of intersections ####
inter_all_quant <- rbind(get_intersect(pcb_tu %>%
                                         dplyr::filter(!operon %in% OperonSet$operon,
                                                       !operon %in% smrt_operons$operon) %>%
                                         select(operon),
                                       "NANO"),
                         get_intersect(OperonSet %>%
                                         dplyr::filter(!operon %in% pcb_tu$operon,
                                                       !operon %in% smrt_operons$operon) %>%
                                         select(operon),
                                       "DB"),
                         get_intersect(smrt_operons %>%
                                         dplyr::filter(!operon %in% pcb_tu$operon,
                                                       !operon %in% OperonSet$operon) %>%
                                         select(operon),
                                       "SMRT"),
                         get_intersect(pcb_tu %>%
                                         dplyr::filter(operon %in% OperonSet$operon,
                                                       !operon %in% smrt_operons$c) %>%
                                         select(operon),
                                       "DB&NANO"),
                         get_intersect(smrt_operons %>%
                                         dplyr::filter(operon %in% OperonSet$operon,
                                                       !operon %in% pcb_tu$operon) %>%
                                         select(operon),
                                       "DB&SMRT"),
                         get_intersect(smrt_operons %>%
                                         dplyr::filter(!operon %in% OperonSet$operon,
                                                       operon %in% pcb_tu$operon) %>%
                                         select(operon),
                                       "NANO&SMRT"),
                         get_intersect(smrt_operons %>%
                                         dplyr::filter(operon %in% OperonSet$operon,
                                                       operon %in% pcb_tu$operon) %>%
                                         select(operon),
                                       "DB&NANO&SMRT"))


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


### Operon comparison ONT methods/SMRT/DB - Supplementary Fig. 21A #### 
upset(data = fromExpression(upset_table_all), 
      sets = c("db", "ONT", "smrt"),
      keep.order = TRUE, order.by = c("degree"),
      sets.x.label = "Number of TUs", mainbar.y.label = "Intersections", 
      line.size = 3, point.size = 11, 
      sets.bar.color = "gray60",
      matrix.color = "gray60")

### Operon size per intersection - Supplementary Fig. 21B #### 
ggplot(data = inter_all_quant_scaled, aes(x = size_operon, y = how_many, color = intersection, fill = intersection)) +
  geom_line() +
  geom_point(shape = 21, color = "black", size = 4, alpha = 1) +
  scale_y_log10() +
  theme_Publication_white() 

