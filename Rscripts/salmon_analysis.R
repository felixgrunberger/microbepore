# >> Salmon quantification << #

# load libraries ----
library(here)
source(here("Rscripts/load_libraries.R"))

# functions & defs ----

modify_salmon_output <- function(input, method){
  suppressMessages(vroom(input, num_threads = 8)) %>%
    arrange(desc(NumReads)) %>%
    mutate(counts = NumReads, gene = Name, salmon_tpm = TPM) %>%
    dplyr::filter(!gene %in% names_rRNA) %>%
    select(counts, gene, salmon_tpm,EffectiveLength) %>%
    mutate(gene = str_split_fixed(gene,"-",2)[,2]) %>%
    arrange(desc(counts)) %>%
    left_join(ecoli_gff, by = "gene") %>%
    mutate(type_fine = ifelse(type == "rRNA", as.character(locus_name), as.character(type)),
           type_fine = ifelse(type == "tRNA", "ncRNA", as.character(type_fine))) %>%
    mutate(sample = method,
           rpk = (counts/EffectiveLength*1000),
           TPM_hand = rpk/(sum(rpk, na.rm = T)/1000000),
           rpkm = (counts/(sum(counts)/1000000))/width*1000) %>%
    dplyr::select(gene, type_fine, counts, sample,rpkm, TPM_hand, salmon_tpm)
}

point_cor <-  function(mydf, myx, myy, myfill, mysize){
  df <- {{mydf}} %>%
    dplyr::filter(!is.na({{myx}}),
                  !is.na({{myy}})) %>%
    mutate(density = get_density(log10({{myx}}), log10({{myy}}))) %>%
    dplyr::filter(type_fine == "CDS")
  
  ggplot(data = df, aes(x = {{myx}}, y = {{myy}}, fill = {{myfill}}, size = {{mysize}})) +
    geom_abline(linetype = "dashed", slope = 1) +
    geom_point(alpha = 1, shape = 21, color = "black") +
    scale_fill_gradientn(colours = brewer.pal(name = "Blues", n = 9)) +
    theme_Publication_white() +
    scale_x_log10(limits = c(0.1, 1000000), expand = c(0.1,0)) +
    scale_y_log10(limits = c(0.1, 1000000), expand = c(0.1,0)) +
    stat_cor(method = "pearson", label.x = 1) +
    coord_equal() 
}



# load & tidy data ----
dir <- here()

## get names for rRNAs ====
names_rRNA <- ecoli_gff$id_name[ecoli_gff$type == "rRNA"]

## mean gc content per gene table ==== 
gff_table_gc <- ecoli_gff %>%
  dplyr::rowwise() %>%
  dplyr::mutate(gc = GC.content(as.DNAbin(ecoli_fasta$chr[start_feature:end_feature])))

## salmon data ====
### quantification from untrimmed files ####
files          <- list.files(paste0(dir,"/data/salmon_data_notrimming/"), recursive = T,full.names = T, pattern = "quant.sf")
salmon_frame   <- pmap_dfr(list(files,str_split_fixed(str_split_fixed(files, "\\/", n = 9)[,8],"_fu",2)[,1]),modify_salmon_output)

### write quantification data ####
fwrite(salmon_frame, paste0(dir,"/tables/salmon_table.tsv"), col.names = T, row.names = F, quote = F, sep = "\t")
# salmon_frame <- vroom(paste0(dir, "/tables/salmon_table.tsv"))

### write to Supplementary Table 3 ####
salmon_Table %>%
  dplyr::rename(TPM = TPM_hand,
                sample = method) %>%
  dplyr::select(gene, TPM, sample) %>%
  dplyr::mutate(sample = ifelse(sample == "illumina", "SRR1927169", sample)) %>%
  dplyr::filter(!is.na(sample)) %>%
  distinct(gene, sample, .keep_all = T)  %>%
  pivot_wider(names_from = sample, values_from = TPM) %>%
  arrange(gene) %>%
  write_xlsx(path = here("tables/Supplementary_Table3.xlsx"))

### 1 dataset per column, 1 row per gene ####
salmon_frame_wide <- salmon_frame %>%
  dplyr::select(gene, type_fine, sample, TPM_hand) %>%
  dplyr::filter(sample %in% c(bc_to_sample$sample, "SRR7533627", "SRR7533626", "illumina")) %>%
  dplyr::filter(TPM_hand > 0 & !is.na(TPM_hand) & is.finite(TPM_hand)) %>%
  dplyr::select(gene, type_fine, sample, TPM_hand) %>%
  pivot_wider(names_from = sample, values_from = TPM_hand, values_fn = {sum}) %>%
  left_join(gff_table_gc) %>%
  dplyr::filter(!is.na(type_fine), type == "CDS") 

### correlation matrix ####
#### pairwise complete Pearson correlation ####
res             <- cor(log10(salmon_frame_wide[c(3,12,13,10,11,4,6,5,7:9, 14:16)]),  
                       method = "pearson", use = "pairwise.complete.obs")
res_gg          <- reshape2::melt(get_upper_tri(res))

#### pairwise complete Pearson observations ####
res_counts      <- pairwiseCount(salmon_frame_wide[c(3,12,13,10,11,4,6,5,7:9, 14:16)], diagonal = F)
res_counts_gg   <- reshape2::melt(get_lower_tri(res_counts))

### GC CORRELATION ANALYSIS ####
# > calc correlation for different gc filtering
set.seed(1)
gc_interesting <- seq(from = 0.3,to = 0.6, by = 0.02)
set_interesting <- c(5,8,10)

gc_tab1 <- data.table()
gc_tab2 <- data.table()
for(i in seq_along(gc_interesting)){
  for(j in seq_along(set_interesting)){
    
    salmon_frame_test <- salmon_frame_wide %>%
      dplyr::select(gc, colnames(salmon_frame_wide)[c(set_interesting[j],14)]) %>%
      remove_missing() %>%
      dplyr::filter(gc >= gc_interesting[i])
    
    salmon_frame_test_all <- salmon_frame_wide %>%
      dplyr::select(colnames(salmon_frame_wide)[c(set_interesting[j],14)]) %>%
      remove_missing() %>%
      sample_n(nrow(salmon_frame_test))
    
    if(nrow(salmon_frame_test) > 0){
      gc_tab1 <- data.table(gc = gc_interesting[i],
                            size = nrow(salmon_frame_test),
                            cor  = cor(log10(salmon_frame_test[c(2,3)]),  method = "pearson", use = "pairwise.complete.obs")[1,2],
                            cor_all = cor(log10(salmon_frame_test_all[c(1,2)]),  method = "pearson", use = "pairwise.complete.obs")[1,2],
                            dataset = colnames(salmon_frame_wide)[set_interesting[j]])
    }
    gc_tab2 <- rbind(gc_tab2, gc_tab1)
  }
}

# > make table
gc_tab_all <- gc_tab2 %>%
  dplyr::filter(size > 100) %>%
  #mutate(value = cor_all - cor) %>%
  pivot_longer(cor:cor_all,names_to = "cor",values_to = "counts") %>%
  mutate(mode = str_sub(dataset,12,17),
         set = ifelse(cor == "cor_all", "random", "gc_set")) %>%
  dplyr::mutate(size = ifelse(set == "random", NA, size))

# PLOTS ----

## Gene expression correlation - Fig. 2A ====
point_cor(salmon_frame_wide, 
          DCS109_Ecoli_NOTEX_replicate2, 
          RNA002_Ecoli_NOTEX_replicate2,
          density,
          width) 

## Gene expression correlation - Fig. 2B ====
point_cor(salmon_frame_wide, 
          DCS109_Ecoli_NOTEX_replicate2, 
          illumina,
          density,
          width) 

## Correlation matrix - Fig. 2C ====

### Part 1 ####
corr_matrix_plot(res_gg, Var2, Var1, value) +
  geom_tile(color = "black", size = 0.3, aes(width = value, height = value)) +
  geom_text(aes(label=round(value, digits = 2)), color = "white", size = 4) 

### Part 2 ####  
corr_matrix_plot(res_counts_gg, Var2, Var1, value) +
  geom_point(shape = 21, color = "black") +
  scale_fill_gradientn(colours = brewer.pal(name = "Blues", n = 9),
                       limit = c(500,4294), space = "Lab", 
                       name="Pearson\nCorrelation")
  
## Salmon quantification/GC - Supplementary Fig. 11A ===
point_cor(salmon_frame_wide %>% dplyr::filter(gc >= 0.52), 
          PCB109_PCR15_Ecoli_TEX_replicate5, 
          illumina,
          density,
          width) 

## Pearson correlation per GC - Supplementary Fig. 11B ===
ggplot(data = gc_tab_all, 
       aes(x = gc,  y = counts, color = dataset, linetype = set)) +
  geom_line(size = 2) +
  geom_point(aes(size = size), shape = 21, fill = "white", stroke = 2) +
  geom_vline(xintercept = GC.content(as.DNAbin(ecoli_fasta$chr))) +
  scale_y_continuous(limits = c(0.5,1), expand = c(0,0)) +
  theme_Publication_white() +
  scale_color_manual(values = c(cbf1, "black"))
