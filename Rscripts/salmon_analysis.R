# >> Salmon quantification << #

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
get_upper_tri <- function(cormat){
  cormat[upper.tri(cormat)]<- NA
  return(cormat)
}

get_lower_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

modify_salmon_output <- function(input, method){
  suppressMessages(vroom(input, num_threads = 8)) %>%
    arrange(desc(NumReads)) %>%
    mutate(counts = NumReads, gene = Name, salmon_tpm = TPM) %>%
    dplyr::filter(!gene %in% names_rRNA) %>%
    select(counts, gene, salmon_tpm,EffectiveLength) %>%
    mutate(gene = str_split_fixed(gene,"-",2)[,2]) %>%
    arrange(desc(counts)) %>%
    left_join(gff_table, by = "gene") %>%
    mutate(type_fine = ifelse(type == "rRNA", as.character(locus_name), as.character(type)),
           type_fine = ifelse(type == "tRNA", "ncRNA", as.character(type_fine))) %>%
    mutate(method = method,
           rpk = (counts/EffectiveLength*1000),
           TPM_hand = rpk/(sum(rpk, na.rm = T)/1000000),
           rpkm = (counts/(sum(counts)/1000000))/width*1000) %>%
    dplyr::select(gene, type_fine, counts, method,rpkm, TPM_hand, salmon_tpm)
}

read_in_gff <- function(input_file){
  read.gff(input_file) %>%
    dplyr::filter(!type %in% c("exon", "gene", "region", "origin of replication")) %>%
    as_tibble() %>%
    dplyr::mutate(start_feature = start, end_feature = end,strand_feature = strand) %>%
    dplyr::mutate(Parent = str_split_fixed(str_split_fixed(attributes, ";Parent=",2)[,2],";Dbxref",2)[,1],
                  ecogene = str_split_fixed(str_split_fixed(attributes, ",GeneID", 2)[,1], "EcoGene:",2)[,2],
                  short_gene = str_split_fixed(str_split_fixed(attributes, ";locus_tag", 2)[,1], "gene=",2)[,2],
                  id_name = ifelse(type %in% "repeat_region", str_split_fixed(str_split_fixed(attributes, ";Note=", 2)[,1], "ID=", 2)[,2],
                                   ifelse(type %in% "pseudogene", str_split_fixed(str_split_fixed(attributes, ";Dbxref=", 2)[,1], "ID=", 2)[,2],
                                          ifelse(type %in% "sequence_feature", str_split_fixed(str_split_fixed(attributes, ";Dbxref=", 2)[,1], "ID=", 2)[,2],
                                                 ifelse(type %in% "mobile_genetic_element", str_split_fixed(str_split_fixed(attributes, ";gbkey=", 2)[,1], "ID=", 2)[,2],
                                                        str_split_fixed(str_split_fixed(attributes, ";Parent=", 2)[,1], "ID=", 2)[,2])))),
                  locus_name = ifelse(type %in% c("CDS","mobile_genetic_element", "ncRNA", "recombination_feature"), str_split_fixed(str_split_fixed(attributes, ";product=", 2)[,2], ";", 2)[,1],
                                      ifelse(type ==  "pseudogene",  str_split_fixed(str_split_fixed(attributes, ";gene_biotype", 2)[,1], "gene=", 2)[,2],
                                             ifelse(type == "repeat_region", str_split_fixed(str_split_fixed(attributes, ";gbkey", 2)[,1], "Note=", 2)[,2],
                                                    ifelse(type %in% "sequence_feature", str_split_fixed(str_split_fixed(attributes, ";locus_tag=", 2)[,1], "gene=", 2)[,2],
                                                           ifelse(type %in% "mobile_genetic_element", str_split_fixed(attributes, "insertion sequence:", 2)[,2],
                                                                  ifelse(type == "rRNA", str_split_fixed(str_split_fixed(attributes, ";product=", 2)[,2], " ", 2)[,1], 
                                                                         ifelse(type == "tRNA", str_split_fixed(attributes, ";product=", 2)[,2], NA ))))))),
                  width = abs(start_feature - end_feature)) %>%
    dplyr::select(seqid, id_name, locus_name, start_feature, end_feature, strand_feature, Parent, type, width, ecogene, short_gene) %>%
    mutate(gene = str_split_fixed(Parent,"-",2)[,2])
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

corr_matrix_plot <- function(mydf, myx, myy, myfill){
  
  ggplot(data = mydf, aes(x = {{myx}}, y = {{myy}}, fill = {{myfill}}, color = {{myfill}},size = {{myfill}})) +
    geom_tile(color = "grey", size = 0.3, fill = "white") +
    theme_void() +
    scale_fill_gradientn(colours = brewer.pal(name = "YlGnBu", n = 9),
                         limit = c(0.5,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 10, hjust = 1)) +
    coord_fixed() +
    scale_y_discrete(limits=rev)
}

# load & tidy data ----
dir <- "/Volumes/EX_SSD"

## genome data ====
gff_table   <- read_in_gff(paste0(dir, "/data/genome_data/NC_000913.3.gff3"))
ecoli_fasta <- readDNAStringSet(paste0(dir, "/data/genome_data/NC_000913.3.fasta"))
names(ecoli_fasta) <- "chr"

## get names for rRNAs ====
names_rRNA <- gff_table$id_name[gff_table$type == "rRNA"]

## mean gc content per gene table ==== 
gff_table_gc <- gff_table %>%
  dplyr::rowwise() %>%
  dplyr::mutate(gc = GC.content(as.DNAbin(ecoli_fasta$chr[start_feature:end_feature])))

## salmon data ====
### quantification from untrimmed files ####
files          <- list.files(paste0(dir,"/data/salmon_data_notrimming/"), recursive = T,full.names = T, pattern = "quant.sf")
salmon_frame   <- pmap_dfr(list(files,str_split_fixed(str_split_fixed(files, "\\/", n = 9)[,8],"_fu",2)[,1]),modify_salmon_output)

### 1 dataset per column, 1 row per gene ####
salmon_frame_wide <- salmon_frame %>%
  left_join(old_new, by = c("method" = "old_name")) %>%
  dplyr::mutate(method = ifelse(method %in% c("SRR7533627", "SRR7533626", "illumina"), method,new_name)) %>%
  dplyr::select(gene, type_fine, method, TPM_hand) %>%
  dplyr::filter(method %in% c(bc_to_sample$sample, "SRR7533627", "SRR7533626", "illumina")) %>%
  dplyr::filter(TPM_hand > 0 & !is.na(TPM_hand) & is.finite(TPM_hand)) %>%
  dplyr::select(gene, type_fine, method, TPM_hand) %>%
  pivot_wider(names_from = method, values_from = TPM_hand, values_fn = {sum}) %>%
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
  




# > overamplification of short reads
salmon_frame_long <- salmon_frame_wide %>%
  dplyr::select(gene, type_fine, gc, width, all_210317_RNA002_Ecoli_NOTEX_replicate1, all_201208_PCB109_Ecoli_NOTEX_replicate1) %>%
  pivot_longer(cols = all_210317_RNA002_Ecoli_NOTEX_replicate1:all_201208_PCB109_Ecoli_NOTEX_replicate1,names_to = "dataset", values_to = "counts")

ggplot(data = salmon_frame_wide %>% dplyr::filter(type_fine == "CDS"), 
       aes(x = all_210317_RNA002_Ecoli_NOTEX_replicate1, 
           y = all_201208_PCB109_Ecoli_NOTEX_replicate2,
           fill = density)) +
  
  geom_abline(linetype = "dashed", slope = 1) +
  geom_point(alpha = 1, aes(size = width), shape = 21, color = "black") +
  #stat_density2d(aes(alpha=..level.., fill = ..level..),color = NA,
  #               bins=10, geom="polygon") +
  facet_grid(cols = vars(gc > 0.4)) +
  #geom_density2d(color = "black",contour_var = "ndensity", bins = 10, aes(alpha = ..level..)) +
  scale_fill_gradientn(colours = brewer.pal(name = "Blues", n = 9)) +
  theme_Publication_white() +
  #scale_size(range = c(2,10)) +
  scale_x_log10(limits = c(0.1, 1000000), expand = c(0.1,0)) +
  scale_y_log10(limits = c(0.1, 1000000), expand = c(0.1,0)) +
  stat_cor(method = "pearson", label.x = 1) +
  coord_equal() 
dev.off()







#> 
mean(ecoli_gff_gc$gc[ecoli_gff$type == "CDS"])

raw_salmon_frame_other <- raw_salmon_frame %>%
  left_join(ecoli_gff) %>%
  dplyr::filter(!is.na(type_fine), type == "CDS")

salmon_frame_variance <- rbind(raw_salmon_frame,full_salmon_frame) %>%
  dplyr::select(gene, type_fine, method, TPM_hand) %>%
  dplyr::filter(TPM_hand > 0) %>%
  group_by(gene) %>%
  mutate(gene_mean = mean(TPM_hand),
         gene_variance = sd(TPM_hand)) %>%
  left_join(ecoli_gff_gc) %>%
  dplyr::filter(!is.na(type_fine), type == "CDS") %>%
  distinct(gene, .keep_all = T)

## correlation matrix analysis ----------------------------------------


pdf(here("figures/R_plots/210408_salmon_correlation_matrix.pdf"), 
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
max(res_counts_gg$value, na.rm = T)
min(res_counts_gg$value, na.rm = T)

pdf(here("figures/R_plots/210413_salmon_correlation_matrix_N.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = res_counts_gg, aes(Var2, Var1, fill = value, color = value, size = value)) +
  geom_tile(color = "grey", size = 0.3, fill = "white") +
  geom_point(shape = 21, color = "black") +
  #geom_tile(color = "black", size = 0.3, aes(width = value, height = value)) +
  theme_void() +
  scale_fill_gradientn(colours = brewer.pal(name = "Blues", n = 9),
                       limit = c(500,4294), space = "Lab", 
                       name="Pearson\nCorrelation") +
  #geom_text(aes(label=round(value, digits = 2)), color = "white", size = 4) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1)) +
  coord_fixed() +
  scale_y_discrete(limits=rev)

dev.off()



## GC CORRELATION ANALYSIS =======================================================
# > calc correlation for different gc filtering
colnames(salmon_frame_wide)
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



# The palette with grey:
cbf1 <- c("#EFEAFF","#F5AAA3","#CFCFCF", "#F6B2FB", "#ABC2DA")

gc_tab_all <- gc_tab2 %>%
  dplyr::filter(size > 100) %>%
  #mutate(value = cor_all - cor) %>%
  pivot_longer(cor:cor_all,names_to = "cor",values_to = "counts") %>%
  mutate(mode = str_sub(dataset,12,17),
         set = ifelse(cor == "cor_all", "random", "gc_set")) %>%
  dplyr::mutate(size = ifelse(set == "random", NA, size))

pdf(here("figures/R_plots/210412_correlation_GC.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = gc_tab_all, 
       aes(x = gc,  y = counts, color = dataset, linetype = set)) +
  geom_line(size = 2) +
  geom_point(aes(size = size), shape = 21, fill = "white", stroke = 2) +
  geom_vline(xintercept = GC.content(as.DNAbin(ecoli_fasta$`NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome`))) +
  scale_y_continuous(limits = c(0.5,1), expand = c(0,0)) +
  theme_Publication_white() +
  scale_color_manual(values = c(cbf1, "black"))
dev.off()





salmon_frame_wide <- rbind(salmon_frame %>% mutate(way = "notrimming"),
                           full_salmon_frame %>% mutate(way = "pychopper")) %>%
  dplyr::select(gene, type_fine, method, TPM_hand, way) %>%
  dplyr::filter(method %in% c(bc_to_sample$sample[c(1,5:12,15:16)], "SRR7533627", "SRR7533626", "illumina")) %>%
  dplyr::filter(TPM_hand > 0 & !is.na(TPM_hand) & is.finite(TPM_hand)) %>%
  mutate(method = ifelse(way == "notrimming", paste0("all_", method),paste0("trim_", method))) %>%
  dplyr::select(gene, type_fine, method, TPM_hand) %>%
  pivot_wider(names_from = method, values_from = TPM_hand, values_fn = {sum}) %>%
  left_join(gff_table_gc) %>%
  dplyr::filter(!is.na(type_fine), type == "CDS") %>%
  dplyr::filter(!is.na(all_201208_PCB109_Ecoli_NOTEX_replicate2),
                !is.na(all_illumina))

## direct RNA raw vs Illumina
pdf(here("figures/R_plots/210412_salmon_ILL_PCB_GC.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
salmon_frame_wide$density <- get_density(log10(salmon_frame_wide$all_201208_PCB109_Ecoli_NOTEX_replicate2), 
                                         log10(salmon_frame_wide$all_illumina))

ggplot(data = salmon_frame_wide %>% dplyr::filter(gc >= 0.52), 
       aes(x = all_illumina, 
           y = all_201208_PCB109_Ecoli_NOTEX_replicate2,
           fill = density)) +
  geom_abline(linetype = "dashed", slope = 1) +
  geom_point(alpha = 1, aes(size = width), shape = 21, color = "black") +
  #stat_density2d(aes(alpha=..level.., fill = ..level..),color = NA,
  #               bins=10, geom="polygon") +
  facet_grid(cols = vars(gc > 0.4)) +
  #geom_density2d(color = "black",contour_var = "ndensity", bins = 10, aes(alpha = ..level..)) +
  scale_fill_gradientn(colours = brewer.pal(name = "Blues", n = 9)) +
  theme_Publication_white() +
  #scale_size(range = c(2,10)) +
  scale_x_log10(limits = c(0.1, 1000000), expand = c(0.1,0)) +
  scale_y_log10(limits = c(0.1, 1000000), expand = c(0.1,0)) +
  stat_cor(method = "pearson", label.x = 1) +
  coord_equal() 
dev.off()


# width effect --------------------------------
width_interesting <- seq(from = 0,to = 1000, by = 100)
set_interesting <- c(5,8,10)

width_tab1 <- data.table()
width_tab2 <- data.table()
for(i in seq_along(width_interesting)){
  for(j in seq_along(set_interesting)){
    
    salmon_frame_test <- salmon_frame_wide %>%
      dplyr::select(width, colnames(salmon_frame_wide)[c(set_interesting[j],15)]) %>%
      remove_missing() %>%
      dplyr::filter(width <= width_interesting[i])
    
    salmon_frame_test_all <- salmon_frame_wide %>%
      dplyr::select(colnames(salmon_frame_wide)[c(set_interesting[j],15)]) %>%
      remove_missing() %>%
      sample_n(nrow(salmon_frame_test))
    
    if(nrow(salmon_frame_test) > 0){
      width_tab1 <- data.table(width = width_interesting[i],
                               size = nrow(salmon_frame_test),
                               cor  = cor(log10(salmon_frame_test[c(2,3)]),  method = "pearson", use = "pairwise.complete.obs")[1,2],
                               cor_all = cor(log10(salmon_frame_test_all[c(1,2)]),  method = "pearson", use = "pairwise.complete.obs")[1,2],
                               dataset = colnames(salmon_frame_wide)[set_interesting[j]])
    }
    width_tab2 <- rbind(width_tab2, width_tab1)
  }
  
}



# The palette with grey:
cbf1 <- c("#EFEAFF","#F5AAA3","#CFCFCF", "#F6B2FB", "#ABC2DA")

width_tab_all <- width_tab2 %>%
  dplyr::filter(size > 100) %>%
  #mutate(value = cor_all - cor) %>%
  pivot_longer(cor:cor_all,names_to = "cor",values_to = "counts") %>%
  mutate(mode = str_sub(dataset,12,17),
         set = ifelse(cor == "cor_all", "random", "width_set")) %>%
  dplyr::mutate(size = ifelse(set == "random", NA, size))

pdf(here("figures/R_plots/210412_correlation_width.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = width_tab_all, 
       aes(x = width,  y = counts, color = dataset, linetype = set)) +
  geom_line(size = 2) +
  geom_point(aes(size = size), shape = 21, fill = "white", stroke = 2) +
  scale_y_continuous(limits = c(0.5,1), expand = c(0,0)) +
  theme_Publication_white() +
  scale_color_manual(values = c(cbf1, "black"))
dev.off()

