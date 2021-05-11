# >> Mapped read analysis << #

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

## read in gff file ====
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

## read in BAM file and deal with multi-mapping reads ====
read_bam_files <- function(inputBAM, method){
  
  # read in files
  init <- readGAlignments(inputBAM, use.names = T, param = ScanBamParam(tag=c("NM"), what="mapq"))
  init_t <- GenomicAlignments::as.data.frame(init) %>%
    dplyr::mutate(minion_read_name = names(init),
                  mapped_gene = seqnames) 
  
  left  <- paste(str_split_fixed(string = init_t$cigar, pattern = "M", n = 2)[,1],"M", sep = "")
  right <- paste(str_split_fixed(string = init_t$cigar, pattern = "M", n = 2)[,2],"1M", sep = "")
  
  #................................calculate cigar tables / SOFT AND HARD CLIPPING!!!
  init_t$soft_l <- as_tibble(cigarOpTable(left))$S
  init_t$hard_l <- as_tibble(cigarOpTable(left))$H
  init_t$soft_r <- as_tibble(cigarOpTable(right))$S
  init_t$hard_r <- as_tibble(cigarOpTable(right))$H
  
  # calculate number of aligned reads based on CIGAR operations (M,I)
  init_t$aligned_reads <- unlist(lapply(explodeCigarOpLengths(init_t$cigar, ops = c("M", "I")), function(x) sum(x)))
  
  # calc read identity..
  init_t_final <- init_t %>%
    dplyr::mutate(identity = (1 - NM/aligned_reads)*100) %>%
    dplyr::group_by(minion_read_name) %>%
    dplyr::filter(identity == max(identity),
                  aligned_reads == max(aligned_reads)) %>%
    dplyr::distinct(minion_read_name, .keep_all = T) %>%
    dplyr::mutate(sample = method,
                  gene = str_split_fixed(mapped_gene,"-",2)[,2])
  
  # return table
  return(init_t_final)
}


## add gene mapping information from remapped genes ====
annotate_bams <- function(input_mapped, input_remapped, dataset){
  mapped_t   <- read_bam_files(input_mapped, dataset)
  remapped_t <- read_bam_files(input_remapped, dataset)
  mapped_t_a <- mapped_t %>%
    dplyr::select(-mapped_gene, -gene) %>%
    left_join(remapped_t %>%
                dplyr::select(minion_read_name, mapped_gene, gene), 
              by = "minion_read_name") %>%
    left_join(gff_table, by = "gene")
}

## calculate number of mapped reads per gene and join with gff file ====
mutate_bam_files <- function(inputBAMtable){
  inputBAMtable %>%
    group_by(gene, sample) %>%
    dplyr::mutate(counts = n(),
                  reads = sum(aligned_reads)) %>%
    dplyr::select(gene, counts, reads, sample) %>%
    distinct(gene, counts, reads,sample, .keep_all = T) %>%
    arrange(desc(counts)) %>%
    left_join(gff_table, by = "gene") %>%
    ungroup()
}

## calculate total counts per feature ====
calc_bam_counts <- function(inputBAMtable, mode){
  inputBAMtable %>%
    group_by(sample) %>%
    mutate(total_counts = sum(counts),
           total_reads = sum(reads)) %>%
    dplyr::mutate(type = ifelse(type == "rRNA" & mode == "rRNA", locus_name, as.character(type)),
                  type = ifelse(type == "tRNA", "ncRNA", as.character(type))) %>%
    group_by(type, sample) %>%
    summarise(number_of_counts = sum(counts),
              number_of_reads = sum(reads),
              number_of_counts_p = number_of_counts/total_counts*100,
              number_of_reads_p = number_of_reads/total_reads*100) %>%
    distinct(type, number_of_counts, number_of_reads,sample, .keep_all = T) %>%
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

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
})

# data ----
dir <- "/Volumes/EX_SSD/"

## genome data ====
gff_table      <-  read_in_gff(paste0(dir, "data/genome_data/NC_000913.3.gff3"))

## mapped files | simply from unfiltered fastq to minimap2/remap ====
files          <- list.files(paste0(dir,"data/mapped_data_notrimming"), recursive = T, full.names = T,pattern = ".sorted.bam$")
mapped_frame   <- pmap_dfr(list(files[!str_detect(files, "remapped")],
                                files[str_detect(files, "remapped")],
                                str_split_fixed(str_split_fixed(files, "\\/", n = 8)[,7],"_fu",2)[,1]),
                          annotate_bams)

## save data frame ====
fwrite(mapped_frame, paste0(dir, "/data/mapped_data_no_trimming.tsv"), sep = "\t",col.names = T, nThread = 8)
mapped_frame <- vroom(paste0(dir, "/data/mapped_data_no_trimming.tsv"))

## combine with sequencing summary file to identify unmapped reads ====
summary_frame_sample  <- vroom(paste0(dir, "/data/summary_data_overview.tsv")) %>%
  group_by(sample) %>%
  summarise(all_reads = n(),
            all_bases = sum(sequence_length_template)) %>%
  mutate(mode = substr(sample, 8,10)) 

## calc total number of mapped reads ==== 
merged_frame <- left_join(mapped_frame %>%
  distinct(minion_read_name, .keep_all = T) %>%
  group_by(sample) %>%
  summarise(mapped_reads = n_distinct(minion_read_name),
            mapped_bases= sum(aligned_reads)),
  summary_frame_sample, by = "sample") %>%
  mutate(percentage_reads = mapped_reads/all_reads*100,
         percentage_bases = mapped_bases/all_bases*100) %>%
  mutate(mode = substr(sample, 1,3))

# calculate stats ----
summary_stats <- mapped_frame %>%
  group_by(sample) %>%
  summarise(number_of_mapped_reads = n(),
            number_of_mapped_bases = sum(aligned_reads),
            median_mapped_length   = median(aligned_reads),
            median_mapped_identity = round(median(identity), digits = 2)) %>%
  mutate(mode = substr(sample, 1,3)) %>%
  arrange(factor(sample, levels = bc_to_sample$sample[c(1,10,11,8,9,2,4,3,5,6,7)]))


# PLOTS ----

## reorder levels ====
merged_frame$sample <- factor(merged_frame$sample,
                                      levels = rev(bc_to_sample$sample[c(1,10,11,8,9,2,4,3,5,6,7)]))
merged_frame$mode <- factor(merged_frame$mode,
                                    levels = c("RNA", "DCS", "PCB"))

## color palette ====
cbf1 <- c("#EFEAFF","#F5AAA3","#CFCFCF", "#F6B2FB", "#ABC2DA")

## plotting ==== 

### Proportion of mapped reads - Supplementary Fig. 5A ####
raw_reads_plotting(merged_frame, percentage_reads, sample, mode, cbf1[c(3,2,5)]) +
  geom_bar(stat = "identity", color = "black") +
  scale_x_continuous(expand = c(0,0), limits = c(0,100)) +
  xlab("Proportion of mapped reads (%)") 

### Proportion of mapped bases - Supplementary Fig. 5B ####
raw_reads_plotting(merged_frame, percentage_bases, sample, mode, cbf1[c(3,2,5)]) +
  geom_bar(stat = "identity", color = "black") +
  scale_x_continuous(expand = c(0,0), limits = c(0,100)) +
  xlab("Proportion of mapped reads (%)") 


# Feature-dependent analysis of mapped reads ----
## Summary data ====
summary_total <- vroom(paste0(dir, "/data/summary_data_overview.tsv"), num_threads = 8) %>%
  left_join(mapped_frame %>% dplyr::select(minion_read_name, type, aligned_reads,identity), 
            by = c("read_id" = "minion_read_name")) %>%
  dplyr::mutate(type = ifelse(is.na(type), "unmapped", 
                              ifelse(type == "CDS", "mRNA", 
                                     ifelse(type == "rRNA", "rRNA","other_ncRNA"))),
                group = ifelse(type == "unmapped", "Unmapped","Mapped"))

# PLOTS ----

## reorder levels ====
summary_total$sample <- factor(summary_total$sample,
                                      levels = rev(bc_to_sample$sample[c(1,10,11,8,9,2,4,3,5,6,7)]))

summary_total$type <- factor(summary_total$type,
                                    levels = c("mRNA", "rRNA", "other_ncRNA", "unmapped"))

summary_total$mode <- factor(summary_total$mode,
                             levels = rev(c("RNA", "DCS", "PCB")))

summary_total$group <- factor(summary_total$group,
                             levels = rev(c("Mapped", "Unmapped")))

### Read length distribution of raw reads per mapping class - Supplementary Fig. 6A ####
raw_reads_plotting(summary_total, sequence_length_template, sample, mode, cbf1[c(2,5,3)]) +
  facet_grid(cols = vars(type)) +
  geom_density_ridges(aes(height =..ndensity..), scale = 0.9, color = "black") +
  scale_x_continuous(limits = c(0,3500), expand = c(0,0)) +
  scale_y_discrete(expand = c(0.01,0.1)) +
  xlab("Read length (bases)") 

### Read quality distribution of raw reads per mapping class - Supplementary Fig. 6B ####
raw_reads_plotting(summary_total, mean_qscore_template, sample, mode, cbf1[c(2,5,3)]) +
  facet_grid(cols = vars(type)) +
  geom_density_ridges(aes(height =..ndensity..), scale = 0.9, color = "black") +
  scale_x_continuous(limits = c(0,20), expand = c(0,0)) +
  scale_y_discrete(expand = c(0.01,0.1)) +
  xlab("Read length (bases)") 

### Read length Mapped/Unmapped per mode - Supplementary Fig. 6C ####
ggplot(data = summary_total %>% sample_n(100000), 
       aes(x = mode, y = sequence_length_template, 
           factor = group, alpha = group)) +
  geom_split_violin(trim = T, scale = "width", aes(fill = mode)) +
  coord_flip() +
  theme_Publication_white() +
  scale_y_continuous(limits = c(0,3500), expand = c(0,0)) +
  scale_fill_manual(values = cbf1[c(5,2,3)]) +
  ylab("") +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linetype = "dashed", color = "black")) 

### Read quality Mapped/Unmapped per mode - Supplementary Fig. 6D ####
ggplot(data = summary_total %>% sample_n(100000), 
       aes(x = mode, y = mean_qscore_template, 
           factor = group, alpha = group)) +
  geom_split_violin(trim = T, scale = "width", aes(fill = mode)) +
  coord_flip() +
  theme_Publication_white() +
  scale_y_continuous(limits = c(0,20), expand = c(0,0)) +
  scale_fill_manual(values = cbf1[c(5,2,3)]) +
  ylab("") +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linetype = "dashed", color = "black")) 




# Aligned read length ----------------------------------------------------------
total_frame_untrimmed_f <- total_frame_untrimmed %>%
  dplyr::rename(sample = method) %>%
  dplyr::filter(sample %in% bc_to_sample$sample[c(1,5:12,15:16)]) %>%
  distinct(minion_read_name, .keep_all = T)

total_frame_untrimmed_f_stats <- total_frame_untrimmed_f %>%
  mutate(mode = substr(sample, 8,10),
         type = replace_na(type, "unspecified"))

total_frame_untrimmed_f_stats$sample <- factor(total_frame_untrimmed_f_stats$sample,
                                               levels = rev(bc_to_sample$sample[c(1,15,16,11,12,5,7,6,8:10)]))

total_frame_untrimmed_f_stats$type <- factor(total_frame_untrimmed_f_stats$type,
                                             levels = rev(c("unspecified", "ncRNA", "rRNA", "CDS")))

# plot ----
 
## set color palette ==== 
cbf1 <- c("#EFEAFF","#F5AAA3","#CFCFCF", "#F6B2FB", "#ABC2DA")

### Proportion of mapped reads - Supplementary Fig. 5A ####

### Proportion of mapped bases - Supplementary Fig. 5B ####

pdf(here("figures/R_plots/210407_aligned_read_length.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = total_frame_untrimmed_f_stats, 
       aes(y = sample, x = aligned_reads, 
           fill = mode)) +
  facet_grid(cols = vars(type)) +
  geom_density_ridges(alpha = 1, 
                      aes(height =..ndensity..), 
                      scale = 0.9) +
  theme_Publication_white() +
  geom_vline(xintercept = c(362,1541, 2903)) +
  scale_x_continuous(limits = c(0,3500), expand = c(0,0)) +
  scale_fill_manual(values = cbf1[c(2,5,3)])
dev.off()

pdf(here("figures/R_plots/210407_aligned_read_identity.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = total_frame_untrimmed_f_stats, 
       aes(y = sample, x = identity, 
           fill = mode)) +
  facet_grid(cols = vars(type)) +
  geom_density_ridges(alpha = 1, 
                      aes(height =..ndensity..), 
                      scale = 0.9) +
  theme_Publication_white() +
  geom_vline(xintercept = c(362,1541, 2903)) +
  scale_x_continuous(limits = c(50,100), expand = c(0,0)) +
  scale_fill_manual(values = cbf1[c(2,5,3)])
dev.off()

# ---> median aligned read length values & read identity
total_frame_untrimmed_f_stats %>%
  dplyr::filter(type == "CDS") %>%
  group_by(mode, substr(sample, 1, 10)) %>%
  summarise(median(aligned_reads), median(identity))

pdf(here("figures/R_plots/210406_perc_mapped.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = merged_frame,
       aes(x = sample, y = percentage_mapped, fill = mode)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  geom_bar(stat = "identity", color = "black") +
  theme_Publication_white() +
  coord_flip() +
  scale_fill_manual(values = cbf1[c(2,5,3)])
dev.off()


# > mapped bases
summary_frame_sample_f2 <- summary_frame_sample %>%
  group_by(sample) %>%
  summarise(all = sum(sequence_length_template)) %>%
  mutate(mode = substr(sample, 8,10)) %>%
  dplyr::filter(sample %in% bc_to_sample$sample[c(1,5:12,15:16)])

total_frame_untrimmed_f2 <- total_frame_untrimmed %>%
  dplyr::rename(sample = method) %>%
  dplyr::filter(sample %in% bc_to_sample$sample[c(1,5:12,15:16)]) %>%
  distinct(minion_read_name, .keep_all = T) %>%
  group_by(sample) %>%
  summarise(mapped = sum(aligned_reads)) 

merged_frame2 <- left_join(summary_frame_sample_f2, total_frame_untrimmed_f2, by = "sample") %>%
  mutate(percentage_mapped = mapped/all*100)
merged_frame2$sample <- factor(merged_frame2$sample,
                               levels = rev(bc_to_sample$sample[c(1,15,16,11,12,5,7,6,8:10)]))

pdf(here("figures/R_plots/210406_perc_mapped_bases.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = merged_frame2,
       aes(x = sample, y = percentage_mapped, fill = mode)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  geom_bar(stat = "identity", color = "black") +
  theme_Publication_white() +
  coord_flip() +
  scale_fill_manual(values = cbf1[c(2,5,3)])
dev.off()

merged_frame2 %>% 
  group_by(mode) %>%
  summarise(mean(percentage_mapped))

## mapped files | pychopper auto cutadapt clipped -------------------------------------
files   <- list.files(path = paste0(dir, "/data/mapped_data_pychopper_auto_cutadapt_clipped"), 
                      pattern = ".sorted.bam$",
                      recursive = T, full.names = T)

mapped_files   <- files[which(1:length(files) %% 2 == 1)] 
remapped_files <- files[which(1:length(files) %% 2 == 0)] 
dataset_names  <- str_split_fixed(str_split_fixed(mapped_files, "\\/", n = 9)[,8],"_fu",2)[,1]

total_frame   <- data.table()
partial_frame <- data.table()
for(i in seq_along(dataset_names)){
  print(paste0("file number ", i, " of ", length(dataset_names)))
  tic("calc bam files")
  partial_frame <- annotate_bams(input_mapped   = mapped_files[i], 
                                 input_remapped = remapped_files[i], 
                                 dataset        = dataset_names[i])
  total_frame <- rbind(total_frame, partial_frame) %>%
    as_tibble()
  toc()
}

total_frame %>%
  group_by(type, method) %>%
  summarise(n = n())

# save data frame ------------------------------------------------------------------
fwrite(total_frame, paste0(dir, "/data/mapped_data_pychopper_auto_cutadapt_clipped.tsv"), sep = "\t",col.names = T)




# read in again with ----------------------------------------------
total_frame <- vroom(here("data/mapped_data_no_trimming"))

total_frame_s <- total_frame %>%
  dplyr::filter(method %in% "201210_PCB109_Ecoli_NOTEX_replicate1.sorted.bam") %>%
  distinct(minion_read_name, .keep_all = T) %>%
  dplyr::filter(aligned_reads >= 1)

# plotting --------------------------------------------------------------------------

## prepare quantification summary --------------------------------------------------
# > mutate bam frames
total_frame_q <- mutate_bam_files(total_frame)

# > calc bam counts without rRNA loci
total_frame_c <- calc_bam_counts(total_frame_q, mode = "other") %>%
  dplyr::mutate(type = factor(type, levels=c("unknown","ncRNA", "rRNA","CDS")))

# > calc bam counts with rRNA loci
total_frame_c_rRNA <- calc_bam_counts(total_frame_q, mode = "rRNA") %>%
  dplyr::mutate(type = factor(type, levels=c("unknown","ncRNA", "5S", "16S","23S","CDS")))

## summary plots of mapped reads ====================================================

### relative number of counts (with rRNA) ###################################
pdf(here("figures/R_plots/210406_relative_counts_rRNA.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = total_frame_c_rRNA, 
       aes(x = method, y = number_of_counts_p,
           fill = type)) +
  geom_bar(stat = "identity", color = "black", size = 0.3) +
  theme(axis.text.x = element_text(hjust = 0,angle = 90)) +
  #scale_fill_manual(values = rev(ibm_colors[1:6])) +
  viridis::scale_fill_viridis(option = "cividis", begin = 1, end = 0, discrete = T) +
  theme_Publication_white() +
  scale_y_continuous(expand = c(0,0)) +
  coord_flip()
dev.off()

### relative number of bases (with rRNA) ###################################
pdf(here("figures/R_plots/210406_relative_bases_rRNA.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = total_frame_c_rRNA, 
       aes(x = method, y = number_of_reads_p,
           fill = type)) +
  geom_bar(stat = "identity", color = "black", size = 0.3) +
  theme(axis.text.x = element_text(hjust = 0,angle = 90)) +
  #scale_fill_manual(values = rev(ibm_colors[1:6])) +
  viridis::scale_fill_viridis(option = "cividis", begin = 1, end = 0, discrete = T) +
  theme_Publication_white() +
  scale_y_continuous(expand = c(0,0)) +
  coord_flip()
dev.off()

### relative number of counts without rRNA loci ###################################
total_frame_c <- total_frame_c %>% 
  dplyr::filter(method %in% bc_to_sample$sample[c(1,5:12,15:16)])
total_frame_c$method <- factor(total_frame_c$method,
                               levels = rev(bc_to_sample$sample[c(1,15,16,11,12,5,7,6,8:10)]))


# The palette with highlight:
cbf1_high <- c("#EFEAFF", "#ABC2DA","#F6B2FB","#648FFF")

levels(as.factor(total_frame_c$type))
pdf(here("figures/R_plots/210406_relative_counts.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = total_frame_c, 
       aes(x = method, y = number_of_counts_p,
           fill = type)) +
  geom_bar(stat = "identity", color = "black", size = 0.3) +
  theme(axis.text.x = element_text(hjust = 0,angle = 90)) +
  scale_fill_manual(values = cbf1_high) +
  theme_Publication_white() +
  scale_y_continuous(expand = c(0,0)) +
  coord_flip()
dev.off()

# > calculate values
total_frame_c %>%
  dplyr::filter(type == "CDS") %>%
  mutate(mode = substr(method, 1,10)) %>%
  group_by(mode) %>%
  summarise(mean(number_of_reads_p))

### absolute number of counts (not necessary for plotting) ###################################
# > good to compare to max number of sequenced reads to show how many map
helper <- summary_frame_b %>% 
  dplyr::filter(!barcode %in% c("barcode07", "barcode08", "barcode09", "barcode10", "barcode11", "barcode12")) %>%
  dplyr::filter(sequence_length_template > 0) %>%
  group_by(seq_run, barcode) %>%
  summarise(n = n()) %>%
  mutate(sample = paste0(seq_run, "_", barcode),
         mode = substr(sample, 8,10)) %>%
  ungroup() %>%
  dplyr::filter(n > 31000, barcode != "unclassified") 

helper$method <- total_frame_c$method[c(1:5,7,6,8,9,10)]

ggplot(data = total_frame_c, 
       aes(x = method, y = number_of_counts,
           fill = type)) +
  geom_bar(stat = "identity", color = "black", size = 0.3) +
  theme(axis.text.x = element_text(hjust = 0,angle = 90)) +
  scale_fill_manual(values = ibm_colors[c(6,5,3,1)]) +
  theme_Publication_white() +
  scale_y_continuous(expand = c(0,0))

### relative number of aligned bases ###################################
pdf(here("figures/R_plots/210406_relative_bases.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(total_frame_c,
       aes(x = method, y = number_of_reads_p,
           fill = type)) +
  geom_bar(stat = "identity", color = "black", size = 0.3) +
  theme(axis.text.x = element_text(hjust = 0,angle = 90)) +
  #scale_fill_manual(values = ibm_colors[c(6,5,3,1)]) +
  scale_fill_manual(values = cbf1_high) +
  theme_Publication_white() +
  scale_y_continuous(expand = c(0,0)) +
  coord_flip()
dev.off()

library(chroma)
cubehelix_colors(10)
cubehelix_palette()


## relative number of aligned bases ###################################
ggplot(data = total_frame_c %>% dplyr::filter(type == "CDS"), 
       aes(x = method, y = number_of_reads,
           fill = type)) +
  geom_bar(stat = "identity", color = "black", size = 0.2) +
  theme(axis.text.x = element_text(hjust = 0,angle = 90)) +
  scale_fill_manual(values = ibm_colors)


## Comparison of other stats =============================================
# > make CDS frame
total_frame_c <- total_frame_c %>% 
  dplyr::filter(method %in% bc_to_sample$sample[c(1,5:12,15:16)])
total_frame_c$method <- factor(total_frame_c$method,
                               levels = rev(bc_to_sample$sample[c(1,15,16,11,12,5,7,6,8:10)]))

total_frame_cds <- total_frame %>%
  dplyr::filter(type == "CDS") %>%
  mutate(mode = substr(method, 8,10))

### aligned read comparison ##########################
pdf(here("figures/R_plots/210215_CDS_aligned_reads.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = total_frame_cds, 
       aes(x = aligned_reads, y = method, fill = factor(..quantile..))) +
  stat_density_ridges(geom = "density_ridges_gradient", 
                      aes(height =..ndensity..),
                      quantile_lines = T, scale = 0.95,
                      quantiles = c(0.25,0.5,0.75), lwd = 0.5) +
  theme_Publication_white() +
  scale_fill_manual(values = alpha(rev(c("#00204D", "#40587A", "#7F8FA6", "#BFC7D2")), 1)) +
  scale_x_continuous(limits = c(0,2200), expand = c(0,0)) 
dev.off()

rcartocolor::display_carto_all(colorblind_friendly = T)

### identity comparison
pdf(here("figures/R_plots/210215_CDS_identity.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = total_frame_cds, 
       aes(x = identity,y = method, 
           fill = factor(..quantile..))) +
  scale_x_continuous(limits = c(65,100), expand = c(0,0)) +
  stat_density_ridges(geom = "density_ridges_gradient", alpha = 0.75, 
                      quantile_lines = T, 
                      aes(height =..ndensity..),scale = 0.95,
                      quantiles = c(0.25,0.5,0.75)) +
  theme_Publication_white() +
  scale_fill_manual(values = alpha(rev(c("#00204D", "#40587A", "#7F8FA6", "#BFC7D2")), 1)) 
dev.off()



# combine sequencing summary with mapped reads ----------------------------------------
total_frame_cds_summary <- total_frame_cds %>%
  left_join(summary_frame, by = c("minion_read_name" = "read_id")) %>%
  dplyr::filter(method %in% c("190123_RNA001_Ecoli_TEX_replicate1.sorted.bam", 
                              "201125_DCS108_Ecoli_NOTEX_replicate1.sorted.bam",
                              "201210_PCB109_Ecoli_NOTEX_replicate1.sorted.bam"))

## aligned_Read length vs sequenced length
pdf(here("figures/R_plots/210215_CDS_sequenced_vs_aligned.pdf"), 
    width = 12, height = 7, paper = "special", onefile=FALSE)
ggplot(data = total_frame_cds_summary, 
       aes(x = sequence_length_template, 
           y = as.numeric(aligned_reads))) +
  geom_abline(linetype = "dashed") +
  facet_grid(cols = vars(mode)) +
  xlab("Sequenced length (%)") +
  ylab("Aligned length (nt)") +
  #geom_hex(bins = 70) +
  stat_binhex(aes(fill= ..density..*100,
                  color= ..density..*100),bins = 75) +
  viridis::scale_fill_viridis(option = "cividis") + 
  viridis::scale_color_viridis(option = "cividis") + 
  #scale_alpha(range = c(0.7,1)) +
  theme_Publication_white() +
  guides(alpha = F) +
  scale_x_continuous(limits = c(50,10000), expand = c(0,0), trans = "log10") +
  scale_y_continuous(limits = c(50,10000), expand = c(0,0), trans = "log10") 
dev.off()


total_frame_cds_summary <- total_frame %>%
  dplyr::filter(type == "CDS") %>%
  mutate(mode = substr(method, 8,10)) %>%
  dplyr::filter(method %in% bc_to_sample$sample[c(1,5:12,15:16)]) %>%
  left_join(summary_frame, by = c("minion_read_name" = "read_id"))

subsample_number <- 5000 #0
subsample_raw <- function(data, s_n, mode_t){
  data  %>%
    dplyr::filter(mode == mode_t) %>%
    sample_n(s_n)
}

total_frame_DCS <- subsample_raw(total_frame_cds_summary, subsample_number, "DCS")
total_frame_PCB <- subsample_raw(total_frame_cds_summary, subsample_number, "PCB")
total_frame_RNA <- subsample_raw(total_frame_cds_summary, subsample_number, "RNA")

total_frame_all <- rbind(total_frame_DCS,
                         total_frame_PCB,
                         total_frame_RNA)

total_frame_all$method <- factor(total_frame_all$method,
                                 levels = rev(bc_to_sample$sample[c(1,15,16,11,12,5,7,6,8:10)]))
total_frame_all$mode <- factor(total_frame_all$mode,
                               levels = c("RNA", "DCS", "PCB"))

# > CDS sequenced vs aligned
pdf(here("figures/R_plots/210406_read_length_vs_quality_group.pdf"), 
    width = 6, height = 7, paper = "special", onefile=FALSE)
ggplot(data = total_frame_all,
       aes(x = sequence_length_template, 
           y = aligned_reads, 
           color = mode)) +
  facet_grid(rows = vars(mode)) +
  geom_point(alpha = 0.2, size = 1) +
  geom_abline(linetype = "dashed", slope = 1) +
  geom_abline(linetype = "dashed", slope = .5) +
  stat_density2d(aes(alpha=..level.., fill = mode),color = NA,
                 bins=10, geom="polygon") +
  geom_density2d(color = "black", contour_var = "ndensity", bins = 10, aes(alpha = ..level..)) +
  xlab("Sequenced length") +
  ylab("Aligned length (nt)") +
  theme_Publication_white() +
  guides(alpha = F) +
  scale_x_continuous(limits = c(0,3500), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,3500), expand = c(0,0)) +
  scale_fill_manual(values =(cbf1[c(3,2,5)])) +
  scale_color_manual(values =(cbf1[c(3,2,5)])) +
  coord_equal()
dev.off()



# > CDS qscore vs aligned identity
pdf(here("figures/R_plots/210406_qscore_vs_identity.pdf"), 
    width = 6, height = 7, paper = "special", onefile=FALSE)
ggplot(data = total_frame_all,
       aes(x = mean_qscore_template, 
           y = identity, 
           color = mode)) +
  facet_grid(rows = vars(mode)) +
  geom_point(alpha = 0.2, size = 1) +
  stat_density2d(aes(alpha=..level.., fill = mode),color = NA,size = 0.5,
                 bins=15, geom="polygon") +
  geom_density2d(color = "black", contour_var = "ndensity", bins = 15, size = 0.5,aes(alpha = ..level..)) +
  xlab("Mean qscore template") +
  ylab("Identity") +
  theme_Publication_white() +
  guides(alpha = F) +
  scale_x_continuous(limits = c(0,20), expand = c(0,0)) +
  scale_y_continuous(limits = c(50,100), expand = c(0,0)) +
  scale_fill_manual(values =(cbf1[c(3,2,5)])) +
  scale_color_manual(values =(cbf1[c(3,2,5)])) 
dev.off()



## plot groups --> qscore mapped vs unmapped
total_frame_cds_summary$mode <- factor(total_frame_cds_summary$mode,
                                       levels = rev(c("RNA", "DCS", "PCB")))
pdf(here("figures/R_plots/210407_split_violin_quality_DCS.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = total_frame_cds_summary, 
       aes(x = mode, 
           y = mean_qscore_template, 
           factor = sequence_length_template <= 1.75*aligned_reads)) +
  geom_split_violin(trim = T, scale = "width", aes(fill = mode)) +
  coord_flip() +
  theme_Publication_white() +
  scale_y_continuous(limits = c(0,20), expand = c(0,0)) +
  scale_fill_manual(values = cbf1[c(5,2,3)]) +
  scale_color_manual(values = cbf1[c(5,2,3)])
dev.off()

pdf(here("figures/R_plots/210407_split_violin_identity_DCS.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = total_frame_cds_summary, 
       aes(x = mode, 
           y = identity, 
           factor = sequence_length_template <= 1.75*aligned_reads)) +
  geom_split_violin(trim = T, scale = "width", aes(fill = mode)) +
  coord_flip() +
  theme_Publication_white() +
  scale_y_continuous(limits = c(50,100), expand = c(0,0)) +
  scale_fill_manual(values = cbf1[c(5,2,3)]) +
  scale_color_manual(values = cbf1[c(5,2,3)])
dev.off()


a <- total_frame_cds_summary %>%
  mutate(strandswitch = ifelse(sequence_length_template <= 1.75*aligned_reads, "correct", "not_correct")) %>%
  group_by(mode, strandswitch) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = strandswitch, values_from = n) %>%
  mutate(c_p = correct/(correct+not_correct),
         n_p = not_correct/(correct+not_correct)) %>%
  dplyr::select(-correct, -not_correct) %>%
  pivot_longer(c_p:n_p, values_to = "perc")


pdf(here("figures/R_plots/210407_quantification_2D.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = a, aes(y = mode, x = perc, fill = mode, group = name == "c_p")) +
  geom_col(color = "black") +
  scale_fill_manual(values = cbf1[c(5,2,3)]) +
  theme_Publication_white() +
  scale_x_continuous(limits = c(0,1), expand = c(0,0)) 
dev.off()




