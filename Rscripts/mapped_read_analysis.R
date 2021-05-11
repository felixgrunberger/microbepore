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
dir <- here()

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
  mutate(mode = substr(sample, 1,3)) 

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
ggplot(data = summary_total, 
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
ggplot(data = summary_total, 
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

