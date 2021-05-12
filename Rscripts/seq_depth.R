# >> sequencing depth << #

# load libraries ----
library(here)
source(here("Rscripts/load_libraries.R"))

# functions & defs ----
genes_seen <- function(mydf, myfilter){
  mydf %>%
    dplyr::filter(type == myfilter) %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(total_number = n()) %>%
    dplyr::summarise(CDS_seen = length(unique(gene)), 
                     nr_of_reads = max(total_number)) 
}

subsample_reads_to_gene <- function(mydf, sample_number){
  genes_seen({{mydf}} %>% 
               sample_n(sample_number), "CDS")
}

# load & tidy data ----
dir <- here()

## genome data ====
gff_table      <-  read_in_gff(paste0(dir, "data/genome_data/NC_000913.3.gff3"))

## load saved mapped data tables ====
total_frame <- vroom(paste0(dir, "/data/mapped_data_no_trimming.tsv"))
other_frame <- vroom(paste0(dir, "/data/mapped_data_sra.tsv"))

## calc mRNA sequencing depth X-fold ====
### Nanopore data ####
total_frame_CDS <- total_frame %>%
  dplyr::filter(type == "CDS") %>%
  group_by(sample) %>%
  summarise(depth = sum(aligned_reads)/sum(gff_table$width)) %>%
  mutate(mode = substr(sample, 1,3)) %>%
  arrange(factor(sample, levels = bc_to_sample$sample[c(1,10,11,8,9,2,4,3,5,6,7)]))

### Comparison data/made in the same way as for the Nanopore reads ####
other_frame_CDS <- other_frame %>%
  dplyr::filter(type == "CDS") %>%
  group_by(sample) %>%
  summarise(depth = sum(aligned_reads)/sum(gff_table$width)) %>%
  mutate(mode = substr(sample, 1,7))

### combine frame ####
all_frame <- rbindlist(list(total_frame_CDS, other_frame_CDS))

### calc number of seen genes ####
genes_seen_frame  <- rbindlist(list(genes_seen(total_frame, "CDS"), 
                                    genes_seen(other_frame%>% dplyr::rename(sample = method), "CDS")))

### number of genes seen - subsampled ####
# > selected example data set 
total_frame_selected <- total_frame %>%
  dplyr::filter(sample == "PCB109_PCR12_Ecoli_NOTEX_replicate4",
                type == "CDS")

# > set subset numbers 
subset_numbers <- c(10,30,100,300,500,1000,5000,
                    10000,20000,30000,40000,50000,
                    60000,70000,80000,90000,100000,
                    110000,120000,160000,
                    200000,250000,300000)

# > subsample and calc number of seen genes
depth_frame <- data.table()
for(i in seq_along(subset_numbers)){
  tic("subset")
  
  subset_frame <- subsample_reads_to_gene(total_frame_selected, subset_numbers[i])
  depth_frame <- rbind(depth_frame, subset_frame)
  toc()
}

# PLOTS ----

## reorder levels ====
all_frame$sample <- factor(all_frame$sample,
                           levels = rev(c(bc_to_sample$sample[c(1,10,11,8,9,2,4,3,5,6,7)],
                           rev(c("SRR7533627", "SRR7533626","SRR1927169")))))
all_frame$mode <- factor(all_frame$mode,
                        levels = c("RNA", "DCS", "PCB", "SRR1927", "SRR7533"))

genes_seen_frame$sample <- factor(genes_seen_frame$sample,
                           levels = rev(c(bc_to_sample$sample[c(1,10,11,8,9,2,4,3,5,6,7)],
                                          rev(c("SRR7533627", "SRR7533626","SRR1927169")))))

## plotting ==== 

### Sequencing depth - Supplementary Fig. 10A ####
raw_reads_plotting(all_frame, depth, sample, mode, cbf1[c(2,5,3,1,4)]) +
  geom_bar(stat = "identity", color = "black") +
  scale_x_continuous(expand = c(0,0), limits = c(0,80)) +
  xlab("Sequencing depth (X-fold)") 

### Genes seen - Supplementary Fig. 10B ####
raw_reads_plotting(genes_seen_frame %>% mutate(mode = str_sub(sample, 1,5)), 
                   CDS_seen/nrow(gff_table %>% dplyr::filter(type == "CDS")) * 100, sample, mode, cbf1[c(2,5,3,1,4)]) +
  geom_bar(stat = "identity", color = "black") +
  scale_x_continuous(expand = c(0,0), limits = c(0,100)) +
  xlab("Sequencing depth (X-fold)") 

### Genes seen per subsampled reads - Supplementary Fig. 10C ####
ggplot(data = depth_frame, aes(x = nr_of_reads, 
                                 y = CDS_seen/nrow(gff_table %>% dplyr::filter(type == "CDS")) * 100)) +
  geom_line() +
  geom_area() +
  geom_point(size = 2, shape = 21) +
  scale_y_continuous(limits = c(0,100), expand = c(0,0)) +  
  #scale_x_log10() +
  geom_vline(xintercept = 70000) +
  theme_Publication_white() +
  scale_color_manual(values = cbf1[5]) +
  scale_fill_manual(values = cbf1[5])

