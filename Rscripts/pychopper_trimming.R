# >> Pychopper trimming << #

# load libraries ----
library(here)
source(here("Rscripts/load_libraries.R"))

# functions & defs ----
calc_n_fastq <- function(myinput){
  pos <- str_count(myinput, pattern = "\\/")
  fq <- ShortRead::readFastq(myinput) %>%
    ShortRead::id() %>%
    as.data.table() %>%
    mutate(minion_read_name = if(sum(str_detect(string = x, pattern = "\\|")) > 0){
      str_split_fixed(str_split_fixed(x, "\\|",2)[,2], " runid", 2)[,1]
    } else {str_split_fixed(x, " runid", 2)[,1]}) %>%
    dplyr::filter(minion_read_name %in% total_frame$minion_read_name)
  
  f <- data.table(n = nrow(fq),
                  sample = str_split_fixed(str_split_fixed(str_split_fixed(myinput, "\\/", n = pos+2)[,pos+1], "_fu", 2)[,1], "\\.", 2)[,1]) %>%
    as_tibble()
  return(f)
}

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
# load & tidy data ----

## load saved mapped data table ====
dir <- here()
total_frame <- vroom(paste0(dir, "/data/mapped_data_no_trimming.tsv"), num_threads = 8) %>%
  dplyr::filter(type == "CDS")

## Calculation: Number of reads in pychopper categories ====

### load fastq files pychopper trimmed/DCS pychopper rescued/ raw ####
fl_files              <- grep(list.files(paste0(dir,"/data/pychopper_data_auto/"), full.names = T,recursive = T, pattern = "full_length_output.fastq"), pattern = "RNA", invert = T, value = T)
fl_rescue_files       <- grep(list.files(paste0(dir,"/data/fastq_fl_data/"), full.names = T,recursive = T, pattern = "_full_length_all.fastq"), pattern = "RNA", invert = T, value = T)
al_files              <- grep(list.files(paste0(dir,"/data/fastq_data/"), full.names = T,recursive = T, pattern = ".fastq"), pattern = "RNA", invert = T, value = T)

### calc number of fl, fl_rescue and all reads ####
pychopper_frame <- bind_rows(
  rbindlist(map(fl_files, calc_n_fastq)) %>%
    mutate(method = "full_length"),
  rbindlist(map(fl_rescue_files, calc_n_fastq)) %>%
    mutate(method = "full_length_rescued"),
  rbindlist(map(al_fq, calc_n_fastq)) %>%
    mutate(method = "all")) %>%
  mutate(mode = str_sub(sample,1,3))

### modify output table, scale values to 100 ####
pychopper_frame_group <- pychopper_frame %>%
  pivot_wider(names_from = method, values_from = n) %>%
  dplyr::mutate(full_length = full_length/all*100,
                full_length_rescued = full_length_rescued/all*100 - full_length,
                all = 100 - full_length - full_length_rescued) %>%
  pivot_longer(full_length:all,names_to = "method", values_to = "n")
  

## Calculation: Distance of aligned untrimmed/pychopper reads ====

### genome data ####
gff_table      <-  read_in_gff(paste0(dir, "data/genome_data/NC_000913.3.gff3"))

### mapped files | after pychopper trimming ####
files   <- grep(list.files(path = paste0(dir,"/data/mapped_data_pychopper_auto"), 
                      pattern = ".sorted.bam$",
                      recursive = T, full.names = T), pattern = "RNA", invert = T, value = T)

mapped_frame   <- pmap_dfr(list(files[!str_detect(files, "remapped")],
                                files[str_detect(files, "remapped")],
                                str_split_fixed(str_split_fixed(files[!str_detect(files, "remapped")], "\\/", n = 9)[,8],"_fu",2)[,1]),
                           annotate_bams)

# > filter for CDS mapping reads 
mapped_frame_pychopper <- mapped_frame %>%
  left_join(old_new, by = c("sample" = "old_name")) %>%
  dplyr::mutate(sample = new_name) %>%
  dplyr::select(-new_name) %>%
  dplyr::filter(!is.na(sample)) %>%
  dplyr::filter(type == "CDS") %>%
  dplyr::mutate(mode = substr(sample, 1,3))
  
### merge untrimmed & pychopper bam frames /group per read id, calc difference #####
total_frame_merged <- rbindlist(list(total_frame %>% 
                                       dplyr::mutate(mode = substr(sample, 1,3),
                                                     group = "total"),
                                     mapped_frame_pychopper %>% 
                                       dplyr::mutate(group = "pychopper"))) %>%
  mutate(minion_read_name = ifelse(group == "pychopper", 
                                   str_split_fixed(minion_read_name, "\\|",2)[,2], minion_read_name)) %>%
  dplyr::filter(sample %in% bc_to_sample$sample[c(8,9,2,4,3,5,6,7)]) %>%
  rowid_to_column("id") %>%
  dplyr::select(sample,mode,minion_read_name,group,aligned_reads) %>%
  pivot_wider(names_from = group, values_from = aligned_reads, values_fn = {max}) %>%
  mutate(diff = total - pychopper) 
  

# PLOTS ----

## reorder levels ====
pychopper_frame_group$sample <- factor(pychopper_frame_group$sample,
                                       levels = rev(bc_to_sample$sample[c(8,9,2,4,3,5,6,7)]))
pychopper_frame_group$method <- factor(pychopper_frame_group$method,
                                       levels = rev(c("full_length", "full_length_rescued", "all")))

total_frame_merged$sample    <- factor(total_frame_merged$sample,
                                       levels = rev(bc_to_sample$sample[c(8,9,2,4,3,5,6,7)]))

## plotting ====
### Pychopper categories - Supplementary Fig. 12B ####
raw_reads_plotting(pychopper_frame_group, n, sample, mode, cbf1[c(2,5)]) +
  geom_bar(stat = "identity", color = "black", aes(alpha = method), position = position_stack()) +
  scale_x_continuous(expand = c(0,0), limits = c(0,100)) +
  xlab("Distribution of CDS-mapping reads (%)") 

### Aligned read distance - Supplementary Fig. 12C ####
raw_reads_plotting(total_frame_merged, diff, sample, mode, cbf1[c(2,5)]) +
  geom_density_ridges(aes(height =..ndensity..),stat = "binline", binwidth = 1,
                      scale = 0.9) +
  scale_x_continuous(expand = c(0,0), limits = c(-20,100)) +
  xlab("Distance untrimmed-full-length aligned reads (bases)") +
  scale_y_discrete(expand = c(0.01,0.1)) 

