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

# load & tidy data ----

## load saved mapped data table ====
dir <- here()
total_frame <- vroom(paste0(dir, "/data/mapped_data_no_trimming.tsv"), num_threads = 8) %>%
  dplyr::filter(type == "CDS")

## load fastq files pychopper trimmed/DCS pychopper rescued/ raw ====
fl_files              <- grep(list.files(paste0(dir,"/data/pychopper_data_auto/"), full.names = T,recursive = T, pattern = "full_length_output.fastq"), pattern = "RNA", invert = T, value = T)
fl_rescue_files       <- grep(list.files(paste0(dir,"/data/fastq_fl_data/"), full.names = T,recursive = T, pattern = "_full_length_all.fastq"), pattern = "RNA", invert = T, value = T)
al_files              <- grep(list.files(paste0(dir,"/data/fastq_data/"), full.names = T,recursive = T, pattern = ".fastq"), pattern = "RNA", invert = T, value = T)

## calc number of fl, fl_rescue and all reads ====
pychopper_frame <- bind_rows(
  rbindlist(map(fl_files, calc_n_fastq)) %>%
    mutate(method = "full_length"),
  rbindlist(map(fl_rescue_files, calc_n_fastq)) %>%
    mutate(method = "full_length_rescued"),
  rbindlist(map(al_fq, calc_n_fastq)) %>%
    mutate(method = "all")) %>%
  mutate(mode = str_sub(sample,1,3))

## modify output table, scale values to 100 ====
pychopper_frame_group <- pychopper_frame %>%
  pivot_wider(names_from = method, values_from = n) %>%
  dplyr::mutate(full_length = full_length/all*100,
                full_length_rescued = full_length_rescued/all*100 - full_length,
                all = 100 - full_length - full_length_rescued) %>%
  pivot_longer(full_length:all,names_to = "method", values_to = "n")
  
# PLOTS ----

## reorder levels ====
pychopper_frame_group$sample <- factor(pychopper_frame_group$sample,
                                       levels = rev(bc_to_sample$sample[c(8,9,2,4,3,5,6,7)]))
pychopper_frame_group$method <- factor(pychopper_frame_group$method,
                                       levels = rev(c("full_length", "full_length_rescued", "all")))

## plotting ====
### Pychopper categories - Supplementary Fig. 12B ####
raw_reads_plotting(pychopper_frame_group, n, sample, mode, cbf1[c(2,5)]) +
  geom_bar(stat = "identity", color = "black", aes(alpha = method), position = position_stack()) +
  scale_x_continuous(expand = c(0,0), limits = c(0,100)) +
  xlab("Distribution of CDS-mapping reads (%)") 



ggplot(data = fl_frame_final) +
  geom_col(aes(y = sample, x = 100), fill = "white", color = "black") +
  geom_col(aes(y = sample, x = perc_cds_rescue, fill = mode), alpha = 0.5, color = "black")  +
  geom_col(aes(y = sample, x = perc_cds, fill = mode),color = "black")  +
  scale_fill_manual(values = cbf1[c(2,5,3)]) +
  scale_x_continuous(limits = c(0,100), expand = c(0,0)) +
  theme_Publication_white() 
dev.off()
