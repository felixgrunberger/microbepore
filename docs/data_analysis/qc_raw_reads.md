---
layout: page
title: Quality control of raw reads
parent: Data analysis
nav_order: 2
---

## Quality control of raw reads 


```bash
# Raw read analysis ----

# load libraries ----
library(here)
source(here("Rscripts/load_libraries.R"))

# functions ----
grep_summary <- function(input, dataset){
  
  # read in sequencing summary files
  # only select relevant columns, check if file was demultiplexed by live-basecalling or demultiplexed seperately for barcode detection
  suppressMessages(vroom(file = paste0(dir,"/data/summary_data/",input), 
                         num_threads = 8)) %>%
    dplyr::mutate(barcode = if("barcode_arrangement" %in% colnames(.)) barcode_arrangement else "no_barcode",
                  seq_run = dataset) %>%
    dplyr::select(seq_run,read_id, run_id, sequence_length_template, mean_qscore_template, barcode) 
}

read_barcode <- function(input){
  
  vroom(paste0(paste0(dir,"/data/barcode_data/"), input), num_threads = 8) %>%
    dplyr::select(read_id, barcode_arrangement) %>%
    as_tibble() 
}

raw_reads_plotting <- function(mydf, myx, myy, myfill, mypalette){
  ggplot(data = mydf, aes(x = {{myx}}, y = {{myy}}, fill = {{myfill}})) +
    #scale_x_continuous(expand = c(0,0), limits = c(0,2.5)) +
    #geom_bar(stat = "identity", color = "black") +
    theme_Publication_white() +
    ylab("") +
    #xlab("Number of reads (in Millions)") +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(linetype = "dashed", color = "black")) +
    scale_fill_manual(values = mypalette)
}


# data ----

## directory ====
dir <- "/Volumes/EX_SSD/"

## sequencing summary data ====
files          <- list.files(paste0(dir,"/data/summary_data/"), recursive = T, pattern = "sequencing_summary.txt.gz")
summary_frame  <- pmap_dfr(list(files, str_split_fixed(files, "_seq", n = 2)[,1]), grep_summary)

## barcode summary tables from guppy output ====
b_files          <- list.files(paste0(dir,"/data/barcode_data/"), recursive = T, pattern = "txt.gz")
barcode_frame    <- pmap_dfr(list(b_files), read_barcode)

### merge summary_files with barcode files and detect barcode if not live-demultiplexed ####
summary_frame_f <- summary_frame %>%
  left_join(as_tibble(barcode_frame), by = c("read_id")) %>%
  mutate(barcode = ifelse(is.na(barcode_arrangement), barcode, barcode_arrangement))

## make alternative data frame, 1 for each replicate/sample -------------------------------
bc_to_sample <- data.table(seq_run = c("190123_RNA001_Ecoli",
                                       rep("201208_PCB109_Ecoli",4),
                                       rep("201210_PCB109_Ecoli",2),
                                       rep("210317_DCS109_Ecoli",2),
                                       "210317_RNA002_Ecoli",
                                       "210318_RNA002_Ecoli"),
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

summary_frame_sample <- summary_frame_f %>%
  left_join(bc_to_sample, by = c("seq_run", "barcode")) %>%
  dplyr::filter(!is.na(sample)) %>%
  mutate(mode = substr(sample, 1,3))

# save data frame ----
fwrite(summary_frame_sample, 
       paste0(dir, "/data/summary_data/merged_sequencing_summary_frame.txt"), 
       sep = "\t",
       col.names = T)

# calculate stats ----
summary_stats <- summary_frame_sample %>%
  group_by(sample) %>%
  summarise(number_of_reads = n(),
            number_of_bases = sum(sequence_length_template),
            median_raw_length = median(sequence_length_template),
            median_raw_qscore = round(median(mean_qscore_template), digits = 2)) %>%
  mutate(mode = substr(sample, 1,3)) %>%
  arrange(factor(sample, levels = bc_to_sample$sample[c(1,10,11,8,9,2,4,3,5,6,7)]))

# plotting ----

## reorder levels ====
summary_frame_sample$sample <- factor(summary_frame_sample$sample,
                               levels = rev(bc_to_sample$sample[c(1,10,11,8,9,2,4,3,5,6,7)]))
summary_frame_sample$mode <- factor(summary_frame_sample$mode,
                                      levels = c("RNA", "DCS", "PCB"))
summary_stats$sample <- factor(summary_stats$sample,
                               levels = rev(bc_to_sample$sample[c(1,10,11,8,9,2,4,3,5,6,7)]))

## color palette ====
cbf1 <- c("#EFEAFF","#F5AAA3","#CFCFCF", "#F6B2FB", "#ABC2DA")

## plotting ====

  
### Total number of reads #### Supplementary Fig. 3A
raw_reads_plotting(summary_stats, number_of_reads/1000000, sample, mode, cbf1[c(2,5,3)]) +
  geom_bar(stat = "identity", color = "black") +
  scale_x_continuous(expand = c(0,0), limits = c(0,2.5)) +
  xlab("Number of reads (in Millions)") 

### Total number of bases #### Supplementary Fig. 3A
raw_reads_plotting(summary_stats, number_of_bases/1000000000, sample, mode, cbf1[c(2,5,3)]) +
  geom_bar(stat = "identity", color = "black") +
  scale_x_continuous(expand = c(0,0), limits = c(0,2)) +
  xlab("Number of bases (in Gb)")

### read quality distribution #### Supplementary Fig. 4A
raw_reads_plotting(summary_frame_sample, mean_qscore_template, sample, mode, cbf1[c(2,5,3)]) +
  geom_density_ridges(aes(height =..ndensity..), scale = 0.9, color = "black") +
  scale_x_continuous(limits = c(0,20), expand = c(0,0)) +
  scale_y_discrete(expand = c(0.01,0.1)) +
  xlab("Mean Qscore template (Phred-like score)") 

### read quality distribution #### Supplementary Fig. 4B
raw_reads_plotting(summary_frame_sample, sequence_length_template, sample, mode, cbf1[c(2,5,3)]) +
  geom_density_ridges(aes(height =..ndensity..), scale = 0.9, color = "black") +
  scale_x_continuous(limits = c(0,3500), expand = c(0,0)) +
  scale_y_discrete(expand = c(0.01,0.1)) +
  xlab("Read length (bases)") 
```
Some R code: 

{% highlight R %}
a  1:10
b <- 2:20

{% endhighlight %}


{% highlight R %}
# 1. step: Load all packages we need ----------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(plotly)
library(viridis)
library(gganimate)
library(webshot)
{% endhighlight %}


What was i doing here.
