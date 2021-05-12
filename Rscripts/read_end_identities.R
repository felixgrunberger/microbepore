# >> Read end positions - trimming effects << #

# load libraries ----
library(here)
source(here("Rscripts/load_libraries.R"))

# functions & defs ----
get_freq_nucleotide <- function(wanted_position, grouping, input_fastq){
  if(grouping == "prime5"){
    alphabetFrequency(subseq(input_fastq,wanted_position,wanted_position), 
                      as.prob = T, collapse = T, baseOnly = T) %>%
      as.data.frame() %>%
      rownames_to_column() %>%
      mutate(position = wanted_position) %>%
      as_tibble() %>%
      dplyr::rename(nucleotide = 1,
                    perc = 2) %>%
      mutate(group = grouping)}
  else{
    alphabetFrequency(subseq(input_fastq,-wanted_position,-wanted_position), 
                      as.prob = T, collapse = T, baseOnly = T) %>%
      as.data.frame() %>%
      rownames_to_column() %>%
      mutate(position = -wanted_position) %>%
      as_tibble() %>%
      dplyr::rename(nucleotide = 1,
                    perc = 2) %>%
      mutate(group = grouping)
  }
}

calculate_end_identities <- function(inputFastq, mygroup){
  nucleotide_frame_prime5 <- data.table()
  nucleotide_frame_prime3 <- data.table()
  final_nucleotide_frame  <- data.table()
  
  for(i in 1:5){
    nucleotide_frame_prime5 <- get_freq_nucleotide(wanted_position = i,
                                                   grouping = "prime5",
                                                   input_fastq = {{inputFastq}}@sread[width({{inputFastq}}@sread)>10])
    nucleotide_frame_prime3 <- get_freq_nucleotide(wanted_position = i,
                                                   grouping = "prime3",
                                                   input_fastq = {{inputFastq}}@sread[width({{inputFastq}}@sread)>10])
    final_nucleotide_frame <- rbind(final_nucleotide_frame, nucleotide_frame_prime5, nucleotide_frame_prime3)
  }
  return(final_nucleotide_frame %>%
           mutate(method = mygroup))
}

calc_As <- function(input_fq, input_n){
  data.table(n   = unlist(sum(str_count(string = {{input_fq}}@sread, pattern = paste0("A{",input_n,",}")))),
             n_A = input_n)
}

# load & tidy data ----

## read in fastq data ====
dir <- here()
fl_files       <- list.files(paste0(dir,"/data/fastq_data/"), recursive = T, pattern = ".fastq")
dataset_names  <- str_split_fixed(str_split_fixed(fl_files, "\\/", n = 2)[,2], ".fastq", 2)[,1]

# > select dataset --> cDNA replicate 2
i <- 11

### genome ####
genome_fq <- ShortRead::readFasta(paste0(dir,"/data/genome_data/NC_000913.3.fasta"))

### raw ####
al_fq <- ShortRead::readFastq(paste0(dir,"/data/fastq_data/",paste(str_split_fixed(dataset_names[i],"_",5)[-c(4,5)],collapse = "_"), "/", dataset_names[i], ".fastq"))

### full-length ####
fl_fq <- ShortRead::readFastq(paste0(dir,"/data/fastq_fl_data/",paste(str_split_fixed(dataset_names[i],"_",5)[-c(4,5)],collapse = "_"), "/", dataset_names[i], "", "/",dataset_names[i], "_full_length_all.fastq"))

### polyA and SSP trimmed #####
cut_fq <- ShortRead::readFastq(paste0(dir,"/data/fastq_fl_cutadapt_SSP/",paste(str_split_fixed(dataset_names[i],"_",5)[-c(4,5)],collapse = "_"), "/", dataset_names[i], "_full_length_all", "/",dataset_names[i], "_full_length_all.cutadapt_SSP.fastq"))



## make end frequency table ====
end_table <- rbindlist(list(calculate_end_identities(al_fq, "raw"),
                            calculate_end_identities(fl_fq, "full-length"),
                            calculate_end_identities(cut_fq, "trimmed")))

## polyA trimming analysis ====

### genomic background ####
genome_frame <-rbindlist(pmap(list(list(genome_fq),1:20),calc_As)) %>%
  mutate(method = "genome_bg",
         total  = sum(width(genome_fq@sread)))

### raw ####
al_frame <-rbindlist(pmap(list(list(al_fq),1:20),calc_As)) %>%
  mutate(method = "raw",
         total  = sum(width(al_fq@sread)))

### full-length ####1:20
fl_frame <-rbindlist(pmap(list(list(fl_fq),1:20),calc_As)) %>%
  mutate(method = "full-length",
         total  = sum(width(fl_fq@sread)))

### polyA and SSP trimmed #####
cut_frame <-rbindlist(pmap(list(list(cut_fq),1:20),calc_As)) %>%
  mutate(method = "trimmed",
         total  = sum(width(cut_fq@sread)))

### bind all data ####
polyA_frame <- rbindlist(list(genome_frame, al_frame, fl_frame, cut_frame)) %>%
  mutate(n = n + 1) %>%
  mutate(n2 = (n/total*100)) %>%
  dplyr::select(method, n2, n_A) %>%
  as_tibble() %>%
  group_by(n_A) %>%
  mutate(bg = n2[method == "genome_bg"],
         scaled_value = n2/bg) 

# PLOTS ----

## reorder levels ====
end_table$method <- factor(end_table$method,
                            levels = (c("raw", "full-length", "trimmed")))

end_table$group <- factor(end_table$group,
                           levels = (c("prime5", "prime3")))

## plotting ==== 

### Read end identities cDNA replicate 2 - Supplementary Fig. 13A ####
raw_reads_plotting(end_table, position, perc,nucleotide, acgt_colors) +
  geom_bar(stat = "identity", col = "black") +
  facet_grid(cols = vars(group), rows = vars(method)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.grid.major.x = element_blank())

### polyA end enrichment cDNA replicate 2 - Supplementary Fig. 13B ####
raw_reads_plotting(polyA_frame, n_A, scaled_value,method, ibm_colors) +
  geom_line(size = 2, alpha = 0.7) +
  geom_point(shape = 21, color = "black", alpha = 0.8, size = 7, stroke = 2) +
  scale_y_log10() +
  coord_cartesian(ylim = c(0.1,10000)) +
  ylab("enrichment over genome in %") +
  xlab("Category A+ (X or more As") +
  theme_Publication_white()




  



