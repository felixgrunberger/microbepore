# >> Coverage/end/unit example << #

# load libraries ----
library(here)
source(here("Rscripts/load_libraries.R"))

# functions & defs ----
get_cov <- function(file, dataset, low, high){
  fileset <- paste0(dir,"/data/coverage_data/")
  suppressMessages(vroom(paste0(fileset,file), col_names = F, num_threads = 8)) %>%
    dplyr::rename(position = X2, coverage = X3) %>%
    dplyr::filter(position > low, position < high) %>%
    mutate(sample = dataset,
           cov_scale = scales::rescale(coverage, to=c(0,100))) %>%
    dplyr::select(sample, position, coverage, cov_scale)
}

get_ends <- function(file, dataset, low, high){
  fileset <- paste0(dir,"/data/tts_data/")
  suppressMessages(vroom(paste0(fileset,file), col_names = F, num_threads = 8)) %>%
    dplyr::rename(position = X2, end = X3, coverage = X4) %>%
    dplyr::filter(position > low, position < high) %>%
    mutate(sample = dataset,
           cov_scale = scales::rescale(coverage, to=c(0,100))) %>%
    dplyr::select(sample, position, coverage, cov_scale) %>%
    mutate(bin = cut(position, seq(min(position), max(position), 20))) %>%
    group_by(bin) %>%
    mutate(coverage_bin = log10(sum(coverage)+1)) %>%
    ungroup() %>%
    distinct(bin, coverage_bin, .keep_all = T) %>%
    mutate(coverage_bin_scale = scales::rescale(coverage_bin, to = c(0,100)))
}

get_starts <- function(file, dataset, low, high){
  fileset <- paste0(dir,"/data/tss_data/")
  suppressMessages(vroom(paste0(fileset,file), col_names = F, num_threads = 8)) %>%
    dplyr::rename(position = X2, end = X3, coverage = X4) %>%
    dplyr::filter(position > low, position < high) %>%
    mutate(sample = dataset,
           cov_scale = scales::rescale(coverage, to=c(0,100))) %>%
    dplyr::select(sample, position, coverage, cov_scale) %>%
    mutate(bin = cut(position, seq(min(position), max(position), 20))) %>%
    group_by(bin) %>%
    mutate(coverage_bin = log10(sum(coverage)+1)) %>%
    ungroup() %>%
    distinct(bin, coverage_bin, .keep_all = T) %>%
    mutate(coverage_bin_scale = scales::rescale(coverage_bin, to = c(0,100)))
}


# load & tidy data ----

## positions of interest ====
start_mode <- 189511
end_mode   <- 193725

## binned 3´end sites ====
end_pcr_cdna_scale <- get_ends("tts_data_pychopper_auto_cutadapt_SSP_clipped/201210_PCB109_Ecoli/201210_PCB109_Ecoli_NOTEX_replicate1_full_length_all/201210_PCB109_Ecoli_NOTEX_replicate1_full_length_all.plus.bedgraph",
                               "201210_PCB109_Ecoli_NOTEX_replicate1_starts", start_mode, end_mode)
end_cdna_scale     <- get_ends("tts_data_pychopper_auto_cutadapt_SSP_clipped/210317_DCS109_Ecoli/210317_DCS109_Ecoli_NOTEX_replicate1_full_length_all/210317_DCS109_Ecoli_NOTEX_replicate1_full_length_all.plus.bedgraph",
                               "210317_DCS109_Ecoli_NOTEX_replicate1_starts", start_mode, end_mode)
end_dRNA_scale     <- get_ends("tts_data_pychopper_auto_cutadapt_SSP_clipped/190123_RNA001_Ecoli/190123_RNA001_Ecoli_TEX_replicate1/190123_RNA001_Ecoli_TEX_replicate1.plus.bedgraph",
                               "190123_RNA001_Ecoli_TEX_replicate1_starts", start_mode, end_mode)

end_set <- rbind(end_pcr_cdna_scale,
                 end_cdna_scale,
                 end_dRNA_scale) 

## binned 5´end sites ====
start_pcr_cdna_scale <- get_starts("tss_data_notrimming/201210_PCB109_Ecoli/201210_PCB109_Ecoli_NOTEX_replicate1/201210_PCB109_Ecoli_NOTEX_replicate1.plus.bedgraph",
                                   "201210_PCB109_Ecoli_NOTEX_replicate1_starts", start_mode, end_mode)
start_cdna_scale     <- get_starts("tss_data_notrimming/210317_DCS109_Ecoli/210317_DCS109_Ecoli_NOTEX_replicate1/210317_DCS109_Ecoli_NOTEX_replicate1.plus.bedgraph",
                                   "210317_DCS109_Ecoli_NOTEX_replicate1_starts", start_mode, end_mode)
start_dRNA_scale     <- get_starts("tss_data_notrimming/190123_RNA001_Ecoli/190123_RNA001_Ecoli_TEX_replicate1/190123_RNA001_Ecoli_TEX_replicate1.plus.bedgraph",
                                   "190123_RNA001_Ecoli_TEX_replicate1_starts", start_mode, end_mode)

start_set <- rbind(start_pcr_cdna_scale,
                   start_cdna_scale,
                   start_dRNA_scale) 

## normalized coverage ====
pcr_cdna_scale <- get_cov("coverage_data_notrimming_R/201210_PCB109_Ecoli/201210_PCB109_Ecoli_NOTEX_replicate1/201210_PCB109_Ecoli_NOTEX_replicate1.plus.bedgraph",
                          "201210_PCB109_Ecoli_NOTEX_replicate1", start_mode, end_mode)
cdna_scale     <- get_cov("coverage_data_notrimming_R/210317_DCS109_Ecoli/210317_DCS109_Ecoli_NOTEX_replicate1/210317_DCS109_Ecoli_NOTEX_replicate1.plus.bedgraph",
                          "210317_DCS109_Ecoli_NOTEX_replicate1", start_mode, end_mode)
dRNA_scale     <- get_cov("coverage_data_notrimming_R/190123_RNA001_Ecoli/190123_RNA001_Ecoli_TEX_replicate1/190123_RNA001_Ecoli_TEX_replicate1.plus.bedgraph",
                          "190123_RNA001_Ecoli_TEX_replicate1", start_mode, end_mode)

coverage_set <- rbind(pcr_cdna_scale,
                      cdna_scale,
                      dRNA_scale)

## gene annotation ====
ecoli_gff_filtered <- ecoli_gff %>%
  dplyr::filter(start_feature > start_mode,
                end_feature < end_mode)

## Primary 5´ends ====
dir <- here()
tss_frame_f <- vroom(paste0(dir, "/tables/tss_tables/tss_data_untrimmed.tsv")) %>%
  dplyr::filter(sample %in% c("RNA001_Ecoli_TEX_replicate1",
                               "PCB109_PCR12_Ecoli_NOTEX_replicate4",
                               "DCS109_Ecoli_NOTEX_replicate2"),
                TSS > start_mode,
                TSS < end_mode,
                TSS_type %in% c("primary"),
                strand == "+") %>%
  dplyr::select(sample, TSS, TSS_type) 

## Primary 3´ends ====
tts_frame_f <- vroom(paste0(dir, "/tables/tts_tables/tts_data_trimmed.tsv")) %>%
  dplyr::filter(sample %in% c("RNA001_Ecoli_TEX_replicate1",
                              "PCB109_PCR12_Ecoli_NOTEX_replicate4",
                              "DCS109_Ecoli_NOTEX_replicate2"),
                TTS > start_mode,
                TTS < end_mode,
                TTS_type %in% c("primary"),
                strand == "+") %>%
  dplyr::select(sample, TTS, TTS_type)

## transcriptional unit data ====
pcb_units <- vroom(paste0(dir, "/tables/operon_tables/PCB109_PCR12_Ecoli_NOTEX_replicate4.operons.tsv"))


## genome data ====
ecoli_gff_cds_selected <- ecoli_gff_cds %>%
  dplyr::filter(start_feature %in% start_mode:end_mode, end_feature %in% start_mode:end_mode) %>% 
  mutate(plot = "plot2")

## read in bam-like files to get single-read plots ====
map_table_f <- perform_analysis(dataset_choice = sets[2]) %>%
  dplyr::filter(strand == strand_feature,
                start %in% start_mode:end_mode, end %in% start_mode:end_mode) %>% 
  arrange(desc(width.x)) %>% 
  rownames_to_column("id_plot") %>%
  mutate(id_plot = as.numeric(id_plot)) %>%
  mutate(plot = "plot1") %>%
  arrange(id_plot)

# PLOTS ----

### Coverage, read end histograms and predicted primary ends - Fig. 3A #### 
ggplot(data = coverage_set, aes(x = position, y = cov_scale)) +
  facet_grid(rows = vars(sample)) +
  geom_segment(data = ecoli_gff_filtered %>% mutate(sample = NA), 
               aes(x = start_feature, xend = end_feature, y = 50, yend = 50,color = short_gene, size = 10)) +
  geom_area(alpha = 1, fill = "#595959",color = "#595959") +
  scale_fill_manual(values = ibm_colors) +
  theme_Publication_white() +
  geom_vline(data = tss_frame_f, aes(xintercept = TSS, color = TSS_type)) +
  geom_vline(data = tts_frame_f, aes(xintercept = TTS, color = TTS_type), linetype = "dashed") +
  geom_rect(data = start_set, aes(xmin = position,xmax = position+20,ymin = 0, ymax = coverage_bin_scale), color = NA, fill = "#33A3C2", alpha = 0.7) +
  geom_rect(data = end_set, aes(xmin = position,xmax = position+20,ymin = 0, ymax = coverage_bin_scale), color = NA, fill = "#0A1E5B", alpha = 0.7) +
  scale_y_continuous(limits = c(0,100), expand = c(0,1)) +
  scale_x_continuous(limits = c(start_mode, end_mode))

### Coverage,predicted primary ends, reads, units - Fig. 5A #### 

#### Single read coverage ####
ggplot(data = map_table_f) +
  geom_segment(aes(y = id_plot, yend = id_plot, x = start, xend = end, alpha = aligned_reads, color = strand), size = .5) +
  scale_x_continuous(limits = c(start_mode, end_mode)) + 
  geom_segment(data = ecoli_gff_cds_selected,
               aes(x = start_feature, xend = end_feature, y = 1, yend = 1, color = gene), size = 4) +
  facet_grid(rows = vars(plot)) +
  theme_Publication_white()  

#### Transcriptional units ####
pcb_units %>%
  dplyr::filter(operon_seen >= 10,
                start %in% start_mode:(end_mode+1000),
                end %in% start_mode:(end_mode+1000)) %>%
  ggplot(data = .) +
    geom_path(aes(x = start, y = operon_name, group = operon_name), size = 1) +
    geom_rect(aes(xmin=(start), xmax=(end), ymin=as.integer(operon_name)-0.33, ymax=as.integer(operon_name)+0.33, fill = perc_seen), 
              color = "black", alpha = 1,size = 0.5, linetype = 1) +
    coord_cartesian(xlim = c(start_mode, end_mode+1000)) +
    scale_fill_gradientn(colours = brewer.pal(name = "YlGnBu", n = 9), limit = c(0,100))  +
    theme_Publication_white() 

# > number of units
pcb_units %>%
  dplyr::filter(operon_seen >= 10,
                start %in% start_mode:(end_mode+1000),
                end %in% start_mode:(end_mode+1000)) %>%
  arrange(desc(start)) %>% 
  group_by(operon_name) %>%
  summarise(calc= min(operon_seen)) %>%
  arrange(desc(as.integer(operon_name))) 
