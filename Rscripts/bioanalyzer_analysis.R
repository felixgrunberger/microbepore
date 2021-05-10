# >> Bioanalyzer analysis << #

# load libraries ----
library(here)
source(here("Rscripts/load_libraries.R"))

# functions ----
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

# data ---- 
dir <- "/Volumes/EX_SSD/"

## genome annotation, rRNA loci ====
gff_rRNA <- read_in_gff(paste0(dir,"data/genome_data/NC_000913.3.gff3")) %>% 
  dplyr::filter(type == "rRNA", locus_name %in% c("16S", "23S")) %>%
  group_by(locus_name) %>%
  summarise(mean_v = mean(width))

## Bioanalyzer data ==== 
### read in data and change dataset names ####
index            <- data.table(sample.index = 1:11, 
                               dataset = c("raw_RNA_replicate1", "raw_RNA_replicate2", "raw_RNA_replicate3", "raw_RNA_replicate4", 
                                           "polyA_RNA_replicate1","polyA_RNA_replicate2","polyA_RNA_replicate3","polyA_RNA_replicate4",
                                           "empty", "empty", "ladder"))

RIN_polyA_sample <- bioanalyzeR::read.electrophoresis(here("data/bioanalyzer_data/2020-11-19/201119_1_rna_Prokaryote Total RNA Nano_DE72902515_2020-11-19_11-43-47.xml")) 

RIN_polyA_sample_table <- RIN_polyA_sample$data %>%
  as_tibble() %>%
  left_join(index) %>%
  dplyr::filter(!dataset %in% c("ladder", "empty")) %>%
  dplyr::filter(!is.na(fluorescence)) %>%
  group_by(dataset) %>%
  dplyr::mutate(f_scale = scales::rescale(fluorescence, to = c(0,100))) %>%
  ungroup() %>%
  dplyr::mutate(group = ifelse(str_sub(dataset,1,3) == "raw", "raw", "polyA"),
                rep = paste0("Replicate ",str_sub(dataset,-1,-1)))

### extract ladder peaks
ladder <- RIN_polyA_sample$peaks %>% 
  filter(sample.index == 11,
         peak.observations == "Ladder Peak") 

# plot ----
## RIN & effect of polyA tailining - Supplementary Fig. 2A ====
ggplot(data = RIN_polyA_sample_table, 
       aes(x = aligned.time, y = f_scale, 
           fill = group)) +
  facet_grid(rows = vars(rep)) +
  geom_area(alpha = 0.75) +
  geom_line(color = "black", size = 0.5) +
  scale_x_log10(limits = c(30,55),breaks = c(22.5, (ladder$aligned.time)), labels = c(25, (ladder$length))/1000, expand = c(0,0)) +
  scale_y_continuous(limits = c(0,100), expand = c(0,0)) +
  scale_fill_manual(values = c("#648FFF", "#595959")) +
  theme_Publication_white()  +
  ylab("Relative fluorescence") +
  xlab("Size (kb, bioanalyzer scale") +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linetype = "dashed", color = "black")) 





bio <- bioanalyzeR::read.electrophoresis(here("data/bioanalyzer_data/2020-12-08/201208_cDNA_PCR_High Sensitivity DNA Assay_DE72902515_2020-12-08_11-53-48.xml"))
bio$samples
bio$samples$sample.name <- c("cDNA_NOTEX_replicate1_PCR15", "cDNA_NOTEX_replicate2_PCR15", 
                             "cDNA_TEX_replicate1_PCR15", "cDNA_TEX_replicate2_PCR15",
                             "cDNA_NOTEX_replicate1", "cDNA_NOTEX_replicate2", 
                             "cDNA_TEX_replicate1", "cDNA_TEX_replicate2",
                             "Ladder")

dilutions_dep <- bio
dilutions_dep$samples$dep.type <- dilutions_dep$samples$sample.name 

dilutions_dep$samples <- dilutions_dep$samples %>%
  mutate(dep.type = as.factor(ifelse(grepl("NOTEX", dep.type), "NOTEX", "TEX")))
dilutions_dep$samples$dep.type <- ordered(dilutions_dep$samples$dep.type, 
                                          levels = c("TEX", "NOTEX"))


dilutions_dep$samples$PCR <- as.factor(ifelse(grepl("PCR", dilutions_dep$samples$sample.name), 
                                              "PCR","noPCR"))

sizes <- (subset(dilutions_dep, grepl("Ladder", sample.name)))$peaks %>% 
  filter(peak.observations == "Ladder Peak") 


dilutions_dep <- subset(dilutions_dep, !grepl("Ladder|_2|term|102|80", sample.name))


# other attempt
sample_names <- data.table(sample.index = 1:9, 
                           dataset = c("cDNA_NOTEX_replicate1_PCR15", "cDNA_NOTEX_replicate2_PCR15", 
                                       "cDNA_TEX_replicate1_PCR15", "cDNA_TEX_replicate2_PCR15",
                                       "cDNA_NOTEX_replicate1", "cDNA_NOTEX_replicate2", 
                                       "cDNA_TEX_replicate1", "cDNA_TEX_replicate2",
                                       "Ladder"))

sample_names2 <- data.table(sample.index = 1:12, 
                            dataset = c("cDNA_NOTEX_replicate1_PCR11", 
                                        "cDNA_NOTEX_replicate1_PCR12", "cDNA_NOTEX_replicate1_PCR12_dil", 
                                        "cDNA_NOTEX_replicate1_PCR13","cDNA_NOTEX_replicate1_PCR13_dil",
                                        "cDNA_NOTEX_replicate1_PCR14","cDNA_NOTEX_replicate1_PCR14_dil",
                                        "cDNA_NOTEX_replicate1_PCR15_dil","cDNA_NOTEX_replicate1_PCR15",
                                        "cDNA_NOTEX_replicate2_PCR13", "cDNA_NOTEX_replicate2_PCR13_dil",
                                        "Ladder"))

glimpse(bio)
bio <- bioanalyzeR::read.electrophoresis(here("data/bioanalyzer_data/2020-12-08/201208_cDNA_PCR_High Sensitivity DNA Assay_DE72902515_2020-12-08_11-53-48.xml"))
peaks <- bio$data %>%
  as_tibble() %>%
  left_join(sample_names) %>%
  dplyr::filter(!dataset %in% c("Ladder"),
                dataset %in% c("cDNA_NOTEX_replicate1_PCR15", "cDNA_NOTEX_replicate1",
                               "cDNA_NOTEX_replicate2_PCR15", "cDNA_NOTEX_replicate2")) %>%
  #left_join(marker, by = sample.index)
  group_by(dataset) %>%
  dplyr::filter(!is.na(fluorescence)) %>%
  #mutate(fluorescence = fluorescence + 10) %>%
  #dplyr::mutate(fluorescence = ifelse(time < 55 | time > 105, 1, fluorescence)) %>%
  dplyr::filter(aligned.time >= 60 & aligned.time <= 112) %>%
  dplyr::mutate(fluorescence2 = scales::rescale(fluorescence, to = c(0,100))) %>%
  ungroup() %>%
  rowwise() %>%
  dplyr::mutate(group = paste0(str_split_fixed(dataset, "_", 4)[,c(1:3)], collapse = "_"),
                group2 = ifelse(str_detect(dataset, "PCR"), "PCR", "raw"))

sizes <- bio$peaks %>% 
  filter(sample.index == 9,
         peak.observations == "Ladder Peak") 

# plot CDNA input samples vs PCR15 cycles
pdf(here("figures/R_plots/210215_cDNA_vs_15.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = peaks, aes(x = aligned.time, y = fluorescence2, 
                         fill = group2)) +
  geom_line() +
  facet_grid(rows = vars(group)) +
  geom_area(alpha = 0.75) +
  scale_x_log10(limits = c(60,112),breaks = c(22.5, (sizes$aligned.time)), labels = c(25, (sizes$length)), expand = c(0,0)) +
  theme_Publication_white() 
  scale_fill_manual(values = ibm_colors[c(4,2)])
dev.off()  



bio2 <- bioanalyzeR::read.electrophoresis(here("data/bioanalyzer_data/2020-12-09/201209_cDNA_PCR2_High Sensitivity DNA Assay_DE72902515_2020-12-09_19-22-09.xml"))

peaks <- bio2$data %>%
  as_tibble() %>%
  left_join(sample_names2) %>%
  dplyr::filter(!dataset %in% c("Ladder", "cDNA_NOTEX_replicate1_PCR11", "cDNA_NOTEX_replicate2_PCR13"),
                str_sub(dataset, -3,-1) != "dil") %>%
  #left_join(marker, by = sample.index)
  group_by(dataset) %>%
  dplyr::filter(!is.na(fluorescence)) %>%
  #mutate(fluorescence = fluorescence + 10) %>%
  #dplyr::mutate(fluorescence = ifelse(time < 55 | time > 105, 1, fluorescence)) %>%
  dplyr::filter(aligned.time >= 60 & aligned.time <= 110) %>%
  dplyr::mutate(fluorescence2 = scales::rescale(fluorescence, to = c(0,100))) %>%
  ungroup() %>%
  dplyr::mutate(group_tex = ifelse(str_sub(dataset,6,7) != "TE", "NOTEX", "TEX"))

sizes <- bio2$peaks %>% 
  filter(sample.index == 12,
         peak.observations == "Ladder Peak") 

# plot effect of more cycles
pdf(here("figures/R_plots/210215_cDNA_12_13_14_15.pdf"), 
    width = 7, height = 7, paper = "special", onefile=FALSE)
ggplot(data = peaks, aes(x = aligned.time, y = fluorescence2, 
                         fill = dataset)) +
  geom_line() +
  facet_grid(rows = vars(dataset)) +
  geom_area(alpha = 0.75) +
  scale_x_log10(limits = c(60,110),breaks = c(22.5, (sizes$aligned.time)), labels = c(25, (sizes$length)), expand = c(0,0)) +
  theme_Publication_white() +
  viridis::scale_color_viridis(option = "cividis", discrete = T) +
  viridis::scale_fill_viridis(option = "cividis", discrete = T) 
dev.off()

# depletion and polyA-tailing



