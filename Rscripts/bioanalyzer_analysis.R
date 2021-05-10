# >> Bioanalyzer analysis << #

# load libraries ----
library(here)
source(here("Rscripts/load_libraries.R"))

# data ---- 

## Bioanalyzer data - polyA tailing ==== 
### read in data and change dataset names ####
index            <- data.table(sample.index = 1:11, 
                               dataset = c("raw_RNA_replicate1", "raw_RNA_replicate2", "raw_RNA_replicate3", "raw_RNA_replicate4", 
                                           "polyA_RNA_replicate1","polyA_RNA_replicate2","polyA_RNA_replicate3","polyA_RNA_replicate4",
                                           "empty", "empty", "ladder"))

RIN_polyA_sample <- bioanalyzeR::read.electrophoresis(here("data/bioanalyzer_data/RIN/RIN.xml")) 

### make sample table ####
RIN_polyA_sample_table <- RIN_polyA_sample$data %>%
  as_tibble() %>%
  left_join(index) %>%
  dplyr::filter(!dataset %in% c("ladder", "empty")) %>%
  dplyr::filter(!is.na(fluorescence)) %>%
  group_by(dataset) %>%
  dplyr::filter(!is.na(fluorescence)) %>%
  dplyr::mutate(f_scale = scales::rescale(fluorescence, to = c(0,100))) %>%
  ungroup() %>%
  dplyr::mutate(group = ifelse(str_sub(dataset,1,3) == "raw", "raw", "polyA"),
                rep = paste0("Replicate ",str_sub(dataset,-1,-1)))

### extract ladder peaks #### 
ladder <- RIN_polyA_sample$peaks %>% 
  filter(sample.index == 11,
         peak.observations == "Ladder Peak") 

## Bioanalyzer data - PCR effect ==== 
sample_names2 <- data.table(sample.index = 1:12, 
                            dataset = c("cDNA_NOTEX_replicate1_PCR11", 
                                        "cDNA_NOTEX_replicate1_PCR12", "cDNA_NOTEX_replicate1_PCR12_dil", 
                                        "cDNA_NOTEX_replicate1_PCR13","cDNA_NOTEX_replicate1_PCR13_dil",
                                        "cDNA_NOTEX_replicate1_PCR14","cDNA_NOTEX_replicate1_PCR14_dil",
                                        "cDNA_NOTEX_replicate1_PCR15_dil","cDNA_NOTEX_replicate1_PCR15",
                                        "cDNA_NOTEX_replicate2_PCR13", "cDNA_NOTEX_replicate2_PCR13_dil",
                                        "Ladder"))

bio2 <- bioanalyzeR::read.electrophoresis(here("data/bioanalyzer_data/PCR_effect/PCR.xml"))

peaks <- bio2$data %>%
  as_tibble() %>%
  left_join(sample_names2) %>%
  dplyr::filter(!dataset %in% c("Ladder", 
                                "cDNA_NOTEX_replicate1_PCR11", 
                                "cDNA_NOTEX_replicate2_PCR13"),
                str_sub(dataset, -3,-1) != "dil") %>%
  #left_join(marker, by = sample.index)
  group_by(dataset) %>%
  dplyr::filter(!is.na(fluorescence)) %>%
  dplyr::filter(aligned.time >= 60 & aligned.time <= 110) %>%
  dplyr::mutate(fluorescence2 = scales::rescale(fluorescence, to = c(0,100))) %>%
  ungroup() %>%
  dplyr::mutate(group_tex = ifelse(str_sub(dataset,6,7) != "TE", "NOTEX", "TEX"),
                dataset = str_sub(dataset, nchar(dataset)-1, nchar(dataset)))

sizes <- bio2$peaks %>% 
  filter(sample.index == 12,
         peak.observations == "Ladder Peak") 

# plot ----
## RIN & effect of polyA tailing - Supplementary Fig. 2A ====
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

## PCR effect - Supplementary Fig. 2B ====
ggplot(data = peaks, aes(x = aligned.time, y = fluorescence2)) +
  geom_line() +
  facet_grid(rows = vars(dataset)) +
  geom_area(alpha = 0.75, fill = "#595959") +
  scale_x_log10(limits = c(60,110),breaks = c(22.5, (sizes$aligned.time)), labels = c(25, (sizes$length))/1000, expand = c(0,0)) +
  theme_Publication_white() +
  ylab("Relative fluorescence") +
  xlab("Size (kb, bioanalyzer scale") +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linetype = "dashed", color = "black")) 
   

