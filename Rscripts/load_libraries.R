# set libraries, datasets and standard plotting functions that are used during the analysis

## libraries ====
packages <- c("ggeconodist", "tidyverse", "here", "ggthemes", "gganimate","writexl","seqinr",
              "colorblindr","rcartocolor","ChIPpeakAnno","tictoc","patchwork","psych","bioanalyzeR",
              "data.table", "ggExtra", "Rsamtools", "GenomicAlignments", "UpSetR","vroom", 
              "Rsubread", "ape", "DT", "ggpubr", "ggridges", "ggsci","LncFinder","RColorBrewer",
              "CoverageView", "gghalves", "pryr", "fst", "R.utils", "readxl", "ggseqlogo")

invisible(lapply(packages, require, character.only = TRUE))

## plotting theme ====
theme_Publication_white <- function(base_size=10) {
  (theme_foundation(base_size=base_size, base_family="Helvetica")
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(3), hjust = 0.5),
           text = element_text(color = "black"),
           panel.background = element_rect(colour = NA, fill = "white"),
           plot.background = element_rect(colour = NA, fill = "white"),
           panel.border = element_rect(colour = "black"),
           axis.title = element_text(face = "bold",size = rel(1.5)),
           axis.title.y = element_text(angle=90,vjust =2, size = rel(1.2)),
           axis.title.x = element_text(vjust = -0.2, size = rel(1.2)),
           axis.text = element_text(size = rel(1.2)), 
           axis.line = element_line(colour="black"),
           axis.ticks = element_line(),
           panel.grid.major = element_line(colour="grey80"),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(color = NA, fill = "white"),
           legend.position = "bottom",
           legend.background= element_rect(color = NA, fill = "white"),
           legend.direction = "horizontal",
           legend.key.size= unit(0.2, "cm"),
           legend.spacing = unit(0, "cm"),
           legend.text = element_text(color = "black", size = rel(1.2)),
           legend.title = element_text(face="italic", color = "black"),
           plot.margin=unit(c(10,5,5,5),"mm"),
           strip.background=element_rect(colour="grey90",fill="grey70"),
           strip.text = element_text(face="bold")
   ))
}

## standard plotting function ====
raw_reads_plotting <- function(mydf, myx, myy, myfill, mypalette){
  ggplot(data = mydf, aes(x = {{myx}}, y = {{myy}}, fill = {{myfill}})) +
    theme_Publication_white() +
    ylab("") +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(linetype = "dashed", color = "black")) +
    scale_fill_manual(values = mypalette)
}

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

## color palettes ====

### standard ####
cbf1 <- c("#EFEAFF","#F5AAA3","#CFCFCF", "#F6B2FB", "#ABC2DA")

### highlight ####
cbf1_high <- c("#EFEAFF", "#ABC2DA","#F6B2FB","#648FFF")

