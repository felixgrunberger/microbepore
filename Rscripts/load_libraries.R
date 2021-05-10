# set libraries and standard plotting functions that are used during the analysis

#...................................libraries
packages <- c("ggeconodist", "tidyverse", "here", "ggthemes", "gganimate","writexl","seqinr",
              "colorblindr","rcartocolor","ChIPpeakAnno","tictoc","patchwork","psych","bioanalyzeR",
              "data.table", "ggExtra", "Rsamtools", "GenomicAlignments", "UpSetR","vroom", 
              "Rsubread", "ape", "DT", "ggpubr", "ggridges", "ggsci","LncFinder","RColorBrewer",
              "CoverageView", "gghalves", "pryr", "fst", "R.utils", "readxl", "ggseqlogo")

invisible(lapply(packages, require, character.only = TRUE))

#...................................plotting theme
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
