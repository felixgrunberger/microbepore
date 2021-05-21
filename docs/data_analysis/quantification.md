---
layout: page
title: Transcript quantification
parent: Data analysis
nav_order: 4
---

## Transcript quantification   
For quantification, we used salmon in mapping-based mode.   

```scala
# load functions ----
modify_salmon_output <- function(input, method){
  suppressMessages(vroom(input, num_threads = 8)) %>%
    arrange(desc(NumReads)) %>%
    mutate(counts = NumReads, gene = Name, salmon_tpm = TPM) %>%
    dplyr::filter(!gene %in% c(name_16S,name_23S,name_5S)) %>%
    select(counts, gene, salmon_tpm,EffectiveLength) %>%
    mutate(gene = str_split_fixed(gene,"-",2)[,2]) %>%
    arrange(desc(counts)) %>%
    left_join(gff_table, by = "gene") %>%
    mutate(type_fine = ifelse(type == "rRNA", as.character(locus_name), as.character(type)),
           type_fine = ifelse(type == "tRNA", "ncRNA", as.character(type_fine))) %>%
    mutate(method = method,
           rpk = (counts/EffectiveLength*1000),
           TPM_hand = rpk/(sum(rpk, na.rm = T)/1000000),
           rpkm = (counts/(sum(counts)/1000000))/width*1000) %>%
    dplyr::select(gene, type_fine, counts, method,rpkm, TPM_hand, salmon_tpm)
}
```

```scala
# load & tidy data ----
## salmon data ====
### full-length, polyA-trimmed, clipped ####
files          <- list.files(paste0(dir,"/data/salmon_data_notrimming/"), recursive = T, pattern = "quant.sf")
dataset_names  <- str_split_fixed(str_split_fixed(files, "\\/", n = 3)[,2],"_fu",2)[,1]
salmon_frame   <- data.table()
partial_frame  <- data.table()

### raw
for(i in seq_along(dataset_names)){
  
  print(paste0("file number ", i, " of ", length(dataset_names)))
  tic("comp salmon")
  partial_frame <- modify_salmon_output(paste0(dir,"/data/salmon_data_notrimming/",files[i]), paste0(dataset_names[i]))
  toc()
  
  salmon_frame <- rbind(salmon_frame,partial_frame)
}

```



{% capture Rock %}
{% highlight R linenos %}
# load libraries ----
library(here)
source(here("Rscripts/load_libraries.R"))
{% endhighlight %}
{% endcapture %}
{% include fix_linenos.html code=Rock %}




