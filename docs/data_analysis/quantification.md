---
layout: page
title: Gene abundance estimation
parent: Data analysis
nav_order: 5
---

## Gene abundance estimation   
- To estimate gene abundances `salmon` (v.1.4.0) was applied in alignment-based mode as described in https://salmon.readthedocs.io/en/latest/salmon.html#quantifying-in-alignment-based-mode  
- For installation you can follow the explanations listed at https://combine-lab.github.io/salmon/getting_started/ (installs salmon in own conda environment)

```bash
input=microbepore/data/mapped/raw # input directory with all remapped files

for file in ${input}/*/*/*remapped.sorted.bam
do
  
  # folder and filenames
  f_ex=${file##*/}
  foldername=$(echo ${f_ex} | cut -d"_" -f 1,2,3)
  filename=${f_ex%%.*}
  
  # create dir for quantification using salmon in alignment-based mode (e.g. used in conda environment)
  mkdir microbepore/data/salmon
  mkdir microbepore/data/salmon/${foldername}
  output=microbepore/data/salmon/${foldername}/${filename}
  mkdir ${output}
  
  conda activate salmon # activate conda environment
  
  # use salmon in alignment-based mode
  salmon quant \
  -t ${transcripts} \
  -l A \
  -a ${file} \
  -o ${output} \
  --threads 8 
  
  conda deactivate
  
done
```

- Transcripts per million (TPM) were re-calculated using the salmon-computed effective transcript length, after dropping reads mapping to rRNAs, that are variable between non-depleted and depleted RNA sets (compare custom [Rscripts/salmon_analysis.R](https://github.com/felixgrunberger/microbepore/blob/master/Rscripts/salmon_analysis.R)).   
- Calculation of custom TPMs is performed in the `modify_salmon_output` function:  

```r
modify_salmon_output <- function(input, method){
  # read in salmon output
  suppressMessages(vroom(input, num_threads = 8)) %>%  
    # sort by Number of reads mapping
    arrange(desc(NumReads)) %>% # read in salmon output
    # rename columns
    mutate(counts = NumReads, gene = Name, salmon_tpm = TPM) %>%
    # drop all rRNAs
    dplyr::filter(!gene %in% names_rRNA) %>%
    # select columns
    select(counts, gene, salmon_tpm,EffectiveLength) %>%
    # rename genes
    mutate(gene = str_split_fixed(gene,"-",2)[,2]) %>%
    # sort by number of counts
    arrange(desc(counts)) %>%
    # add gene information (including with to table)
    left_join(ecoli_gff, by = "gene") %>%
    # make RNA type categories
    mutate(type_fine = ifelse(type == "rRNA", as.character(locus_name), as.character(type)),
           type_fine = ifelse(type == "tRNA", "ncRNA", as.character(type_fine))) %>%
    # calculate RPK from salmon-calculated effectiveLength, TPM(hand) = RPK/sum(RPK)/1000000, additionally calc RPKM       
    mutate(sample = method,
           rpk = (counts/EffectiveLength*1000),
           TPM_hand = rpk/(sum(rpk, na.rm = T)/1000000),
           rpkm = (counts/(sum(counts)/1000000))/width*1000) %>%
    # TPM_hand is the value used        
    dplyr::select(gene, type_fine, counts, sample,rpkm, TPM_hand, salmon_tpm)
}
```



