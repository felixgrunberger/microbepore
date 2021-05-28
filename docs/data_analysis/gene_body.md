---
layout: page
title: Gene body coverage
parent: Data analysis
nav_order: 8
---

## Gene body coverage analysis  
{: .no_toc }
____
<details open markdown="block">
  <summary>
    Table of contents
  </summary>
  {: .text-delta }
1. TOC
{:toc}
</details> 
____

To assess the impact of trimmings on gene body coverage, a coverage meta-analysis was performed   

### Preparation of ONT-annotated transcript files    
- First, a transcript file was created for all genes with an ONT-annotated primary 5´ and 3´ end (see previous section) 

### Calculation of coverage using `bedtools coverage`  
- Based on this, strand-specific coverage files were created from the BAM files   

```bash
input=microbepore/data/mapped

# calculate coverage over transcripts with TSS and TTS | for pychopper auto > cutadapt > clipped or RAW 
for file in ${input}/trimmed/*/*/*clipped.sorted.bam # ||  for file in ${input}/raw/*/*/*.sorted.bam
do 
  filename_extended=${file##*/}
  keyword=$(echo $filename_extended | cut -d"." -f 2)
  foldername=$(echo $filename_extended | cut -d"_" -f 1,2,3)
  filename=${filename_extended%%.*}

  # mk dirs
  mkdir microbepore/data/coverage/trimmed
  mkdir microbepore/data/coverage/trimmed/${foldername}
  output=microbepore/data/coverage/trimmed/${foldername}/${filename}
  mkdir ${output}

  # calc coverage
  samtools view -F 16 -o temp.sorted.bam ${file} 
  bedtools coverage \
  -d \
  -a ${dir}/data/bed/transcripts.plus.bedgraph \ # bed file of genes with annotated 5´and 3´end
  -b temp.sorted.bam \
  > ${output}/${filename}.plus.coverage
  
  samtools view -f 16 -o temp.sorted.bam ${file} 
  bedtools coverage \
  -d \
  -a ${dir}/data/bed/transcripts.minus.bedgraph \ # bed file of genes with annotated 5´and 3´end
  -b temp.sorted.bam \
  > ${output}/${filename}.minus.coverage
done
```

### Downstream R analysis   
Coverage analysis performed using a custom R script.  