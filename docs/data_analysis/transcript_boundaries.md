---
layout: page
title: Detection of transcript boundaries
parent: Data analysis
nav_order: 7
---

## Detection of transcript boundaries   
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

The determination of enriched 5´and 3´ends was carried out in the same way, but independently of each other, and is briefly explained in the following:  
- First, strand-specific read ends in bedgraph format were created from BAM files using [`bedtools genomecov`](https://bedtools.readthedocs.io/en/latest/) (-5 or -3 option, -bga)  
- Next, the previously published [`Termseq_peaks`](https://pypi.org/project/termseq-peaks/) script was used to call peaks for each sample individually without including replicates (https://github.com/NICHD-BSPC/termseq-peaks) 
- This script is based on `scipy.signal.find_peaks`, which is running in the background of `Termseq_peaks` with lenient parameters (prominence=(None,None), width=(1,None), rel_height=0.75)
- However, we deliberately used `Termseq_peaks` since its ability to include replicates by applying an Irreproducible Discovery Rate method which can be applied to future studies 
- For end detection, only the leniently called peaks in the narrowPeak file were used after adding the number of counts for each position using `bedtools intersect`. 

### 5´end detection  
5´end peak calling was performed in the following way:   

```bash
input=microbepore/data/mapped

# perform tss detection for pychopper auto > cutadapt_polyA > SSP-cutadapt > clipped  or for raw mapped reads
for file in ${input}/trimmed/*/*/*clipped.sorted.bam # ||  for file in ${input}/raw/*/*/*.sorted.bam
do 
  # file and folder names
  filename_extended=${file##*/}
  keyword=$(echo $filename_extended | cut -d"." -f 2)
  foldername=$(echo $filename_extended | cut -d"_" -f 1,2,3)
  filename=${filename_extended%%.*}
  
  # make directories
  mkdir microbepore/data/tss/trimmed
  mkdir microbepore/data/tss/trimmed/${foldername}
  output=microbepore/data/tss/trimmed/${foldername}/${filename}
  mkdir ${output}

  # step 1: calculate 5´positions for plus and minus strand
  bedtools genomecov \
    -ibam ${file} \
    -bga \
    -5 \
    -strand + > ${output}/${filename}.plus.bedgraph
  
  bedtools genomecov \
    -ibam ${file} \
    -bga \
    -5 \
    -strand - > ${output}/${filename}.minus.bedgraph
    
  # step 2: termseq peaks
  termseq_peaks ${output}/${filename}.plus.bedgraph ${output}/${filename}.plus.bedgraph --peaks ${output}/${filename}.plus.peaks --strand +
  termseq_peaks ${output}/${filename}.minus.bedgraph ${output}/${filename}.minus.bedgraph --peaks ${output}/${filename}.minus.peaks --strand -
    
  # step 3: add coverage information
  bedtools intersect \
    -wao \
    -a ${output}/${filename}.plus.peaks.oracle.narrowPeak \
    -b ${output}/${filename}.plus.bedgraph \
    > ${output}/${filename}.plus.peaks.oracle.narrowPeak.counts
    
  bedtools intersect \
    -wao \
    -a ${output}/${filename}.minus.peaks.oracle.narrowPeak \
    -b ${output}/${filename}.minus.bedgraph \
    > ${output}/${filename}.minus.peaks.oracle.narrowPeak.counts

done
```


### 3´end detection  
3´end peak calling was performed in the following way:   

```bash
input=microbepore/data/mapped

# perform tts detection for pychopper auto > cutadapt_polyA > SSP-cutadapt > clipped  or for raw mapped reads
for file in ${input}/trimmed/*/*/*clipped.sorted.bam # ||  for file in ${input}/raw/*/*/*.sorted.bam
do 
  filename_extended=${file##*/}
  keyword=$(echo $filename_extended | cut -d"." -f 2)
  foldername=$(echo $filename_extended | cut -d"_" -f 1,2,3)
  filename=${filename_extended%%.*}
    
  echo ${filename}

  mkdir microbepore/data/tts/trimmed
  mkdir microbepore/data/tts/trimmed
  mkdir microbepore/data/tts/trimmed/${foldername}
  output=microbepore/data/tts/trimmed/${foldername}/${filename}
  mkdir ${output}

  # step 1: calculate 3´positions for plus and minus strand
  bedtools genomecov \
    -ibam ${file} \
    -bga \
    -3 \
    -strand + > ${output}/${filename}.plus.bedgraph
  
  bedtools genomecov \
    -ibam ${file} \
    -bga \
    -3 \
    -strand - > ${output}/${filename}.minus.bedgraph
    
  # step 2: termseq peaks
  termseq_peaks ${output}/${filename}.plus.bedgraph ${output}/${filename}.plus.bedgraph --peaks ${output}/${filename}.plus.peaks --strand +
    termseq_peaks ${output}/${filename}.minus.bedgraph ${output}/${filename}.minus.bedgraph --peaks ${output}/${filename}.minus.peaks --strand -
    
  # step 3: add coverage information
  bedtools intersect \
    -wao \
    -a ${output}/${filename}.plus.peaks.oracle.narrowPeak \
    -b ${output}/${filename}.plus.bedgraph \
    > ${output}/${filename}.plus.peaks.oracle.narrowPeak.counts
    
  bedtools intersect \
    -wao \
    -a ${output}/${filename}.minus.peaks.oracle.narrowPeak \
    -b ${output}/${filename}.minus.bedgraph \
    > ${output}/${filename}.minus.peaks.oracle.narrowPeak.counts

done
```
