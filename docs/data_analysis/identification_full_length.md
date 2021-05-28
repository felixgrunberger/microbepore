---
layout: page
title: Identification and trimming of full-length reads 
parent: Data analysis
nav_order: 4
---

## Trimming of reads using `pychopper`, `cutadapt` & `samclip`   
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

### Identification of full-length reads using [`pychopper`](https://github.com/nanoporetech/pychopper)  
- Full-length cDNA reads containing SSP and VNP primers in the correct orientation were identified using `pychopper` (v.2.5.0) with standard parameters using the default pHMM backend and autotuned cutoff parameters estimated from subsampled data  
- Save output in *pychopper folder*.

```bash
# files
input=microbepore/data/FASTQ/normal # input directory with all merged FASTQ files, 1 for each barcode or single DRS run

# perform pychopper for all cDNA and (PCR)-cDNA files
for file in ${input}/*/*.fastq
do 

  # folder and filenames
  f_ex=${file##*/}
  foldername=$(echo $f_ex | cut -d"_" -f 1,2,3)
  filename=${f_ex%%.*}
  
  # make directories
  mkdir microbepore/data/pychopper/normal
  mkdir microbepore/data/pychopper/normal/${foldername}
  output=microbepore/data/pychopper/normal/${foldername}/${filename}
  mkdir ${output}

  # perform pychopper using precomputed q
  cdna_classifier.py \
  -r ${output}/${filename}_report.pdf \
  -t 8 \
  -u ${output}/${filename}_unclassified.fastq \
  -w ${output}/${filename}_rescued.fastq \
  -S ${output}/${filename}_stats.txt \
  $file \
  ${output}/${filename}_full_length_output.fastq
done
```

After a first round, a second round of `pychopper` was applied to the unclassified direct cDNA reads with DCS-specific read rescue enabled.  

```bash
# files
input=microbepore/data/pychopper/normal # input directory with all merged FASTQ files, 1 for each barcode or single DRS run

# perform pychopper using the -x rescue option for DCS files
for file in ${input}/*unclassified.fastq # only use unclassified reads from first round as input
do 

  # folder and filenames
  filename_extended=${file##*/}
  foldername=$(echo $filename_extended | cut -d"_" -f 1,2,3)
  filename=${filename_extended%%.*}
  
  # make directories
  mkdir ${dir}/data/pychopper/rescued
  mkdir ${dir}/data/pychopper/rescued/${foldername}
  output=microbepore/data/pychopper/rescued/${foldername}/${filename}
  mkdir ${output}
  
  # perfrom pychopper using -X option for native cDNA datasets
  cdna_classifier.py \
  -r ${output}/${filename}_report.pdf \
  -t 8 \
  -x rescue \
  -u ${output}/${filename}_unclassified.fastq \
  -w ${output}/${filename}_rescued.fastq \
  -S ${output}/${filename}_stats.txt \
  $file \
  ${output}/${filename}_full_length_output.fastq
done
```

Reads from rescued and normal folders were merged and used for subsequent steps. 

```bash
# files
input=microbepore/data/pychopper/

# merge all full-length and rescued reads as full-length
for file in ${input}/normal/*/*/*full_length_output.fastq # both normal and rescued folders
do 
  filename_extended=${file##*/}
  foldername=$(echo $filename_extended | cut -d"_" -f 1,2,3)
  filename=$(echo $filename_extended | cut -d"_" -f 1,2,3,4,5)

  keyword=$(echo $foldername | cut -d"_" -f 2) # get libary kit ID
  
  mkdir microbepore/data/FASTQ/full_length
  mkdir microbepore/data/FASTQ/full_length/${foldername}
  output=microbepore/data/FASTQ/full_length/${foldername}/${filename}
  mkdir ${output}
  
  if [[ $keyword =~ "PCB109" ]]; then
    cat $file ${input}/normal/${foldername}/${filename}/${filename}_rescued.fastq > ${output}/${filename}_full_length_all.fastq
  elif [[ $keyword =~ "DCS109" ]]; then
    cat $file ${input}/normal/${foldername}/${filename}/${filename}_rescued.fastq
    ${input}/rescued/${foldername}/${filename}_unclassified/${filename}_unclassified_full_length_output.fastq
    ${input}/rescued/${foldername}/${filename}_unclassified/${filename}_unclassified_rescued.fastq > ${output}/${filename}_full_length_all.fastq
  fi
done
```

For easier handling in the subsequent steps, DRS FASTQ files are also moved to the microbepore/data/FASTQ/full_length folder and adding *_full_length_all* to the filename.  

### Remove polyA-tails using [`cutadapt`](https://cutadapt.readthedocs.io/en/stable/)  
- To evaluate the influence of different trimming approaches on the accuracy of transcript boundary analysis, we applied additional 5´ and 3´ trimming steps using `cutadapt` v3.2   
- To this end, polyA sequences were removed from the 3´ends:   

```bash
# files
input=microbepore/data/FASTQ/full_length # input directory with all merged FASTQ files, 1 for each barcode or single DRS run

for file in ${input}/*/*/*_full_length_all.fastq
do 

  # folder and filenames
  filename_extended=${file##*/}
  keyword=$(echo $filename_extended | cut -d"." -f 2)
  foldername=$(echo $filename_extended | cut -d"_" -f 1,2,3)
  filename=${filename_extended%%.*}
  
  mkdir microbepore/data/FASTQ/cutadapt
  mkdir microbepore/data/FASTQ/cutadapt/${foldername}
  output=microbepore/data/FASTQ/cutadapt/${foldername}/${filename}
  mkdir ${output}
  
  # cutadapt
  cutadapt \
    -a "A{10}" \ # trim polyAs longer than 10 bases from the 3´end
    -e 1 \ # allowed error rate
    -j 0 \ # auto-detect cores
    -o ${output}/${filename}.cutadapt.fastq \
    ${file}
done
````

### Remove remaining SSP adapter using [`cutadapt`](https://cutadapt.readthedocs.io/en/stable/)  
- Remove remaining SSP sequences from the 5´ends of the cDNA reads using:   

```bash
input=microbepore/data/FASTQ/cutadapt

# >  SSP adapter
for file in ${input}/*/*/*cutadapt.fastq
do 
  filename_extended=${file##*/}
  keyword=$(echo $filename_extended | cut -d"." -f 2)
  foldername=$(echo $filename_extended | cut -d"_" -f 1,2,3)
  filename=${filename_extended%%.*}
  
  mkdir microbepore/data/FASTQ/cutadapt_SSP
  mkdir microbepore/data/FASTQ/cutadapt_SSP/${foldername}
  output=microbepore/data/FASTQ/cutadapt_SSP/${foldername}/${filename}
  mkdir ${output}

  cutadapt \
    -g "TTTCTGTTGGTGCTGATATTGCTGGG" \
    -e 1 \
    -j 0 \
    -o ${output}/${filename}.cutadapt_SSP.fastq \
    ${file}
done
```

### Mapping of trimmed reads, removing clips using `samclip`    
- Finally, trimmed reads were mapped using `minimap2` as described before  
- Reads with more than 10 clipped bases on either side were removed from the alignments using [`samclip`](https://github.com/tseemann/samclip) (v.0.4.0)   

1. Step: Align   

```bash
input=microbepore/data/FASTQ/cutadapt_SSP
fasta=microbepore/data/genome/NC_000913.3.fasta # downloaded from GenBank

# map (pychopper) > polyA_trimmed > SSP trimmed fastqs
for file in ${input}/*/*/*fastq
do 
  filename_extended=${file##*/}
  foldername=$(echo ${filename_extended} | cut -d"_" -f 1,2,3)
  filename=${filename_extended%%.*}

  mkdir microbepore/data/mapped/adapter_trimmed
  mkdir microbepore/data/mapped/adapter_trimmed/${foldername}
  output=microbepore/data/mapped/adapter_trimmed/${foldername}/${filename}
  mkdir ${output}

  ## align using minimap2
  if [[ $filename =~ "RNA" ]]; 
  then
  # align using minimap2
  minimap2 -ax splice -p 0.99 -uf -k14 --MD -t 8 ${fasta} ${file} > ${output}/${filename}.sam
  else
    minimap2 -ax splice -p 0.99 -k14 --MD -t 8 ${fasta} ${file} > ${output}/${filename}.sam
  fi
done
```

2. Step: Remove clipping > 10 bases

```bash
input=microbepore/data/mapped/adapter_trimmed
fasta=microbepore/data/genome/NC_000913.3.fasta # downloaded from GenBank
transcripts=microbepore/data/genomeNC_000913.3.transcripts.fasta # transcripts file made using gffread

# remove reads with more than 10 bases that are clipped on either side. 
for file in ${input}/*/*/*.sam
do 
  filename_extended=${file##*/}
  keyword=$(echo $filename_extended | cut -d"." -f 2)
  foldername=$(echo $filename_extended | cut -d"_" -f 1,2,3)
  filename=${filename_extended%%.*}
  
  if [[ $keyword =~ "sam" ]]; then
    echo ${foldername}
    echo ${filename}
    echo ${keyword}
    
    mkdir microbepore/data/mapped/trimmed
    mkdir microbepore/data/mapped/trimmed/${foldername}
    output=microbepore/data/mapped/trimmed/${foldername}/${filename}
    mkdir ${output}
  
    # remove mapped reads with a Maximum clip length to allow (10, 5 is default)
    samclip --max 10 --ref ${fasta} < ${file} > ${output}/${filename}.clipped.sam
    
    # convert to sorted.bam file
    samtools flagstat ${output}/${filename}.clipped.sam > ${output}/${filename}.clipped.stats.txt
    samtools view -bS ${output}/${filename}.clipped.sam -o ${output}/${filename}.clipped.bam
    samtools sort ${output}/${filename}.clipped.bam -o ${output}/${filename}.clipped.sorted.bam
    samtools index ${output}/${filename}.clipped.sorted.bam
    
    ## remap fastq converted reads
  bedtools bamtofastq -i ${output}/${filename}.clipped.sorted.bam -fq ${output}/${filename}.remapped.fastq
  
  ## map again
  if [[ $filename =~ "RNA" ]]; 
  then
  # align using minimap2
  minimap2 -ax splice -p 0.99 -uf -k14 --MD -t 8 ${transcripts} ${file} > ${output}/${filename}.remapped.sam
  else
    minimap2 -ax splice -p 0.99 -k14 --MD -t 8 ${transcripts} ${file} > ${output}/${filename}.remapped.sam
  fi
  
  # convert to sorted.bam file
  samtools view -bS ${output}/${filename}.remapped.sam -o ${output}/${filename}.remapped.bam
  samtools sort ${output}/${filename}.remapped.bam -o ${output}/${filename}.remapped.sorted.bam
  samtools index ${output}/${filename}.remapped.sorted.bam
  fi
done
```


