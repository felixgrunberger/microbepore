---
layout: page
title: Read alignment
parent: Data analysis
nav_order: 3
---

## Read alignment   
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
### Genome files   
- Genome FASTA and GFF3 files have been downloaded from [GenBank](https://www.ncbi.nlm.nih.gov/nuccore/545778205)  
- Transcript file was made using [`gffread`](https://github.com/gpertea/gffread) with `gffread microbepore/data/genome/NC_000913.3.gff -g microbepore/data/genome/NC_000913.3.fasta -w microbepore/data/genome/NC_000913.3.transcripts.fasta`  


### Mapping of reads to the genome using [`minimap2`](https://github.com/lh3/minimap2)      
- Files were mapped to the reference genome from *Escherichia coli* K-12 MG1655 ([GenBank](https://www.ncbi.nlm.nih.gov/nuccore/545778205): U00096.3) using `minimap2` (Release 2.18-r1015)  
- Output alignments in the SAM format were generated with `-ax splice -k14` for **Nanopore cDNA-seq** and `-ax splice, -uf, -k14` for **DRS** with i) `-p 0.99`, to return primary and secondary mappings and ii) with `--MD`, to include the MD tag for calculating mapping identities 
- Alignment files were further converted to BAM files, sorted and indexed using [`SAMtools`(https://github.com/samtools/)     
- To analyse single reads in more detail with respect to the RNA type (mRNA, rRNA, other ncRNA, unspecified) they map to, BAM files were first converted back to FASTQ using [`bedtools`](https://bedtools.readthedocs.io/en/latest/) v2.29.2   
- Next FASTQ files were remapped to a transcriptome file using `minimap2` with the previously mentioned parameters to assign single read names with feature IDs       

```bash
# files
input=microbepore/data/FASTQ/normal # input directory with all merged FASTQ files, 1 for each barcode or single DRS run
fasta=microbepore/data/genome/NC_000913.3.fasta # downloaded from GenBank
transcripts=microbepore/data/genomeNC_000913.3.transcripts.fasta # transcripts file made using gffread

# Mapping & Remapping - loop through all FASTQs
for file in ${input}/*/*.fastq
do
  
  # folder and filenames
  f_ex=${file##*/}
  foldername=$(echo ${f_ex} | cut -d"_" -f 1,2,3) # depending on how you name your files 
  filename=${f_ex%%.*}
  
  # make directories
  mkdir microbepore/data/mapped/raw # direct output to mapped folder for raw reads
  mkdir microbepore/data/mapped/raw/${foldername} # run_id
  output=microbepore/data/mapped/raw/${foldername}/${filename} # run_id/barcode_id
  mkdir ${output}

  if [[ $filename =~ "RNA" ]]; 
  then
  # align using minimap2
  minimap2 -ax splice -p 0.99 -uf -k14 --MD -t 8 ${fasta} ${file} > ${output}/${filename}.sam # DRS
  else
    minimap2 -ax splice -p 0.99 -k14 --MD -t 8 ${fasta} ${file} > ${output}/${filename}.sam # (PCR-)cDNA
  fi
 
  # convert to sorted.bam file
  samtools view -bS ${output}/${filename}.sam -o ${output}/${filename}.bam
  samtools sort ${output}/${filename}.bam -o ${output}/${filename}.sorted.bam
  samtools index ${output}/${filename}.sorted.bam
  
  # bam to fastq for remapping of mapped reads
  bedtools bamtofastq -i ${output}/${filename}.sorted.bam -fq ${output}/${filename}.remapped.fastq
  
  # map again
  if [[ $filename =~ "RNA" ]]; 
  then
  minimap2 -ax splice -p 0.99 -uf -k14 --MD -t 8 ${transcripts} ${output}/${filename}.remapped.fastq > ${output}/${filename}.remapped.sam
  else
    minimap2 -ax splice -p 0.99 -k14 --MD -t 8 ${transcripts} ${output}/${filename}.remapped.fastq > ${output}/${filename}.remapped.sam
 fi
 
  # convert to sorted.bam file
  samtools view -bS ${output}/${filename}.remapped.sam -o ${output}/${filename}.remapped.bam
  samtools sort ${output}/${filename}.remapped.bam -o ${output}/${filename}.remapped.sorted.bam
  samtools index ${output}/${filename}.remapped.sorted.bam
done
```