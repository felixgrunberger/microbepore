Pipeline
================

-   [Data generation](#data-generation)
-   [Navigation](#navigation)
-   [Bash pipe](#bash-pipe)
    -   [Basecalling using
        guppy\_basecaller](#basecalling-using-guppy_basecaller)
    -   [Demultiplexing using
        guppy\_barcoder](#demultiplexing-using-guppy_barcoder)
    -   [Mapping using minimap2](#mapping-using-minimap2)
        -   [Mapping using `minimap2`](#mapping-using-minimap2-1)
            -   [Mapping of reads to reference
                genomes](#mapping-of-reads-to-reference-genomes)
        -   [Summarise metadata information on single-read
            level](#summarise-metadata-information-on-single-read-level)
    -   [Gene abundance estimation using `salmon` in alignment-based
        mode](#gene-abundance-estimation-using-salmon-in-alignment-based-mode)
    -   [Trimming of reads](#trimming-of-reads)
        -   [Identification of full-length reads using
            pychopper](#identification-of-full-length-reads-using-pychopper)
        -   [Remove polyA-tails using
            cutadapt](#remove-polya-tails-using-cutadapt)
        -   [Remove remaining SSP adapter using
            cutadapt](#remove-remaining-ssp-adapter-using-cutadapt)
        -   [Mapping of trimmed reads](#mapping-of-trimmed-reads)
        -   [Remove clipping &gt; 10 bases](#remove-clipping--10-bases)
    -   [Detection of transcript
        boundaries](#detection-of-transcript-boundaries)
        -   [5´end detection](#5end-detection)
        -   [3´end detection](#3end-detection)
    -   [Gene body coverage analysis](#gene-body-coverage-analysis)
    -   [Transcriptional unit analysis](#transcriptional-unit-analysis)
        -   [Quality control](#quality-control)
            -   [Analysis of raw reads](#analysis-of-raw-reads)
            -   [Analysis of mapped reads](#analysis-of-mapped-reads)
            -   [Run statistics](#run-statistics)
        -   [Detection of transcriptional units
            (TU)](#detection-of-transcriptional-units-tu)
        -   [Annotation of transcription start sites
            (TSS)](#annotation-of-transcription-start-sites-tss)
        -   [Annotation of transcription termination sites
            (TTS)](#annotation-of-transcription-termination-sites-tts)

## Data generation

Libraries for Nanopore sequencing were prepared from poly(A)-tailed RNAs
according to the SQK-RNA001 Kit protocol (Oxford Nanopore, Version:
DRS\_9026\_v1\_revP\_15Dec2016) with minor modifications for barcoded
libraries. In this case, Agencourt AMPure XP magnetic beads (Beckman
Coulter) in combination with 1 µl of RiboGuard RNase Inhibitor (Lucigen)
were used instead of the recommended Agencourt RNAclean XP beads to
purify samples after enzymatic reactions. For the barcoded libraries,
the RTA adapter was replaced by custom adapters described in
<https://github.com/hyeshik/poreplex> and reverse transcription (RT) was
performed in individual tubes for each library. After RT reactions, cDNA
was quantified using the Qubit DNA HS assay kit (Thermo Fisher
Scientific) and equimolar amounts of DNA for the multiplexed samples
were used in the next step for ligation of the RNA Adapter (RMX) in a
single tube. Subsequent reactions were performed according to the
protocols recommended by ONT. The libraries were sequenced on a MinION
using R9.4 flow cells and subsequently, FAST5 files were generated using
the recommended script in MinKNOW.

## Navigation

We managed our folders in the following way:

``` bash
Native_RNAseq_Microbes/
├── data/
|   ├── genome_data
|   ├── tidy_data
|   ├── summary_data
|   ├── mapped_data
|   ├── guppy_data
|   ├── fastq_data
|   ├── coverage_data
|   ├── meme_data
|   ├── enolase_data
|   ├── poly_data
|   ├── tombo_data
|   └── operon_data
├── Rscrips
├── figures
├── tables/
|   ├── tss_tables
|   ├── tts_tables
|   ├── tu_tables
|   └── counts_tables
├── LICENSE
└── README
```

Relative paths in the custom `R` and `bash` scripts are included for the
complete analysis.

# Bash pipe

## Basecalling using guppy\_basecaller

Demultiplexed raw FAST5 files and other raw FAST5 files that have not
been barcoded (raw MinKNOW output) can be basecalled (*translated* in
FASTQ data) using `guppy`, the ONT-developed basecaller (available in
the ONT Community). We used version 3.0.3 for basecalling of all of our
reads:

``` bash
#!/bin/bash

## Basecalling of RNA samples
vol=/data/devices/SanDisk-Extreme_SSD-p1
input=${vol}/210317_RNA002_Ecoli/210317_RNA002_Ecoli/20210317_1630_MN24615_FAP64575_e847ffe7
output=${vol}/210317_RNA002_Ecoli/basecalling
nohup guppy_basecaller \
--input_path $input \
--save_path $output \
-c rna_r9.4.1_70bps_hac.cfg  \
--calib_detect \
--reverse_sequence true \
--u_substitution true \
--compress_fastq \
--fast5_out \
--recursive \
--progress_stats_frequency 60 \
--chunks_per_runner 256 \
--gpu_runners_per_device 4 \
--num_callers 1 \
-x 'auto' &
```

Demultiplexed files from *fast5\_failed* and *fast5\_passed* folders can
be basecalled seperately. Final output files from `guppy` are merged
with:

``` bash
## Basecalling of cDNA samples
ssh minit@mc-110727.local
vol=/data/devices/SAMSUNG-HM320JI-p2
input=${vol}/FGseq021_201208_PCB109_Ecoli/201208_PCB109_Ecoli/20201208_1459_MN24615_FAL36510_8287229a
output=/data/201208_PCB109_Ecoli

nohup guppy_basecaller \
--input_path $input \
--save_path $output \
-c dna_r9.4.1_450bps_hac.cfg  \
--compress_fastq \
--fast5_out \
--recursive \
--progress_stats_frequency 60 \
--chunks_per_runner 256 \
--gpu_runners_per_device 4 \
--num_callers 1 \
-x 'auto' &
```

## Demultiplexing using guppy\_barcoder

``` bash
#!/bin/bash
guppy_barcoder \
--input_path $input \
--save_path $output \
--config configuration.cfg \
--barcode_kits SQK-PCB109 \
--progress_stats_frequency 60
```

## Mapping using minimap2

### Mapping using [`minimap2`](https://github.com/lh3/minimap2)

#### Mapping of reads to reference genomes

FAST5 passed and FAST5 failed reads that have been demultiplexed using
`poreplex` and basecalled using `guppy` can now be mapped to the
reference genomes using
[`minimap2`](https://github.com/lh3/minimap2)([Li 2018](#ref-Li2018)).
Release 2.17-r941 was used for our analysis. Output alignments in the
SAM format were generated with the recommended options for noisy
Nanopore Direct RNA-seq (-ax splice, -uf, -k14) and also with (1) -p set
to 0.99, to return primary and secondary mappings and (2) with –MD
turned on, to include the MD tag for calculating mapping identities.
Alignment files were further converted to bam files, sorted and indexed
using `SAMtools` ([Li et al. 2009](#ref-Li2009)). Strand-specific wig
and bigwig files were finally created using `bam2wig` (Version 1.5,
<https://github.com/MikeAxtell/bam2wig>).

``` bash
#!/bin/bash
dir=/Volumes/EX_SSD
for file in ${dir}/data/fastq_data/*/*.fastq
do
  
  # folder and filenames
  f_ex=${file##*/}
  foldername=$(echo ${f_ex} | cut -d"_" -f 1,2,3)
  filename=${f_ex%%.*}
  
  # make directories
  mkdir ${dir}/data/mapped_data_notrimming
  mkdir ${dir}/data/mapped_data_notrimming/${foldername}
  output=${dir}/data/mapped_data_notrimming/${foldername}/${filename}
  mkdir ${output}

  if [[ $filename =~ "RNA" ]]; 
  then
  # align using minimap2
  minimap2 -ax splice -p 0.99 -uf -k14 --MD -t 8 ${dir}/data/genome_data/NC_000913.3.fasta ${file} > ${output}/${filename}.sam
  else
    minimap2 -ax splice -p 0.99 -k14 --MD -t 8 ${dir}/data/genome_data/NC_000913.3.fasta ${file} > ${output}/${filename}.sam
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
  minimap2 -ax splice -p 0.99 -uf -k14 --MD -t 8 ${dir}/data/genome_data/NC_000913.3.transcripts.fasta ${output}/${filename}.remapped.fastq > ${output}/${filename}.remapped.sam
  else
    minimap2 -ax splice -p 0.99 -k14 --MD -t 8 ${dir}/data/genome_data/NC_000913.3.transcripts.fasta ${output}/${filename}.remapped.fastq > ${output}/${filename}.remapped.sam
 fi
 
  # convert to sorted.bam file
  samtools view -bS ${output}/${filename}.remapped.sam -o ${output}/${filename}.remapped.bam
  samtools sort ${output}/${filename}.remapped.bam -o ${output}/${filename}.remapped.sorted.bam
  samtools index ${output}/${filename}.remapped.sorted.bam
```

### Summarise metadata information on single-read level

For each data set one matrix including raw-read, mapped-read and count
information was prepared using the [`tidy_data`](Rscripts/tidy_data.R)
script. Following files are needed as input:

-   sequencing\_summary.txt from `guppy`  
-   genome fasta file
-   genome gff file  
-   Sorted .bam file from `minimap2` output

Files are stored as `R` files, but can also be exported as tables
(commented in the code). As each line represents one single read or
sometimes one read refers to multiple lines (multimapping of rRNA reads
to multiple loci) these tables can become quite large.

## Gene abundance estimation using `salmon` in alignment-based mode

We calculated transcript abundances for long-read Nanopore native RNA
reads using `featurecounts`([Liao, Smyth, and Shi 2019](#ref-Liao2019a))
with the following command in `R`:

``` bash
#!/bin/bash
dir=/Volumes/EX_SSD

for file in ${dir}/data/mapped_data_notrimming/*/*/*remapped.sorted.bam
do
  
  # folder and filenames
  f_ex=${file##*/}
  foldername=$(echo ${f_ex} | cut -d"_" -f 1,2,3)
  filename=${f_ex%%.*}
  
  # create dir for quantification using salmon in alignment-based mode (used in conda environment)
  mkdir ${dir}/data/salmon_data_notrimming
  mkdir ${dir}/data/salmon_data_notrimming/${foldername}
  output_salmon=${dir}/data/salmon_data_notrimming/${foldername}/${filename}
  mkdir ${output_salmon}
  
  conda activate salmon
  
  # make transcripts files
  # gffread $genome_gff -g $genome_fasta -w $genome_transcripts
  
  salmon quant \
  -t ${dir}/data/genome_data/NC_000913.3.transcripts.fasta \
  -l A \
  -a ${file} \
  -o ${output_salmon} \
  --threads 8 
  
  conda deactivate
  
done
```

## Trimming of reads

### Identification of full-length reads using pychopper

``` bash
#!/bin/bash

# perform pychopper for all DCS and PCB files
for file in $dir/data/fastq_data/*/*.fastq
do 

  # folder and filenames
  f_ex=${file##*/}
  foldername=$(echo $f_ex | cut -d"_" -f 1,2,3)
  filename=${f_ex%%.*}
  
  # make directories
  mkdir ${dir}/data/pychopper_data_auto
  output1=${dir}/data/pychopper_data_auto/${foldername}
  mkdir ${output1}
  output=${output1}/${filename}
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

# perform pychopper using the -x rescue option for DCS files
for file in $dir/data/pychopper_data_auto/*/*/*unclassified.fastq
do 

  # folder and filenames
  filename_extended=${file##*/}
  foldername=$(echo $filename_extended | cut -d"_" -f 1,2,3)
  filename=${filename_extended%%.*}
  
  # make directories
  mkdir ${dir}/data/pychopper_data_rescue
  output1=${dir}/data/pychopper_data_rescue/${foldername}
  mkdir ${output1}
  output=${output1}/${filename}
  mkdir ${output}
  
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

# merge all full-length and rescued reads as full-length
for file in $dir/data/pychopper_data_auto/*/*/*full_length_output.fastq
do 
  filename_extended=${file##*/}
  foldername=$(echo $filename_extended | cut -d"_" -f 1,2,3)
  filename=$(echo $filename_extended | cut -d"_" -f 1,2,3,4,5)

  keyword=$(echo $foldername | cut -d"_" -f 2)
  
  mkdir ${dir}/data/fastq_fl_data
  output1=${dir}/data/fastq_fl_data/${foldername}
  mkdir ${output1}
  output=${output1}/${filename}
  mkdir ${output}
  
  if [[ $keyword =~ "PCB109" ]]; then
    cat $file $dir/data/pychopper_data_auto/${foldername}/${filename}/${filename}_rescued.fastq > ${output}/${filename}_full_length_all.fastq
  elif [[ $keyword =~ "DCS109" ]]; then
    cat $file $dir/data/pychopper_data_auto/${foldername}/${filename}/${filename}_rescued.fastq $dir/data/pychopper_data_rescue/${foldername}/${filename}_unclassified/${filename}_unclassified_full_length_output.fastq $dir/data/pychopper_data_rescue/${foldername}/${filename}_unclassified/${filename}_unclassified_rescued.fastq > ${output}/${filename}_full_length_all.fastq
  fi
done
```

### Remove polyA-tails using cutadapt

``` bash
#!/bin/bash

for file in ${dir}/data/fastq_fl_data/*/*/*_full_length_all.fastq
do 

  # folder and filenames
  filename_extended=${file##*/}
  keyword=$(echo $filename_extended | cut -d"." -f 2)
  foldername=$(echo $filename_extended | cut -d"_" -f 1,2,3)
  filename=${filename_extended%%.*}
  
    
  mkdir ${dir}/data/fastq_fl_cutadapt
  mkdir ${dir}/data/fastq_fl_cutadapt/${foldername}
  output=${dir}/data/fastq_fl_cutadapt/${foldername}/${filename}
  mkdir ${output}
  
  cutadapt \
    -a "A{10}" \
    -e 1 \
    -j 0 \
    -o ${output}/${filename}.cutadapt.fastq \
    ${file}
done
```

### Remove remaining SSP adapter using cutadapt

``` bash
#!/bin/bash

# >  SSP adapter
for file in ${dir}/data/fastq_fl_cutadapt/*/*/*cutadapt.fastq
do 
  filename_extended=${file##*/}
  keyword=$(echo $filename_extended | cut -d"." -f 2)
  foldername=$(echo $filename_extended | cut -d"_" -f 1,2,3)
  filename=${filename_extended%%.*}
  
    
    mkdir ${dir}/data/fastq_fl_cutadapt_SSP
    mkdir ${dir}/data/fastq_fl_cutadapt_SSP/${foldername}
    output=${dir}/data/fastq_fl_cutadapt_SSP/${foldername}/${filename}
    mkdir ${output}
  
   cutadapt \
    -g "TTTCTGTTGGTGCTGATATTGCTGGG" \
    -e 1 \
    -j 0 \
    -o ${output}/${filename}.cutadapt_SSP.fastq \
    ${file}
done
```

### Mapping of trimmed reads

``` bash
#!/bin/bash

# map (pychopper) > polyA_trimmed > SSP trimmed fastqs
for file in ${dir}/data/fastq_fl_cutadapt_SSP/*/*/*fastq
do 
  filename_extended=${file##*/}
  foldername=$(echo ${filename_extended} | cut -d"_" -f 1,2,3)
  filename=${filename_extended%%.*}

  mkdir ${dir}/data/mapped_data_pychopper_auto_cutadapt_SSP
  output=${dir}/data/mapped_data_pychopper_auto_cutadapt_SSP/${foldername}/${filename}
  mkdir ${dir}/data/mapped_data_pychopper_auto_cutadapt_SSP/${foldername}
  mkdir ${output}

  ## align using minimap2
  if [[ $filename =~ "RNA" ]]; 
  then
  # align using minimap2
  minimap2 -ax splice -p 0.99 -uf -k14 --MD -t 8 ${dir}/data/genome_data/NC_000913.3.fasta ${file} > ${output}/${filename}.sam
  else
    minimap2 -ax splice -p 0.99 -k14 --MD -t 8 ${dir}/data/genome_data/NC_000913.3.fasta ${file} > ${output}/${filename}.sam
  fi
done
```

### Remove clipping &gt; 10 bases

``` bash
#!/bin/bash

# remove reads with more than 10 bases that are clipped on either side. 
for file in ${dir}/data/mapped_data_pychopper_auto_cutadapt_SSP/*/*/*.sam
do 
  filename_extended=${file##*/}
  keyword=$(echo $filename_extended | cut -d"." -f 2)
  foldername=$(echo $filename_extended | cut -d"_" -f 1,2,3)
  filename=${filename_extended%%.*}
  
  if [[ $keyword =~ "sam" ]]; then
    echo ${foldername}
    echo ${filename}
    echo ${keyword}
    
    mkdir ${dir}/data/mapped_data_pychopper_auto_cutadapt_SSP_clipped
    mkdir ${dir}/data/mapped_data_pychopper_auto_cutadapt_SSP_clipped/${foldername}
    output=${dir}/data/mapped_data_pychopper_auto_cutadapt_SSP_clipped/${foldername}/${filename}
    mkdir ${output}
  
    # remove mapped reads with a Maximum clip length to allow (10, 5 is default)
    samclip --max 10 --ref ${dir}/data/genome_data/NC_000913.3.fasta < ${file} > ${output}/${filename}.clipped.sam
    
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
  minimap2 -ax splice -p 0.99 -uf -k14 --MD -t 8 ${dir}/data/genome_data/NC_000913.3.transcripts.fasta ${file} > ${output}/${filename}.remapped.sam
  else
    minimap2 -ax splice -p 0.99 -k14 --MD -t 8 ${dir}/data/genome_data/NC_000913.3.transcripts.fasta ${file} > ${output}/${filename}.remapped.sam
  fi
  
  # convert to sorted.bam file
  samtools view -bS ${output}/${filename}.remapped.sam -o ${output}/${filename}.remapped.bam
  samtools sort ${output}/${filename}.remapped.bam -o ${output}/${filename}.remapped.sorted.bam
  samtools index ${output}/${filename}.remapped.sorted.bam
  fi
done
```

## Detection of transcript boundaries

### 5´end detection

``` bash
#!/bin/bash

# perform tss detection for pychopper auto > cutadapt_polyA > SSP-cutadapt > clipped
for file in ${dir}/data/mapped_data_pychopper_auto_cutadapt_SSP_clipped/*/*/*clipped.sorted.bam
do 
  # file and folder names
  filename_extended=${file##*/}
  keyword=$(echo $filename_extended | cut -d"." -f 2)
  foldername=$(echo $filename_extended | cut -d"_" -f 1,2,3)
  filename=${filename_extended%%.*}
  
  # make directories
  mkdir ${dir}/data/tss_data
  mkdir ${dir}/data/tss_data/tss_data_pychopper_auto_cutadapt_SSP_clipped
  mkdir ${dir}/data/tss_data/tss_data_pychopper_auto_cutadapt_SSP_clipped/${foldername}
  output=${dir}/data/tss_data/tss_data_pychopper_auto_cutadapt_SSP_clipped/${foldername}/${filename}
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
    -a ${output}/${filename}.plus.peaks \
    -b ${output}/${filename}.plus.bedgraph \
    > ${output}/${filename}.plus.peaks.counts
    
  bedtools intersect \
    -wao \
    -a ${output}/${filename}.minus.peaks \
    -b ${output}/${filename}.minus.bedgraph \
    > ${output}/${filename}.minus.peaks.counts
    
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

``` bash
# perform tts detection for pychopper auto > cutadapt polyA > cutadapt SSP > clipped
for file in ${dir}/data/mapped_data_pychopper_auto_cutadapt_SSP_clipped/*/*/*clipped.sorted.bam
do 
  filename_extended=${file##*/}
  keyword=$(echo $filename_extended | cut -d"." -f 2)
  foldername=$(echo $filename_extended | cut -d"_" -f 1,2,3)
  filename=${filename_extended%%.*}
    
  echo ${filename}

  mkdir ${dir}/data/tts_data
  mkdir ${dir}/data/tts_data/tts_data_pychopper_auto_cutadapt_SSP_clipped
  mkdir ${dir}/data/tts_data/tts_data_pychopper_auto_cutadapt_SSP_clipped/${foldername}
  output=${dir}/data/tts_data/tts_data_pychopper_auto_cutadapt_SSP_clipped/${foldername}/${filename}
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
    -a ${output}/${filename}.plus.peaks \
    -b ${output}/${filename}.plus.bedgraph \
    > ${output}/${filename}.plus.peaks.counts
    
  bedtools intersect \
    -wao \
    -a ${output}/${filename}.minus.peaks \
    -b ${output}/${filename}.minus.bedgraph \
    > ${output}/${filename}.minus.peaks.counts
    
    
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

## Gene body coverage analysis

``` bash
# calculate coverage over transcripts with TSS and TTS | for pychopper auto > cutadapt > clipped 
for file in ${dir}/data/mapped_data_pychopper_auto_cutadapt_SSP_clipped/*/*/*clipped.sorted.bam
do 
  filename_extended=${file##*/}
  keyword=$(echo $filename_extended | cut -d"." -f 2)
  foldername=$(echo $filename_extended | cut -d"_" -f 1,2,3)
  filename=${filename_extended%%.*}
    
  echo ${filename}
  echo ${foldername}

  # mk dirs
  mkdir ${dir}/data/coverage_data/coverage_data_pychopper_auto_cutadapt_SSP_clipped_stranded
  mkdir ${dir}/data/coverage_data/coverage_data_pychopper_auto_cutadapt_SSP_clipped_stranded/${foldername}
  output=${dir}/data/coverage_data/coverage_data_pychopper_auto_cutadapt_SSP_clipped_stranded/${foldername}/${filename}
  mkdir ${output}

  # calc coverage
  samtools view -F 16 -o temp.sorted.bam ${file} 
  bedtools coverage \
  -d \
  -a ${dir}/data/bed_data/transcripts.plus.bedgraph \
  -b temp.sorted.bam \
  > ${output}/${filename}.plus.coverage
  
  samtools view -f 16 -o temp.sorted.bam ${file} 
  bedtools coverage \
  -d \
  -a ${dir}/data/bed_data/transcripts.minus.bedgraph \
  -b temp.sorted.bam \
  > ${output}/${filename}.minus.coverage
done
```

## Transcriptional unit analysis

### Quality control

Quality control of sequencing runs was performed at the raw-read level
(unmapped reads from `guppy` output) and at the mapped-read level (after
`minimap2` mapping and `tidy_data` conversion).

#### Analysis of raw reads

Raw read plots and statistics were generated from `guppy`s
*sequencing\_summary.txt* files and are shown in the
[`raw_read_plots`](Rscripts/raw_read_plots.R) script.  
To control for problems during loading (air bubbles, …) an interactive
heatmap output of MinKNOW can be rebuild using following script:
<!--[`animated_qc`](Rscripts/animated_qc.R).  ![](figures/animated_throughput.gif) -->

In a similar way the total cumulative throughput over time can be
animated using the [’animted\_qc\`](Rscripts/animated_qc.R) script:  
<!-- ![](figures/animated_hours_yield.gif)  -->

#### Analysis of mapped reads

Mapped read analysis was performed using the
[`mapped_read_plots`](Rscripts/mapped_read_plots.R) script. The number
of reads mapping to different genomic features was calculated with
`featurecounts`, visualized using the
[`featurecounts_categories`](Rscripts/featurecounts_categories.R) script
and are stored in the [`counts_table`](tables/counts_tables/) files.

#### Run statistics

Supplementary Table 1 in the manuscript based on raw and mapped features
was calculated using the following script:
[`calculate_run_statistics`](Rscripts/calculate_run_statistics.R).

### Detection of transcriptional units (TU)

The detection of transcriptional units is a two-step process:

-   Collapsing of overlapping reads to TU-clusters is described in
    [`tu_cluster_generation.R`](Rscripts/tu_cluster_generation.R)  
-   Splitting of of clusters based on sequencing depth on the 3´end is
    described in
    [`tu_subcluster_annotation`](Rscripts/tu_subcluster_annotation.R)

The clusters for each genome are saved in the
[*tu\_tables*](tables/tu_tables/) section.  
Comparison to databases (manuscript Supplementary Figures 10a-c) were
plotted using the
[`tu_comparison_database`](Rscripts/tu_comparison_database.R) script.

Prerequisites for the TU annotation, like sequencing of full-lenght
transcripts and coverage bias towards 3´ end was evaluated based on the
scripts [`fraction_full_length`](Rscripts/fraction_full_length.R) and
[`coverage_drop_3prime`](Rscripts/coverage_drop_3prime.R).

### Annotation of transcription start sites (TSS)

The detection of transcription start sites (TSS) was based on the
detection of transcriptional units and is described in
[`transcription_start_sites`](Rscripts/transcription_start_sites.R).  
A transcriptional start site table for each organism can be found here:
[*tss\_tables*](tables/tss_tables/)

### Annotation of transcription termination sites (TTS)

Transcription termination sites mapping was performed in a similar way
as TSS detection and is described in the
[`transcription_termination_sites`](Rscripts/transcription_termination_sites.R)
script. Additionally, single-gene tracks were analysed using the
[`tts_single_gene`](Rscripts/tts_single_gene.R) script.  
A transcriptional termination site table for each organism can be found
here: [*tts\_tables*](tables/tts_tables/)

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Li2018" class="csl-entry">

Li, Heng. 2018. “<span class="nocase">Minimap2: Pairwise alignment for
nucleotide sequences</span>.” *Bioinformatics*.
<https://doi.org/10.1093/bioinformatics/bty191>.

</div>

<div id="ref-Li2009" class="csl-entry">

Li, Heng, Bob Handsaker, Alec Wysoker, Tim Fennell, Jue Ruan, Nils
Homer, Gabor Marth, Goncalo Abecasis, and Richard Durbin. 2009. “<span
class="nocase">The Sequence Alignment/Map format and SAMtools</span>.”
*Bioinformatics* 25 (16): 2078–79.
<https://doi.org/10.1093/bioinformatics/btp352>.

</div>

<div id="ref-Liao2019a" class="csl-entry">

Liao, Yang, Gordon K. Smyth, and Wei Shi. 2019. “<span
class="nocase">The R package Rsubread is easier, faster, cheaper and
better for alignment and quantification of RNA sequencing reads</span>.”
*Nucleic Acids Research*. <https://doi.org/10.1093/nar/gkz114>.

</div>

</div>