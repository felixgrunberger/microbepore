Expanding the transcriptomic toolbox in prokaryotes by Nanopore
sequencing of RNA and cDNA molecules
================
<a href="https://orcid.org/0000-0001-7444-2408">Felix
Grünberger<sup>1</sup></a>,
<a href="https://orcid.org/0000-0002-0522-843X">Sébastien
Ferreira-Cerca<sup>2</sup></a>, and
<a href="https://orcid.org/0000-0002-0570-2517">Dina
Grohmann<sup>1°</sup></a>  

<sup>1</sup> Department of Biochemistry, Genetics and Microbiology,
Institute of Microbiology, Single-Molecule Biochemistry Lab &
Biochemistry Centre Regensburg, University of Regensburg,
Universitätsstraße 31, 93053 Regensburg, Germany

<sup>2</sup> Biochemistry III – Institute for Biochemistry, Genetics and
Microbiology, University of Regensburg, Universitätsstraße 31, 93053
Regensburg, Germany.

<sup>°</sup> Corresponding authors

<!-- README.md is generated from README.Rmd. Please edit that file -->

------------------------------------------------------------------------

![](docs/assets/images/microbepore_logo.png)

## About this repository

This is the repository for the manuscript “Exploring prokaryotic
transcription, operon structures, rRNA maturation and modifications
using Nanopore-based native RNA sequencing,” which can be found on
<a href = "https://www.biorxiv.org/content/10.1101/2019.12.18.880849v2">bioRxiv</a>.
It contains a description of the bioinformatical tools used to process
native RNA sequencing data and the downstream analysis mostly based on
custom [Rscripts](Rscripts).

The repository is currently actively developed.

[![Active
Development](https://img.shields.io/badge/Maintenance%20Level-Actively%20Developed-brightgreen.svg)](https://gist.github.com/cheerfulstoic/d107229326a01ff0f333a1d3476e068d)

## Full documentation here

<https://felixgrunberger.github.io/microbepore/>

## What can you find here

Here you can find all of the Rscripts used during the analysis (compare
)

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

## Figures

``` r
library(here)
#> here() starts at /Users/felix/Documents/R/GITHUB/microbepore
library(readxl)
devtools::install_github("glin/reactable")
#> Downloading GitHub repo glin/reactable@HEAD
#> 
#>      checking for file ‘/private/var/folders/m1/j6m28r317pngj1mw_l3v5zxw0000gn/T/RtmpqzTrr2/remotes1b6f62ccec15/glin-reactable-8197c00/DESCRIPTION’ ...  ✓  checking for file ‘/private/var/folders/m1/j6m28r317pngj1mw_l3v5zxw0000gn/T/RtmpqzTrr2/remotes1b6f62ccec15/glin-reactable-8197c00/DESCRIPTION’
#>   ─  preparing ‘reactable’:
#>      checking DESCRIPTION meta-information ...  ✓  checking DESCRIPTION meta-information
#>   ─  checking for LF line-endings in source and make files and shell scripts
#>   ─  checking for empty or unneeded directories
#>   ─  building ‘reactable_0.2.3.9000.tar.gz’
#>      
#> 
library(reactable)
figure_tables <- read_xlsx(here("tables/Supplementary_Table2.xlsx"))
reactable(figure_tables)
```

![](README-unnamed-chunk-2-1.png)<!-- -->

## Data analysis

Before starting a sequencing run, different running options can be
selected in MinKNOW. In case reads are stored in multi-FAST5-containing
files (e.g. 4000 reads per file), files can be converted to single-read
FAST5 files using the
[ont\_fast5\_api](https://github.com/nanoporetech/ont_fast5_api), as
some workflows
(e.g. [`nanopolish`](https://nanopolish.readthedocs.io/en/latest/) and
[`tombo`](https://nanoporetech.github.io/tombo/)) rely on single-FAST5
files for further analysis.  
After a run, reads are stored in two folders (*fast5\_failed*,
*fast5\_passed*). To prevent actual good reads from beeing discarded we
**included all reads from both folders** in the following steps of the
analysis.  
First, we converted multi-FAST5-files with the `multi_to_single_fast5`
command from the
[ont\_fast5\_api](https://github.com/nanoporetech/ont_fast5_api):

``` bash
#!/bin/bash

# set input_path, save_path, search in all folders for fast5 files, set number of threads
multi_to_single_fast5 \
    --input_path <path folder containing multi_read_fast5 files> \
    --save_path <path to folder where single_read fast5 files will be output> \
    --recursive <recursively search sub-directories> \
    --threads <number of CPU threads to use>
```

The output will be single-read FAST5 files in the *save\_path* folder
with one subfolder per multi-read input file.

### Navigation

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

### Demultiplexing using [`poreplex`](https://github.com/hyeshik/poreplex)

Multiplexed libraries (*how to* is described here:
<https://github.com/hyeshik/poreplex>) can be demultiplexed using
[`poreplex`](https://github.com/hyeshik/poreplex). Following this
approach four direct RNA sequencing libraries can be barcoded, pooled
and sequenced together. `Poreplex` can demultiplex the libraries into
separate folders with:

``` bash
#!/bin/bash

# trim adapters, basecall using albacore, de-multiplex, create symbolic fast5 link, increase number of working processes, sort reads to folders according to barcodes
poreplex \
    -i <path/to/fast5> \
    -o <path/to/output> \
    --trim-adapter <trim 3′ adapter sequences from FASTQ outputs> \
    --barcoding <sort barcoded reads into separate outputs> \
    --basecall <call the ONT albacore for basecalling on-the-fly> \
    --symlink-fast5 <create symbolic links to FAST5 files in output directories even when hard linking is possible> \
    --parallel <number of worker processes>
```

Reads can be basecalled automatically during demultiplexing using
`albacore`, the outdated basecaller of ONT. As we observed dramatic
differences in the detected read qualities and number of reads that can
be mapped, we chose to only sort the reads based on `poreplex` and used
`guppy` (Version 3.0.3) for basecalling. In this case the `poreplex` can
be shortened to:
`poreplex -i <input> -o <output> --barcoding --parallel`.

### Basecalling using `guppy`

Demultiplexed raw FAST5 files and other raw FAST5 files that have not
been barcoded (raw MinKNOW output) can be basecalled (*translated* in
FASTQ data) using `guppy`, the ONT-developed basecaller (available in
the ONT Community). We used version 3.0.3 for basecalling of all of our
reads:

``` bash
#!/bin/bash

guppy_basecaller \
    --flowcell FLO-MIN106 <flowcell-version> \
    --kit SQK-RNA001 <sequencing-kit> \
    --input $input <path to input FAST5 files> \
    --save_path $output <path to output> \
    --recursive <search input folders recursively> \
    --reverse_sequence yes <as RNA is sequenced from 3´to 5´> \
    --hp_correct 1 <enable homopolymer correction> \
    --enable_trimming <trim RNA reads> \
    --cpu_threads_per_caller <set number of threads> \
    --calib_detect <detect RCS spike in>
```

Demultiplexed files from *fast5\_failed* and *fast5\_passed* folders can
be basecalled seperately. Final output files from `guppy` are merged
with:

``` bash
#!/bin/bash

# FASTQ
### e.g. for BC1 failed and passed folders
cat $dir/guppy_data/BC1_fail/*.fastq $dir/guppy_data/BC1_pass/*.fastq > $dir/fastq_data/BC1_combined.fastq

# Sequencing summary
### e.g. for BC1 failed and passed folders
cat $dir/guppy_data/BC1_fail/sequencing_summary.txt $dir/guppy_data/BC1_pass/sequencing_summary.txt > $dir/summary_data/BC1_sequencing_summary.txt
```

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

# set file paths
dir=<path_to_basedir>
genome_fasta=<path_to_fasta_file>(downloaded from NCBI)
out_dir=$dir/mapped_data

# map using minimap2 | convert using samtools and bam2wig
for file in $dir/fastq_data/*_combined.fastq # path to all of the merged .fastq files
do 
  filename_extended=$file##*/
  filename=$filename_extended%%.*
  minimap2 -p 0.99 -ax splice -k14 --MD -uf $genome_fasta $file > $out_dir/$filename".sam" # map using minimaps2
  samtools flagstat $out_dir/$filename".sam" > $out_dir/$filename"_stats.txt" # calculate mapping statistics
  samtools view -bS $out_dir/$filename".sam" -o $out_dir/$filename".bam" # sam to bam
  samtools sort $out_dir/$filename".bam" -o $out_dir/$filename"_sorted.bam" # bam to sorted bam 
  samtools index $out_dir/$filename"_sorted.bam" # index sorted bam file 
  samtools view -h -F 16 $out_dir/$filename".bam" > $out_dir/$filename"_forward.bam" # create strand specific bam files
  samtools view -h -f 16 $out_dir/$filename".bam" > $out_dir/$filename"_reverse.bam" # strand specific bam files
  samtools sort $out_dir/$filename"_forward.bam" -o $out_dir/$filename"_forward_sorted.bam" # sort 
  samtools sort $out_dir/$filename"_reverse.bam" -o $out_dir/$filename"_reverse_sorted.bam" # sort
  samtools depth -a -d 0 $out_dir/$filename"_forward_sorted.bam" >  $out_dir/$filename"_forward_sorted_depth.txt" # calculate number of reads for each position in a strand specific way (no threshold set) and write to bed-like txt file
  samtools depth -a -d 0 $out_dir/$filename"_reverse_sorted.bam" > $out_dir/$filename"_reverse_sorted_depth.txt" # calculate number of reads for each position in a strand specific way (no threshold set) and write to bed-like txt file
  ./bam2wig $out_dir/$filename"_sorted.bam" # create wig files from bam files, all reads
  ./bam2wig -s top $out_dir/$filename"_sorted.bam" # create wig files from bam files, top strand
  ./bam2wig -s bottom $out_dir/$filename"_sorted.bam" # create wig files from bam files, bottom strand
  echo $filename "mapping finished" 
done
```

#### Mapping of reads to spike-in control

`Guppy` filters out the calibration reads (make sure to enable
`--calib_detect`) that can be mapped to the enolase gene to perform
quality control analysis.

``` bash
#!/bin/bash

# set file paths
dir=<path_to_basedir>
genome_fasta=<path_to_enolase_file>
out_dir=$dir/mapped_data

# before mapping merge failed and passed files
cat $dir/guppy_data/BC1_fail/calibration_strands/*.fastq $dir/guppy_data/BC1_pass/calibration_strands/*.fastq > $dir/enolase_data/BC1_calibration.fastq

# map using minimap2
for file in $dir/enolase_data/*calibration.fastq
do 
  filename_extended=$file##*/
  filename=$filename_extended%%.*
  minimap2 -ax splice -k14  --MD -uf $genome_fasta $file > $out_dir/$filename".sam"
  samtools flagstat $out_dir/$filename".sam" > $out_dir/$filename"_stats.txt"
  samtools view -bS $out_dir/$filename".sam" | samtools sort -o $out_dir/$filename"_sorted.bam"
  samtools index $out_dir/$filename"_sorted.bam"
  samtools depth -a -d 0 $out_dir/$filename"_sorted.bam" > $out_dir/$filename"_depth.txt"
  echo $filename "mapping finished"!  
done
```

### Poly(A)-tail analysis using [`nanopolish`](https://nanopolish.readthedocs.io/en/latest/quickstart_polya.html)

Poly(A) tail length was estimated by nanopolish following the
recommended workflow (Version 0.10.2,
<https://nanopolish.readthedocs.io/en/latest/quickstart_polya.html>)
([Loman, Quick, and Simpson 2015](#ref-Loman2015)).  
The workflow includes indexing of the reads using the `nanopolish index`
module and mapping of the reads using `minimap2`. The poly(A)-tail
length estimator can then be run on mapped files with
`nanopolish polya`, with genome fasta, fastq read file and mapped bam
file as input. The output is stored in a .tsv file that was analyzed
using a [custom R script](Rscripts/poly_a_tail_analysis.R).

### Gene expression analysis using `featurecounts`

We calculated transcript abundances for long-read Nanopore native RNA
reads using `featurecounts`([Liao, Smyth, and Shi 2019](#ref-Liao2019a))
with the following command in `R`:

``` r
featureCounts(allowMultiOverlap = T, 
              files = input_bam_file, 
              annot.ext = input_gff_file, 
              isGTFAnnotationFile = T, 
              GTF.featureType = name, <name argument was used multiple times for CDS, rRNA and tRNAs>
              GTF.attrType = "ID", 
              isLongRead = T)
```

Number of counts for each gene were added to the metadata information,
the `Rscript` for the complete analysis can be found in the next
section.

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

### Single-read analysis

Single-read analysis was used to look at long 3´UTRs and processing
events on the 16S and 23S rRNAs during ribosomal maturation.

#### mRNA

Examples are shown in Supplementary Fig. 8 of the manuscript and were
plotted using the
[`single_reads_pilin_hvo`](Rscripts/single_reads_pilin_hvo.R) and
[`single_reads_alba_pfu`](Rscripts/single_reads_alba_pfu.R) scripts.

#### rRNA

Reads for plotting of Figure 4 were plotted using the
[`single_reads_rrnac`](Rscripts/single_reads_rrnc.R) script.

### Detection of rRNA processing sites and classification of rRNA intermediates

To detect processing sites in archaeal and bacterial rDNA operons we
performed enrichment analysis of read start and end positions. Code is
provided for the *E. coli* analysis in
[`rRNA_read_boundaries_ecoli`](Rscripts/rRNA_read_boundaries_ecoli.R).  
Next, co-occurence analysis was performed by (i) categorizing reads
according to enriched and literature-expected 5´ positions, (ii)
selecting all reads that start within +/-1 from the relevant 5´ position
and (iii) analysing the respective read ends
([`rRNA_read_coocurence_ecoli`](Rscripts/rRNA_read_cooccurence_ecoli.R)).  
Exemplary reads of selected categories with enriched connected terminal
positions were visualised in a genome browser-like view utlizing scripts
similar to the [`single_reads_rrnac`](Rscripts/single_reads_rrnc.R)
Rscript. Circular rRNA precursor in archaea were initially observed in a
subset of reads, which end near/at the 5´cleavage site of the
bulge-helix-bulge (BHB), but are extensively left-clipped. To
investigate circ rRNA in more detail, we mapped reads to a permuted
linear sequence that contained joined 3´BHB/5´-BHB ends.

``` r
# in R

#...................................haloferax genome/fasta information
gff <- read.gff(here("data/genome_data/hvo.gff")) %>%
  dplyr::filter(type == "rRNA")

fasta <- readDNAStringSet(here("data/genome_data/hvo.fasta"))

# > 1598084 5´end BHB
# > 1599741 3´end BHB

#...................................write permuted 16S circ-rRNA sequence
circ16_junc_open <- paste(fasta$chr[(gff$end[1]-500):(1599741)],
                          fasta$chr[(1598084):(gff$start[1]+500)],
                          sep = "")

write.fasta(sequences = circ16_junc_open, 
            names = "circ_16S_open",
            file.out = here("data/nanopore_reanalysis/hvo/circ16S.fasta"))
```

Reads were mapped using `minimap2`, sorted and indexed using `samtools`.
BAM files were imported again in R and read start and end positions
visualised using the `geom_linerange` option in ggplot2
([`circ_read_plotting`](Rscripts/circ_read_plotting.R)).

### Modified base detection

We estimated the amount of modified bases using two different
approaches:

#### [`Tombo`](https://nanoporetech.github.io/tombo/)

We used `Tombo` (Version 1.5, <https://nanoporetech.github.io/tombo>) to
identify modified bases based on a comparison to a theoretical
distribution (*de novo* model) and based on the comparison to a
reference data set (sample-compare model)(compare [Stoiber et al.
2016](#ref-Stoiber2016)). The different steps include `preprocessing`,
`resquiggling` and `plotting` of the raw signals at specific genomic
coordinates and are well-described in the documentation.  
The probability of modified bases was calculated using the
`detect_modification` *de\_novo* command, written to *.wig* files using
the `text_output browser_files` command and added to the
`plot genome_locations` plot using the
[`tombo_fractions_positions`](Rscripts/tombo_fractions_positions.R)
script. For Fig. 6g the signals were calculated for both samples
(wildtype and deletion mutant) and compared using the
`control-fast5-basedirs` and `overplot Boxplot` option. For Fig 7b in
the manuscript a reference data set was created by sorting the reads
mapping to the 16S rRNA based on the position of the read start (before
or after gene start), thereby dividing the data set in reads that belong
to mature and unprocessed 16S rRNAs (script see here:
[`filter_rrna_reads_for_tombo`](Rscripts/filter_rrna_reads_for_tombo.R)).
The selected reads can be extracted using the custom Rscript or using
the `fast5_subset` module of the
[ont\_fast5\_api](https://github.com/nanoporetech/ont_fast5_api). Minion
read names of reads that fulfill certain criteria therefore have to be
written to a `read_id_list`. Probabilities were calculated for the
sample-compare model for all read categories and plotted using custom
R-scripts ([`tombo_intermediates`](Rscripts/tombo_intermediates.R)).

#### Pileup mapping assignments

For calculating the frequency of correct, deleted, inserted and wrong
nucleotides at a genomic position `pysamstats`
(<https://github.com/alimanfoo/pysamstats>) was used.

``` bash
#!/bin/bash

#...................................PILEUP BASES FREQUENCY
dir=/data

#....NOTEX WT
pysamstats -D 8000000 -S nofilter --window-size=1 -t variation_strand -f $dir/genome_data/hvo.fasta $dir/reanalysis_tombo/hvo/mapped_data/primary_pre_rRNA_wt_sorted.bam > $dir/pileups/primary_pre_rRNA_wt_pileup.tsv

pysamstats -D 8000000 -S nofilter --window-size=1 -t variation_strand -f $dir/genome_data/hvo.fasta $dir/reanalysis_tombo/hvo/mapped_data/closed_circ_rRNA_wt_sorted.bam > $dir/pileups/closed_circ_rRNA_wt_pileup.tsv

pysamstats -D 8000000 -S nofilter --window-size=1 -t variation_strand -f $dir/genome_data/hvo.fasta $dir/reanalysis_tombo/hvo/mapped_data/open_circ_rRNA_wt_sorted.bam > $dir/pileups/open_circ_rRNA_wt_pileup.tsv

pysamstats -D 8000000 -S nofilter --window-size=1 -t variation_strand -f $dir/genome_data/hvo.fasta $dir/reanalysis_tombo/hvo/mapped_data/mature_rRNA_wt_sorted.bam > $dir/pileups/mature_rRNA_wt_pileup.tsv
```

Plots were generated using custom R scripts
([`pileup_mod_bases`](Rscripts/pileup_mod_bases.R)). The results were
compared to known modification sites in 16S rRNA for *H. volcanii*
([Grosjean et al. 2008](#ref-Grosjean2008)) and *P. furiosus*. Note that
the positions of modified RNA base modifications for *P. furiosus* are
derived from a recently published study in *P. abyssi*.

## Data availability

### Raw sequencing files

The raw sequencing data in fast5 format have been submitted to the NCBI
sequence read archive
(<a href="https://www.ncbi.nlm.nih.gov/sra">SRA</a>) under BioProject
accession number PRJNA632538.

### Fastq files

<https://drive.google.com/drive/folders/195xIehurQuKQVgCu92TRREmFiv6oEdov?usp=sharing>

### Mapped files

------------------------------------------------------------------------

## License

This project is under the general MIT License - see the
[LICENSE](LICENSE) file for details

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Grosjean2008" class="csl-entry">

Grosjean, Henri, Christine Gaspin, Christian Marck, Wayne A. Decatur,
and Valérie de Crécy-Lagard. 2008. “<span class="nocase">RNomics and
Modomics in the halophilic archaea Haloferax volcanii: Identification of
RNA modification genes</span>.” *BMC Genomics* 9: 1–26.
<https://doi.org/10.1186/1471-2164-9-470>.

</div>

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

<div id="ref-Loman2015" class="csl-entry">

Loman, Nicholas J., Joshua Quick, and Jared T. Simpson. 2015. “<span
class="nocase">A complete bacterial genome assembled de novo using only
nanopore sequencing data</span>.” *Nature Methods* 12 (8): 733–35.
<https://doi.org/10.1038/nmeth.3444>.

</div>

<div id="ref-Stoiber2016" class="csl-entry">

Stoiber, Marcus H, Joshua Quick, Rob Egan, Ji Eun Lee, Susan E Celniker,
Robert Neely, Nicholas Loman, Len Pennacchio, and James B Brown. 2016.
“<span class="nocase">De novo Identification of DNA Modifications
Enabled by Genome-Guided Nanopore Signal Processing</span>.” *bioRxiv*,
094672. <https://doi.org/10.1101/094672>.

</div>

</div>
