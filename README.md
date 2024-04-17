Expanding the transcriptomic toolbox in prokaryotes by Nanopore
sequencing of RNA and cDNA molecules
================
<a href="https://orcid.org/0000-0001-7444-2408">Felix
Grünberger<sup>1°</sup></a>,
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

This is the repository for the manuscript “Expanding the transcriptomic
toolbox in prokaryotes by Nanopore sequencing of RNA and cDNA
molecules.” In this study, we applied and benchmarked all currently
available RNA-seq kits from Oxford Nanopore technologies to analyse RNAs
in the prokaryotic model organism *Escherichia coli* K-12. These
include:  
- Direct sequencing of native RNAs (DRS) using RNA001 & RNA002
chemistry  
- Native cDNA sequencing (cDNA) using DCS109 chemistry  
- PCR-cDNA sequencing (PCR-cDNA) using PCB109 chemistry

The repository is currently actively developed.

[![Active
Development](https://img.shields.io/badge/Maintenance%20Level-Actively%20Developed-brightgreen.svg)](https://gist.github.com/cheerfulstoic/d107229326a01ff0f333a1d3476e068d)

## Full documentation here

<https://felixgrunberger.github.io/microbepore/>

## Preprint

Preprint will soon be available at bioRxiv.  
In case you are interested, have a look at our previous work:[Exploring
prokaryotic transcription, operon structures, rRNA maturation and
modifications using Nanopore-based native RNA
sequencing](%22https://www.biorxiv.org/content/10.1101/2019.12.18.880849v2%22).

## What can you find here

A description of the workflow using publicly available tools used to
basecall, demultiplex, trim, map and count data can be found in the
[pipeline](pipeline) section. Downstream analysis, including quality
control, annotation of transcript boundaries, gene body coverage
analysis, transcriptional unit annotation are based on custom
[Rscripts](Rscripts).

## Figures

Here you can find links to the scripts used to make all of the figures
based on numeric data.

<details>
<summary>
click to expand
</summary>

|               |     |     |                                                                          |
|---------------|-----|-----|--------------------------------------------------------------------------|
| Main          | 1   | A   | NA                                                                       |
| Main          | 1   | B   | NA                                                                       |
| Main          | 1   | C   | NA                                                                       |
| Main          | 1   | D   | NA                                                                       |
| Main          | 2   | A   | [Rscripts/salmon\_analysis.R](Rscripts/salmon_analysis.R)                |
| Main          | 2   | B   | [Rscripts/salmon\_analysis.R](Rscripts/salmon_analysis.R)                |
| Main          | 2   | C   | [Rscripts/salmon\_analysis.R](Rscripts/salmon_analysis.R)                |
| Main          | 3   | A   | [Rscripts/example\_coverage\_plots.R](Rscripts/example_coverage_plots.R) |
| Main          | 3   | B   | [Rscripts/end5\_detection.R](Rscripts/end5_detection.R)                  |
| Main          | 3   | C   | [Rscripts/end3\_detection.R](Rscripts/end3_detection.R)                  |
| Main          | 4   | A   | [Rscripts/gene\_body\_coverage.R](Rscripts/gene_body_coverage.R)         |
| Main          | 4   | B   | NA                                                                       |
| Main          | 4   | C   | [Rscripts/gene\_body\_coverage.R](Rscripts/gene_body_coverage.R)         |
| Main          | 4   | D   | [Rscripts/gene\_body\_coverage.R](Rscripts/gene_body_coverage.R)         |
| Main          | 4   | E   | [Rscripts/gene\_body\_coverage.R](Rscripts/gene_body_coverage.R)         |
| Main          | 4   | F   | [Rscripts/gene\_body\_coverage.R](Rscripts/gene_body_coverage.R)         |
| Main          | 5   | A   | [Rscripts/example\_coverage\_plots.R](Rscripts/example_coverage_plots.R) |
| Main          | 5   | B   | [Rscripts/operon\_analysis.R](Rscripts/operon_analysis.R)                |
| Main          | 5   | C   | [Rscripts/operon\_analysis.R](Rscripts/operon_analysis.R)                |
| Supplementary | 1   | NA  | NA                                                                       |
| Supplementary | 2   | A   | [Rscripts/bioanalyzer\_analysis.R](Rscripts/bioanalyzer_analysis.R)      |
| Supplementary | 2   | B   | [Rscripts/bioanalyzer\_analysis.R](Rscripts/bioanalyzer_analysis.R)      |
| Supplementary | 3   | A   | [Rscripts/raw\_read\_analysis.R](Rscripts/raw_read_analysis.R)           |
| Supplementary | 3   | B   | [Rscripts/raw\_read\_analysis.R](Rscripts/raw_read_analysis.R)           |
| Supplementary | 4   | A   | [Rscripts/raw\_read\_analysis.R](Rscripts/raw_read_analysis.R)           |
| Supplementary | 4   | B   | [Rscripts/raw\_read\_analysis.R](Rscripts/raw_read_analysis.R)           |
| Supplementary | 4   | C   | [Rscripts/raw\_read\_analysis.R](Rscripts/raw_read_analysis.R)           |
| Supplementary | 5   | A   | [Rscripts/mapped\_read\_analysis.R](Rscripts/mapped_read_analysis.R)     |
| Supplementary | 5   | B   | [Rscripts/mapped\_read\_analysis.R](Rscripts/mapped_read_analysis.R)     |
| Supplementary | 6   | A   | [Rscripts/mapped\_read\_analysis.R](Rscripts/mapped_read_analysis.R)     |
| Supplementary | 6   | B   | [Rscripts/mapped\_read\_analysis.R](Rscripts/mapped_read_analysis.R)     |
| Supplementary | 6   | C   | [Rscripts/mapped\_read\_analysis.R](Rscripts/mapped_read_analysis.R)     |
| Supplementary | 6   | D   | [Rscripts/mapped\_read\_analysis.R](Rscripts/mapped_read_analysis.R)     |
| Supplementary | 7   | A   | [Rscripts/mapped\_read\_analysis2.R](Rscripts/mapped_read_analysis2.R)   |
| Supplementary | 7   | B   | [Rscripts/mapped\_read\_analysis2.R](Rscripts/mapped_read_analysis2.R)   |
| Supplementary | 8   | A   | [Rscripts/mapped\_read\_analysis2.R](Rscripts/mapped_read_analysis2.R)   |
| Supplementary | 8   | B   | [Rscripts/mapped\_read\_analysis2.R](Rscripts/mapped_read_analysis2.R)   |
| Supplementary | 9   | A   | [Rscripts/mapped\_read\_analysis2.R](Rscripts/mapped_read_analysis2.R)   |
| Supplementary | 9   | B   | [Rscripts/mapped\_read\_analysis2.R](Rscripts/mapped_read_analysis2.R)   |
| Supplementary | 9   | C   | [Rscripts/mapped\_read\_analysis2.R](Rscripts/mapped_read_analysis2.R)   |
| Supplementary | 9   | D   | [Rscripts/mapped\_read\_analysis2.R](Rscripts/mapped_read_analysis2.R)   |
| Supplementary | 9   | E   | [Rscripts/mapped\_read\_analysis2.R](Rscripts/mapped_read_analysis2.R)   |
| Supplementary | 10  | A   | [Rscripts/seq\_depth.R](Rscripts/seq_depth.R)                            |
| Supplementary | 10  | B   | [Rscripts/seq\_depth.R](Rscripts/seq_depth.R)                            |
| Supplementary | 10  | C   | [Rscripts/seq\_depth.R](Rscripts/seq_depth.R)                            |
| Supplementary | 11  | A   | [Rscripts/salmon\_analysis.R](Rscripts/salmon_analysis.R)                |
| Supplementary | 11  | B   | [Rscripts/salmon\_analysis.R](Rscripts/salmon_analysis.R)                |
| Supplementary | 12  | A   | NA                                                                       |
| Supplementary | 12  | B   | [Rscripts/pychopper\_trimming.R](Rscripts/pychopper_trimming.R)          |
| Supplementary | 12  | C   | [Rscripts/pychopper\_trimming.R](Rscripts/pychopper_trimming.R)          |
| Supplementary | 13  | A   | [Rscripts/read\_end\_identities.R](Rscripts/read_end_identities.R)       |
| Supplementary | 13  | B   | [Rscripts/read\_end\_identities.R](Rscripts/read_end_identities.R)       |
| Supplementary | 14  | A   | [Rscripts/end5\_detection.R](Rscripts/end5_detection.R)                  |
| Supplementary | 14  | B   | [Rscripts/end5\_detection.R](Rscripts/end5_detection.R)                  |
| Supplementary | 15  | A   | [Rscripts/end5\_detection.R](Rscripts/end5_detection.R)                  |
| Supplementary | 15  | B   | [Rscripts/end5\_detection.R](Rscripts/end5_detection.R)                  |
| Supplementary | 15  | C   | [Rscripts/end5\_detection.R](Rscripts/end5_detection.R)                  |
| Supplementary | 16  | A   | [Rscripts/end5\_detection.R](Rscripts/end5_detection.R)                  |
| Supplementary | 16  | B   | [Rscripts/end5\_detection.R](Rscripts/end5_detection.R)                  |
| Supplementary | 17  | A   | [Rscripts/end5\_detection.R](Rscripts/end5_detection.R)                  |
| Supplementary | 17  | B   | [Rscripts/end5\_detection.R](Rscripts/end5_detection.R)                  |
| Supplementary | 18  | A   | [Rscripts/end3\_detection.R](Rscripts/end3_detection.R)                  |
| Supplementary | 18  | B   | [Rscripts/end3\_detection.R](Rscripts/end3_detection.R)                  |
| Supplementary | 19  |     | [Rscripts/end3\_detection.R](Rscripts/end3_detection.R)                  |
| Supplementary | 20  | A   | [Rscripts/end3\_detection.R](Rscripts/end3_detection.R)                  |
| Supplementary | 20  | B   | [Rscripts/end3\_detection.R](Rscripts/end3_detection.R)                  |
| Supplementary | 21  | A   | [Rscripts/operon\_analysis.R](Rscripts/operon_analysis.R)                |
| Supplementary | 21  | B   | [Rscripts/operon\_analysis.R](Rscripts/operon_analysis.R)                |
| Supplementary | 22  | A   | NA                                                                       |
| Supplementary | 22  | B   | [Rscripts/example\_coverage\_plots.R](Rscripts/example_coverage_plots.R) |
| Supplementary | 22  | C   | [Rscripts/example\_coverage\_plots.R](Rscripts/example_coverage_plots.R) |
| Supplementary | 22  | D   | [Rscripts/example\_coverage\_plots.R](Rscripts/example_coverage_plots.R) |
| Supplementary | 23  | A   | NA                                                                       |
| Supplementary | 23  | B   | [Rscripts/example\_coverage\_plots.R](Rscripts/example_coverage_plots.R) |
| Supplementary | 23  | C   | [Rscripts/example\_coverage\_plots.R](Rscripts/example_coverage_plots.R) |
| Supplementary | 23  | D   | [Rscripts/example\_coverage\_plots.R](Rscripts/example_coverage_plots.R) |

</details>

## Data availability

### Sequencing files in original FAST5 format

Sequencing files in original FAST5 format are publicly available in the
Sequence Read Archive [SRA](https://www.ncbi.nlm.nih.gov/sra) (RNA001:
PRJNA632538, all other datasets: PRJNA731531).

### Basecalled and demultiplexed FASTQ files

For easier access, basecalled & demultiplexed FASTQ files are available
in a [Google Drive
Folder](https://drive.google.com/drive/folders/1RO1yIAXWSnKgfe4-XuYFUCioHnHHLrpC?usp=sharing) and on [Zenodo](zenodo.org/record/4879174#.YLSkjy221pQ).

### Mapped BAM files

Minimap2-mapped untrimmed reads are also available in the [Google Drive
Folder](https://drive.google.com/drive/folders/1RO1yIAXWSnKgfe4-XuYFUCioHnHHLrpC?usp=sharing) and on [Zenodo](zenodo.org/record/4879174#.YLSkjy221pQ).

------------------------------------------------------------------------

## License

This project is under the general MIT License - see the
[LICENSE](LICENSE) file for details

## References
