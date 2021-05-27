---
layout: page
title: Sample preparation
permalink: /sample_preparation/
nav_order: 3
has_children: false
---

## Sample preparation  
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

### Cell growth and RNA extraction  
- Grow *Escherichia coli* K-12 MG1655 cells in rich medium (10 g tryptone, 5 g yeast extract, 5 g NaCl per liter, pH 7.2) to an OD<sub>600nm</sub> of 0.5-0.6 (or your organism of choice in the conditions you´re interested in)   
- Stabilize RNAs by adding two volumes of RNAlater (Thermo Fisher Scientific)  
- Store cultures at -20°C or harvest cells by centrifugation at 4°C  
- Extract RNA with the method of your choice. We (and others) had good experience with the RNeasy Kit from Qiagen  
> **Note:** Silica-membrane columns have a cut-off size of about 200 nucleotides 
- Check integrity of total RNA (Bioanalyzer)   
> **Note:** RNA degradation has a massive influence on the results of your experiment  

### Poly(A) tailing, rRNA depletion and additional RNA treatment   
- Perform **poly(A)–tailing** using the *E. coli* poly(A) polymerase (New England Biolabs):  
  - Incubate RNAs at 70°C for 2 min and snap cool it on a pre-chilled freezer block  
  - Incubate 5 µg RNA, 20 units poly(A) polymerase, 5 µl reaction buffer, 1 mM ATP for 15 min at 37°C in a total reaction volume of 50 µl. 
  - Stop and clean up the reaction following the RNeasy Micro clean-up protocol (Qiagen)  
  - Evaluate efficiency of poly(A)-tailing (Bioanalyzer, peaks of rRNAs)  
- Perform **rRNA depletion**, e.g. using the [Pan-Prokaryote riboPOOL by siTOOLs](https://www.sitoolsbiotech.com/ribopools.php) (Clean-up can also be performed using the RNeasy protocol)    
> **Note:** After depletion, only **~5%** of the initial RNA input quantity are left.     
- Perform additional treatment of your choice and clean-up, e.g. using the Terminator 5´-Phosphate-Dependent Exonuclease (TEX, Lucigen). 

> **Important**: Before library preparation, check: 
1. extent of remaining buffer and DNA contamination that could have a negative impact on your library and 
2. RNA size and quantity (Qubit) to determine molarities as accurately as possible  


### Library preparation  
Libraries for Nanopore sequencing were prepared from poly(A)-tailed RNAs according to the protocols provided by Oxford Nanopore (Oxford Nanopore Technologies Ltd, Oxford, UK) for direct sequencing of native RNAs (SQK-RNA001, SQK-RNA002), direct cDNA native barcoding (SQK-DCS109 with EXP-NBD104) and PCR-cDNA barcoding (SQK-PCB109) with the following minor modifications: Agencourt AMPure XP magnetic beads (Beckman Coulter) in combination with 1 µl of RiboGuard RNase Inhibitor (Lucigen) were used instead of the recommended Agencourt RNAclean XP beads to clean up samples. For reverse transcription, Maxima H Minus Reverse Transcriptase (Thermo Fisher Scientific) was not only used for all cDNA samples, but also for the RNA002 samples (SuperScript III Reverse Transcriptase from Thermo Fisher Scientific used for RNA001 sample).  
The amount of input RNA, barcoding strategy, number of PCR cycles, extension times and library preparation kit used can be found in Supplementary Table 1 and are also summarized in part in the following Figure.





