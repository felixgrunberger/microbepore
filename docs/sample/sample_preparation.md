---
layout: page
title: Sample preparation
permalink: /sample_preparation/
nav_order: 3
has_children: false
---

## Sample preparation  
{: .no_toc }

1. TOC
{:toc}

### Experimental design  

### RNA extraction and treatment  
#### Cell growth and RNA extraction   
Escherichia coli K-12 MG1655 cells were grown in rich medium (10 g tryptone, 5 g yeast extract, 5 g NaCl per liter, pH 7.2) to an OD600nm of 0.5-0.6. To stabilize RNAs, two volumes of RNAlater (Thermo Fisher Scientific) were immediately added to the cultures and stored at -20°C until cells were harvested by centrifugation at 4°C.
Total RNA of all samples except RNA001 was extracted using RNeasy Mini Kit (Qiagen) according to the manufacturer´s instructions. RNA001 RNA was purified using the Monarch® Total RNA Miniprep Kit (New England Biolabs). The integrity of total RNA from E. coli was assessed via a Bioanalyzer (Agilent) run using the RNA 6000 Pico Kit (Agilent), and only RNAs with RNA integrity numbers (RIN) above 9.5 were used for subsequent treatments and sequencing. 


#### Poly(A) tailing, rRNA depletion and additional RNA treatment   
RNAs were heat incubated at 70°C for 2 min and snap cooled on a pre-chilled freezer block before poly(A)–tailing using the E. coli poly(A) polymerase (New England Biolabs). Briefly, 5 µg RNA, 20 units poly(A) polymerase, 5 µl reaction buffer and 1 mM ATP were incubated for 15 min at 37°C in a total reaction volume of 50 µl. To stop and clean up the reaction, poly(A)-tailed RNAs were purified following the RNeasy Micro clean-up protocol (Qiagen), which was used for all subsequent RNA clean-ups. The efficiency of poly(A)-tailing was evaluated via a Bioanalyzer run. rRNA depletion was performed using the Pan-Prokaryote riboPOOL by siTOOLs, which effectively removes rRNAs from E. coli. For TEX-treated samples, partial digestion of RNAs that are not 5´-triphosphorylated (e.g. tRNAs, rRNAs) was achieved by incubation of the RNA with the Terminator 5´-Phosphate-Dependent Exonuclease (TEX, Lucigen). 
For the RNA001 sample, 10 µg of RNA were incubated with 1 unit TEX, 2 µl TEX reaction buffer (Lucigen) and 0.5 µl RiboGuard RNase Inhibitor (Lucigen) in a total volume of 20 µl for 60 minutes at 30°C. Additionally, 20 ng of rRNA-depleted samples subsequently used in the PCR-cDNA workflow were TEX-treated using the same enzyme and buffer concentrations but reducing the reaction time to 15 minutes. All reactions were terminated by adding EDTA and the RNAs cleaned up.
Before library preparation, the extent of remaining buffer and DNA contamination were tested by performing standard spectroscopic measurements (Nanodrop One) and using the Qubit 1X dsDNA HS assay kit (Thermo Fisher Scientific). Input RNAs were finally quantified using the Qubit RNA HS assay kit. 


### Library preparation  
Libraries for Nanopore sequencing were prepared from poly(A)-tailed RNAs according to the protocols provided by Oxford Nanopore (Oxford Nanopore Technologies Ltd, Oxford, UK) for direct sequencing of native RNAs (SQK-RNA001, SQK-RNA002), direct cDNA native barcoding (SQK-DCS109 with EXP-NBD104) and PCR-cDNA barcoding (SQK-PCB109) with the following minor modifications: Agencourt AMPure XP magnetic beads (Beckman Coulter) in combination with 1 µl of RiboGuard RNase Inhibitor (Lucigen) were used instead of the recommended Agencourt RNAclean XP beads to clean up samples. For reverse transcription, Maxima H Minus Reverse Transcriptase (Thermo Fisher Scientific) was not only used for all cDNA samples, but also for the RNA002 samples (SuperScript III Reverse Transcriptase from Thermo Fisher Scientific used for RNA001 sample).  
The amount of input RNA, barcoding strategy, number of PCR cycles, extension times and library preparation kit used can be found in Supplementary Table 1 and are also summarized in part in the following Figure.





