---
layout: page
title: Project
nav_order: 2
permalink: /project
---

## Project 
At the beginning of this project, there was no comprehensive analysis of Nanopore sequencing of RNA and cDNA molecules using 3<sup>rd</sup> generation Nanopore technology in **prokaryotes** available.
Therefore, the aim of this study was to benchmark current RNA-seq protocols provided by [Oxford Nanopore Technologies](http://nanoporetech.com) (ONT) using *Escherichia coli* as one of the most popular bacterial model organisms. We were especially interested, how well Nanopore sequencing captures multiple transcriptomic features at once and how the data compare to similar protocols, like [SMRT-Cappable-seq](https://www.nature.com/articles/s41467-018-05997-6).  
<br> 
In this documentation, we would like give a more detailed description of things to consider for [sample preparation](../sample_preparation) and [sequencing](../sequencing) and guide you through the individual steps of the [data analysis](../data_analysis) using tools and scripts developed by ONT, many other developers and custom R workflows.  

## Experimental design  
We evaluated the performance of all RNA-seq protocols currently available from ONT, namely:  
- **DRS**: Direct sequencing of native RNAs (using SQK-RNA001 & SQK-RNA002)  
- **cDNA**: Direct sequencing of cDNAs (using SQK-DCS109)   
- **PCR-cDNA**: Sequencing of PCR-amplifed cDNAs (using SQK-PCB109)  

  