---
layout: page
title: Basecalling and demultiplexing of raw reads 
parent: Data analysis
nav_order: 2
---

## Basecalling and demultiplexing of raw reads  
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
- For better comparability we re-basecalled the FAST5 file from all runs using the same version on guppy on a Mk1C (ont-guppy-for-mk1c v4.3.4)   
- Once the Mk1C is set up correctly and connected to the internet, updates are displayed automatically and can be completed by click   
- Other software downloads are available from the [ONT downloads page](https://community.nanoporetech.com/downloads)    

### Basecalling of raw reads using `guppy_basecaller`   
- After sequencing (and despite live-basecalling) all datasets in the *raw_FAST5 folder* were re-basecalled using `guppy` (ont-guppy-for-mk1c v4.3.4) in high-accuracy mode (rna_r9.4.1_70bps_hac.cfg, dna_r9.4.1_450bps_hac.cfg) without quality filtering   
-  The output files in FASTQ format were written to the *basecalled folder*.   

**Note**: Mk1C   
- [We (and others)](https://community.nanoporetech.com/posts/fastest-way-to-re-basecall) noticed the apparently acces via GUI is not the fastest way to perform basecalling  
- In case you want to optimise GPU-accelerated basecalling on the Mk1C (settings for NVIDIA Jetson TX2) using have a look at [Miles Benton`s GitHub](https://github.com/sirselim/jetson_nanopore_sequencing)   

According this his suggestions (which really speed up the process), we performed basecalling like this (using ssh access to the Mk1C):     

```bash
# files
input=microbepore/data/raw_FAST5/run_id # add run id
output_DRS=microbepore/data/FASTQ/normal/run_id # add run id
output_cDNA=microbepore/data/basecalled/run_id # add run id

# Basecalling of DRS files
guppy_basecaller \
--input_path ${input} \ # input path
--save_path ${output_DRS} \ # output path
-c rna_r9.4.1_70bps_hac.cfg  \ # config file: high accuracy RNA
--calib_detect \ # detect calibration spike-in
--reverse_sequence true \ # reverse since sequenced 3´-->5´
--u_substitution true \ # replace U´s with T´s
--compress_fastq \ # compress output
--fast5_out \ # output FAST5
--recursive \ # look for FAST5 recursively in path
--progress_stats_frequency 60 \ # output progress every minute
--chunks_per_runner 256 \ # options for Mk1C
--gpu_runners_per_device 4 \ # options for Mk1C
--num_callers 1 \ # options for Mk1C
-x auto # options for Mk1C

# Basecalling of cDNA files 
guppy_basecaller \
--input_path ${input} \
--save_path ${output_cDNA} \
-c dna_r9.4.1_450bps_hac.cfg \ # config file: high accuracy cDNA 
--compress_fastq \
--fast5_out \
--recursive \
--progress_stats_frequency 60 \
--chunks_per_runner 256 \
--gpu_runners_per_device 4 \
--num_callers 1 \
-x auto
```

**Note**: 
- DRS & (PCR-)cDNA runs require different options.    
- Config file selection based on selected accuracy, flowcell version, library preparation kit are listed with `guppy_basecaller --print_workflows`    

Since basecalling is probably taking quite a while, but you want to keep the process running without interruption you can use the [`nohup`](https://linuxhint.com/how_to_use_nohup_linux/) command. You can use it like `nohup guppy_basecaller ..... &`. This starts basecalling using `guppy` and puts to process in the background (that´s what the `&` does). Now, you can log out from the Mk1C and the process keeps running. To monitor the process you can check the shell output that is written to `nohup.out`.    

With the selected options `guppy` produces fast5_pass, fast5_fail, fastq, summary and report files that are written to the *FASTQ folder*. 
FASTQ are not grouped in pass and fail groups since `--min_qscore` is not enabled. Multiple FASTQs can be merged using `cat microbepore/data/basecalled/run_id/*.fastq > microbepore/data/basecalled/run_id/run_id.fastq`.  

Sequencing summary files are also  written to the *FASTQ folder* and are used during the quality control of the runs and reads. For better viewing they can be moved to the *summary folder* using
`mv microbepore/data/FASTQ/run_id/sequencing_summary.txt microbepore/data/summary/run_id.txt`

### Demultiplexing of basecalled reads using `guppy_barcoder`   
Next, multiplexed cDNA libraries are demultiplexed in a separate step using `guppy_barcoder`.   
```bash
# files
input=microbepore/data/basecalled/run_id # add run id
output=microbepore/data/FASTQ/normal/run_id # add run id

# Demultiplexing of (PCR-)cDNA files
guppy_barcoder \
--input_path ${input} \
--save_path ${output} \
--config configuration.cfg \
--barcode_kits SQK-PCB109 \
--progress_stats_frequency 60
```

Multiple FASTQs are written to the *FASTQ folder* and can be merged with e.g. `cat microbepore/data/FASTQ/run_id/barcode01/*.fastq > microbepore/data/FASTQ/run_id/run_id_barcode01.fastq`. Barcode summary files are written to the *FASTQ folder* and can be moved to the *barcode folder*. 




