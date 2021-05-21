---
layout: page
title: Project organization
permalink: /project_organization/
nav_order: 5
has_children: false
---


## Project organization   
The easiest way to find your way around your project is to 

```R
library(here)
#> here() starts at /Users/jenny/rrr/here_here
here()
#> [1] "/Users/jenny/rrr/here_here"

```

```bash
microbepore/
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