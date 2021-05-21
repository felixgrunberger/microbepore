---
layout: page
title: Basecalling and demultiplexing of raw reads 
parent: Data analysis
nav_order: 2
---

## Basecalling and demultiplexing of raw reads 
All fast5 reads were basecalled and demultiplexed using guppy (ont-guppy-for-mk1c v4.3.4).  

### Direct RNA sequencing   

```bash
## > basecalling 210317_RNA002_Ecoli
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

### cDNA sequencing   


{% highlight ruby %}
def print_hi(name)
  puts "Hi, #{name}"
end
print_hi('Tom')
#=> prints 'Hi, Tom' to STDOUT.
{% endhighlight %}

This is the base Jekyll theme. You can find out more info about customizing your Jekyll theme, as well as basic Jekyll usage documentation at [jekyllrb.com](https://jekyllrb.com/)

You can find the source code for Minima at GitHub:
[jekyll][jekyll-organization] /
[minima](https://github.com/jekyll/minima)

You can find the source code for Jekyll at GitHub:
[jekyll][jekyll-organization] /
[jekyll](https://github.com/jekyll/jekyll)


[jekyll-organization]: https://github.com/jekyll
