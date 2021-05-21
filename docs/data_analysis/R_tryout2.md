---
layout: page
title: R tryout
parent: Data analysis
nav_order: 10
---
  

To get straight to the point: This post is about <strong>data manipulation and visualization in R</strong>.  
Obviously, that is something completely different from my first post and with more content coming up, I have to think about giving the blog a clearer structure 🤔  
But for now, it should be fine. As you will definitely see some code during your read, please do not be intimidated! You will not miss the final statement, even if you are not familiar with programming. As some kind of a reminiscence to the previous topic, this post again deals with movies and the development of their duration over the years.  

<h2>
Should I start learning how to code?
</h2>
Learning how to code from scratch is awfully painful.  
In the beginning, I realized, that you first have to know which language fits your needs. But without previous knowledge, it is nearly impossible to choose one. Luckily, there are many <a target="_blank" href="https://www.codementor.io/codementorteam/beginner-programming-language-job-salary-community-7s26wmbm6">articles</a> trying to help you to come to a decision.  
For those of you working with Microsoft Excel on a daily basis, but not being fully satisfied and wanting to try something new, <a target="_blank" href="https://www.python.org">Python</a> or <a target="_blank" href="https://www.r-project.org">R</a> might be worth a look. I personally started to make full use of R half a year ago and I can highly recommend it to anybody dealing with data of whatever nature.  

<h2>
Why use R?
</h2>
The benefits of using R are easy to identify and I will only list some of them.  
<strong>R</strong>:
<ul>
<li>
is actively developed
</li>
<li>
is available for every system
</li>
<li>
allows you to perform a reproducible, high-quality analysis
</li>
<li>
has powerful graphical capabilities
</li>
<li>
is <strong>100% free</strong>
</li>
</ul>
<h2>
4 easy steps to look at your data
</h2>
This blog post will not explain how to <a target="_blank" href="https://www.r-project.org">install R</a> or <a target="_blank" href="https://www.rstudio.com/products/rstudio/download/">RStudio</a> (use RStudio!).  
What I want to show is a typical analysis pipeline, that can easily be adopted to any kind of spreadsheet data. I will go step-by-step from data import to drawing a final conclusion.  
So let´s start by loading the packages we need. These packages (<a target="_blank" href="https://www.datacamp.com/community/tutorials/r-packages-guide#gs.oJviIjA">install them first</a>) add new functions and therefore increase the power of R dramatically. <br> <br>


``` r

# 1. step: Load all packages we need ----------------------------------------------------------------------------

library(dplyr)          # dplyr is a grammar of efficient data manipulation developed by Hadley Wickham
library(ggplot2)        # ggplot2 is a powerful plotting system in R
library(plotly)         # plotly is a graphing library to do interactive, publication-quality graphs
library(viridis)        # viridis comes with different color scales to make pretty plots 
library(ggjoy)          # joyplots provides a way of visualizing changes in distributions over time or space
library(gridExtra)      # gridExtra enables you to arrange multiple grid-based plots on a page

```

