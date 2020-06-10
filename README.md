
---
title: "Workflow for Complete Network Analysis of MDM using mdmnets package"
author: Tatyana Zamkovaya
output: github_document 
---


## Introduction and Overview
Networks of the co-occurrence relationships present within a microbial community can detect otherwise overlooked associations between species. Additionally, network metrics can quantitatively determine the most important species of a given environment. While many advances have been made in microbial ecology thanks to a network-based approach, most analyses to date have neglected the role of uncharacterized, unknown species called Microbial Dark Matter (MDM). The presence and dominance of MDM across all non-human biomes suggests that this exclusion is a mistake and that the ecological role of MDM must be understood to gain complete knowledge of all microbial interactions encompassing an environment. We need to better understand why MDM persist and with whom they interact on a global scale. More importantly, we need to evaluate whether certain MDM taxa confer a value to the environment that cannot be provided by other, well-characterized microbes. 

**Here, we present a network-based approach to determine the ecological relevance of MDM**. First, all taxa designated as unassigned, ambiguous, or uncultured are defined as _**MDM**_. The proportion and prevalence of MDM taxa is then evaluated and compared to other, well-characterized (_**Known**_) taxa. Next, to assess the global impact of MDM on the microbial environment, three networks are constructed - including all taxa (called _**"Original"**_), excluding MDM (called _**"Without Unknown"**_), and excluding an equal number of random knowns (called _**"Bootstrap"**_), for each taxonomic classification level from Phylum to Genus. Network centrality measures _**degree**_, _**betweenness**_, and _**closeness**_ are calculated and compared between the three network types, to quantitatively determine if MDM removal significantly impacts network shape and structure. Changes in network metrics are statistically evaluated and visualized as boxplots. Finally, to assess which MDM or known taxa most influence their respective microbial communities, _**hub score**_ is calculated for all taxa present in the _**Original**_ networks. All results can be validated by checking for early filtering, correlation metric, and sample biases using a range of sample prevalence thresholds, correlation metric tools, and constructing an OTU abundance heatmap for all samples used.    

## MDM Workflow Steps and accompanying scripts/packages
### MDM Workflow Overview

|Workflow Step                          | Script/Package                                | 
|:------------------------------------- |:---------------------------------------------:|
|1. Import Data into R                  | `met_min_occ.sh` & `met_min_occ.R`            |
|2. Analyze MDM and Known OTUs          |`getsharedOTUs_known_v_unk.R` (& `mdmnets`)|
|3. Create Networks with MDM            | `mdmnets`                                 |
|4. Create Networks without MDM         | `mdmnets`                                 |
|5. Create Bootstrap Networks           | `mdmnets`                                 |
|6. Compare Network Measure Changes     | `mdmnets`                                 | 
|7. Create Hub Networks                 | `mdmnets`                                 |
|8. Compare Node Neighbor Interactions  | `mdmnets`                                 |
|9. Validate Results                    | `mdmnets`                                 |

![Workflow](https://github.com/tatyanazam/mdmnets/blob/master/MDMWorkflowpic.jpg)

## Necessary Packages for MDM Analysis 
Before starting the MDM Workflow analysis with mdmnets, make sure to have the following installed in R:
```{r nec_librairies, warning=FALSE, message=FALSE}
library(data.table)
library(reshape2)
library(phyloseq)
library(plyr)
library(dplyr)
library(SpiecEasi)
library(ggplot2)
library(igraph)
library(gridExtra)
library(ggpubr)
library(ggsignif)
library(gplots)
library(devtools)
library(stats)
library(parallel)
library(RColorBrewer)
library(Hmisc)
library(viridis)
library(Matrix)
```

Now you are ready to begin the MDM Analysis!

## Installation of mdmnets package 

![logo](https://github.com/ConesaLab/MDM/blob/master/mdmnetslogo_v2.jpg)

You can install mdmnets with:
 
 * `devtools::install_github("ConesaLab/MDM")`
 
 To see how to use the mdmnets package, please ![go here](https://github.com/ConesaLab/MDM/blob/master/MDMnets_vignette.md)
