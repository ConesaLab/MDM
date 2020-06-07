# How to use `mdmnets` package to analyze data
## Table of Contents
1. [Load Data into R and convert to phyloseq objects](#intro)
    * [Loading QIIME/OTU data](#qiimeimport) 
    * [Loading dada2 ASV data](#dada2import)
    * [Loading QIIME2 ASV data](#qiime2import)
2. [Calculate Known vs Unknown (MDM) Proportions](#prop)
3. [Create Networks with MDM](#netswmdm)
4. [Create Networks without MDM](#netswomdm)
5. [Calculate Network Measure Changes](#meas)
6. [Create Bootstrap Networks](#boot)
7. [Visualize and Statistically Evaluate Network Measure Changes](#boxplot)
8. [Create Hub Networks](#hub)
9. [Compare Node Neighbor Interactions](#neigh)
10. [Validate Results](#val)

## 1. Load data into phyloseq format <a name="intro"> </a>

#### Load and convert QIIME/ OTU format data to phyloseq <a name="qiimeimport"> </a>

```r
library(phyloseq)
#OTU table format is typically taxa as rows, samples as columns 
otutable <- import_biom(BIOMfilename = "/path/to/yourbiom/yourbiom.biom")
mapping <- import_qiime_sample_data(mapfilename = "/path/to/your/metadata/mapping_file.txt")

#Merge mapping file with otu table into one phyloseq object
phylo <- merge_phyloseq(otutable, mapping)
```
#### Load/convert ASV data generated from QIIME2 or dada2 to phyloseq
##### dada2 formatted data conversion <a name="dada2import"> </a>
If data was made with dada2 and in ASV format, data can be converted like so, using [the dada2 tutorial](https://benjjneb.github.io/dada2/tutorial.html) as a guide

```r
library(dada2)
library(phyloseq)
library(biostrings)
library(RDPutils)
#dada2 output format is typically taxa as columns and samples as rows
#using the dada2 tutorial as a guide, and the seqtab.nochim as the input

# 1. read in ASV dada2 output file:
seqtab.nochim <- readRDS("seqtab_nochim.rds")

#2. format fasta sequences of ASVs into R with biostrings package
rep.seqs <- colnames(seqtab.nochim)
rep.seqs <- Biostrings::DNAStringSet(rep.seqs)

#3. generate taxa names for ASVs
otu.names <- RDPutils::make_otu_names(1:length(rep.seqs))
colnames(seqtab.nochim) <- otu.names
names(rep.seqs) <- otu.names

#4. create phyloseq object
asv.table <- otu_table(seqtab.nochim, taxa_are_rows=FALSE)
phylo <- phyloseq(asv.table, rep.seqs) 

#sample data and tax table information can also added to the phyloseq command above
```
##### QIIME2 ASV data conversion  <a name="qiime2import"> </a>
 ###### 1. If ASVs were formed in QIIME2 directly, use the following to import/convert data
```r
#QIIME2 data (in .qza format) can be imported using the 
#qiime2R package
if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")
library(qiime2R)

ASV_table <- read_qza("table.qza")
#just like with dada2, ASVs are columns and samples are rows

tax_table <- read_qza("taxonomy.qza")
#taxonomy table is imported separately, using same function 

#separate out taxonomy information by rank with 
#parse_taxonomy function
tax_table <- parse_taxonomy(taxonomy$data)
# the output will be ASVs as rows and 7 columns, 
    #from Kingdom to Species, designating the taxonomy of each ASV

#merge together with phyloseq function from phyloseq package
phylo <- phyloseq(otu_table(ASV_table$data, taxa_are_rows=T), 
tax_table(as.data.frame(tax_table)) %>% column_to_rownames("Feature ID") %>% as.matrix())

#or merge ASV and taxonomy info together, along with sample data into phyloseq with qza_to_phyloseq function
phylo <- qza_to_phyloseq(features = "table.qza", tree = "tree.qza", "taxonomy.qza", metadata = "metadata.tsv")
``` 
######  2. If you prefer to import without using qiime2R package, data can also be converted directly in qiime2
1. First export asv table

```
qiime tools export \
	--input-path table.qza \
	--output-path phyloseq
```	

2. convert output biom file to tab-separated text 

```
biom convert \
	-i phyloseq/feature-table.biom \
	-o phyloseq/otu_table.tsv \
	--to-tsv
```	
3. Then in R do the following:
```r
library(phyloseq)
library(Biostrings)
asv_table <- read.table("otu_table.tsv", row.names=1, header=TRUE, sep="\t")

#if dealing with large data (many GB), read in the file using 
#data.table package and fread function
library(data.table)
asv_table <- data.table::fread("otu_table.tsv")

phylo <- otu_table(asv_table, taxa_are_rows=TRUE)
```
## 2. Calculate Known vs Unknown (MDM) Proportions <a name="prop"> </a>
To calculate prevalence and abundance of taxa present, if you have a biom file, use the function `make_prev_df` and your biom file as input, or if you have dada2/QIIME2 data imported as a phyloseq object already, you may do the following to replicate the results of this function:
```r
prevdf = apply(X = otu_table(phylo),MARGIN = ifelse(taxa_are_rows(ps), yes=1, no=2), FUN=function(x){sum(x>0)}) #find prevalence of each taxa per sample
  #bind together taxonomic information and prevalence and total abundance of each taxa into one dataframe
  prevdf = data.frame(Prevalence = prevdf, TotalAbundance = taxa_sums(phylo), tax_table(phylo))
  prevdf <- prevdf[prevdf$Prevalence > 0,] #remove taxa not present in samples
  return(prevdf)
```
Now, using the resulting _prevdf_ dataframe as input, a lineplot of the prevalence of known vs unknown taxa at genus level can be created
using the function `get_back_counts_for_line_plotf`. 
```r
prev_lineplot <- get_back_counts_for_line_plotf(prevdf, "Environment name here")
```
The above functions can be very computationally intensive so, if necessary, calculate known vs unknown proportions using an HPC environemnt or R (not RStudio.)  

## 3. Create Networks with MDM <a name="netswmdm"> </a>
To create networks, first decide on the best sample prevalence for your data. We recommend using as high a percent sample prevalence as possible (ideally 40-60 %), or a minimum of 30 percent sample prevalence for extremely sparse data. 
Use the function `get_back_res_meeting_min_occ` on the `phylo` object created when loading data to create the phyloseq and graph objects needed to build networks. 
The default percent sample prevalence is currently 40 % (written as 0.4) but this parameter (`filter_val_percent`) may be modified to best fit each user's data. 
```r
phylo_for_net <- get_back_res_meeting_min_occ(phylo, filter_val_percent=0.4)
#resulting data will be list of 2- phyloseq, graph
#phyloseq with only taxa meeting sample prevalence kept
#graph created using SpiecEasi neighborhood algorithm on phyloseq
```
Make sure to rename any ambiguous or unassigned taxa to NA to easily visualize unknowns within networks using provided `unassignedOTUs_MDM` list within mdmnets package.
```r
colnames(tax_table(phylo_for_net[[1]])) = c("Kingdom", "Phylum", "Class","Order", "Family", "Genus", "Species")
tax_info = tax_table(phylo_for_net[[1]])
taxa_info$Phylum[taxa_info$Phylum %in% unassignedOTUs_MDM] <- "<NA>"
    taxa_info$Class[taxa_info$Class %in% unassignedOTUs_MDM] <- "<NA>"
    taxa_info$Order[taxa_info$Order %in% unassignedOTUs_MDM] <- "<NA>"
    taxa_info$Family[taxa_info$Family %in% unassignedOTUs_MDM] <- "<NA>"
    taxa_info$Genus[taxa_info$Genus %in% unassignedOTUs_MDM] <- "<NA>"
    taxa_info$Species[taxa_info$Species %in% unassignedOTUs_MDM] <- "<NA>"
    taxa_info[] <- lapply(taxa_info, as.character) #convert to character
    taxa_info[] <- lapply(taxa_info, function(x) { str_replace_all(x, "uncultured", NA_character_)}) #replace all uncultured to NA
    taxa_info <- as.matrix(taxa_info) 
    taxa2 <- tax_table(taxa_info) #convert to phyloseq format
    final_phylo <- merge_phyloseq(otu_table(phylo_for_net[[1]]), taxa2) 
#final_phylo = merged modified taxonomic info with otu info 
```

Now visualize networks using the `get_net_plots_all_ranks` function using  _phylo_for_net_ and _final_phylo_ as input. 
```r
#First double check that taxonomy names are "Kingdom" -> "Species" and rename if necessary
colnames(tax_table(final_phylo)) = c("Kingdom", "Phylum", "Class","Order", "Family", "Genus", "Species")
#now create a network at each level from Phylum-Genus
phylo_nets <- get_net_plots_all_ranks(orig_graph = phylo_for_net[[2]], orig_phylo = final_phylo, met_name = "Environment/Data name here")
```
_Phylo_nets_ will be a list of 5 network plots, from Phylum to Genus, with nodes colored by taxonomic rank level. 
To illustrate this example, we will use the provided _hs_graph_ and _new_hs_phyloseqobj_final_ files from hot springs.
```r 
hs_nets <- get_net_plots_all_ranks(hs_graph, new_hs_phyloseqobj_final)
```
All plots can then be visualized at once using the `plot_grid` function from `cowplot` 
```r
library(cowplot)
plot_grid(plotlist = hs_nets)
```
![net_plots](https://github.com/ConesaLab/MDM/blob/master/hs_all_net_plots.png)

## 4. Create Networks without MDM <a name="netswomdm"> </a>
To create networks without MDM, first find and remove any MDM at taxonomic rank. To do so, use the functions  `met_wo_unk` and `get_graph_wo_unk`, in that order. 
```r 
phylo_wo_unk <- met_wo_unk(final_phylo)
phylo_graph_wo_unk <- get_graph_wo_unk(orig_graph = phylo_for_net[[2]], met_wo_unk = phylo_wo_unk)
```
_phylo_graph_wo_unk_ is a list of 7 igraph objects, from Kingdom to Species, each with no MDM/MDM-related edges. 

Here, we illustrate this example with hot springs data by using the provided _new_hs_phyloseq_final_ and _hs_graph_ files.
```r 
hs_wo_unk <- met_wo_unk(new_hs_phyloseqobj_final)
hs_graph_wo_unk <- get_graph_wo_unk(hs_graph, hs_wo_unk)
```
Visualize resultant networks without MDM using the same function (`get_net_plots_all_ranks`) as previously. 
Using the hot springs data as an example:
```r 
hs_nets_wo_mdm <- lapply(hs_graph_wo_unk, function(orig_graph){ 
new_nets <- get_net_plots_all_ranks(orig_graph, new_hs_phyloseqobj_final, "HS")
})

#now, to visually compare networks at the same level, use grid.arrange 
library(gridExtra)
grid.arrange(hs_nets$Class, hs_nets_wo_mdm$Class$Class)
```
![class_plots](https://github.com/ConesaLab/MDM/blob/master/hs_net_class_plots_w_and_wo_mdm.png)

## 5. Calculate Network Measure Changes <a name="meas"> </a>
Now that networks with and without MDM have been created, we can now evaluate changes in the degree, betweenness, and closeness scores of each node. To calculate these network centrality scores, we use the _`igraph`_ package functions _`degree`_, _`betweenness`_, and _`closeness`_ within our function `degree_calc_f`. For this function, a graph object is input and as output, we have a dataframe of all three centrality measure values for each node present in the network. To show examples for networks with and without MDM, we will use as input _hs_graph_ and _hs_graph_wo_unk_. 
```r
#Calculate degree, closeness and betweenness for original graph including MDM
hs_degree_df <- degree_calc_f(hs_graph)
#hs_degree_df is dataframe, where rows are taxonomic id's and columns are degree, names, betweenness, and closeness. 

#Calculate degree, closeness, and betweenness for graphs excluding MDM (for all levels at once)
hs_degree_df_wo_unk <- lapply(hs_graph_wo_unk, function(orig_graph){
    df <- degree_calc_f(orig_graph)
})
#result will now be 7 dataframes, for each taxonomic level and the number of rows will be the number of taxa remaining (that were not MDM) at that level
```
## 6. Create Bootstrap Networks <a name="boot"> </a>
To make sure that changes are occurring due to the removal of MDM themselves and not due to the number of nodes removed, we randomly subsample 100 times, removing random knowns from the original network,and re-evaluate degree,closeness, and betweenness for each network without random knowns using the function `comp_by_deleting_random_knowns_t_v3`. 
As input, this function calls for a graph object of the network with MDM (_orig_graph_), the list of graphs without MDM (_new_wo_unk_graph_), the phyloseq object (orig_phylo), the dataframe of centrality scores for the network with MDM (orig_df), and the list of dataframes of the networks without MDM (degree_df_wo_unk). The number of iterations is currently default at 100 but may be increased or decreased if necessary. 
 
This step is computationally intensive and as default will use **2** cores (mc.coreval=2) to work in parallel. Up to **4** cores may be specified when working on a typical MacOSX RStudio environment but if necessary, for large datasets, an HPC environment and 4+ cores may be specified. 

As output, we have a dataframe of all centrality scores for the original network with MDM, without MDM, and without random knowns that we can statistically evaluate. The dataframe has 4 columns - **rank_level** (Kingdom -> Species evaluated), **type** (Original, Without_Unknown, and Bootstrap- the network types evaluated), **measure** (degree, bw (betweenness), and closeness), and **data** (numeric value) which reflect the taxonomic rank, network type, network centrality metric,and value of the network centrality metric respectively. 

We will demonstrate this using the hot springs data (hs_graph, hs_graph_wo_unk, new_hs_phyloseqobj_final, hs_degree_df, hs_degree_df_wo_unk) as input, keeping the default iteration/core values.
```r
hs_comp_net_df <- comp_by_deleting_random_knowns_t_v3(hs_graph, hs_graph_wo_unk, new_hs_phyloseqobj_final, hs_degree_df, hs_degree_df_wo_unk)
head(hs_comp_net_df)
 rank_level     type measure data
1    Kingdom Original  degree    5
2    Kingdom Original  degree    2
3    Kingdom Original  degree    3
4    Kingdom Original  degree    3
5    Kingdom Original  degree    1
6    Kingdom Original  degree    3
```
## 7. Visualize and Statistically Evaluate Network Measure Changes <a name="boxplot"> </a>
Now that we have calculated degree, betweenness, and closeness for all three network types, we can accurately statistically evaluate whether or not removing MDM had a significant effect on any network measure. To do so, we use the function `vis_comp_net_meas_boxplots2`. 
As input, we use the dataframe produced from the function `comp_by_deleting_random_knowns_t_v3`. and specify the environment/data name we would like. As default, all measures are plotted but if you want to examine just one or two measures at a time, you may write `meas_name_options = "degree"` or `meas_name_options = c("degree", "betweenness")`. As default, all rank levels from Phylum to Genus are plotted as well. Again, to view just one or two ranks at a time, you may write in the rank levels you would like. Use `p.signif` or `p.format` for the **`p_label`** parameter to see the _stars of significance_ or _exact p-value_ respectively. The default size of the p-value/star is set to 2 but may be decreased or increased to improve figure readability. 
```r
hs_net_boxplots <- vis_comp_net_meas_boxplots2(hs_comp_net_df, "HS", p_label = "p.signif")
```
![hs_boxplots](https://github.com/ConesaLab/MDM/blob/master/hs_comp_net_meas_boxplots.png)

## 8. Create Hub Networks <a name="hub"> </a>
To evaluate hub networks at each level, we use the function `get_hub_plots_all_ranks` and use as input a graph and phyloseq object.
5 network plots are produced, for Phylum to Genus, with nodes colored by the respective taxonomic rank level. 
```r
hs_hub_plots <- get_hub_plots_all_ranks(hs_graph, new_hs_phyloseqobj_final, "HS")
```
To visualize all hub networks at once, use the `plot_grid` function from **`cowplot`** as previously.
```r
plot_grid(plotlist=hs_hub_plots)
```
![hub_plots](https://github.com/ConesaLab/MDM/blob/master/hs_all_hub_plots.png)

## 9. Compare Node Neighbor Interactions <a name="neigh"> </a>
To find with which taxa each node most frequently shares edges (which taxa co-occur with each other most), we can use the function `get_bar_and_dens_class_interactions`. 
With this function, we can produce barplots and pie plots of the most frequent neighbors of each node, for both known and unknown taxa at Class level. We can also plot the density distribution comparing the self-self interactions for known and unknown taxa and self-self interactions of unknown taxa compared to all other interactions present in the network. 


## 10. Validate Results <a name="val"> </a>
To make sure that the observations found (resultant networks, changes in network measures, neighbors, etc) are not due to biases like samples chosen, network tool chosen, or percent sample prevalence chosen, users can validate and compare results using different network correlation/regression tools and a range of percent sample prevalences. We provide the function `vis_comp_net_meas_checks` to statistically evaluate and visualize these changes using boxplots and also provide the function `create_abund_heatmap` to create heatmaps of OTU abundance across samples. We hope that these checks help with your data analysis. 

