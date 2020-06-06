# How to use `mdmnets` package to analyze data
## Table of Contents
1. [Load Data into R and convert to phyloseq objects](#intro)
    * [Loading QIIME/OTU data](#qiimeimport) 
    * [Loading dada2 ASV data](#dada2import)
    * [Loading QIIME2 ASV data](#qiime2import)
2. [Calculate Known vs Unknown (MDM) Proportions](#prop)
3. [Create Networks with MDM](#netswmdm)
4. [Create Networks without MDM](#netswomdm)
5. [Create Bootstrap Networks](#boot)
6. [Compare Network Measure Changes](#meas)
7. [Create Hub Networks](#hub)
8. [Compare Node Neighbor Interactions](#neigh)
9. [Validate Results](#val)

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
If data was made with dada2 and in ASV format, data can be converted like so, using [dada2 tutorial][dada2guide] as a guide: 

[dada2guide]: (https://benjjneb.github.io/dada2/tutorial.html) 
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

## 3. Create Networks with MDM <a name="netswmdm"> </a>

## 4. Create Networks without MDM <a name="netswomdm"> </a>
## 5. Create Bootstrap Networks <a name="boot"> </a>
## 6. Compare Network Measure Changes <a name="meas"> </a>
## 7. Create Hub Networks <a name="hub"> </a>
## 8. Compare Node Neighbor Interactions <a name="neigh"> </a>
## 9. Validate Results <a name="val"> </a>


