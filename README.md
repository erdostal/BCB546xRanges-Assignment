---
title: "BCB546x Ranges Assignment"
author: "Emma Miller"
date: "April 17, 2017"
output: pdf_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

#Load required packages for this pipeline
```{r}
library(rtracklayer)
library(IRanges)
library(GenomicRanges)
library("GenomicFeatures")
```

#take mouse annotation library and turn it into an object in my environment (shortcut)
```{r}
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene
```

#Import file containing variant information for chromosome 1 of mouse
```{r}
dbsnp137 <- import("mm10_snp137_chr1_trunc.bed.gz")
```

#Subset and reduce mouse library to exons found only in Chromosome 1
```{r}
collapsed_exons <- reduce(exons(txdb), ignore.strand=TRUE)
chr1_collapsed_exons <-collapsed_exons
```


#inspect data, it is not trustworthy. Summary shows that some values have a width of 0
```{r}
summary(width(dbsnp137))
```

#Remove ranges with a width of 0 since they cannot have any overlap
```{r}
dbsnp137_resized <- dbsnp137
zw_i <- width(dbsnp137_resized) == 0
dbsnp137_resized[zw_i] <- resize(dbsnp137_resized[zw_i], width=1)
```

#With adjust ranges, pull out variants that overlap chr1 regions
```{r}
hits <- findOverlaps(dbsnp137_resized, chr1_collapsed_exons, 
                       ignore.strand=TRUE)
#examine hits file
hits
```
#number of variants that overlap with exons
```{r}
length(unique(queryHits(hits)))
#percentage of variants from chromosome 1 represented 
length(unique(queryHits(hits)))/length(dbsnp137_resized)
```
#Get a closer look at the variants which overlap in these regions
```{r}
subsetByOverlaps(dbsnp137_resized, chr1_collapsed_exons, ignore.strand=TRUE)
```

#count the number of variants which had overlap in each exonic region
```{r}
var_counts <- countOverlaps(chr1_collapsed_exons, dbsnp137_resized, ignore.strand=TRUE)
#Add the counts of variants to our exon file
chr1_collapsed_exons$num_vars <- var_counts
#examine file 
chr1_collapsed_exons
```
