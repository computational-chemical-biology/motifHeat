motif\_heatmap
================

## motifHeat

This is a short tutorial to introduce motifHeat, associating
hierarchical clustering to motif detection.

First we need to install the dependencies

``` r
packageList <- c("devtools", "dendextend")
newPackages <- packageList[!(packageList %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages)

library("devtools")
library("dendextend")
```

    ## 
    ## ---------------------
    ## Welcome to dendextend version 1.12.0
    ## Type citation('dendextend') for how to cite the package.
    ## 
    ## Type browseVignettes(package = 'dendextend') for the package vignette.
    ## The github page is: https://github.com/talgalili/dendextend/
    ## 
    ## Suggestions and bug-reports can be submitted at: https://github.com/talgalili/dendextend/issues
    ## Or contact: <tal.galili@gmail.com>
    ## 
    ##  To suppress this message use:  suppressPackageStartupMessages(library(dendextend))
    ## ---------------------

    ## 
    ## Attaching package: 'dendextend'

    ## The following object is masked from 'package:stats':
    ## 
    ##     cutree

``` r
if (!require("motifHeat")) {
    install_github("computational-chemical-biology/motifHeat")
} else {
    library("motifHeat")
}
```

    ## Loading required package: motifHeat

## Load data from GNPS

You need your GNPS task id and optionally
ms2lda:

``` r
dlist <- access_gnps('60e8f7096e994d5c9717d471445f8bfe',            '56f089a1772340a0aa2978d41478e3e4')
```

    ## Downloading MZmine feature table...
    ## Downloading metadata table...
    ## Downloading ms2lda nodes table...
    ## Downloading ms2lda motifs table...
    ## Downloading ms2lda edges table...

## Recover most frequent motifs

Display motifs
    annotated:

``` r
head(dlist$motifs_in_scans)
```

    ##   scan precursor.mass retention.time                   motif probability
    ## 1  593       670.4670              0 euphorbia_motif_158.m2m   0.1237446
    ## 2 1200       548.3069              0               motif_333   0.4879933
    ## 3 1200       548.3069              0               motif_374   0.1437222
    ## 4 4024       586.9921              0               motif_521   0.8134345
    ## 5 4025       599.0103              0               motif_521   0.1510719
    ## 6 4020       341.3160              0       gnps_motif_43.m2m   0.2011264
    ##     overlap                            motifdb_url
    ## 1 0.9915291 http://ms2lda.org/motifdb/motif/151115
    ## 2 0.3328333                                       
    ## 3 0.8781329                                       
    ## 4 0.7017773                                       
    ## 5 0.5447280                                       
    ## 6 0.9999885 http://ms2lda.org/motifdb/motif/151021
    ##                                                                             motifdb_annotation
    ## 1                                                                   Loss of NH3 adducts in DSF
    ## 2                                                                                             
    ## 3                                                                                             
    ## 4                                                                                             
    ## 5                                                                                             
    ## 6 Water loss - indicative of a free hydroxyl group – (in beer often seen in sugary structures)

Count the motifs and display the most frequent

``` r
motif_count <- table(dlist$motifs_in_scans[,'motif'])
motif_count[order(motif_count, decreasing=TRUE)][1:5]
```

    ## 
    ##        motif_480 gnps_motif_7.m2m        motif_295 gnps_motif_4.m2m 
    ##              330              273              258              200 
    ##        motif_337 
    ##              149

Select cluster indexes associated to most
frequent

``` r
smotif <- dlist$motifs_in_scans[dlist$motifs_in_scans[,'motif'] %in% c('motif_480'), c('scan', 'motif')]
head(smotif)
```

    ##     scan     motif
    ## 92  2442 motif_480
    ## 97  2446 motif_480
    ## 106 2449 motif_480
    ## 163 3389 motif_480
    ## 202 1199 motif_480
    ## 206 1191 motif_480

## Format color for sample classes

factorCorLis is a nested list with one color for each level of a
factor:

``` r
factorColList <- list(list(colors=c("darkorchid","darkred"), factor='Extraction'), list(colors=c("green", "darkgreen"), factor='Trimethoprim'))
```

## Plot a basic heatmap

You can also embed plots, for example:

``` r
tab <- dlist$features
meta <- dlist$metadata
h <- format_heatmap(tab, meta, selectField='StrainName', selectValue='Burkholderia dolosa AU0645  Genomovar type VI', factorColList=factorColList)
```

![](motif_heatmap_files/figure-gfm/unnamed-chunk-7-1.png)<!-- --> We can
also add motif labels to the features in the
heatmap:

``` r
h <- format_heatmap(tab, meta, selectField='StrainName', labCol=smotif, selectValue='Burkholderia dolosa AU0645  Genomovar type VI', factorColList=factorColList)
```

![](motif_heatmap_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->
Subselect metadata and repeat
heatmap:

``` r
meta2 <- meta[meta$TimePoint=='48h' & grepl('EA', meta$Extraction) & grepl('Burkholderia', meta$StrainName),]
factorColList <- list(list(colors=rainbow(length(unique(meta2$StrainName))), factor='StrainName'))
h <- format_heatmap(tab, meta2, labCol=smotif, factorColList=factorColList)
```

![](motif_heatmap_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->
Changing color scale:

``` r
library("gplots")
```

    ## 
    ## Attaching package: 'gplots'

    ## The following object is masked from 'package:stats':
    ## 
    ##     lowess

``` r
h <- format_heatmap(tab, meta2, labCol=smotif, factorColList=factorColList, colorScale = redgreen(75))
```

![](motif_heatmap_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->
Adding
normalization:

``` r
h <- format_heatmap(tab, meta2, labCol=smotif, norm=TRUE, factorColList=factorColList, colorScale = redgreen(75))
```

![](motif_heatmap_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->
Creating a pdf:

    ## png 
    ##   2

Repeating the process for other condition:

    ## png 
    ##   2

## Recover data from heatmap’s hierarchical clustering

The heatmap result contains hierarchical clustering from samples and
features. We can associate the sample grouping to colors

``` r
# needs preprocessing
head(cbind(t(h$colors), colnames(tab)[grep('Peak area', colnames(tab))])[h$heatmap$rowInd,])
```

    ## Warning in cbind(t(h$colors), colnames(tab)[grep("Peak area",
    ## colnames(tab))]): number of rows of result is not a multiple of vector
    ## length (arg 2)

    ##                                 StrainName 
    ## Burkholderia thailandensis E264 "#CC00FFFF"
    ## Burkholderia thailandensis E264 "#CC00FFFF"
    ## Burkholderia cenocepacia MC0-3  "#00FF66FF"
    ## Burkholderia cenocepacia MC0-3  "#00FF66FF"
    ## Burkholderia cenocepacia MC0-3  "#00FF66FF"
    ## Burkholderia cenocepacia MC0-3  "#00FF66FF"
    ##                                                                                                
    ## Burkholderia thailandensis E264 "TestMix18_P1-F-2_01_2856.mzXML filtered Peak area"            
    ## Burkholderia thailandensis E264 "TestMix7_P1-F-2_01_2704.mzXML filtered Peak area"             
    ## Burkholderia cenocepacia MC0-3  "BccA40645_EA_24_4_P1-B-1_01_2835.mzXML filtered Peak area"    
    ## Burkholderia cenocepacia MC0-3  "BmvCF2_EA_24_TMP_P1-B-5_01_2733.mzXML filtered Peak area"     
    ## Burkholderia cenocepacia MC0-3  "BtE265_EA_48_TMP_P1-E-6_01_2756.mzXML filtered Peak area"     
    ## Burkholderia cenocepacia MC0-3  "BccA40645_EA_24_TMP_2_P1-B-5_01_2776.mzXML filtered Peak area"

Create hclut object from dendrogram and plot the cluters with arbitrary
number of groups

``` r
hc <- as.hclust(h2$heatmap$colDendrogram)
plot(hc)
```

![](motif_heatmap_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
plot(hc, labels=FALSE)
rect.hclust(hc, 6)
```

![](motif_heatmap_files/figure-gfm/unnamed-chunk-15-2.png)<!-- -->

Change the criteria to hight and record grouping

``` r
plot(hc, labels=FALSE)
rect.hclust(hc, h=1500)
```

![](motif_heatmap_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
hcgrp <- cutree(hc, h=1500)
```

Use grouping to create colored bars

``` r
dend <- h2$heatmap$colDendrogram
dend %>% plot 
colored_bars(colors = hcgrp, dend = dend) 
```

![](motif_heatmap_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

## Associate hierachical clustering to network component and motif frequency

Create a table with all component indexes and associated cluster
ids:

``` r
ms2lda <- cbind(c(dlist$ms2lda_edges[, 'CLUSTERID1'], dlist$ms2lda_edges[, 'CLUSTERID2']), rep(dlist$ms2lda_edges[, 'ComponentIndex'], 2))
ms2lda <- unique(ms2lda)
colnames(ms2lda) <- c('CLUSTERID','ComponentIndex')
ms2lda <- ms2lda[ms2lda[,2]!=-1,]
head(ms2lda)
```

    ##      CLUSTERID ComponentIndex
    ## [1,]       725              1
    ## [2,]      3461              2
    ## [3,]      1670              2
    ## [4,]      3718              2
    ## [5,]      1897              2
    ## [6,]      1141              2

Merge component index and moftif
tables

``` r
ms2lda <- merge(ms2lda, dlist$ms2lda_nodes[,c('scans', 'motif')], by.x='CLUSTERID', by.y='scans')
tgrp <- data.frame(CLUSTERID=as.numeric(names(hcgrp)), HCA_GRP=hcgrp) 
ms2lda <- merge(ms2lda, tgrp, by.x='CLUSTERID', by.y='CLUSTERID')
head(ms2lda)
```

    ##   CLUSTERID ComponentIndex                       motif HCA_GRP
    ## 1         1            387                   motif_451       1
    ## 2         2            387                   motif_451       1
    ## 3         3            387                   motif_451       2
    ## 4         5            982                                   1
    ## 5         6              7                                   2
    ## 6         7            167 gnps_motif_43.m2m,motif_451       1

Obtain the proportion of nodes from each connected component that are
associated to the same group detectep by hierarchical
clustering

``` r
comp2grp <- tapply(ms2lda[,4], ms2lda[,2], function(x) c(names(table(x))[which.max(table(x))], max(table(x))/length(x), length(x)))

tcomp2grp <- do.call(rbind, comp2grp)
tcomp2grp <- cbind(names(comp2grp), tcomp2grp)
colnames(tcomp2grp) <- c('ComponentIndex', 'HCA_GRP', 'Max_Proportion', 'Total_Nodes')
tcomp2grp <- tcomp2grp[order(as.numeric(tcomp2grp[,4]), decreasing=TRUE),]
mfmotif <- tapply(ms2lda[,3], ms2lda[,2], function(x) names(which.max(table(unlist(sapply(x, strsplit, ',')))))) 
mfmotif[unlist(lapply(mfmotif, is.null))] <- '' 
tcomp2grp <- cbind(tcomp2grp, '') 
colnames(tcomp2grp)[5] <- 'Most_Freq_Motif' 
tcomp2grp[,5] <- unlist(mfmotif[tcomp2grp[,1]]) 
ulinks <- unique(dlist$motifs_in_scans[,c('motif', 'motifdb_url')]) 
ulinks <- as.matrix(ulinks) 
tcomp2grp <- cbind(tcomp2grp, '')
ids <- match(tcomp2grp[,5], ulinks[,1])
tcomp2grp[!is.na(ids), 6] <- ulinks[ids[!is.na(ids)],2]
head(tcomp2grp)
```

    ##     ComponentIndex HCA_GRP Max_Proportion      Total_Nodes Most_Freq_Motif
    ## 2   "2"            "3"     "0.519607843137255" "102"       "motif_295"    
    ## 7   "7"            "4"     "0.3625"            "80"        "motif_391"    
    ## 60  "60"           "2"     "0.513157894736842" "76"        "motif_486"    
    ## 12  "12"           "1"     "0.4"               "70"        "motif_295"    
    ## 51  "51"           "6"     "0.883720930232558" "43"        "motif_477"    
    ## 297 "297"          "1"     "0.395348837209302" "43"        "motif_448"    
    ##       
    ## 2   ""
    ## 7   ""
    ## 60  ""
    ## 12  ""
    ## 51  ""
    ## 297 ""

Save table

``` r
colnames(tcomp2grp)[6] <- 'Link'
write.table(tcomp2grp, 'component_group_association.tsv', sep='\t', row.names=FALSE)
```
