# differential methylation ignoring strand info
# differential analysis on the level of bases
# then differential analysis at the level of cpg islands

library(methylKit)
library(tidyverse)
library(writexl)
library(GenomicRanges)

# loading in files

file.list=list( "df_tumor.zip",
                "df_normal.zip" )

# read the files to a methylRawList object: myobj_full

# we limit coverage
myobj=methRead(file.list,
               pipeline=list(fraction=FALSE, chr.col=1,start.col=2,end.col=2,coverage.col=4,
                             strand.col=3,freqC.col=5),
               sample.id=list("test", "ctrl"),
               assembly="hg19",
               treatment=c(1, 0),
               context="CpG",
               sep=" ",
               mincov = 10
)

##### inspect object

myobj

# methylRawList object with 2 methylRaw objects
# 
# methylRaw object with 1761485 rows
# --------------
#   chr start   end strand coverage numCs numTs
# 1 chr2 10270 10270      -       11    11     0
# 2 chr2 10360 10360      -       27    27     0
# 3 chr2 10379 10379      -       28    26     2
# 4 chr2 10385 10385      -       25    25     0
# 5 chr2 10391 10391      -       24    23     1
# 6 chr2 10398 10398      -       24    24     0
# --------------
#   sample.id: test 
# assembly: hg19 
# context: CpG 
# resolution: base 
# 
# methylRaw object with 1674449 rows
# --------------
#   chr start   end strand coverage numCs numTs
# 1 chr2 10360 10360      -       19    19     0
# 2 chr2 10379 10379      -       20    19     1
# 3 chr2 10385 10385      -       17    17     0
# 4 chr2 10391 10391      -       17    17     0
# 5 chr2 10398 10398      -       15    15     0
# 6 chr2 10405 10405      -       11    11     0
# --------------
#   sample.id: ctrl 
# assembly: hg19 
# context: CpG 
# resolution: base 
# 
# treatment: 1 0 

# # inspecting the content per object
# 
# myobj@.Data[[1]]
# myobj@.Data[[2]]

# methylation statistics

# the distribution of methylation percentages
# ignoring strand info
par(mfrow=c(1,2))
getMethylationStats(myobj[[1]],plot=TRUE,both.strands=F)
getMethylationStats(myobj[[2]],plot=TRUE,both.strands=F)

# read coverage
# numbers on bars denote what percentage of locations are contained in that bin
par(mfrow=c(1,2))
getCoverageStats(myobj[[1]],plot=TRUE,both.strands=F)
getCoverageStats(myobj[[2]],plot=TRUE,both.strands=F)

# we might want to filter too high and too low coverage
# we 

filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL,
                                hi.count=NULL,hi.perc=99.9)

# check the difference between the unfiltered and filtered datasets

for (i in 1:2){
  print(i)
  print(length(myobj@.Data[[i]]$chr))
  print(length(filtered.myobj@.Data[[i]]$chr))
}

# [1] 1
# [1] 1761485
# [1] 1759659
# [1] 2
# [1] 1674449
# [1] 1672724

################################################################################

# Comparative analysis

################################################################################

###################### dataprep

# we have our samples (tumor, normal) apart -> we want them in one file
# we want to have info of methylation per base

# unite function
# output: methylBase object
# it keeps bases which are present in both samples, you can change this setting

# meth_withstrands=methylKit::unite(myobj, destrand=FALSE)
meth_destranded=methylKit::unite(myobj, destrand=TRUE)

class(meth_destranded)

# [1] "methylBase"
# attr(,"package")
# [1] "methylKit"

# inspect data

head(meth_destranded)

# methylBase object with 6 rows
# --------------
#   chr start   end strand coverage1 numCs1 numTs1 coverage2 numCs2 numTs2
# 1 chr2 10359 10359      +        27     27      0        19     19      0
# 2 chr2 10378 10378      +        28     26      2        20     19      1
# 3 chr2 10384 10384      +        25     25      0        17     17      0
# 4 chr2 10390 10390      +        24     23      1        17     17      0
# 5 chr2 10397 10397      +        24     24      0        15     15      0
# 6 chr2 10404 10404      +        20     19      1        11     11      0
# --------------
#   sample.ids: test ctrl 
# destranded TRUE 
# assembly: hg19 
# context: CpG 
# treament: 1 0 
# resolution: base 

getSampleID(meth_destranded)
# column names
names(meth_destranded)
# what type of methylation do we have
getContext(meth_destranded)
# column id's of the number of methylated C's per base
meth_destranded@numCs.index

# [1] 6 9


###################### sample correlation

getCorrelation(meth_destranded,plot=TRUE)

# test      ctrl
# test 1.0000000 0.8455818
# ctrl 0.8455818 1.0000000


################################################################################

# differentially methylated bases

################################################################################

# analysis

# in this section we get differentially methylated base pairs keeping strand awareness
# thesholds: difference 25, q-value: 0.01
# we can perform only fisher exact test because we have one replicate per group
# we can use regression based methods if we have multiple replicates per group

# fisher exact test
new.meth=reorganize(meth_destranded,sample.ids=c("test","ctrl"),treatment=c(1,0))
myDiff_fisher_base <- calculateDiffMeth(new.meth)


# more specific results: the function serves only to filter the results

# get all differentially methylated bases
# where the difference is at least 25 and q-value is is lower than 0.01
myDiff25_all_bases <- getMethylDiff(myDiff_fisher_base,difference=25,qvalue=0.01)

# # we can get hyper and hypomethylated bases apart, if we want
# myDiff25_hyper_bases <- getMethylDiff(myDiff_fisher_base,difference=25,qvalue=0.01, type="hyper")
# myDiff25_hypo_bases <- getMethylDiff(myDiff_fisher_base,difference=25,qvalue=0.01, type="hypo")

# # we can save these results as RData, if we want
# save(meth_destranded, myDiff_fisher_base, myDiff25_all_bases, myDiff25_hyper_bases,
#      myDiff25_hypo_bases, file="methylation_fisher_base_destranded.RData")

################################################################################

# differentially methylated regions - cpg islands

################################################################################

# next, we ask if there are differentially methylated regions between the tumor and normal sample
# for this we need to aggragte methylated information for cpg islands and perform the same analysis

### read in annotation data

load("annotation_data.RData")

# inspect cpg data

annot_cpg

# GRanges object with 144965 ranges and 5 metadata columns:
#   seqnames        ranges strand |          id     tx_id   gene_id    symbol             type
# <Rle>     <IRanges>  <Rle> | <character> <logical> <logical> <logical>      <character>
#   [1]           chr1   28736-29810      * |    island:1      <NA>      <NA>      <NA> hg19_cpg_islands
# [2]           chr1 135125-135563      * |    island:2      <NA>      <NA>      <NA> hg19_cpg_islands
# [3]           chr1 327791-328229      * |    island:3      <NA>      <NA>      <NA> hg19_cpg_islands
# [4]           chr1 437152-438164      * |    island:4      <NA>      <NA>      <NA> hg19_cpg_islands
# [5]           chr1 449274-450544      * |    island:5      <NA>      <NA>      <NA> hg19_cpg_islands
# ...            ...           ...    ... .         ...       ...       ...       ...              ...
# [144961] chrUn_gl000245       1-36651      * | inter:20609      <NA>      <NA>      <NA>   hg19_cpg_inter
# [144962] chrUn_gl000246       1-38154      * | inter:20610      <NA>      <NA>      <NA>   hg19_cpg_inter
# [144963] chrUn_gl000247       1-36422      * | inter:20611      <NA>      <NA>      <NA>   hg19_cpg_inter
# [144964] chrUn_gl000248       1-39786      * | inter:20612      <NA>      <NA>      <NA>   hg19_cpg_inter
# [144965] chrUn_gl000249       1-38502      * | inter:20613      <NA>      <NA>      <NA>   hg19_cpg_inter
# -------
#   seqinfo: 93 sequences (1 circular) from hg19 genome

seqnames(annot_cpg)

unique(annot_cpg$type)
# "hg19_cpg_islands" "hg19_cpg_shores"  "hg19_cpg_shelves" "hg19_cpg_inter" 

# this data includes all chromosomes and all methylation related features
# we need only hg19_cpg_islands

# we can transform this GRanges object into a GRanges List object based on type
# this results in a list with four GRanges obhjects in it, one per type

# making grangeslist object

cpg_list <- split(annot_cpg, annot_cpg$type)

# the names of elements - 1 per type

names(cpg_list)

# [1] "hg19_cpg_inter"   "hg19_cpg_islands" "hg19_cpg_shelves" "hg19_cpg_shores" 

# inspect first list element, inter
cpg_list[[1]]

# inspect second list element, islands - this is what we need
cpg_list[[2]]



# we make a table with methylation information per cpg island

# region counts over cpg islands

cpg_islands_counts <- regionCounts(myobj,cpg_list$hg19_cpg_islands)

cpg_islands_counts

# methylRawList object with 2 methylRaw objects
# 
# methylRaw object with 1688 rows
# --------------
#   chr  start    end strand coverage numCs numTs
# 1 chr2  45511  46559      *     5040   187  4853
# 2 chr2 197438 198535      *     3032  2921   111
# 3 chr2 214012 214245      *     1533  1406   127
# 4 chr2 263401 265238      *     6094   282  5812
# 5 chr2 287354 290099      *     8236   730  7506
# 6 chr2 306676 307167      *     1453  1268   185
# --------------
#   sample.id: test 
# assembly: hg19 
# context: CpG 
# resolution: region 
# 
# methylRaw object with 1688 rows
# --------------
#   chr  start    end strand coverage numCs numTs
# 1 chr2  45511  46559      *     5625  1419  4206
# 2 chr2 197438 198535      *     3140  3032   108
# 3 chr2 214012 214245      *     1422  1301   121
# 4 chr2 263401 265238      *     5594   930  4664
# 5 chr2 287354 290099      *     7878   575  7303
# 6 chr2 306676 307167      *     1499  1350   149
# --------------
#   sample.id: ctrl 
# assembly: hg19 
# context: CpG 
# resolution: region 
# 
# treatment: 1 0 

# we have one object per sample
# we need to use the unite function again

cpg_islands_counts_flat <- methylKit::unite(cpg_islands_counts, destrand=TRUE)

# perform fisher exact test

myDiff_fisher_cpg_islands <- calculateDiffMeth(cpg_islands_counts_flat)

# get all differentially methylated cpg islands
# where the difference is at least 25 and q-value is is lower than 0.01
myDiff25_all_cpg_islands <- getMethylDiff(myDiff_fisher_cpg_islands,difference=25,qvalue=0.01)

# # we can save these results as RData, if we want
# save(cpg_islands_counts_flat, myDiff_fisher_cpg_islands, myDiff25_all_cpg_islands, 
#  file="methylation_fisher_cpg_islands_destranded.RData")

# for the visualization we need only the following data from the Environment:
# 1. annot_cpg
# 2. myDiff_25_all_bases - differentially methylated bases (diff>=25, q<=0.01)

rm(list=ls()[! ls() %in% c("annot_cpg","myDiff25_all_bases")])