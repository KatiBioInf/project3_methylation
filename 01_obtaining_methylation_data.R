# import libraries
# if you do not have these installed, install these first

library(GenomicRanges)
library(rtracklayer)
library(tidyverse)

# you need to have the original bed files (Colon_Primary_Tumor.BiSeq.bed, 
# Colon_Adjacent_Normal.BiSeq.bed) downloaded and saved in the same map
# as this project

################################################################################
#                                 Obtaining data
################################################################################    

# from WGBS experiment
# bioproject id: PRJNA201480
# https://www.ncbi.nlm.nih.gov/bioproject/PRJNA201480

# instrument: Illumina HiSeq 2000
# Genome_build: hg19
# Supplementary_files_format_and_content: bed file containing all seen CpGs 
# within this library. The number of methylated reads/number of total reads is listed in the score column

# all data
# https://www.ncbi.nlm.nih.gov/gds/?term=GSE46644[ACCN]%20AND%20gsm[ETYP]

# primary colon tumor
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1204465
# tumor adjacent normal
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1204466

# import bed files using the import.bed function from rtracklayer package

tumor_bed <- import.bed(con="Colon_Primary_Tumor.BiSeq.bed")
normal_bed <- import.bed(con="Colon_Adjacent_Normal.BiSeq.bed")

# inspecting data

tumor_bed

# GRanges object with 27631247 ranges and 2 metadata columns:
#   seqnames            ranges strand |        name     score
# <Rle>         <IRanges>  <Rle> | <character> <numeric>
#   [1]     chr1       10469-10470      + |       '7/7'      1000
# [2]     chr1       10471-10472      + |       '8/8'      1000
# [3]     chr1       10484-10485      + |       '7/8'       875
# [4]     chr1       10489-10490      + |       '9/9'      1000
# [5]     chr1       10493-10494      + |      '9/11'       818
# ...      ...               ...    ... .         ...       ...
# [27631243]     chrY 59033368-59033369      - |       '1/1'      1000
# [27631244]     chrY 59033547-59033548      + |       '8/8'      1000
# [27631245]     chrY 59033727-59033728      + |       '8/9'       889
# [27631246]     chrY 59033729-59033730      + |       '6/9'       667
# [27631247]     chrY 59033848-59033849      + |       '4/8'       500
# -------
#   seqinfo: 24 sequences from an unspecified genome; no seqlengths

normal_bed

# GRanges object with 27381577 ranges and 2 metadata columns:
#   seqnames            ranges strand |        name     score
# <Rle>         <IRanges>  <Rle> | <character> <numeric>
#   [1]     chr1       10469-10470      - |       '2/2'      1000
# [2]     chr1       10471-10472      - |       '1/2'       500
# [3]     chr1       10484-10485      + |       '5/5'      1000
# [4]     chr1       10489-10490      - |       '4/4'      1000
# [5]     chr1       10493-10494      + |       '3/5'       600
# ...      ...               ...    ... .         ...       ...
# [27381573]     chrY 59033368-59033369      + |       '3/3'      1000
# [27381574]     chrY 59033547-59033548      + |       '4/5'       800
# [27381575]     chrY 59033727-59033728      + |     '16/17'       941
# [27381576]     chrY 59033729-59033730      + |     '17/18'       944
# [27381577]     chrY 59033848-59033849      + |       '3/6'       500
# -------
#   seqinfo: 24 sequences from an unspecified genome; no seqlengths

################################################################################
#                                 data preprocessing
################################################################################

# how should it look like?
# the structure of the data in the MethylKit tutorial looks a bit different from our imported data
# goal: transform our data into a similar format as in the MethylKit tutorial

# this is from the methylkit exdata site: https://github.com/al2na/methylKit/blob/master/inst/extdata/test2.myCpG.txt

# chrBase	chr	base	strand	coverage	freqC	freqT
# chr21.9764513	chr21	9764513	R	20	75.00	25.00
# chr21.9764539	chr21	9764539	R	20	95.00	5.00
# chr21.9767701	chr21	9767701	F	10	10.00	90.00
# chr21.9804584	chr21	9804584	F	10	90.00	10.00
# chr21.9804595	chr21	9804595	F	10	100.00	0.00
# chr21.9802620	chr21	9802620	F	13	61.54	38.46
# chr21.9802590	chr21	9802590	F	13	76.92	23.08


# we turn our data into a data frame
# we focus only on chromosome 2, otherwise the data are too large
df_tumor <- data.frame(tumor_bed) %>% dplyr::filter(seqnames=="chr2")
df_normal <- data.frame(normal_bed) %>% dplyr::filter(seqnames=="chr2")

# we can remove the original files from the Environment
rm(normal_bed, tumor_bed)

# is the width always 2? yes, it means that we have info on the level of bases

df_tumor %>% dplyr::group_by(width) %>% dplyr::summarise(N=n())

# # A tibble: 1 × 2
# width       N
# <int>   <int>
#   1     2 2123600

df_normal %>% dplyr::group_by(width) %>% dplyr::summarise(N=n())

# # A tibble: 1 × 2
# width       N
# <int>   <int>
#   1     2 2100842

############################## fixing positional info

# rename chromosome column
# in the final table we will need only the start coordinate because 
# we analyse bases (start and end coordinate are the same)
# we remove score column because we will use our own measure

df_tumor <- df_tumor %>% rename(chr=seqnames) %>% dplyr::select(-end, -width, -score)
df_normal <- df_normal %>% rename(chr=seqnames) %>% dplyr::select(-end, -width, -score)

# strand info is good as it is

############################# fixing methylation info

# splitting values in name column
# first value is methylated bases (C's)
# second value is total reads (coverage)

df_tumor <- df_tumor %>% 
  rowwise() %>% 
  dplyr::mutate(countC=strsplit(name, "[/]")[[1]][1],
                countTotal=strsplit(name, "[/]")[[1]][2]) %>% 
  ungroup()


df_normal <- df_normal %>% 
  rowwise() %>% 
  dplyr::mutate(countC=strsplit(name, "[/]")[[1]][1],
                countTotal=strsplit(name, "[/]")[[1]][2]) %>% 
  ungroup()

# cleaning up and turning into numeric

df_tumor <- df_tumor %>% 
  dplyr::mutate(countC=str_remove(countC, pattern = "'"),
                coverage=str_remove(countTotal, pattern = "'")) %>% 
  dplyr::mutate(countC=as.numeric(countC),
                coverage=as.numeric(coverage))

df_normal <- df_normal %>% 
  dplyr::mutate(countC=str_remove(countC, pattern = "'"),
                coverage=str_remove(countTotal, pattern = "'")) %>% 
  dplyr::mutate(countC=as.numeric(countC),
                coverage=as.numeric(coverage))

########################## final data

# percentages

df_tumor <- df_tumor %>% 
  select(-name, -countTotal) %>% 
  dplyr::mutate(freqC=round(countC*100/coverage, 2),
                freqT=100-freqC) %>% 
  select(-countC) %>% 
  select( chr, start, strand, coverage, freqC, freqT)

df_normal <- df_normal %>% 
  select(-name, -countTotal) %>% 
  dplyr::mutate(freqC=round(countC*100/coverage, 2),
                freqT=100-freqC) %>% 
  select(-countC) %>% 
  select( chr, start, strand, coverage, freqC, freqT)

# # saving data as txt
# 
# write.table(df_tumor, quote = F, row.names = F, "df_tumor.txt")
# write.table(df_normal, quote = F, row.names = F, "df_normal.txt")