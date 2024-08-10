library(EnrichedHeatmap)
library(tidyverse)
library(GenomicRanges)

################################################################################
# I want to visualize differenrially methylated bases per cpg island on chr 2
################################################################################

# what I need:

# a file with cpg islands on chr2

# a file with differentially methylated bases

# remarks:
# we make use of destranded data

################################################################################

# data import

######################################## example data

######################################## cpg islands

# we create a GRanges list with 4 Granges object by type
# we did the same in scrip 03
cpg <- split(annot_cpg, annot_cpg@elementMetadata$type)

# cpg islands is the second element of the list

# we transform the GRnges object into a data frame to be able to edit is easily
# we limit the file to chr 2 only
# we keep only relevant columns
# we then transfrom the data back to a GRanges object because that is
# what the package EnrichedHeatmap requires

# extracting cpg islands

cpg_islands <- cpg[[2]]

# limit to chr and keep only required variables

df_cpg_islands_chr2 <- data.frame(cpg$hg19_cpg_islands) %>% 
  # keep only chr2
  dplyr::filter(seqnames=="chr2") %>% 
  # keep only necessary columns
  dplyr::select(seqnames, start, end, strand, id, type) %>% 
  arrange(id)

# transfrom to GRanges object again
cpg_islands_GRanges <- makeGRangesFromDataFrame(df_cpg_islands_chr2,
                                                keep.extra.columns = TRUE,
                                                seqnames.field = c("seqnames", "chr"),
                                                start.field = "start",
                                                end.field = "end")
cpg_islands_GRanges

################################### differentially methylated bases

# we could keep the differences and use those as scores or values to be visualized
# instead we assign the value of +1 to differentially hypermethylated bases
# -1 to differentially hypomethylated bases
# tumor vs. normal

# again, we transform the differentially methylated data which is a methylKit object to datafram

class(myDiff25_all_bases)

# [1] "methylDiff"
# attr(,"package")
# [1] "methylKit"

diff_methylated_bases <- data.frame(myDiff25_all_bases) %>% 
  dplyr::mutate(value=case_when(meth.diff<0 ~ -1,
                                meth.diff>0 ~ 1)) %>% 
  dplyr::select(chr, start, end, value) 


# then transform it to a GRnages onject
diff_meth_bases_GRanges <- makeGRangesFromDataFrame(diff_methylated_bases,
                           keep.extra.columns = TRUE,
                           seqnames.field = c("seqnames", "chr"),
                           start.field = "start",
                           end.field = "end")
diff_meth_bases_GRanges

###############################################################################
# heatmap
###############################################################################


# step 1: prepare base matrix with smoothing
# instead of NA we use 0 (background argument) for base coordinates for which we do not have info
# focus on cpg islands
# and surrounding of 1000 bases upstream and 1000 downstream
# window width=50

mat <- normalizeToMatrix(signal=diff_meth_bases_GRanges,
                         target = cpg_islands_GRanges,
                         value_column = "value", mean_mode = "absolute",
                         extend = 1000, w = 50, background = 0, smooth = T)

EnrichedHeatmap(mat, col = c("blue", "white", "red"), name = "methylation", axis_name_rot = 90,
                column_title = "methylation near CpG islands")

# we see different types of methylation patterns
# 1. cpg island hypermethylated, surrounding also
# 2. no differentally methylated C's (neither in cpg islands, nor in surrounding)
# 3. cpg island not differentially methylated, surrounding hypomethylated
# 4. cpg island and surrounding hypomethylated

# we could try to identify clusters, maybe these 4 initially

set.seed(123)
EnrichedHeatmap(mat, col = c("blue", "white", "red"), name = "methylation", row_km = 4,
                top_annotation = HeatmapAnnotation(enriched = anno_enriched(gp = gpar(col = 2:4, lty = 1:3))),
                column_title = "methylation in and near CpG islands", row_title_rot = 0)

set.seed(123)
EnrichedHeatmap(mat, col = c("blue", "white", "red"), name = "methylation", row_km = 3,
                top_annotation = HeatmapAnnotation(enriched = anno_enriched(gp = gpar(col = 2:4, lty = 1:3))),
                column_title = "methylation in and near CpG islands", row_title_rot = 0)