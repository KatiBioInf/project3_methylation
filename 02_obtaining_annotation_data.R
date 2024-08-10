# we need information on where the cpg islands are located on the human genome
# we use the annotatr package

# https://www.bioconductor.org/packages/release/bioc/manuals/annotatr/man/annotatr.pdf
# https://www.bioconductor.org/packages/release/bioc/vignettes/annotatr/inst/doc/annotatr-vignette.html#cpg-annotations

# library(AnnotationHub)
# BiocManager::install("annotatr")
library(annotatr)

# getting cpg island

annotation_cpg = c("hg19_cpgs")
annot_cpg = build_annotations(genome="hg19", annotations=annotation_cpg)
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

# # if you need gene annotations as well
# # getting genes
# annotation_genes = c("hg19_basicgenes")
# annot_genes = build_annotations(genome="hg19", annotations=annotation_genes)

# # save annotations data
# 
# save(annot_cpg, file = "annotation_data.RData")