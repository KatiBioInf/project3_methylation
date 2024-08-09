# project3_methylation
Differential methylation analysis and visualization

We perform differential methylation analysis on two samples ([colon tumor](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1204465) (sample id: GSM1204465) and [adjacent normal sample](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1204466) (sample id: GSM1204466)) from [this bioproject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA201480) (accession: PRJNA201480; GEO: GSE46644). Although the scripts demonstrate how to process the original data (assuming the data are downloaded in the same map as the R project) we also provide the processed data for analysis and visualization. Thus, the analysis can be performed by directly running script 03 and 04.

We used the Bioconductor package [MethylKit](https://www.bioconductor.org/packages/release/bioc/html/methylKit.html) for analysis and the [EnrichedHeatmap](https://bioconductor.org/packages/release/bioc/html/EnrichedHeatmap.html) package for visualization. We rely heavily on the [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) package as well. Annotation data (cpg islands) were obtained from the ??? package.

The repository consists of the following scripts:

01_obtaining_methylation_data.R

02_obtaining_annotation_data.R

03_Diff_Meth_destranded.R

04_visualization_heatmap.R


Resources:
MethylKit tutorial: https://www.bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html



