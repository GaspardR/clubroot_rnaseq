

##############################################################
##### Differential gene expression analysis using DESeq2 #####
##############################################################

## generate the results from the dds objects

##############################################################
##### LOAD LIBRARIES
##############################################################


library(data.table)
suppressMessages(library(DESeq2))
library("BiocParallel")
register(MulticoreParam(snakemake@threads))

##############################################################
##### LOAD DATA
##############################################################

## load the DDS data
DESeq2_dds <- readRDS(
    snakemake@input[['DESeq2_dds_init']]
)

##############################################################
##### DATA FORMATING FOR DESEQ2 
##############################################################

## DESeq2 function for gene expression analysis
deseq2_analysis <- function(
    dds
) {

    ## filtered genes that do not have a lot of counting
    keep <- rowSums(counts(dds)) >= 10
    dds2 <- dds[keep, ]

    ## differential expression analysis
    dds2 <- DESeq(
        dds2,
        parallel = TRUE
    )

    ## return dds2
    return(dds2)
}


##############################################################
##### WRITE THE OUTPUTS
##############################################################

dds2 <- deseq2_analysis(
    DESeq2_dds
)

## save the dds2 object
saveRDS(
    dds2,
    file = snakemake@output[["DESeq2_dds_execute"]]
)

