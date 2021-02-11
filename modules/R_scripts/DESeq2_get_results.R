

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
DESeq2_dds_execute <- readRDS(
    snakemake@input[['dds_execute']]
)

##############################################################
##### DATA FORMATING FOR DESEQ2 
##############################################################

## DESeq2 function for getting results
deseq2_get_results <- function(
    dds, # dds that is generated in the DESeq2_execute rules
    contrast, # vector that contain : the variable names, condition 1 and condition 2 (ie : c("condition", "C", "I"))
    padj_treshold, # pvalue adjusted threshold used for filtered the genes
    fc_threshold
) {
    ## extract the results
    DE_gene_results <- as.data.frame(
        results(
            dds,
            contrast = contrast,
            independentFiltering = F,
            pAdjustMethod = "bonferroni",
            parallel = TRUE
        )
    )

    ## transform the data frame to data table and reorder the genes with the pvalues
    DE_gene_results <- as.data.table(
        DE_gene_results,
        keep.rownames = T
    )[order(padj, decreasing = FALSE),]

    ## filter by the padj
    DE_gene_results <- DE_gene_results[padj < padj_treshold, ]

    ## filter by the log2FoldChange
    DE_gene_results <- DE_gene_results[(log2FoldChange > fc_threshold) | (log2FoldChange < -(fc_threshold)), ]

    return(DE_gene_results)
}


condition_results <- deseq2_get_results(
    dds = DESeq2_dds_execute,
    contrast = c(
        "condition",
        "C",
        "I"
    ),
    padj_treshold = 0.001,
    fc_threshold = 2
)

##############################################################
##### WRITE THE OUTPUTS
##############################################################
## create the directory
#dir.create(snakemake@output[["DESeq2_results_directory"]])

fwrite(
    condition_results,
    snakemake@output[['results']]   
)
