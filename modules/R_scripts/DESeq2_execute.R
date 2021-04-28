

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
    snakemake@input[['dds_init']]
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


## extract expression count and rlog normalized counts
if (snakemake@params[['condition']] == 'all') {



    ## expression count
    normalized_count <- as.data.table(counts(dds2, normalized = T), keep.rownames = 'gene_name')

    fwrite(
        normalized_count,
        #'data/DESeq2/DESeq2_sizefactor_normalized_count.csv',
        snakemake@output[['DESeq2_size_factor_normalization']],
        sep = ','
    )

    ## rlog expression count
    #vst_normalized_count <- as.data.table(vst(DESeq2_dds), keep.rownames = 'gene_name')
    vst_normalized_count <- as.data.table((assay(vst(DESeq2_dds))), keep.rownames = 'gene_name')

    fwrite(
        vst_normalized_count,
        #'data/DESeq2/DESeq2_vst_normalized_count.csv',
        snakemake@output[['DESeq2_vst_normalization']],
        sep = ','
    )
} 

if (snakemake@params[['condition']] != 'all') {

    empty_output <- data.table(
        empty = character(),
        output = character()
    )

    fwrite(
        empty_output,
        snakemake@output[['DESeq2_size_factor_normalization']],
        sep = ','
    )

    fwrite(
        empty_output,
        snakemake@output[['DESeq2_vst_normalization']],
        sep = ','
    )

}


## create the directory that contain the dds_excute output
#dir.create(snakemake@output[["DESeq2_dds_execute"]])

## save the dds2 object
saveRDS(
    dds2,
    file = snakemake@output[['dds_execute']],

)

