
##################################################
## Enrichment analysis
##################################################

##################################################
## Load the libraries
##################################################

library(data.table)
library(tidyverse)
library(topGO)

##################################################
## Load the data
##################################################

## load the gene ontology id related to the genes of interest
go_id <- readMappings(snakemake@input[["go_id"]])

## load the data generated with the heatmap
heatmap_directory <- (
    snakemake@input[["heatmap_dir"]]
)

## retreive the path of each files into the heatmap directory
file_path_vector <- list.files(
    heatmap_directory,
    full.names = T,
)

## read the file that contain the genes for each cluster that was generated
cluster_genes_dt <- fread(
    grep(file_path_vector, pattern = 'gene_cluster', value = T),
    sep = ','
)

## extract the genes of interest
cluster_genes_of_interest <- cluster_genes_dt[specie == "clubroot", gene_name]


##################################################
## ANALYSIS 
##################################################

enrichment_analysis <- function(
    formatted_go_id, # output generated with the readMappings function
    gene_of_interest_vector, # vector that contain the gene of interest to use for the enrichment analysis
    go_term_type # chose between BP, MF or CC
) {

    ## extract the gene name from the go id data
    gene_names_vector <- names(formatted_go_id)

    ## mark the genes of interest
    geneList <- factor(as.integer(gene_names_vector %in% gene_of_interest_vector))
    names(geneList) <- gene_names_vector

    ## create a new topGodata
    GOdata <- new(
        "topGOdata",
        ontology = go_term_type,
        allGenes = geneList,
        annot = annFUN.gene2GO,
        gene2GO = formatted_go_id
    )

    ## set the statistic for the enrichment
    fisher_test <- new(
        "classicCount",
        testStatistic = GOFisherTest,
        name = "Fisher test"
    )
    KS_test <- new(
        "elimScore",
        testStatistic = GOKSTest,
        name = "KS tests"
    )

    ## perform the enrichment analysis
    results_fisher <- getSigGroups(GOdata, fisher_test)
    results_KS <- getSigGroups(GOdata, KS_test)

    ## extract the pval
    pval <- sort(score(results_KS))

    ## correct it with a fdr
    adj_pval <- p.adjust(
        p = pval,
        method = 'fdr'
    )

    ## generate a table that contain the results
    allRes <- as.data.table(
        GenTable(
            GOdata,
            classic = results_fisher,
            KS = results_KS,
            #weight = resultWeight,
            orderBy = "classic"
        )
    )

    ## put the go term branch into the data table
    allRes[,GO := go_term_type]

    ## return the results
    return(allRes)
}


##################################################
## Call the function for enrichment analysis
##################################################

## set the three ontology branch
go_type <- c(
    "BP",
    "MF",
    "CC"
)

## do the go enrichment for the three branches
results_list <- lapply(
    go_type,
    function(x) enrichment_analysis(
        formatted_go_id = go_id, 
        gene_of_interest_vector = cluster_genes_of_interest,
        go_term_type = x
    )
)

## put into a same data table all the enrichment results
results_dt <- Reduce(
    x = results_list,
    function(x, y) rbind(x, y)
)

##################################################
## Write the data
##################################################

## write the data table that contain all the results
fwrite(
    results_dt,
    snakemake@output[["enrichment_results"]],
    sep = ","
)
