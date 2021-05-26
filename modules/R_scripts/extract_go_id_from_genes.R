
##################################################
## Extract all the go id for each gene of interest
##################################################

##################################################
## Load the libraries
##################################################

library(data.table)
library(tidyverse)

##################################################
## Load the data
##################################################

### load the data generated with the heatmap
#heatmap_directory <- (
#    snakemake@input[["heatmap_dir"]]
#)

### retreive the path of each files into the heatmap directory
#file_path_vector <- list.files(
#    heatmap_directory,
#    full.names = T,
#)

### read the file that contain the genes for each cluster that was generated
#cluster_genes_dt <- fread(
#    grep(file_path_vector, pattern = 'gene_cluster', value = T),
#    sep = ','
#)

## load the gene ontology annotation related to p brassicae
goa_p_brassicae <- fread(
    snakemake@input[["goa_plasmodiophora_brassicae"]],
    sep = "\t"
)

##################################################
## ANALYSIS 
##################################################

## extract only the genes related to p brassicae
#p_brassicae_genes_vector <- unlist(cluster_genes_dt[specie == "clubroot", 'gene_name'])

p_brassicae_genes_vector <- unique(unlist(goa_p_brassicae[, V3]))



## merge the go id for a same gene
get_go_id <- function(
    gene_name, # gene name of interest
    goa_datatable # data table that contain the go annotations
) {
    print(gene_name)
    ## extract the go annotation related to the gene name
    goa_of_interest <- goa_datatable[V3 == gene_name, ]

    ### if there is not information related to the gene name of interest
    #if (nrow(goa_of_interest) == 0) {
    #    ## create a empty data table
    #    go_id_dt <- data.table(
    #        gene = character(),
    #        go_id = character()
    #    )

    #    return(go_id_dt)
    #}

    ## extract the go id related to the gene name
    go_id_of_interest_vector <- unlist(goa_of_interest[,V5]) 

    ## formate the go id
    formated_go_id <- paste(
        go_id_of_interest_vector,
        collapse = ", "
    )
    
    ## put the formated go id into a data table
    go_id_dt <- data.table(
        gene = gene_name,
        go_id = formated_go_id
    )

    return(go_id_dt)

}

## extract all the goa term into a list
formated_goa_of_interest <- lapply(
    p_brassicae_genes_vector,
    function(x) get_go_id(
        gene_name = x,
        goa_datatable = goa_p_brassicae
    )
)

## merge all the information into a data table
formated_goa_of_interest <- Reduce(
    x = formated_goa_of_interest,
    function(x,y) rbind(as.data.table(x),as.data.table(y))
)

##################################################
## WRITE OUTPUT
##################################################

## write the file that contain the go id
fwrite(
    formated_goa_of_interest,
    snakemake@output[['go_id']],
    sep = "\t",
    col.names = F
)
