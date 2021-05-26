
######################################
## Create the bed file of the genes of interest
######################################

######################################
## Load the libraries
######################################

library(data.table)
library(tidyverse)
options(width=system("tput cols", intern=TRUE))


######################################
## load the data
######################################

## load the gtf of clubroot and take only gene information
p_brassicae_gtf <- fread(
    snakemake@input[["gtf_plasmodiophora_brassicae"]],
    sep = '\t'
)[V3 == 'gene']

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

## extract the names of each cluster into a vector
cluster_name_vector <- unique(unlist(cluster_genes_dt[,cluster_name]))

## extract the species 
specie_name_vector <- unique(unlist(cluster_genes_dt[, specie]))

## extract genes
clubroot_gene_vector <- unlist(unlist(cluster_genes_dt[specie == 'clubroot',gene_name]))

######################################
## Load the data
######################################

## initialize the list that will contain gene information for every cluster
gene_information_list <- list()

## initialize the data table that will contain gene information for the generation of the bed file of a cluster
gene_information <- data.table(
    'chromosome' = character(),
    'from' = numeric(),
    'to' = numeric(),
    'gene_name' = character()
)

## for each gene
for (i_gene in seq(1, length(clubroot_gene_vector), 1)) {
    
    ## extract gene name
    gene_name <- clubroot_gene_vector[i_gene]

    ## extract the infromation associated with the gene
    clubroot_of_interest_information <- p_brassicae_gtf[grep(p_brassicae_gtf[, V9], pattern = gene_name),]

    ## put all the relevant information into a data table
    clubroot_of_interest_information_dt <- data.table(
        'chromosome' = unlist(clubroot_of_interest_information[,1]),
        'from' = unlist(clubroot_of_interest_information[,V4]),
        'to' = unlist(clubroot_of_interest_information[,V5]),
        'gene_name' = gene_name
    )

    ## put the information into the gene_information data table
    gene_information <- rbind(
        gene_information,
        clubroot_of_interest_information_dt
    )
}

######################################
## WRITE OUTPUT
######################################

fwrite(
    gene_information,
    snakemake@output[['p_brassicae_bed_file']],
    sep = '\t',
    col.names = F
)
