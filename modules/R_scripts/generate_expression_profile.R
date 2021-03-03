
#################################################
## Generate the expression profile of each groups for the genes
#################################################

#################################################
## Load the libraries
#################################################

library(data.table)
library(tidyverse)
suppressMessages(library(DESeq2))

#################################################
##  Load the data
#################################################

## load the sample information
sample_information <- fread(
    snakemake@input[['sample_information']],
    sep = ","
)

## load the dds
dds_all <- readRDS(
    snakemake@input[['dds_all']]
)

## extract expression count
normalized_count <- as.data.table(counts(dds_all, normalized = T), keep.rownames = 'gene_name')

## load the deseq results for dai comparisons
results_dai_7vs14 <- fread(
    snakemake@input[["results_dai_7vs14"]],
    sep = ','
)

results_dai_7vs21 <- fread(
    snakemake@input[["results_dai_7vs21"]],
    sep = ','
)

results_dai_14vs21 <- fread(
    snakemake@input[["results_dai_14vs21"]],
    sep = ','
)


#################################################
## Functions
#################################################

## extract the genes that are significant differentially expressed
extract_DEgenes <- function(
    deseq_result,
    padj_treshold,
    fc_threshold
) {

    ## copy the input
    DE_gene_names <- copy(deseq_result)

    ## filter by the p values
    DE_gene_names <- deseq_result[padj < padj_treshold, ]

    ## filter by the fold change
    DE_gene_names <- DE_gene_names[(log2FoldChange > fc_threshold) | (log2FoldChange < -(fc_threshold)), ]

    ## extract only gene names
    DE_gene_names <- unlist(DE_gene_names[,gene_name])

    return(DE_gene_names)
}

#################################################
## Extract genes
#################################################

padj_cutoff <- 0.001
fc_cutoff <- 2

## Extract the DEgenes
DEgenes_dai_7vs14 <- extract_DEgenes(
    results_dai_7vs14,
    padj_treshold = padj_cutoff,
    fc_threshold = fc_cutoff
)
DEgenes_dai_7vs21 <- extract_DEgenes(
    results_dai_7vs21,
    padj_treshold = padj_cutoff,
    fc_threshold = fc_cutoff
)
DEgenes_dai_14vs21 <- extract_DEgenes(
    results_dai_14vs21,
    padj_treshold = padj_cutoff,
    fc_threshold = fc_cutoff
)

## extract the genes that are in common in the three comparison
DEgenes_in_common <- Reduce(
    intersect,
    list(
        DEgenes_dai_7vs14,
        DEgenes_dai_7vs21,
        DEgenes_dai_14vs21
    ),
)

## regex spectific to canola or clubroot
canola_id_prefix_regex <- 'ENSRNA|GSBRNA2T'
clubroot_id_prefix_regex <- 'ENSRNAG|PBRA'

## create a table that contain the association of each gene with the host specie
specie_distribution <- data.table(
    gene_name = DEgenes_in_common
)

specie_distribution[
    grep(
        x = gene_name,
        pattern = canola_id_prefix_regex
    ),
    specie_name := 'canola'
]  

specie_distribution[
    grep(
        x = gene_name,
        pattern = clubroot_id_prefix_regex
    ),
    specie_name := 'clubroot'
]

## retreive the distribution information of canola and clubroot genes
gene_distribution <- as.data.table(table(specie_distribution[,specie_name]))[, setnames(.SD, 'V1', 'specie')]
gene_distribution <- gene_distribution[, percentage := N/sum(N)]

#################################################
## function for data formating and plot generation
#################################################

## Function for extracted DEgene of a specific data and formate the data for gene expression profile generation
formate_expression_profil <- function(
    specie_input
) {

    ## extract genes associated with the specie
    DEgenes_in_common <- unlist(specie_distribution[specie_name %in% specie_input, gene_name])

    ## extract the expression data associated with the DEgenes in common
    normalized_count_filtered <- copy(normalized_count[gene_name %in% DEgenes_in_common])

    ## do the transposition
    normalized_count_filtered <- as.data.table(t(normalized_count_filtered), keep.rownames = TRUE)

    ## rename column names
    normalized_count_filtered <- normalized_count_filtered[
        ,
        setnames(
            .SD,
            colnames(normalized_count_filtered),
            unlist(normalized_count_filtered[rn == 'gene_name',])
        )
    ][
        gene_name != 'gene_name',
    ][
        ,
        setnames(
            .SD,
            'gene_name',
            'sample_id'
        )
    ]

    ## formate the expression_data datatable structure
    normalized_count_filtered <- melt.data.table(
        normalized_count_filtered,
        id.vars = "sample_id",
        measure.vars = DEgenes_in_common
    )

    ## merge with the sample information
    normalized_count_filtered <- merge(
        sample_information[, c("sample_id", "dai")],
        normalized_count_filtered,
        by = 'sample_id',
        sort = F
    )

    ## do the mean by genes and by dai
    normalized_count_filtered <- normalized_count_filtered[
        ,
        mean_expression := mean(as.numeric(value)),
        by = c('variable', 'dai')
    ][, value := NULL][order(variable),]

    normalized_count_filtered <- unique(normalized_count_filtered[order(mean_expression),][, c('dai', 'variable', 'mean_expression')])

    #normalized_count_filtered <- normalized_count_filtered[variable != 'ENSRNAG00050137108'][variable != 'ENSRNAG00050137230']

    ## return the normalized count filtered
    return(normalized_count_filtered)
}

## function for expression profile generation
generate_expressionprofile <- function(
    formated_data, # take as input the date generated by the function formate_expression_profil
    file_path
) {
    ## generate the expression profil figures
    p <- ggplot(
            data = formated_data,
            aes(
                x = as.factor(dai),
                y = as.numeric(mean_expression),
                color = as.factor(variable)
            )
        ) + 
            theme_minimal(base_size = 28) + 
            geom_point(alpha = 0) +
            geom_line(aes(group = variable), alpha = 0.7, lwd = 3) +
            labs(
                x = "Day After Innoculation",
                y = "Normalized Expression"
            ) +
            theme(legend.position="none")

    ggsave(
        plot = p,
        filename = file_path,
        device = "png",
        height = 8,
        width = 14,
        limitsize = F
    )
}




#################################################
## Generate figures
#################################################

## create the directory that will contain all the figures of this rules
dir.create(
    snakemake@output[["expression_profile"]]
)

#####
## Generate figures that describes the gene distribution amongst species
#####

## create the table that contain the informations
fwrite(
    x = gene_distribution,
    file = paste(
        snakemake@output[["expression_profile"]],
        '/',
        'gene_distribution.csv',
        sep = ''
    ),
    sep = ','
)

## pie chart of the gene distribution
distribution_bar_chart <- ggplot(
    data = gene_distribution,
    aes(x = '', y = percentage, fill = specie)
) +
    geom_bar(stat="identity") +
    coord_polar("y", start=0) +
    theme_minimal(base_size = 28) + 
    theme(
        axis.text.x=element_blank(),
        legend.key.size = unit(2, 'cm')
    ) +
    labs(x = '', y = '')

ggsave(
    plot = distribution_bar_chart,
    paste(
        snakemake@output[["expression_profile"]],
        '/',
        'gene_distribution_piechart.png',
        sep = ''
    ),
    device = "png",
    height = 8,
    width = 8
)

#####
## generate the expression profil figures
#####

## extract and formate data
data_canola <- formate_expression_profil(specie_input = 'canola')
data_clubroot_full <- formate_expression_profil(specie_input = 'clubroot')
data_clubroot_filtered <- formate_expression_profil(specie_input = 'clubroot')[variable != 'ENSRNAG00050137108'][variable != 'ENSRNAG00050137230']

## generate expression profile plots
generate_expressionprofile(
    formated_data = data_canola,
    file_path = paste(
       snakemake@output[["expression_profile"]],
        '/',
        'canola_gene_expressionprofile.png',
        sep = ''
    )
)
generate_expressionprofile(
    formated_data = data_clubroot_full,
    file_path = paste(
       snakemake@output[["expression_profile"]],
        '/',
        'clubroot_gene_expressionprofile_full.png',
        sep = ''
    )
)
generate_expressionprofile(
    formated_data = data_clubroot_filtered,
    file_path = paste(
       snakemake@output[["expression_profile"]],
        '/',
        'clubroot_gene_expressionprofile_filtered.png',
        sep = ''
    )
)
