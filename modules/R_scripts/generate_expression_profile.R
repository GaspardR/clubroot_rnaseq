
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
    DE_gene_names <- deseq_result[(log2FoldChange > fc_threshold) | (log2FoldChange < -(fc_threshold)), ]

    ## extract only gene names
    DE_gene_names <- unlist(DE_gene_names[,gene_name])

    return(DE_gene_names)
}

#################################################
## Analysis
#################################################


## Extract the DEgenes
DEgenes_dai_7vs14 <- extract_DEgenes(
    results_dai_7vs14,
    padj_treshold = 0.05,
    fc_threshold = 0
)
DEgenes_dai_7vs21 <- extract_DEgenes(
    results_dai_7vs21,
    padj_treshold = 0.05,
    fc_threshold = 0
)
DEgenes_dai_14vs21 <- extract_DEgenes(
    results_dai_14vs21,
    padj_treshold = 0.001,
    fc_threshold = 2
)

## extract the genes that are in common in the three comparison
DEgenes_in_common <- Reduce(
    intersect,
    list(
        DEgenes_dai_7vs14,
        DEgenes_dai_7vs21,
        DEgenes_dai_14vs21
    ),
)[1:5]

## extract the expression data associated with the DEgenes in common
normalized_count_filtered <- normalized_count[gene_name %in% DEgenes_in_common]

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












print(normalized_count_filtered[,])

normalized_count_filtered[
    ,
    mean_expression := mean(as.numeric(value)),
    by = c('dai', 'variable')
]

print(normalized_count_filtered)
break



p <- ggplot(
        data = normalized_count_filtered,
        aes(
            x = as.factor(dai),
            y = as.numeric(value)
            #color = as.factor(dai)
        )
    ) + 
        theme_minimal(base_size = 28) + 
        geom_point(alpha = 0) +
        geom_path(aes(group = dai), alpha = 1, lwd = 1) +

        #geom_point(data = expression_data_input[cluster == "mean_profile",],  aes(x = variable, y = value, color = cluster), color = "black", size = 5, alpha = 0) +
        #geom_path(data = expression_data_input[cluster == "mean_profile",],  aes(x = variable, y = value, color = cluster, group = sample_id), lwd = 3, color = "black") +
        labs(
            x = "Day After Innoculation",
            y = "Normalized Expression"
        )
        #scale_y_continuous(limits=c(0, 2.5)) + # set the scale of the y axe
        #theme(axis.text.x=element_blank())
        #scale_color_manual(
            #values = color_table_input[cluster %in% cluster_name_vector_input,color] # set the color
        #)

print(p)

ggsave(
    plot = p,
    filename = 'data/test.png',
    device = "png",
    height = 15,
    width = 20,
    limitsize = F
)


break