

##############################################################
##### Differential gene expression analysis using DESeq2 #####
##############################################################

## generate the dds object from the raw count matrix

##############################################################
##### LOAD LIBRARIES
##############################################################


library(data.table)
# suppressMessages(library(DESeq2))
library(tidyverse)

##############################################################
##### LOAD DATA
##############################################################


## load the raw combined counts data
dt_combined_counts <- fread(
    file = snakemake@input[['combined']]
)


##############################################################
##### DATA FORMATING FOR DESEQ2 
##############################################################


## extract the colnames that correspond to the sample id (removing the 'gene' names)
sample_id <- colnames(dt_combined_counts)[colnames(dt_combined_counts) != 'gene']

## extract the sample information
sample_information <- lapply(
    X = sample_id,
    function(x) str_split(
        x,
        pattern = '[.]',
    )[[1]]
)

## put the sample information into a data table
sample_information <- cbind(
    data.table('sample_id' = sample_id),
    as.data.table(
        transpose(sample_information)
    )
)

## rename column
sample_information <- sample_information [
    ,
    setnames(
        .SD,
        c('V1', 'V2', 'V3'),
        c('treat', 'condtion', 'dai')
    )
]

## transform to data frame to have rownames
df_sample_information <- data.frame(
    sample_information,
    row.names = sample_information[, sample_id]
)


##############################################################
##### EXECUTE DESEQ2
##############################################################


dds <- DESeqDataSetFromMatrix(
    countData = dt_combined_counts,
    colData = df_sample_information,
    design = ~ dai + condition
)


##############################################################
##### WRITE THE OUTPUTS
##############################################################

## create the directory
dir.create(snakemake@output[["DESeq2_dds_directory"]])

## save the dds object
saveRDS(
    dds,
    file = paste(
        snakemake@output[["DESeq2_dds_directory"]],
        '/',
        'FILENAMES',
        '.dds',
        sep = ''
    )
)