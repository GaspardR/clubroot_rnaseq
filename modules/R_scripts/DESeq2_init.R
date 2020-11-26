

##############################################################
##### Differential gene expression analysis using DESeq2 #####
##############################################################

## generate the dds object from the raw count matrix

##############################################################
##### LOAD LIBRARIES
##############################################################


library(data.table)
suppressMessages(library(DESeq2))
break


##############################################################
##### LOAD DATA
##############################################################


## load the raw count matrix
dt_raw_count <- fread(
    file = snakemake@input[['rawcount_matrix']]
)


##############################################################
##### DATA FORMATING FOR DESEQ2 
##############################################################


# Il nous les donnees d'expression (dt_raw_count), les donnees d'information sur les echantillons (dt_sample_data)




##############################################################
##### EXECUTE DESEQ2
##############################################################


dds <- DESeqDataSetFromMatrix(
    countData = dt_raw_count,
    colData = dt_sample_data,
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