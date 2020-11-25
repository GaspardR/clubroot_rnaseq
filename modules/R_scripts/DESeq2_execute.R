

##############################################################
##### Differential gene expression analysis using DESeq2 #####
##############################################################

## generate the results from the dds objects

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




##############################################################
##### WRITE THE OUTPUTS
##############################################################

## create the directory
dir.create(snakemake@output[["DESeq2_results_directory"]])

## save the dds object