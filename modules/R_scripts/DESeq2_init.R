

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
        c('treat', 'condition', 'dai')
    )
]

## extract the information for each condition
treats <- unique(unlist(sample_information[, 'treat']))
conditions <- unique(unlist(sample_information[, 'condition']))
dai <- unique(unlist(sample_information[, 'dai']))

extract_colData_and_expressionData <- function(
    sample_information_datatable_input, # data table containing the sample_id, treat, condition and dai column names
    comparaison_name
) {

    ## extract the sample of interest contained in the sample information input
    sample_of_interest <- unlist(sample_information_datatable_input[, sample_id])

    ## extract expression data of the sample of interest
    dt_combined_counts_subset <- copy(
        dt_combined_counts[, c('gene', sample_of_interest), with = F]
    )
    
    ## transform the sample information data table input into a data frame
    sample_information_dataframe <- data.frame(
        sample_information_datatable_input,
        row.names = sample_information_datatable_input[, sample_id]
    )

    ## put the expression data and the count into a list and rename it
    output <- list()
    output[[comparaison_name]]['count'] <- list(dt_combined_counts_subset)
    output[[comparaison_name]]['information'] <- list(sample_information_dataframe)

    ## return the output
    return(output)
}


#############################################
## create DESeq2 inputs 
#############################################

## initialize the list that will contain the deseq2 input
deseq2_data_input_list <- list()

## for each dai
for (i in seq(1, length(dai), 1)) {

    # cat('\n\n\n')

    ## copy the data information
    dt_sample_information_subset <- as.data.table(copy(sample_information))

    ## extract the value of the dai of interest
    dai_values <- dai[i]

    ## extract the data associated with the dai of interest
    dt_sample_information_subset <- dt_sample_information_subset[dai == dai_values,]

    ## create the name of the comparison
    comparaison_name <- paste(
        'dai_',
        dai_values,
        sep = ''
    )

    ## extract the expression data associated with the sample id of interest and parse the data for DESeq2
    data_output <- extract_colData_and_expressionData(
        dt_sample_information_subset,
        comparaison_name
    )

    ## put the data into the deseq2_data_input_list
    deseq2_data_input_list[comparaison_name] <- data_output
}

## for each dai combinations

## do all the two combiination of the dai
dai_combination <- combn(
    dai,
    m = 2
)

## for each condition and each combination of the dai

for (condition_index in seq(1, length(conditions), 1)) {

    ## extract the condition value
    condition_value <- conditions[condition_index]

    for (dai_combination_index in seq(1, ncol(dai_combination), 1)) {

        # cat('\n\n\n')

        ## copy the data information
        dt_sample_information_subset <- as.data.table(copy(sample_information))

        ## extract the combination of dai of interest
        dai_values <- dai_combination[,dai_combination_index]

        ## extract the data associated with the condition the dai combination of interest
        dt_sample_information_subset <- dt_sample_information_subset[dai %in% dai_values,]
        dt_sample_information_subset <- dt_sample_information_subset[condition == condition_value, ]

        ## create the name of the comparison
        comparison_name <- paste(
            condition_value,
            '.',
            'dai_',
            dai_values[1],
            '_vs_',
            'dai_',
            dai_values[2],
            sep = ''
        )

        ## extract the expression data associated with the sample id of interest and parse the data for DESeq2
        data_output <- extract_colData_and_expressionData(
            dt_sample_information_subset,
            comparison_name
        )

        ## put the data into the deseq2_data_input_list
        deseq2_data_input_list[comparison_name] <- data_output
    }
}

print(deseq2_data_input_list)

# break













##############################################################
##### EXECUTE DESEQ2 for each input
##############################################################

for (input_index in seq(1, length(deseq2_data_input_list), 1)) {

    ## extract the comparison name
    comparison_name <- names(deseq2_data_input_list[input_index])

    print(comparison_name)
    next

    ## create the dds object
    dds <- DESeqDataSetFromMatrix(
        countData = deseq2_data_input_list[[comparison_name]]['count'],
        colData = deseq2_data_input_list[[comparison_name]]['information'],
        design = ~ dai + condition + treat
    )

    ## save the dds object
    saveRDS(
        dds,
        file = paste(
            snakemake@output[["DESeq2_dds_directory"]],
            '/',
            'comparison_name',
            '.dds',
            sep = ''
        )
    )

}


break
