
#################################################
## Generate the heatmap of the genes that are differentially expressed
#################################################

#################################################
## Load libraries
#################################################

library(data.table)
library(tidyverse)
library(ComplexHeatmap)
library(colorDF)
library(RColorBrewer)


options(width=system("tput cols", intern=TRUE))
print.data.table <- colorDF:::print_colorDF

#################################################
## Load the data
#################################################

## Load the normalized expression data
normalized_count <- fread(
    snakemake@input[['vst_normalized_count']]
)

test <- (grep(colnames(normalized_count), pattern = 'I', value = TRUE))
normalized_count <- normalized_count[, c('gene_name', test), with = F]

#print(normalized_count)
#break


## extract file names of thhe expression profile directory
files <- list.files(snakemake@input[['generate_expression_dir']], full.names = T)

## load the de genes
DE_genes_canola <- scan(
    grep(files, pattern = 'DEgenes_canola', value = T),
    sep = '\n',
    what = character()
)
## load the de genes
DE_genes_clubroot <- scan(
    grep(files, pattern = 'DEgenes_clubroot', value = T),
    sep = '\n',
    what = character()
) 

## extract expression data associated with the canola or the clubroot
normalized_count_canola <- normalized_count[gene_name %in% DE_genes_canola,] 
normalized_count_clubroot <- normalized_count[gene_name %in% DE_genes_clubroot,] 

## merge the canola and clubroot gene expression data
normalized_count_all <- rbind(
    normalized_count_canola,
    normalized_count_clubroot
)


#################################################
## Function for the generation of the heatmap
#################################################

## function for creating the annotation for the heatmap
create_annotation <- function(
    variable,
    clinicalData,
    column_annotation = TRUE,
    color_vector
) {
    ## create the colors
    #color_vector <- c(
    #    brewer.pal(9, "Set1"),
    #    brewer.pal(8, "Dark2"),
    #    brewer.pal(8, "Set2"),
    #    brewer.pal(12, "Paired")
    #)

    ## extract data table with the specific variables
    clinicalData <- clinicalData[, variable, with = F]

    ## initialise the color list
    colors_for_variables <- list()

    ## for each data associated with each variables, do a counting of occurence and extract colors
    for (i in (1:length(variable)))
    {

        ## extract the number of catagories in the variable i-th
        categories <- names(table(clinicalData[, variable[i], with = F]))

        ## take the i-th first colors from the "color_vector" variable
        colors_for_categories <- color_vector[1:length(categories)]

        ## remove the colors used from the initial color vector
        color_vector <- color_vector[!(color_vector %in% colors_for_categories)]
 

        ## name each color by the category
        names(colors_for_categories) <- categories

        ## append into the list the colors categories
        colors_for_variables[[i]] <- colors_for_categories
    }

    ## give variable names for each vector in the list
    names(colors_for_variables) <- variable


    if (column_annotation == TRUE)
    {
        ## create the annotations
        annotation <- HeatmapAnnotation(
            df = clinicalData,
            col = colors_for_variables
        )
        return(annotation)

    }

    if (column_annotation == FALSE)
    {
        ## create the annotations
        annotation <- rowAnnotation(
            df = clinicalData,
            col = colors_for_variables
        )
        return(annotation)

    }
}


## function for creating the heatmap
create_heatmap <- function(
    expression_data,
    figure_path,
    show_row_names = T,
    column_split,
    row_split,
    row_annotation = T,
    width = 15,
    height = 30
) {

    ## extract the gene names
    row_names <- unlist(expression_data[,gene_name])

    print(length(row_names))

    ## transform the data table into data freame
    expression_data <- data.frame(expression_data[, !c('gene_name')])

    print(dim(expression_data))


    ## rename the names of the row
    rownames(expression_data) <- row_names
    
    ## transform all the values to numeric
    expression_data <- apply(
        expression_data,
        c(1,2),
        function(x) as.numeric(x)
    )

    colors_vector <- c(
        brewer.pal(9, "Set1"),
        brewer.pal(8, "Dark2"),
        brewer.pal(8, "Set2"),
        brewer.pal(12, "Paired")
    )

    ##########
    ## For the column annotation of the heatmap
    ##########

    ## retreive the information for each samples
    sampple_information <- as.data.table(sapply(
        colnames(expression_data),
        function(x) str_split(
            x,
            pattern = '[.]'
        )[[1]]
    ))

    ## add all the sample information in a data table
    dt_sample_information <- data.table(
        sample_id = colnames(sampple_information),
        Strain = unlist(sampple_information[1,]),
        DAI = unlist(sampple_information[3,])

    )

    ## create the annotation for the heatmap
    column_annotations <- create_annotation(
        variable = c(
            'Strain',
            'DAI'
        ),
        clinicalData = dt_sample_information,
        column_annotation = TRUE,
        color_vector = colors_vector[1:6]
    )

    ##########
    ## For the column annotation of the heatmap
    ##########

    ## extract the gene name
    dt_gene_information <- data.table(
        gene_name = rownames(expression_data)
    )

    canola_id_prefix_regex <- 'ENSRNA|GSBRNA2T'
    clubroot_id_prefix_regex <- 'ENSRNAG|PBRA'

    ## add the gene type for each gene
    dt_gene_information[
        grep(
            x = gene_name,
            pattern = canola_id_prefix_regex
        ),
        Specie := 'canola'
    ]  

    dt_gene_information[
        grep(
            x = gene_name,
            pattern = clubroot_id_prefix_regex
        ),
        Specie := 'clubroot'
    ]
    

    ## create the row annotation
    row_annotations <- create_annotation(
        variable = c(
            'Specie'
        ),
        clinicalData = dt_gene_information,
        column_annotation = FALSE,
        color_vector = colors_vector[6:7]
    )

    ##########
    ## generate the heatmap
    ##########

    if (row_annotation == TRUE) {
        heatmap <- Heatmap(
            expression_data,
            border = TRUE,

            ## for the columns
            column_title = 'Samples',
            column_title_gp = gpar(fontsize = 22, fontface = "bold"),
            column_split = column_split,
            column_gap = unit(6, "mm"),
            clustering_distance_columns = "minkowski",
            top_annotation = column_annotations,
            left_annotation = row_annotations,

            ## for thr rows
            row_title = 'Genes',
            row_title_gp = gpar(fontsize = 22, fontface = "bold"),
            show_row_names = FALSE,
            row_split = row_split,
            row_gap = unit(3, "mm"),
            clustering_distance_rows = "minkowski",
            clustering_method_rows = "ward.D",
            show_row_dend = T,

            heatmap_legend_param = list(
                title = "EXPRESSION LEVELS",
                title_position = "leftcenter-rot",
                grid_width = unit(0.75, "cm"),
                legend_height = unit(7, "cm"),
                border = "black"
            )
        )
    }

    if (row_annotation == FALSE) {
        heatmap <- Heatmap(
            expression_data,
            border = TRUE,

            ## for the columns
            column_title = 'Samples',
            column_title_gp = gpar(fontsize = 22, fontface = "bold"),
            column_split = column_split,
            column_gap = unit(6, "mm"),
            clustering_distance_columns = "minkowski",
            top_annotation = column_annotations,

            ## for thr rows
            row_title = 'Genes',
            row_title_gp = gpar(fontsize = 22, fontface = "bold"),
            show_row_names = FALSE,
            row_split = row_split,
            row_gap = unit(3, "mm"),
            clustering_distance_rows = "minkowski",
            clustering_method_rows = "ward.D",
            show_row_dend = T,

            heatmap_legend_param = list(
                title = "EXPRESSION LEVELS",
                title_position = "leftcenter-rot",
                grid_width = unit(0.75, "cm"),
                legend_height = unit(7, "cm"),
                border = "black"
            )
        )
    }




    ## create the heatmap and save into a file
    svg(
        file = figure_path,
        width = width,
        height = height,

    )
    print(heatmap)
    dev.off()    

}

#################################################
## Create the heatmao using the create_heatmap function
#################################################

## create the directory that will contain all the generated heatmaps
dir.create(snakemake@output[['heatmap_dir']])

## heatmap using all the genes (P brassicae and B napus)
create_heatmap(
    normalized_count_all,
    width = 6,
    height = 8,
    column_split = 3,
    row_split = 2,
    figure_path = paste(
        snakemake@output[['heatmap_dir']],
        '/heatmap_all.svg',
        sep = ''
    )

)

## heatmap using only canola genes
create_heatmap(
    normalized_count_canola,
    width = 6,
    height = 8,
    column_split = 3,
    row_split = 3,
    row_annotation = F,
    figure_path = paste(
        snakemake@output[['heatmap_dir']],
        '/heatmap_canola.svg',
        sep = ''
    )
)

## heatnap using only P brassicae genes
create_heatmap(
    normalized_count_clubroot,
    width = 6,
    height = 8,
    column_split = 3,
    row_split = 2,
    row_annotation = F,
    figure_path = paste(
        snakemake@output[['heatmap_dir']],
        '/heatmap_clubroot.svg',
        sep = ''
    )

)


#################################################
## END
#################################################
