## Author: Aakanksha Singh
## Date: 12 January 2024
## Version: 1.0


####
## Aim: Loading the libraries required for the project
library("TENxPBMCData")
library("SingleCellExperiment")
library("scater")
library("scran")
library("scuttle")
library("celldex")
library("SingleR")
library("tibble")
library("dplyr")
library("biomaRt")
library("tidyr")
library("tibble")
library("BiocSingular")
library("ensembldb")
library("ggplot2")
library("SPOTlight")
library("SpatialExperiment")
library("NMF")
library("SpatialDecon")
library("GeomxTools")
library("immunedeconv")
library("scRNAseq")
library(scrapper)
library(data.table)
library(stringr)
library(RColorBrewer)
library("ggpubr")
library(preprocessCore)
library(tidyverse)
library(matrixStats)
library(readxl)
library(ggpubr)
library(RColorBrewer)
library(broom)
library(rstatix)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(ggthemes)
library(forcats)
library(singleCellTK)
library("spacexr")

# Sys.setenv(CUDA_VISIBLE_DEVICES = "")
# Sys.setenv(PYTHONNOUSERSITE = "1")
# 
# 
# library(reticulate)
# use_condaenv("cell2loc_env", required = TRUE)
# use_python("/home/aakanksha/miniforge3/envs/cell2loc_env/bin/python", required = TRUE)
# 
# py_run_string("import cell2location; print('ok')")
# source_python("run_cell2location.py")


####
## Aim: To process single cell data and return a SingleCellExperiment Object 
##     with cell type annotation
# For the function processing_sce
#   Input: (a) A SingleCellExperiment dataset
#          (b) The title of the datset (string)
#   Output: (a) Plot highlight the removed mitochondrial gens
#           (b) Plot of mean log expression vs variance of log expression with a treand line
#           (c) A tSNE plot colored by predicted celltype
#           (d) A table with the number of the predicted cell types
#           (e) The output variable stores the analysed SingleCellExperiment Object with cell type annotations. 


processing_sce <- function(dataset,title) 
{   
  # Step 1: Quality Control
  stats <- high.mito <- list()
  unfiltered<-dataset
  
  current <- dataset
  #is.mito <- grep("MT", rowData(dataset)$Symbol_TENx)
  is.mito <- grep("MT", rowData(dataset)$Symbol)
  stats <- perCellQCMetrics(dataset, subsets=list(Mito=is.mito))
  high.mito <- isOutlier(stats$subsets_Mito_percent, type="higher")
  dataset <- dataset[,!high.mito]
  
  current <- unfiltered 
  colData(current) <- cbind(colData(current), stats)
  current$discard <- high.mito
  
  # Note: The code line below generated a plot highlighting the mito genes as outliers
  qcplots<- plotColData(current, x="sum", y="subsets_Mito_percent",
                        colour_by="discard") + scale_x_log10()
  print("Processing sce : Step 1/7 Done")
  
  # Step 2: Log Normalisation
  dataset <- logNormCounts(dataset)
  print("Processing sce : Step 2/7 Done")
  
  # Step 3: Variance Modelling
  dec <- modelGeneVar(dataset)
  hvgs <-getTopHVGs(dataset, prop=0.1)
  
  curdec <- dec
  curfit <- metadata(curdec)
  curve(curfit$trend(x), col='dodgerblue', lwd=2)
  print("Processing sce : Step 3/7 Done")
  
  # Step 4: Dimensionality Reduction (PCA, tSNE, UMAP)
  
  set.seed(10000)
  dataset <- runPCA(dataset)
  
  set.seed(100000)
  dataset <- runTSNE(dataset)
  
  set.seed(1000000)
  dataset <- runUMAP(dataset)
  print("Processing sce : Step 4/7 Done")
  
  # Step 5: Clustering
  g <- buildSNNGraph(dataset, k=10)
  clust <- igraph::cluster_walktrap(g)$membership
  colLabels(dataset)  <- factor(clust)
  print("Processing sce : Step 5/7 Done")
 
  # Step 6: Cell type annotation
  ref<-celldex::MonacoImmuneData(ensembl=TRUE, cell.ont="all")
  predicted_celltype <- SingleR(dataset, ref=ref, labels=ref$label.main)
  
  colData(dataset)$predicted_celltype <- predicted_celltype$labels
  tsne <- scater::plotTSNE(dataset, colour_by="predicted_celltype") + ggtitle(title)
  print("Processing sce : Step 6/7 Done")
  
  # Step 7: Adding an addition column that annotes all T cells the same
  # Adding an additional column with consolidated cell types, 
  # removing rows with no Hugo Symbol,
  # Setting the Hugo symbol as rownames instead of ensmbl ids
  
  dataset$cell<-dataset$predicted_celltype
  colData(dataset)$cell <- ifelse(colData(dataset)$cell %in% c("CD8+ T cells", "CD4+ T cells"), "T cells", colData(dataset)$cell)
  rownames(dataset) <- rowData(dataset)$Symbol
  dataset<-dataset[!(is.na(row.names(dataset))),]
  print("Processing sce : Step 7/7 Done")
  print(table(dataset$cell))
  
  return(dataset)
  
}

####
## Aim: To generate pseudobulk samples with known numbers of cell types from Single
##      SingleCellExperiment Object

# For the function generating_matrix_for_deconv
#   Input: (a) A processed SingleCellExperiment dataset (output from processing_sce())
#          (b) The title of the datset (string)
#          (c) Two vectors containing name and number of each celltype respectively (atm 4 types of cells) 
#   Output: (a) A matrix summarising the SingleCellExperiment that can be used as a bulk sample for further analysis

generating_a_pseudobulk_sample<-function(sce,
                                         title,
                                         cells_to_be_put_in_pseudobulk,
                                         number_of_cells_in_each_celltype, 
                                         perform_TPM_normalisation=c(TRUE,FALSE))
{ 
  # Create a S4 object with the required number of cells for each celltype
  known_sample=data.frame()
  for (i in 1:length(cells_to_be_put_in_pseudobulk))
  {
    cellcol= sce[,sce$cell==cells_to_be_put_in_pseudobulk[i]]
    if (i==1)
      known_sample= cellcol[,sample(ncol(cellcol),number_of_cells_in_each_celltype[i])]
    else
    known_sample=cbind(known_sample,cellcol[,sample(ncol(cellcol),number_of_cells_in_each_celltype[i])])
  }
  
  # Summing the counts of all the cell rowwise to get a pseudobulk sample
  samp<-as.data.frame(Matrix::rowSums(assays(known_sample)$counts))
  samp<-cbind(rownames(known_sample),samp)
  
  # Aggregating rows with common gene symbol
  samp<-aggregate(samp[,-1],by=samp[1],sum)
  
  # Getting TPM ( where TPM= (count/total count)* 1000000)
  count_sum<-sum(samp$x)
  if(perform_TPM_normalisation==TRUE)
    {
      samp[,2]<-sapply(samp[,2],function(i) (i/count_sum)*1000000)
      colnames(samp)=c("HUGO",title)
    }
  else
    {
      colnames(samp)=c("HUGO",title)
    }
  return(samp)
}

####
## Aim: To dconvolute with SpatialDeconv

# For the function spatial_decon
#   Input: (a) A dataframe with the peseudobulk samples (GeneSymbols as rows and colnames as sample)
#          (b) A Processed SingleCellExperiment dataset if a custom signature is to be used, else the defalut safeTME is used as signature
#   Output: (a) A dataframe with proportions of all the cell types

spatial_decon <- function(spatial_data, single_cell_data, matrix_type = c("safeTME", "custom"))
  
{
  # Data wrangling for Spatial Deconv
  
  ## Coverting col to rownames and converting dataframe to numerical matrix
  rownames(spatial_data)<-spatial_data$HUGO
  spatial_data<-data.matrix(spatial_data[-1])
  
  
  ## Convert the deconv_table to a sce which is then converted to seurat with a spatial component
  sc_deconv_table <- SpatialExperiment(
    assays = list(counts = spatial_data, Spatial=spatial_data),
    rowData = rownames(spatial_data),
    colData = colnames(spatial_data))
  
  logcounts(sc_deconv_table)<-log2(counts(sc_deconv_table))
  sc_deconv_table_seurat<- as.Seurat(sc_deconv_table)
  
  ## Adding the Spatial Assay to the seurat object to be able to run spatial decon
  sc_deconv_table_seurat@assays[["Spatial"]]<-CreateAssayObject(data = spatial_data)
  
  # Running spatialdeconv according to the type of signature matrix
  
  if (matrix_type=="custom")
  {
    # Custom profile matrix 
    single_cell_data$cell<-single_cell_data$predicted_celltype
    colData(single_cell_data)$cell <- ifelse(colData(single_cell_data)$cell %in% c("CD8+ T cells", "CD4+ T cells"), "T cells", colData(single_cell_data)$cell)
    rownames(single_cell_data) <- rowData(single_cell_data)$Symbol
    single_cell_data<-single_cell_data[!(is.na(row.names(single_cell_data))),]
    
    annot<-as.data.frame(cbind(single_cell_data$Barcode,single_cell_data$cell))
    
    x<-as.data.frame(assay(single_cell_data))
    colnames(x)<-c(single_cell_data$Barcode)
    
    # The function does not accept a S4 object
    custom_matrix<-create_profile_matrix(mtx = x, 
                                         cellAnnots = annot, 
                                         cellTypeCol = "V2", 
                                         cellNameCol = "V1", 
                                         matrixName = "custom_mini_colon",
                                         outDir = NULL, 
                                         normalize = FALSE, 
                                         minCellNum = 5, 
                                         minGenes = 10,
                                         discardCellTypes = FALSE)
    
    res = runspatialdecon(sc_deconv_table_seurat,
                          bg = 0.01,
                          X = custom_matrix,
                          align_genes = TRUE)
    
    result_df<-as.data.frame(res$prop_of_nontumor)
    result_df$cell_type<-row.names(result_df)
    row.names(result_df)=NULL
    result_df= result_df %>% dplyr::mutate(cell_type= sub("s$", "", cell_type))
  }
  
  if ((matrix_type=="safeTME") | (missing(single_cell_data)))
  {
    
    res = runspatialdecon(sc_deconv_table_seurat,
                          bg = 0.01,
                          X = safeTME,
                          align_genes = TRUE)
    
    result_df<-as.data.frame(res$prop_of_nontumor)
    result_df$cell_type<-row.names(result_df)
    row.names(result_df)=NULL
    
    
    result_df$cell_type<- mgsub::mgsub((result_df$cell_type), c( "macrophages","mast","B.naive", "B.memory" ,
                                                                 "plasma" ,"T.CD4.naive","T.CD4.memory","T.CD8.naive" ,  "T.CD8.memory" ,"NK", "pDCs",
                                                                 "mDCs","monocytes.C" , "monocytes.NC.I" ,   "neutrophils","Treg", "endothelial.cells", 
                                                                 "fibroblasts"), c("Macrophage","Mast cell","B cell","B cell","Plasma cell","T cell",
                                                                                   "T cell","T cell","T cell","NK cell","Dendritic cell","Dendritic cell","Monocyte",
                                                                                   "Monocyte","Neutrophil","T cell","Endothelial cell","Fibroblast cell"))
    
  }
  
  return(result_df)
  
}

####
## Aim: Processing a datset using spatial deconv which is designed specifically for data from Nanostring GeoMx spatial profiler.
##      In this updated method, I process the data like I would process a normal nanostring dataset
##      How is it different from the spatial_decon()
##      1. Use the matrix (cut out the commands to convert it to a seurat object)
##      2. Use the cell_counts feature of the spatialdecon to get cellsper100

# There are two functions in this file, spatial_decon_updated() when you have the number of cells in the bulk sample and a
# second function called spatial_decon_updated_without_cell_number() for when we do not have the number of cells in the 
# the pseudobulk sample
#   Input: (a) The deconvolution table in the form of a data frame
#   Output: (a) A dataframe with cell deconvolutions from spatial_decon with tme and spatial decon with your signature matrix

spatial_decon_updated <- function(spatial_data, single_cell_data, matrix_type = c("safeTME", "custom","LM22","TIL10","Immunostates"), actual_cell_proportions)
  
 {  
    #spatial_data=deconv_table
    #single_cell_data=signature_reference
   

  ## Coverting col to rownames and converting dataframe to numerical matrix
  spatial_data = spatial_data %>% 
                  column_to_rownames(var="HUGO")
  
  a <- as.matrix(spatial_data)
  aoi_table <- data.frame("Sample"=c(colnames(actual_cell_proportions %>% dplyr::select(-cell_type, -ends_with(".actual_frac")))),
                          "nuclei"=c(colSums(actual_cell_proportions %>% dplyr::select(-cell_type,-ends_with(".actual_frac"))))) %>% 
               mutate(Sample_name= str_replace(Sample, "\\.actual_values", ""))
  
  # Running spatialdeconv according to the type of signature matrix
  
  if ((matrix_type=="custom") & (!missing(single_cell_data)))
  {
    # Custom profile matrix 
    single_cell_data$cell<-single_cell_data$predicted_celltype
    colData(single_cell_data)$cell <- ifelse(colData(single_cell_data)$cell %in% c("CD8+ T cells", "CD4+ T cells"), "T cells", colData(single_cell_data)$cell)
    rownames(single_cell_data) <- rowData(single_cell_data)$Symbol
    single_cell_data<-single_cell_data[!(is.na(row.names(single_cell_data))),]
    
    annot<-as.data.frame(cbind(single_cell_data$Barcode,single_cell_data$cell))
    
    #x<-as.data.frame(assay(single_cell_data))
    x <- as.matrix(assay(single_cell_data))
    colnames(x)<-c(single_cell_data$Barcode)
    
    # The function does not accept a S4 object
    custom_matrix<-create_profile_matrix(mtx = x, 
                                         cellAnnots = annot, 
                                         cellTypeCol = "V2", 
                                         cellNameCol = "V1", 
                                         matrixName = "custom_matrix",
                                         outDir = NULL, 
                                         normalize = FALSE, 
                                         minCellNum = 5, 
                                         minGenes = 10,
                                         discardCellTypes = FALSE)

    resu = spatialdecon(a,
                        bg = 0.1,
                        X = custom_matrix,
                        cell_counts = aoi_table$nuclei
    )
    
    cellsper100<-as.data.frame(resu$cell.counts$cells.per.100)
    sum_of_cols <- colSums(cellsper100 [,])
    
    cellsper100<-rbind(cellsper100, c(100-sum_of_cols))
    rownames(cellsper100)[nrow(cellsper100)] <-"uncharacterized cell"
    cellsper100$cell_type<-rownames(cellsper100)
    cellsper100 = cellsper100 %>% dplyr::mutate(cell_type= sub("s$", "", cell_type)) # Change the naming if the cell annotation differs in the single cell data
    
  }
  
  else if (matrix_type=="safeTME")
  { 
    resu = spatialdecon(a,
                        bg = 0.1,
                        X = safeTME,
                        cell_counts = aoi_table$nuclei
    )
    
    cellsper100<-as.data.frame(resu$cell.counts$cells.per.100)
    sum_of_cols <- colSums(cellsper100 [,])
    cellsper100<-rbind(cellsper100, c(100-sum_of_cols))
    rownames(cellsper100)[nrow(cellsper100)] <-"uncharacterized cell"
    cellsper100$cell_type<-rownames(cellsper100)
    cellsper100$cell_type<- mgsub::mgsub((cellsper100$cell_type), c( "macrophages","mast","B.naive", "B.memory" ,
                                                                     "plasma" ,"T.CD4.naive","T.CD4.memory","T.CD8.naive" ,  "T.CD8.memory" ,"NK", "pDCs",
                                                                     "mDCs","monocytes.C" , "monocytes.NC.I" ,   "neutrophils","Treg", "endothelial.cells", 
                                                                     "fibroblasts"), c("Macrophage","Mast cell","B cell","B cell","Plasma cell","T cell",
                                                                                       "T cell","T cell","T cell","NK cell","Dendritic cell","Dendritic cell","Monocyte",
                                                                                       "Monocyte","Neutrophil","T cell","Endothelial cell","Fibroblast cell"))
    
  }
  else if (matrix_type=="TIL10")
  { 
    resu = spatialdecon(a,
                        bg = 0.1,
                        X = til10,
                        cell_counts = aoi_table$nuclei
    )
    
    cellsper100<-as.data.frame(resu$cell.counts$cells.per.100)
    sum_of_cols <- colSums(cellsper100 [,])
    cellsper100<-rbind(cellsper100, c(100-sum_of_cols))
    rownames(cellsper100)[nrow(cellsper100)] <-"uncharacterized cell"
    cellsper100$cell_type<-rownames(cellsper100)
    cellsper100$cell_type<- mgsub::mgsub((cellsper100$cell_type), c("B.cells" ,"Macrophages.M1" , "Macrophages.M2" , "Monocytes", "Neutrophils" ,"NK.cells" ,"T.cells.CD4","T.cells.CD8","Tregs","Dendritic.cells"), 
                                         c("B cell","Macrophage","Macrophage","Monocyte","Neutrophil","NK cell","T cell","T cell","T cell","Dendritic cell"))
    
  }
  
  else if (matrix_type=="LM22")
  { 
    resu = spatialdecon(a,
                        bg = 0.1,
                        X = lm22,
                        cell_counts = aoi_table$nuclei
    )
    
    cellsper100<-as.data.frame(resu$cell.counts$cells.per.100)
    sum_of_cols <- colSums(cellsper100 [,])
    cellsper100<-rbind(cellsper100, c(100-sum_of_cols))
    rownames(cellsper100)[nrow(cellsper100)] <-"uncharacterized cell"
    cellsper100$cell_type<-rownames(cellsper100)
    cellsper100$cell_type<- mgsub::mgsub((cellsper100$cell_type), c("B cells naive","B cells memory" ,"Plasma cells","T cells CD8","T cells CD4 naive","T cells CD4 memory resting" ,"T cells CD4 memory activated",
                                                                    "T cells follicular helper","T cells regulatory (Tregs)","T cells gamma delta","NK cells resting","NK cells activated","Monocytes","Macrophages M0",              
                                                                    "Macrophages M1","Macrophages M2","Dendritic cells resting","Dendritic cells activated","Mast cells resting","Mast cells activated","Eosinophils",                 
                                                                    "Neutrophils","T cells regulatory (Tregs)"), 
                                         c("B cell","B cell" ,"Plasma cell","T cell","T cell","T cell" ,"T cell",
                                           "T cell","T cell","T cell","NK cell","NK cell","Monocyte","Macrophage",              
                                           "Macrophage","Macrophage","Dendritic cell","Dendritic cell","Mast cell","Mast cell","Eosinophil",                 
                                           "Neutrophil","T cell"))
  }
  
  else if (matrix_type=="Immunostates")
  {
    resu = spatialdecon(a,
                        bg = 0.1,
                        X = immunostates,
                        cell_counts = aoi_table$nuclei
    )
    
    cellsper100<-as.data.frame(resu$cell.counts$cells.per.100)
    sum_of_cols <- colSums(cellsper100 [,])
    cellsper100<-rbind(cellsper100, c(100-sum_of_cols))
    rownames(cellsper100)[nrow(cellsper100)] <-"uncharacterized cell"
    cellsper100$cell_type<-rownames(cellsper100)
    cellsper100$cell_type<- mgsub::mgsub((cellsper100$cell_type), c("CD14_positive_monocyte","CD16_positive_monocyte","CD4_positive_alpha_beta_T_cell","CD56bright_natural_killer_cell","CD56dim_natural_killer_cell",   
                                                                    "CD8_positive_alpha_beta_T_cell","MAST_cell","basophil","eosinophil","gamma_delta_T_cell","hematopoietic_progenitor",      
                                                                    "macrophage_m0","macrophage_m1","macrophage_m2","memory_B_cell","myeloid_dendritic_cell","naive_B_cell","neutrophil","plasma_cell","plasmacytoid_dendritic_cell"), 
                                         c("Monocyte","Monocyte","T cell","NK cell","NK cell",   
                                           "T cell","Mast cell","Basophil","Eosinophil","T cell","Hematopoietic progenitor",      
                                           "Macrophage","Macrophage","Macrophage","B cell","Dendritic cell","B cell","Neutrophil","Plasma cell","Plasmacytoid Dendritic cell"))
  }
  
 
  else
  { print("Check your input")
  }
  
  result_of_function= cellsper100 %>% 
    remove_rownames() %>% 
    group_by(cell_type) %>% 
    summarise_all(sum, na.rm=TRUE) 
  
  return(result_of_function)
  
}

# spatial_data=deconv_table
# single_cell_data=signature_reference


spatial_deconv_updated_without_cell_number <- function(spatial_data, single_cell_data, matrix_type = c("safeTME", "custom","LM22","TIL10","Immunostates","Glioma_custom_matrix"))
  
{  
  ## Coverting col to rownames and converting dataframe to numerical matrix
  spatial_data = spatial_data %>% 
  column_to_rownames(var="HUGO")
  a <- as.matrix(spatial_data)
  
  # Running spatialdeconv according to the type of signature matrix
  
  if ((matrix_type=="custom") & (!missing(single_cell_data)))
  {
    # Custom profile matrix 
    single_cell_data$cell<-single_cell_data$predicted_celltype
    colData(single_cell_data)$cell <- ifelse(colData(single_cell_data)$cell %in% c("CD8+ T cells", "CD4+ T cells"), "T cells", colData(single_cell_data)$cell)
    rownames(single_cell_data) <- rowData(single_cell_data)$Symbol
    single_cell_data<-single_cell_data[!(is.na(row.names(single_cell_data))),]
    
    annot<-as.data.frame(cbind(single_cell_data$Barcode,single_cell_data$cell))
    
    #x<-as.data.frame(assay(single_cell_data))
    x <- as.data.frame(as.matrix(assay(single_cell_data, "counts")))
    colnames(x)<-c(single_cell_data$Barcode)
    
    # The function does not accept a S4 object
    custom_matrix<-create_profile_matrix(mtx = x, 
                                         cellAnnots = annot, 
                                         cellTypeCol = "V2", 
                                         cellNameCol = "V1", 
                                         matrixName = "custom_mini_colon",
                                         outDir = NULL, 
                                         normalize = FALSE, 
                                         minCellNum = 5, 
                                         minGenes = 10,
                                         discardCellTypes = FALSE)
    
    resu = spatialdecon(a,
                        bg = 0.1,
                        X = custom_matrix
    )
    
    res_df <-as.data.frame(resu$prop_of_all)
    res_df$cell_type<-row.names(res_df)
    row.names(res_df)=NULL
    
    res_df$cell_type<- mgsub::mgsub((res_df$cell_type), c("B cells","Monocytes","T cells","NK cells","Dendritic cells","Progenitors" ), 
                                    c("B cell","Monocyte","T cell","NK cell","Dendritic cell","Progenitor cell"))
    res_df = res_df %>% dplyr::mutate(cell_type= sub("s$", "", cell_type))
    # Change the naming if the cell annotation differs in the single cell data
    
  }
  
  else if (matrix_type=="safeTME")
  { 
    resu = spatialdecon(a,
                        bg = 0.1,
                        X = safeTME)
    
    res_df <-as.data.frame(resu$prop_of_all)
    res_df$cell_type<-row.names(res_df)
    row.names(res_df)=NULL
    
    res_df$cell_type<- mgsub::mgsub((res_df$cell_type), c( "macrophages","mast","B.naive", "B.memory" ,
                                                           "plasma" ,"T.CD4.naive","T.CD4.memory","T.CD8.naive" ,  "T.CD8.memory" ,"NK", "pDCs",
                                                           "mDCs","monocytes.C" , "monocytes.NC.I" ,   "neutrophils","Treg", "endothelial.cells", 
                                                           "fibroblasts"), c("Macrophage","Mast cell","B cell","B cell","Plasma cell","T cell",
                                                                             "T cell","T cell","T cell","NK cell","Dendritic cell","Dendritic cell","Monocyte",
                                                                             "Monocyte","Neutrophil","T cell","Endothelial cell","Fibroblast cell"))
    
  }
  else if (matrix_type=="TIL10")
  { 
    resu = spatialdecon(a,
                        bg = 0.1,
                        X = til10)
    
    res_df <-as.data.frame(resu$prop_of_all)
    res_df$cell_type<-row.names(res_df)
    row.names(res_df)=NULL
    res_df$cell_type<- mgsub::mgsub((res_df$cell_type), c("B.cells" ,"Macrophages.M1" , "Macrophages.M2" , "Monocytes", "Neutrophils" ,"NK.cells" ,"T.cells.CD4","T.cells.CD8","Tregs","Dendritic.cells"), 
                                    c("B cell","Macrophage","Macrophage","Monocyte","Neutrophil","NK cell","T cell","T cell","T cell","Dendritic cell"))
    
  }
  
  else if (matrix_type=="LM22")
  { 
    resu = spatialdecon(a,
                        bg = 0.1,
                        X = lm22)
    
    res_df <-as.data.frame(resu$prop_of_all)
    res_df$cell_type<-row.names(res_df)
    row.names(res_df)=NULL
    res_df$cell_type<- mgsub::mgsub((res_df$cell_type), c("B cells naive","B cells memory" ,"Plasma cells","T cells CD8","T cells CD4 naive","T cells CD4 memory resting" ,"T cells CD4 memory activated",
                                                          "T cells follicular helper","T cells regulatory (Tregs)","T cells gamma delta","NK cells resting","NK cells activated","Monocytes","Macrophages M0",              
                                                          "Macrophages M1","Macrophages M2","Dendritic cells resting","Dendritic cells activated","Mast cells resting","Mast cells activated","Eosinophils",                 
                                                          "Neutrophils","T cells regulatory (Tregs)"), 
                                    c("B cell","B cell" ,"Plasma cell","T cell","T cell","T cell" ,"T cell",
                                      "T cell","T cell","T cell","NK cell","NK cell","Monocyte","Macrophage",              
                                      "Macrophage","Macrophage","Dendritic cell","Dendritic cell","Mast cell","Mast cell","Eosinophil",                 
                                      "Neutrophil","T cell"))
  }
  
  else if (matrix_type=="Immunostates")
  { 
    resu = spatialdecon(a,
                        bg = 0.1,
                        X = immunostates)
    
    res_df <-as.data.frame(resu$prop_of_all)
    res_df$cell_type<-row.names(res_df)
    row.names(res_df)=NULL
    res_df$cell_type<- mgsub::mgsub((res_df$cell_type), c("CD14_positive_monocyte","CD16_positive_monocyte","CD4_positive_alpha_beta_T_cell","CD56bright_natural_killer_cell","CD56dim_natural_killer_cell",   
                                                          "CD8_positive_alpha_beta_T_cell","MAST_cell","basophil","eosinophil","gamma_delta_T_cell","hematopoietic_progenitor",      
                                                          "macrophage_m0","macrophage_m1","macrophage_m2","memory_B_cell","myeloid_dendritic_cell","naive_B_cell","neutrophil","plasma_cell","plasmacytoid_dendritic_cell"), 
                                    c("Monocyte","Monocyte","T cell","NK cell","NK cell",   
                                      "T cell","Mast cell","Basophil","Eosinophil","T cell","Hematopoietic progenitor",      
                                      "Macrophage","Macrophage","Macrophage","B cell","Dendritic cell","B cell","Neutrophil","Plasma cell","Plasmacytoid Dendritic cell"))
  }
 
  else
  { print("Check your input")
  }
  
  result_of_function= res_df %>% 
    remove_rownames() %>% 
    group_by(cell_type) %>% 
    summarise_all(sum, na.rm=TRUE)
  
  return(result_of_function)
  
}


# For the function spotlight
#   Input: (a) A dataframe with the peseudobulk samples (GeneSymbols as rows and colnames as sample)
#          (b) A Processed SingleCellExperiment dataset
#   Output: (a) A dataframe with proportions of all the cell types

# sce_NS=signature_reference
# spatial_data_NS=deconv_table

# spatial_data_NS = wta %>% as.data.frame() %>% rownames_to_column(var="HUGO")
# sce_NS = sce

spotlight_decon<-function(spatial_data_NS, sce_NS, mean_auc_cutoff=0.8)
  
{   if (is.data.frame(spatial_data_NS)) {
      rownames(spatial_data_NS) <- spatial_data_NS$HUGO
      spatial_data_NS <- spatial_data_NS[, -1] %>% as.matrix()
  
      spe_NS <- SpatialExperiment(
      assays  = list(counts = spatial_data_NS),
      rowData = DataFrame(gene_name = rownames(spatial_data_NS)),
      colData = DataFrame(sample_id = colnames(spatial_data_NS))
      )
      } else {
      spe_NS <- spatial_data_NS
      }
    
    # Fixing the duplicated row names
    sce_NS=singleCellTK::dedupRowNames(sce_NS)
    
    ## Variance modelling
    dec <-scrapper::modelGeneVariances(assay(sce_NS,"logcounts"))
    dec <- scran::modelGeneVar(sce_NS)
    #plot(dec$mean, dec$total, xlab = "Mean log-expression", ylab = "Variance")
    #curve(metadata(dec)$trend(x), col = "blue", add = TRUE)
    
    ## Getting the top 3000 genes
    hvg <- getTopHVGs(dec, n = 3000)
    colLabels(sce_NS) <- colData(sce_NS)$cell
    
    # Get vector indicating which genes are neither ribosomal or mitochondrial
    genes <- !grepl(pattern = "^Rp[l|s]|Mt", x = rownames(sce_NS))
    
    # Compute marker genes
    mgs <- scran::scoreMarkers(sce_NS,subset.row = genes)
    # mgs_fil <- lapply(names(mgs),  function(i) {
    #   x <- mgs[[i]]
    #   # Filter and keep relevant marker genes, those with AUC > 0.8
    #   x <- x[x$mean.AUC > 0.8, ] ## Changed for testing the hodgkins lymphoma case
    #   # Sort the genes from highest to lowest weight
    #   x <- x[order(x$mean.AUC, decreasing = TRUE), ]
    #   # Add gene and cluster id to the dataframe
    #   x$gene <- rownames(x)
    #   x$cluster <- i
    #   data.frame(x)
    # })
    
    mgs_fil <- lapply(names(mgs), function(i) {
      x <- mgs[[i]]
      x <- x[x$mean.AUC > mean_auc_cutoff, ]
      if (nrow(x) == 0) return(NULL)  # skip clusters with no markers passing threshold
      x <- x[order(x$mean.AUC, decreasing = TRUE), ]
      x$gene <- rownames(x)
      x$cluster <- i
      data.frame(x)
    })
    
    # Remove NULL entries
    mgs_fil <- Filter(Negate(is.null), mgs_fil)
    mgs_df <- do.call(rbind, mgs_fil)
    
    idx <- split(seq(ncol(sce_NS)), sce_NS$cell)
    
    # downsample to at most 20 per identity & subset
    n_cells <- 20
    cs_keep <- lapply(idx, function(i) {
      n <- length(i)
      if (n < n_cells)
        n_cells <- n
      sample(i, n_cells)
    })
    sce_NS <- sce_NS[, unlist(cs_keep)]
    
    ## Running the SPOTLight function (Training the NMF model)
    colnames(sce_NS)<-sce_NS$cell
    # if some cell didn't pass the filter, it is not present, hence we subset to only keep the cells that pass the filter
    keep <- sce_NS$cell %in% unique(mgs_df$cluster)
    sce_NS_fil <- sce_NS[, keep]
    
    #mgs_df
    res <- SPOTlight(
      x = sce_NS_fil,
      y = spe_NS,
      groups = as.character(sce_NS_fil$cell),
      mgs = mgs_df,
      hvg = hvg,
      weight_id = "mean.AUC",
      group_id = "cluster",
      gene_id = "gene")
    
    mat<-res$mat
    mod<-res$NMF
    
    # plotTopicProfiles(
    #   x = mod,
    #   y = sce_NS$cell,
    #   facet = FALSE,
    #   min_prop = 0.01,
    #   ncol = 1) +
    #   theme(aspect.ratio = 1)
    # 
    # 
    # plotTopicProfiles(
    #   x = mod,
    #   y = sce_NS$cell,
    #   facet = TRUE,
    #   min_prop = 0.01,
    #   ncol = 6)
    
    
    # sign <- basis(mod)
    # colnames(sign) <- paste0("Topic", seq_len(ncol(sign)))
    
    mat_plot<-t(mat)
    mat_plot<-as.data.frame(mat_plot)
    mat_plot <- tibble::rownames_to_column(mat_plot, "cell_type")
    
    # Making the cell type singular
    mat_plot$cell_type=gsub("*s","",mat_plot$cell_type)
    return(mat_plot)

}

# For the function rctd
#   Input: (a) A dataframe with the peseudobulk samples (GeneSymbols as rows and colnames as sample)
#          (b) A Processed SingleCellExperiment dataset
#   Output: (a) A dataframe with proportions of all the cell types



# sp_data=bootstraped_data
# signature=signature_reference

# bootstraped_data= wta %>% as.data.frame() %>% rownames_to_column(var="HUGO")
# signature=sce

rctd_decon<-function(deconv_table_no_norm, signature)
{
  # Replacing NA with 0
  sp= deconv_table_no_norm %>% column_to_rownames(var="HUGO") %>% as.matrix
  sp[is.na(sp)] <- 0
  
  
  # Creating a valid SpatialExperiment object for RCTD
  spe=SpatialExperiment(assay = list(counts = sp))
  spe <- spe[!duplicated(rownames(spe)), ] # remove duplicated rows
  
  # Removing duplicates from single cell reference
  signature<- signature[!duplicated(rownames(signature)),]
  colnames(signature)=signature$Barcode
  
  # Running RCTD is "multi" mode that contstrains the regression to add to 1, and expects the cell-type
  # all of which are defined in the single cell reference
  rctd_data <- createRctd(spe, signature,cell_type_col = "cell")
  results_spe <- runRctd(rctd_data, rctd_mode = "multi", max_cores = 4)
  res= as.data.frame(as.matrix(assay(results_spe, "weights"))) %>% 
          rownames_to_column(var="cell_type") %>% 
          mutate(cell_type= sub("s$", "", cell_type))
  return(res)
}


# For the function to run cell2location from R
#   Input: (a) A matrix with the peseudobulk samples (GeneSymbols as rows and colnames as sample)
#          (b) tissue= pbmc or aml
#          (c) parameter = experimental_factor+subfactor
#          (d) min cells = no of cells in the sample/spot
#          (e) regularisation_constant = 100 (default)
#
#   Output: (a) A dataframe with proportions of all the cell types

# run_cell2location_from_r <- function(count_matrix, tissue, parameter, min_cells, regularisation_constant) {
#   
#   count_mat_py <- r_to_py(t(count_matrix))  # cell2location expects spots x genes
#   gene_names   <- r_to_py(rownames(count_matrix))
#   spot_names   <- r_to_py(colnames(count_matrix))
#   
#   # call python function
#   props <- run_cell2location(
#     count_matrix   = count_mat_py,
#     gene_names     = gene_names,
#     spot_names     = spot_names,
#     tissue         = tissue,
#     parameter      = parameter,
#     min_cells      = min_cells,
#     regularisation_constant = regularisation_constant
#   )
#   
#   return(py_to_r(props))
# }
# 
# test=run_cell2location_from_r(bootstraped_data_matrix,"pbmc","deafult",2000,100)

####
## Aim: Perform all deconvolutions for one bootstrap value 
# For the function deconvolution_for_bootstrap
#  Input: (a) a dataframe with the pseudobulk samples with the known numbers of each cell type
#          (b) character vector of the cell types in the datasets
#          (c) character vector the datasets included
#          (d) a data frame of the actual values of the cell types in the datasets included
#          
#  Output: (a) A dataframe with cell deconvolutions from Quantiseq,epic,mcp,spotlight, spatial_deco with tme and spatial decon with your signature matrix

# Modifying the function on 30 Jan 2025:
# Which trying to calculate aitchison distances, as our data it compositional, both the actual_fractions and the Deconvolution factions should add upto 1 at the end.
# However in the result of this function, I notices that when certain celltypes aren't predicted by the deconvolution methodolody for a particular Gene_panel, it doesn't
# show in the table after joining. And hence the factions don't add up to 1. So, modifying the function to have all the cells in the actual proportions added to the table
# irrespective of the status whether they are found in the prediction or not.

deconvolution_for_bootstrap <- function(deconv_table, cell_types,signature_reference, actual_cell_proportions)
  
{
  
  ## Calling the deconvolution methods
  
  # Spotlight
  spl <- tryCatch({
    spotlight_decon(deconv_table, signature_reference)
  }, error = function(e) {
    cat("Error occurred during Spotlight deconvolution. Skipping...\n")
    return(NULL)  # Return NULL if there was an error
  })
  
  if (!is.null(spl)) {
    spl$deconv_type <- "SPOTlight"
  }
  
  # Spatial deconv updated with known cell numbers (with TME)
  # sd_tme_known <- spatial_decon_updated(deconv_table,matrix_type = "safeTME",actual_cell_proportions=actual_cell_proportions) %>% 
  #                 mutate(across(-cell_type, ~ . / sum(.))) %>% 
  #                 mutate(deconv_type = "SpatialDecon (TME)\n cells known") 

  # Spatial deconv updated with known cell numbers(with custom matrix)
  # sd_cus_known <- spatial_decon_updated(deconv_table,signature_reference,matrix_type = "custom",actual_cell_proportions=actual_cell_proportions) %>% 
  #                 mutate(across(-cell_type, ~ . / sum(.))) %>% 
  #                 mutate(deconv_type="SpatialDecon (Custom)\n cells known")  
  
  # Spatial deconv updated with unknown cell numbers (with TME)
  sd_tme <- spatial_deconv_updated_without_cell_number(deconv_table,matrix_type = "safeTME") %>% 
                    mutate(deconv_type="SpatialDecon (TME)\n cells unknown")
  
  # Spatial deconv updated with unknown cell numbers(with custom matrix)
  # sd_cus_unknown <- spatial_deconv_updated_without_cell_number(deconv_table,single_cell_data = signature_reference,matrix_type = "custom") %>% 
  #                   mutate(deconv_type="SpatialDecon (Custom)\n cells unknown")
  
  # To be put in the generating matrix function, used in all of them
  deconv_mat <- deconv_table %>% 
                  column_to_rownames(var = "HUGO") %>% 
                  drop_na()
                  
  
  # Quantiseq
  quant <- deconvolute(deconv_mat, "quantiseq",tumor = FALSE) %>% 
           mutate(deconv_type="quanTIseq")
  
  
  # Epic ( Using mRNA proportion as we are using in silico samples)
  epic_res <- EPIC(deconv_mat)
  epic=t(epic_res$mRNAProportions) %>% 
            as.data.frame() %>% 
            rownames_to_column(var = "cell_type") %>% 
            add_row(cell_type = "Monocyte") %>% # Adding empty row for monocytes
            mutate(deconv_type="EPIC")
  
  # MCP
  mcp <- deconvolute(deconv_mat, "mcp_counter")
  
  # Calculating fractions for the mcp values
  mcp= mcp %>%
       mutate(across(-cell_type, ~ . / sum(.)))  %>% 
       complete(cell_type = cell_types) %>% 
       mutate(deconv_type="MCPcounter")
 
  
  print("All the deconvolutions completed successfully")
  
  # Making a singular quantitative table
  
  # Joining all the tables
  quantitative_tables <- rbind(quant,epic,mcp,if (!is.null(spl)) spl,sd_tme)
  
  
  print("All quantitative tables joined succesfully")
  
  # Grouping one of cells
  quantitative_tables$cell_type<-str_replace_all(quantitative_tables$cell_type, fixed(c("T cell CD4+ (non-regulatory)" = "T cell" ,
                                                                                        "T cell CD8+" = "T cell",
                                                                                        "T cell regulatory (Tregs)" = "T cell",
                                                                                        "T cell CD4+" = "T cell", 
                                                                                        "Macrophage M1" = "Macrophage",    
                                                                                        "Macrophage M2" = "Macrophage",
                                                                                        "CD4_Tcells"="T cell",
                                                                                        "CD8_Tcells"="T cell",
                                                                                        "Macrophages"="Macrophage", 
                                                                                        "NKcells"="NK cell",
                                                                                        "Bcells"="B cell",
                                                                                        "otherCells"="Uncharacterized cell",
                                                                                        "CAFs"="Cancer associated fibroblast",
                                                                                        "Endothelial"="Endothelial cell")))
  
  
  # Summarising all the same type of cells together
  summarised_quantitative_table = quantitative_tables %>% 
                                  dplyr::group_by(deconv_type, cell_type) %>% 
                                  dplyr::summarise_if(is.numeric, sum)
  
  print("Summarised quantitative table is being processed")
  
  # For Qualitative Plots
  
  # Adding actual proportions to the tables
   summarised_quantitative_table=left_join(summarised_quantitative_table, 
             actual_cell_proportions %>% mutate(cell_type= sub("s$", "", cell_type)), 
             by=c("cell_type"="cell_type"))
  
  return(summarised_quantitative_table)
}




## Author: Aakanksha Singh
## Date: 17 August 2023
## Version: 1.0
## Aim: Bootstrapping the deconvolution  for the  different signature matrices

# For the function processing_sce
#   Input: (a) a dataframe with the pseudobulk samples with the known numbers of each cell type
#          (b) Signature reference to be used when deconvoluting while using sce
#          (c) character vector of the cell types in the datasets
#          (d) character vector the datasets included
#          (e) a data frame of the actual values of the cell types in the datasets included
#          
#   Output: (a) A dataframe with cell deconvolutions from SpatialDecon with differnt signature matrices
#                   1. safeTME
#                   2. Custom pbmc_8k matrix
#                   3. TIL10
#                   4. LM22
#                   5. Immunostates
#                   6. Glioma_sce_matrix

signature_matrix_bootstrapping <- function(deconv_table, cell_types,dataset_included, actual_cell_proportions)
{ 
  
  # Spatial deconv with the safeTME matrix
  sd_tme <- spatial_deconv_updated_without_cell_number(deconv_table,matrix_type = "safeTME") %>% 
            dplyr::mutate(deconv_type = "SafeTME")
  
  # Spatial deconv (with custom matrix)
  sd_cus <- spatial_deconv_updated_without_cell_number(deconv_table,signature_reference,matrix_type = "custom") %>% 
             dplyr::mutate(deconv_type = "Custom")
  
  #Spatial decon with TIL10 signature
  # sd_til10 <- spatial_deconv_updated_without_cell_number(deconv_table,matrix_type = "TIL10") %>% 
  #             dplyr::mutate(deconv_type = "TIL10")
  
  sd_til10 <- tryCatch({
    spatial_deconv_updated_without_cell_number(deconv_table, matrix_type = "TIL10") %>%
      dplyr::mutate(deconv_type = "TIL10")
  }, error = function(e) {
    message("Error in TIL10 deconvolution: ", e$message)
    return(NULL)
  })
  
  # Spatial decon with LM22 signature
  sd_lm22 <- spatial_deconv_updated_without_cell_number(deconv_table,matrix_type = "LM22") %>% 
              dplyr::mutate(deconv_type = "LM22")
  
  # Spatial decon with Immunostates matrix
  # sd_immunostates <- spatial_deconv_updated_without_cell_number(deconv_table,matrix_type = "Immunostates") %>% 
  #                    dplyr::mutate(deconv_type = "Immunostates")
  sd_immunostates <- tryCatch({
    spatial_deconv_updated_without_cell_number(deconv_table, matrix_type = "Immunostates") %>%
      dplyr::mutate(deconv_type = "Immunostates")
  }, error = function(e) {
    message("Error in Immunostates deconvolution: ", e$message)
    return(NULL)
  })
  
  print("All the deconvolutions with SpatialDecon with different signature matrices completed successfully")
  
  # Making a singular qunatitative table
  
  # Joining all the tables
  quantitative_tables <- rbind(sd_tme,sd_cus,sd_til10,sd_lm22,sd_immunostates)
  
  print("All quantitative tables joined succesfully")
  
  # Grouping one of cells
  quantitative_tables$cell_type<-str_replace_all(quantitative_tables$cell_type, fixed(c("T cell CD4+ (non-regulatory)" = "T cell" ,
                                                                                        "T cell CD8+" = "T cell",
                                                                                        "T cell regulatory (Tregs)" = "T cell",
                                                                                        "T cell CD4+" = "T cell", 
                                                                                        "Macrophage M1" = "Macrophage",    
                                                                                        "Macrophage M2" = "Macrophage")))
  # Summarising all the same type of cells together
  summarised_quantitative_table = quantitative_tables %>% 
    dplyr::group_by(deconv_type, cell_type) %>% 
    dplyr::summarise_if(is.numeric, sum)
  
  print("Summarised quantitative table is being processed")
  
  # For Qualitative Plots
  
  # Adding actual proportions to the tables
  summarised_quantitative_table=left_join(summarised_quantitative_table, 
                                          actual_cell_proportions %>% mutate(cell_type= sub("s$", "", cell_type)), 
                                          by=c("cell_type"="cell_type"))
  
  return(summarised_quantitative_table)
}


# Updating mcp for the scaling function
mcp_deconvolution_summary <- function(deconv_table,actual_cell_proportions)
  
{ #Step 0: TPM normalisation
  deconv_mat <- deconv_table %>% 
    column_to_rownames(var = "HUGO") %>% 
    drop_na()
  
  
  # Step 1: Deconvolution with mcp
  mcp <- deconvolute(deconv_mat, "mcp_counter")
  
  # Step 2: Normalising the MCP results on a scale of 0 to 1
  min_max_df <- tibble( Min_val = rowMins(as.matrix(mcp[, -1])),
                        Max_val = rowMaxs(as.matrix(mcp[, -1])))
  
  
  mcp_res <- mcp
  mcp_res[, -1] <- NA
  
  for (i in 1:nrow(mcp)) 
  {
    for (j in 2:ncol(mcp)) 
    { mcp_res[i, j] <- (mcp[i, j] - min_max_df$Min_val[i]) /
      (min_max_df$Max_val[i] - min_max_df$Min_val[i])
    }
  }
  
  mcp_res$deconv_type="MCPcounter*"
  
  # Step 3: Normalising the actual cell values in the same way for comparison  
  # Applying the same normalisation!? in the data table with the actual cell counts. 
  # Performed it on the columns "actual_values", as didn't want to   
  # perfrom dual normalisation on actual_proportions  
  
  actual_cell_proportions_for_mcp=actual_cell_proportions %>% 
    dplyr::select("cell_type", ends_with(".actual_values"))
  
  min_max_df_actual_cell <- tibble( Min_val = rowMins(as.matrix(actual_cell_proportions_for_mcp[, -1])),
                                    Max_val = rowMaxs(as.matrix(actual_cell_proportions_for_mcp[, -1])))
  
  actual_cell_proportions_for_mcp_res <- actual_cell_proportions_for_mcp
  actual_cell_proportions_for_mcp_res[, -1] <- NA
  
  for (i in 1:nrow(actual_cell_proportions_for_mcp)) 
  {
    for (j in 2:ncol(actual_cell_proportions_for_mcp)) 
    { actual_cell_proportions_for_mcp_res[i, j] <- (actual_cell_proportions_for_mcp[i, j] - min_max_df_actual_cell$Min_val[i]) /
      (min_max_df_actual_cell$Max_val[i] - min_max_df_actual_cell$Min_val[i])
    }
  }
  
  # Changing colnames from values to fractions and merging with the true cellular values  
  colnames(actual_cell_proportions_for_mcp_res) <- gsub("actual_values", "actual_frac", colnames(actual_cell_proportions_for_mcp_res))
  actual_cell_proportions_for_mcp_res=left_join(actual_cell_proportions_for_mcp_res,
                                                bootstraped_data_cell_counts %>% dplyr::select("cell_type", ends_with(".actual_values")),
                                                by=c("cell_type"="cell_type")) %>% 
    dplyr::select(colnames(bootstraped_data_cell_counts))
  
  # Step 4: Putting the normalised mcp results with the normalised actual values 
  summarised_quantitative_table_mcp=left_join(mcp_res, 
                                              actual_cell_proportions_for_mcp_res %>% mutate(cell_type= sub("s$", "", cell_type)), 
                                              by=c("cell_type"="cell_type"))
  
  return(summarised_quantitative_table_mcp)
}


# Updated function for deconvolution on the bootstrapped samples

deconvolution_for_bootstrap_v2 <- function(deconv_table, cell_types,signature_reference, actual_cell_proportions)
  
{
  
  ## Calling the deconvolution methods
  
  # Spotlight
  spl <- tryCatch({
    spotlight_decon(deconv_table, signature_reference)
  }, error = function(e) {
    cat("Error occurred during Spotlight deconvolution. Skipping...\n")
    return(NULL)  # Return NULL if there was an error
  })
  
  if (!is.null(spl)) {
    spl$deconv_type <- "SPOTlight"
  }
  
  
  # Spatial deconv updated with unknown cell numbers (with TME)
  sd_tme <- spatial_deconv_updated_without_cell_number(deconv_table,matrix_type = "safeTME") %>% 
    mutate(deconv_type="SpatialDecon (TME)\n cells unknown")
  
  
  # To be put in the generating matrix function, used in all of them
  deconv_mat <- deconv_table %>% 
    column_to_rownames(var = "HUGO") %>% 
    drop_na()
  
  
  # Quantiseq
  quant <- deconvolute(deconv_mat, "quantiseq",tumor = FALSE) %>% 
    mutate(deconv_type="quanTIseq")
  
  # Epic ( Using mRNA proportion as we are using in silico samples)
  epic_res <- EPIC(deconv_mat)
  epic=t(epic_res$mRNAProportions) %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "cell_type") %>% 
    mutate(deconv_type="EPIC")
  
  #Summarised mcp results
  mcp_res  =mcp_deconvolution_summary(deconv_table,actual_cell_proportions)
  
  print("All the deconvolutions completed successfully")
  
  # Making a singular quantitative table
  
  # Joining all the tables
  quantitative_tables <- rbind(quant,epic,if (!is.null(spl)) spl,sd_tme)
  
  
  print("All quantitative tables joined succesfully")
  
  # Grouping one of cells
  quantitative_tables$cell_type<-str_replace_all(quantitative_tables$cell_type, fixed(c("T cell CD4+ (non-regulatory)" = "T cell" ,
                                                                                        "T cell CD8+" = "T cell",
                                                                                        "T cell regulatory (Tregs)" = "T cell",
                                                                                        "T cell CD4+" = "T cell", 
                                                                                        "Macrophage M1" = "Macrophage",    
                                                                                        "Macrophage M2" = "Macrophage",
                                                                                        "CD4_Tcells"="T cell",
                                                                                        "CD8_Tcells"="T cell",
                                                                                        "Macrophages"="Macrophage", 
                                                                                        "NKcells"="NK cell",
                                                                                        "Bcells"="B cell",
                                                                                        "otherCells"="uncharacterized cell",
                                                                                        "CAFs"="Cancer associated fibroblast",
                                                                                        "Endothelial"="Endothelial cell",
                                                                                        )))
  
  
  # Summarising all the same type of cells together
  summarised_quantitative_table = quantitative_tables %>% 
    dplyr::group_by(deconv_type, cell_type) %>% 
    dplyr::summarise_if(is.numeric, sum)
  
  print("Summarised quantitative table is being processed")
  
  # For Qualitative Plots
  
  # Adding actual proportions to the tables
  summarised_quantitative_table2=left_join(summarised_quantitative_table, 
                                           actual_cell_proportions %>% mutate(cell_type= sub("s$", "", cell_type)), 
                                           by=c("cell_type"="cell_type")) 
  
  res=rbind(rbind(summarised_quantitative_table2, mcp_res %>% dplyr::select(colnames(summarised_quantitative_table2))))
  
  return(res)
}

# Deconvolution function with updated normalisation on the mcp samples (min max normalisation per subfactor)

#Min Max Normalisation
# Function parameters
# colname_selection= contains the substring in the colnames that needs to be 
# selected eg. ("actual_values","Raw")
# column_character_specifications_for_mcp= conatins info about the substring 
# that is needed to further perform the min_max_normalisation within specific 
# groups. This is particulary relevanet when in the bootstraped datasets which 
# conatins data with vastly different cell distributions that need to processed 
# separately
min_max_normalisation_for_mcp <- function(df, colname_selection, column_character_specifications_for_mcp = NA) {
  
  df_all <- data.frame(col = 1:nrow(df))
  
  # If column_character_specifications_for_mcp is NA or NULL, process all columns matching colname_selection
  if (is.na(column_character_specifications_for_mcp[1]) || length(column_character_specifications_for_mcp) == 0) {
    column_character_specifications_for_mcp <- ""
  }
  
  for (k in column_character_specifications_for_mcp) {
    
    df_test <- df %>%
      dplyr::select(cell_type, contains(colname_selection))
    
    # Only filter further if k is not empty
    if (k != "") {
      df_test <- df_test %>%
        dplyr::select(cell_type, contains(k))
    }
    
    min_max_df <- tibble(
      Min_val = rowMins(as.matrix(df_test[, -1])),
      Max_val = rowMaxs(as.matrix(df_test[, -1]))
    )
    
    df_res <- df_test
    df_res[, -1] <- NA
    
    for (i in 1:nrow(df_test)) {
      for (j in 2:ncol(df_test)) {
        df_res[i, j] <- (df_test[i, j] - min_max_df$Min_val[i]) /
          (min_max_df$Max_val[i] - min_max_df$Min_val[i])
      }
    }
    
    df_all <- cbind(df_all, df_res)
  }
  
  result <- df_all[, !duplicated(colnames(df_all))]
  return(result %>% select(-col))
}

# Mcp results after min max normalisation

mcp_deconvolution_summary_v2 <- function(deconv_table,actual_cell_proportions,column_character_specifications_for_mcp)
  
{ #Step 0: TPM normalisation
  deconv_mat <- deconv_table %>% 
    column_to_rownames(var = "HUGO") %>% 
    drop_na()
  
  
  # Step 1: Deconvolution with mcp
  mcp <- deconvolute(deconv_mat, "mcp_counter") %>% 
     mutate(cell_type= case_when(cell_type=="T cell CD8+"~"T cell",
                                 cell_type=="cytotoxicity score"~"Cytotoxicity Score",
                                 .default = cell_type))
  
  # Step 2: Normalising the MCP results on a scale of 0 to 1
  mcp_res <- min_max_normalisation_for_mcp(mcp,".Raw", column_character_specifications_for_mcp)
  mcp_res$deconv_type="MCPcounter*"
  
  # Step 3: Normalising the actual cell values in the same way for comparison  
  # Applying the same normalisation!? in the data table with the actual cell counts. 
  # Performed it on the columns "actual_values", as didn't want to   
  # perfrom dual normalisation on actual_proportions  
  
  actual_cell_proportions_for_mcp_res=min_max_normalisation_for_mcp(actual_cell_proportions,"actual_values", column_character_specifications_for_mcp)
  
  # Changing colnames from values to fractions and merging with the true cellular values  
  colnames(actual_cell_proportions_for_mcp_res) <- gsub("actual_values", "actual_frac", colnames(actual_cell_proportions_for_mcp_res))
  actual_cell_proportions_for_mcp_res=left_join(actual_cell_proportions_for_mcp_res,
                                                actual_cell_proportions %>% dplyr::select("cell_type", ends_with(".actual_values")),
                                                by=c("cell_type"="cell_type")) %>% 
    dplyr::select(colnames(actual_cell_proportions))
  
  # Step 4: Putting the normalised mcp results with the normalised actual values 
  summarised_quantitative_table_mcp=left_join(mcp_res, 
                                              actual_cell_proportions_for_mcp_res %>% mutate(cell_type= sub("s$", "", cell_type)), 
                                              by=c("cell_type"="cell_type"))
                                              
  
  return(summarised_quantitative_table_mcp)
}

# Undated function for deconvolution with the new min max mcp results (06/10/2025)

deconvolution_for_bootstrap_v3 <- function(deconv_table, cell_types,signature_reference, actual_cell_proportions,column_character_specifications_for_mcp)
  
{
  
  ## Calling the deconvolution methods
  
  # Spotlight
  spl <- tryCatch({
    spotlight_decon(deconv_table, signature_reference)
  }, error = function(e) {
    cat("Error occurred during Spotlight deconvolution. Skipping...\n")
    return(NULL)  # Return NULL if there was an error
  })
  
  if (!is.null(spl)) {
    spl$deconv_type <- "SPOTlight"
  }
  
  
  # Spatial deconv updated with unknown cell numbers (with TME)
  sd_tme <- spatial_deconv_updated_without_cell_number(deconv_table,matrix_type = "safeTME") %>% 
    mutate(deconv_type="SpatialDecon (TME)\n cells unknown")
  
  
  # To be put in the generating matrix function, used in all of them
  deconv_mat <- deconv_table %>% 
    column_to_rownames(var = "HUGO") %>% 
    drop_na()
  
  
  # Quantiseq
  quant <- deconvolute(deconv_mat, "quantiseq",tumor = FALSE) %>% 
    mutate(deconv_type="quanTIseq")
  
  # Epic ( Using mRNA proportion as we are using in silico samples)
  epic_res <- EPIC(deconv_mat)
  epic_res <- EPIC(deconv_mat)
  epic=t(epic_res$mRNAProportions) %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "cell_type") %>% 
    mutate(deconv_type="EPIC")
  
  #Summarised mcp results
  mcp_res  =mcp_deconvolution_summary_v2(deconv_table,actual_cell_proportions,column_character_specifications_for_mcp)
  
  print("All the deconvolutions completed successfully")
  
  # Making a singular quantitative table
  
  # Joining all the tables
  quantitative_tables <- rbind(quant,epic,if (!is.null(spl)) spl,sd_tme)
  
  
  print("All quantitative tables joined succesfully")
  
  # Grouping one of cells
  quantitative_tables$cell_type<-str_replace_all(quantitative_tables$cell_type, fixed(c("T cell CD4+ (non-regulatory)" = "T cell" ,
                                                                                        "T cell CD8+" = "T cell",
                                                                                        "T cell regulatory (Tregs)" = "T cell",
                                                                                        "T cell CD4+" = "T cell", 
                                                                                        "Macrophage M1" = "Macrophage",    
                                                                                        "Macrophage M2" = "Macrophage",
                                                                                        "CD4_Tcells"="T cell",
                                                                                        "CD8_Tcells"="T cell",
                                                                                        "Macrophages"="Macrophage", 
                                                                                        "NKcells"="NK cell",
                                                                                        "Bcells"="B cell",
                                                                                        "otherCells"="Uncharacterized cell",
                                                                                        "CAFs"="Cancer associated fibroblast",
                                                                                        "Endothelial"="Endothelial cell",
                                                                                        "T cell CD8+"="T cell",
                                                                                        "cytotoxicity score"="Cytotoxicity score")))
  
  
  # Summarising all the same type of cells together
  summarised_quantitative_table = quantitative_tables %>% 
    dplyr::group_by(deconv_type, cell_type) %>% 
    dplyr::summarise_if(is.numeric, sum)
  
  print("Summarised quantitative table is being processed")
  
  # For Qualitative Plots
  
  # Adding actual proportions to the tables
  summarised_quantitative_table2=left_join(summarised_quantitative_table, 
                                           actual_cell_proportions %>% mutate(cell_type= sub("s$", "", cell_type)), 
                                           by=c("cell_type"="cell_type")) 
  
  res=rbind(rbind(summarised_quantitative_table2, mcp_res %>% dplyr::select(colnames(summarised_quantitative_table2))))
  
  return(res)
}


# Undated function for deconvolution a different bootstrapped sample for EPIC with three cell types (13/10/2025) (Latest)
# deconv_table=bootstraped_data
# signature_reference

deconvolution_for_bootstrap_v4 <- function(deconv_table,
                                           deconv_tabe_epic, 
                                           signature_reference, 
                                           actual_cell_proportions,
                                           actual_cell_proportions_epic,
                                           column_character_specifications_for_mcp,
                                           deconv_table_no_norm)
  
{
  
  ## Calling the deconvolution methods
  
  # Spotlight
  spl <- tryCatch({
    spotlight_decon(deconv_table, signature_reference)
  }, error = function(e) {
    cat("Error occurred during Spotlight deconvolution. Skipping...\n")
    return(NULL)  # Return NULL if there was an error
  })
  
  if (!is.null(spl)) {
    spl$deconv_type <- "SPOTlight"
  }
  
  
  # Spatial deconv updated with unknown cell numbers (with TME)
  sd_tme <- spatial_deconv_updated_without_cell_number(deconv_table,matrix_type = "safeTME") %>% 
    mutate(deconv_type="SpatialDecon (TME)\n cells unknown")
  
  
  # To be put in the generating matrix function, used in all of them
  deconv_mat <- deconv_table %>% 
    column_to_rownames(var = "HUGO") %>% 
    drop_na()
  
  deconv_mat_epic <- deconv_table_epic %>% 
    column_to_rownames(var = "HUGO") %>% 
    drop_na()
  
  # Quantiseq
  quant <- deconvolute(deconv_mat, "quantiseq",tumor = FALSE) %>% 
    mutate(deconv_type="quanTIseq")
  
  # RCTD
  rctd <- rctd_decon(deconv_table_no_norm, signature_reference) %>% # require non-normalised values
          mutate(deconv_type="RCTD")
   
  missing_spots_rctd <- setdiff(colnames(deconv_mat), colnames(rctd)[2:ncol(rctd)])
  
  if (length(missing_spots_rctd) > 0) {
    missing_df <- matrix(0, 
                         nrow = nrow(rctd), 
                         ncol = length(missing_spots_rctd),
                         dimnames = list(rownames(rctd), missing_spots_rctd))
    rctd_results <- cbind(rctd, missing_df)
  }else{
    rctd_results=rctd
  }
  
  # Reorder columns to match original spot order
  rctd_results <- rctd_results[, colnames(quant)]
  
  # Epic ( Using mRNA proportion as we are using in silico samples)
  epic <- EPIC(deconv_mat_epic)
  epic_res=t(epic$mRNAProportions) %>% 
  #epic_res=t(epic$cellFractions) %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "cell_type") %>% 
    mutate(deconv_type="EPIC") %>% 
    mutate(cell_type=case_when(cell_type=="CD4_Tcells"~"T cell",
                               cell_type=="CD8_Tcells"~"T cell",
                               cell_type=="Bcells"~"B cell",
                               cell_type=="Macrophages"~"Macrophage",
                               cell_type=="NKcells"~"NK cell",
                               cell_type=="otherCells"~"uncharacterized cell",
                               cell_type=="CAFs"~"Cancer associated fibroblast",
                               cell_type=="Endothelial"~"Endothelial cell")) %>% 
    dplyr::group_by(deconv_type, cell_type) %>% 
    dplyr::summarise_if(is.numeric, sum) %>% 
    left_join(actual_cell_proportions_epic %>% mutate(cell_type= sub("s$", "", cell_type)), 
              by=c("cell_type"="cell_type")) 
    
  
  #Summarised mcp results
  mcp_res  =mcp_deconvolution_summary_v2(deconv_table,actual_cell_proportions,column_character_specifications_for_mcp)
  
  print("All the deconvolutions completed successfully")
  
  # Making a singular quantitative table
  
  # Joining all the tables
  quantitative_tables <- rbind(quant,if (!is.null(spl)) spl,sd_tme, rctd_results)
  
  
  print("All quantitative tables joined succesfully")
  
  # Grouping one of cells
  quantitative_tables$cell_type<-str_replace_all(quantitative_tables$cell_type, fixed(c("T cell CD4+ (non-regulatory)" = "T cell" ,
                                                                                        "T cell CD8+" = "T cell",
                                                                                        "T cell regulatory (Tregs)" = "T cell",
                                                                                        "T cell CD4+" = "T cell", 
                                                                                        "Macrophage M1" = "Macrophage",    
                                                                                        "Macrophage M2" = "Macrophage",
                                                                                        "CD4_Tcells"="T cell",
                                                                                        "CD8_Tcells"="T cell",
                                                                                        "Macrophages"="Macrophage", 
                                                                                        "NKcells"="NK cell",
                                                                                        "Bcells"="B cell",
                                                                                        "otherCells"="Uncharacterized cell",
                                                                                        "CAFs"="Cancer associated fibroblast")))
  
  
  # Summarising all the same type of cells together
  summarised_quantitative_table = quantitative_tables %>% 
    dplyr::group_by(deconv_type, cell_type) %>% 
    dplyr::summarise_if(is.numeric, sum)
  
  print("Summarised quantitative table is being processed")
  
  # For Qualitative Plots
  
  # Adding actual proportions to the tables
  summarised_quantitative_table2=left_join(summarised_quantitative_table, 
                                           actual_cell_proportions %>% mutate(cell_type= sub("s$", "", cell_type)), 
                                           by=c("cell_type"="cell_type")) 
  # Binding the mcp results
  res=rbind(rbind(summarised_quantitative_table2, mcp_res %>% dplyr::select(colnames(summarised_quantitative_table2))))
  
  # Binding the EPIC results
  res=rbind(rbind(res, epic_res %>% dplyr::select(colnames(res))))
  
  return(res)
}


#General theme for plots
mytheme=list(theme_foundation(base_size=14, base_family="sans")
             + theme(plot.title = element_text(face = "bold",
                                               size = rel(1.2), hjust = 0.5),
                     text = element_text(),
                     panel.background = element_rect(fill = "white", colour = NA),
                     plot.background = element_rect(fill = "white", colour = NA),
                     panel.border = element_rect(colour = NA),
                     axis.title = element_text(face = "bold",size = rel(1)),
                     axis.title.y = element_text(angle=90,vjust =2),
                     axis.title.x = element_text(vjust = -0.2),
                     axis.text = element_text(), 
                     axis.line.x = element_line(colour="black"),
                     axis.line.y = element_line(colour="black"),
                     axis.ticks = element_line(),
                     panel.grid.major = element_line(colour="#f0f0f0"),
                     panel.grid.minor = element_blank(),
                     legend.key = element_rect(colour = NA),
                     legend.position = "bottom",
                     legend.direction = "horizontal",
                     legend.key.size= unit(0.2, "cm"),
                     legend.spacing = unit(0, "cm"),
                     legend.title = element_text(face="italic"),
                     plot.margin=unit(c(10,5,5,5),"mm"),
                     strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                     strip.text = element_text(face="bold")
             ))
