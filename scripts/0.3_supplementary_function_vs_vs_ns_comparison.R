#title: "Getting_the_relevant_nanostring_data"
#author: "Aakanksha Singh"
#date: "2022-10-12"
#aim: to get a aggregated nanostring counts table that has corresponding samples to the visium table

library(dplyr)
library(tidyr)
library(SpatialExperiment)
library(Seurat)
library(ggplot2)
library(SpatialDecon)
library(RColorBrewer)
library(ggpubr)
library(DropletUtils)
library(Matrix)
library(tidyverse)
library(DESeq2)
library(ggrepel)
library(PCAtools)
library(ggthemes)

# For Differential gene analysis
library(stringi)
library("EnhancedVolcano")
library(pheatmap)

visium_counts_tls_V_tumor <-function(file_folder,annotated_file,Samplename)
{   colnames(annotated_file)=c("Barcode","TLS")
tumor_only<-annotated_file[annotated_file$TLS=="tumor_only",]
tls <- annotated_file[!annotated_file$TLS=="tumor_only",]


output_matrix <-read10xCounts(paste0(file_folder,"filtered_feature_bc_matrix" ))
tls_matrix <- output_matrix[,colData(output_matrix)$Barcode %in% tls$Barcode ]
tumor_matrix <- output_matrix[,colData(output_matrix)$Barcode %in% tumor_only$Barcode]

tumor_samp<-as.data.frame(Matrix::rowSums(assays(tumor_matrix)$counts))
tumor_samp<-cbind(rowData(tumor_matrix)$Symbol,tumor_samp)

# Aggregating rows with common gene symbol
tumor_samp<-aggregate(tumor_samp[,-1],by=tumor_samp[1],sum)
colnames(tumor_samp)=c("HUGO","Tumor")


tls_samp<-as.data.frame(Matrix::rowSums(assays(tls_matrix)$counts))
tls_samp<-cbind(rowData(tls_matrix)$Symbol,tls_samp)

# Aggregating rows with common gene symbol
tls_samp<-aggregate(tls_samp[,-1],by=tls_samp[1],sum)
colnames(tls_samp)=c("HUGO","TLS")

# Binding the tls_samp and the tumor samp together  
res<- full_join(tumor_samp,tls_samp, by="HUGO")

# Normalisation
res$Tumor<-(res$Tumor*1000000)/sum(res[,"Tumor"])
res$TLS <-(res$TLS*1000000)/sum(res[,"TLS"])


new_name <- c(paste0(Samplename,"_TLS"),paste0(Samplename,"_Tumor"))
old_name <- c("TLS","Tumor")

res=res %>% dplyr::rename_with(~ new_name, all_of(old_name))
return(res) 
}

get_visium_counts<- function()
{
  sample1<-"~/mnt_TLS_project/Visium/Run2-V11M15-283/Visium02_S1/outs/"
  sample1_annotated_file <- read.csv(paste0(sample1,"Visium02_S1_TLSannotation_loupe.csv"), header = TRUE)
  sample1_counts <- visium_counts_tls_V_tumor(sample1,sample1_annotated_file,"Sample1")
  
  sample2<-"~/mnt_TLS_project/Visium/Run2-V11M15-283/Visium02_S2/outs/"
  sample2_annotated_file <-read.csv(paste0(sample2,"Visium02_S2_TLSannotation_loupe.csv"), header = TRUE)
  sample2_counts <- visium_counts_tls_V_tumor(sample2,sample2_annotated_file,"Sample2")
  
  sample3<-"~/mnt_TLS_project/Visium/Run2-V11M15-283/Visium02_S3/outs/"
  sample3_annotated_file <-read.csv(paste0(sample3,"Visium02_S3_TLSannotation_loupe.csv"), header = TRUE)
  sample3_counts <- visium_counts_tls_V_tumor(sample3,sample3_annotated_file,"Sample3")
  
  sample4<-"~/mnt_TLS_project/Visium/Run2-V11M15-283/Visium02_S4/outs/"
  sample4_annotated_file <-read.csv(paste0(sample4,"Visium02_S4_TLSannotation_loupe.csv"), header = TRUE)
  sample4_counts <- visium_counts_tls_V_tumor(sample3,sample3_annotated_file,"Sample4")
  
  visium_counts_table <- merge(sample1_counts,sample2_counts, by="HUGO") %>% merge(sample3_counts, by="HUGO")%>%
    merge(sample4_counts,by = "HUGO")
}

get_visium_metadata <- function()
{
  sample1<-"~/mnt_TLS_project/Visium/Run2-V11M15-283/Visium02_S1/outs/"
  sample1_annotated_file <- read.csv(paste0(sample1,"Visium02_S1_TLSannotation_loupe.csv"), header = TRUE)
  colnames(sample1_annotated_file)=c("Barcode","TLS")
  
  sample2<-"~/mnt_TLS_project/Visium/Run2-V11M15-283/Visium02_S2/outs/"
  sample2_annotated_file <-read.csv(paste0(sample2,"Visium02_S2_TLSannotation_loupe.csv"), header = TRUE)
  colnames(sample2_annotated_file)=c("Barcode","TLS")
  
  sample3<-"~/mnt_TLS_project/Visium/Run2-V11M15-283/Visium02_S3/outs/"
  sample3_annotated_file <-read.csv(paste0(sample3,"Visium02_S3_TLSannotation_loupe.csv"), header = TRUE)
  colnames(sample3_annotated_file)=c("Barcode","TLS")
  
  sample4<-"~/mnt_TLS_project/Visium/Run2-V11M15-283/Visium02_S4/outs/"
  sample4_annotated_file <-read.csv(paste0(sample4,"Visium02_S4_TLSannotation_loupe.csv"), header = TRUE)
  colnames(sample4_annotated_file)=c("Barcode","TLS")
  
  res <- data.frame("Sample_name"=c("Sample1_Tumor_V","Sample1_TLS_V","Sample2_Tumor_V","Sample2_TLS_V",
                                    "Sample3_Tumor_V","Sample3_TLS_V","Sample4_Tumor_V","Sample4_TLS_V"),
                    "N_spots"=c(sum(sample(6:18,nrow(sample1_annotated_file[sample1_annotated_file$TLS=="tumor_only",]),replace = TRUE )),
                                sum(sample(6:18,nrow(sample1_annotated_file[!sample1_annotated_file$TLS=="tumor_only",]),replace = TRUE )),
                                sum(sample(6:18,nrow(sample2_annotated_file[sample2_annotated_file$TLS=="tumor_only",]),replace = TRUE )),
                                sum(sample(6:18,nrow(sample2_annotated_file[!sample2_annotated_file$TLS=="tumor_only",]),replace = TRUE )),
                                sum(sample(6:18,nrow(sample3_annotated_file[sample3_annotated_file$TLS=="tumor_only",]),replace = TRUE )),
                                sum(sample(6:18,nrow(sample3_annotated_file[!sample3_annotated_file$TLS=="tumor_only",]),replace = TRUE )),
                                sum(sample(6:18,nrow(sample4_annotated_file[sample4_annotated_file$TLS=="tumor_only",]),replace = TRUE )),
                                sum(sample(6:18,nrow(sample4_annotated_file[!sample4_annotated_file$TLS=="tumor_only",]),replace = TRUE )))
  )
  return(res)
}


get_nanostring_counts<- function()
{
  nanostring<- read.csv("~/mnt_TLS_project/Nanostring TAP/Results/raw_data/All_Data_WTA_q3_norm.csv",header = TRUE)
  
  ## wrangling the relevant nanostring data
  
  # Patient: 286/22
  p286.22 <- nanostring[,c("TargetName","GBM.003.Full.ROI","GBM.002.Full.ROI")] %>%
    dplyr::rename( "Sample1_Tumor"="GBM.003.Full.ROI","Sample1_TLS" = "GBM.002.Full.ROI"  ) %>% mutate("Sample1_Tumor" = (Sample1_Tumor/sum(Sample1_Tumor))*1000000, "Sample1_TLS"=(Sample1_TLS/sum(Sample1_TLS))*1000000)
  
  # Patient 328/21
  p328.21 <-nanostring[,c("TargetName","GBM.007.Full.ROI","GBM.004.Full.ROI","GBM.006.Full.ROI","GBM.014.Full.ROI")] %>% 
    dplyr::rename(  "Sample2_Tumor" = "GBM.007.Full.ROI") %>%  rowwise() %>% 
    mutate(Sample2_TLS = sum(GBM.004.Full.ROI,GBM.006.Full.ROI,GBM.014.Full.ROI)) %>% dplyr::select(TargetName,Sample2_Tumor,Sample2_TLS) %>% 
    as.data.frame()%>% mutate("Sample2_Tumor" = (Sample2_Tumor/sum(Sample2_Tumor))*1000000, "Sample2_TLS"=(Sample2_TLS/sum(Sample2_TLS))*1000000)
  
  # Patient 1803/21
  p1803.21 <-nanostring[,c("TargetName","BM.006.Full.ROI","BM.004.Full.ROI","BM.005.Full.ROI")] %>% 
    dplyr::rename("Sample3_Tumor" = "BM.006.Full.ROI") %>% rowwise() %>% 
    mutate(Sample3_TLS = sum(BM.004.Full.ROI,BM.005.Full.ROI)) %>% dplyr::select(TargetName,Sample3_Tumor,Sample3_TLS) %>% as.data.frame() %>% mutate("Sample3_Tumor" = (Sample3_Tumor/sum(Sample3_Tumor))*1000000, "Sample3_TLS"=(Sample3_TLS/sum(Sample3_TLS))*1000000)
  
  # Patient 1670/21
  p1670.21 <- nanostring[,c("TargetName","BM.002.Full.ROI","BM.001.Full.ROI","BM.003.Full.ROI")] %>% 
    dplyr::rename("Sample4_Tumor"= "BM.002.Full.ROI") %>% rowwise() %>% 
    mutate(Sample4_TLS = sum(BM.001.Full.ROI,BM.003.Full.ROI)) %>% dplyr::select(TargetName,Sample4_Tumor,Sample4_TLS) %>% as.data.frame() %>% mutate("Sample4_Tumor" = (Sample4_Tumor/sum(Sample4_Tumor))*1000000, "Sample4_TLS"=(Sample4_TLS/sum(Sample4_TLS))*1000000)
  
  ns_count_data <- merge(p286.22,p328.21, by="TargetName") %>% merge(p1803.21, by="TargetName") %>%  merge(p1670.21, by="TargetName")
  return(ns_count_data)
}
