#script to manipulate gene expression data
#setwd("~/bioinfromatics")

#load libraries
library(dplyr)
library(tidyverse)
library(GEOquery)

#read in the data
dat <- read.csv(file = "../documents/GSE183947_fpkm.csv")
dim(dat)
#get metadata
 gse <- getGEO(GEO = "GSE183947", GSEMatrix = TRUE)

 gse
 
metadata <- pData(phenoData(gse[[1]]))
head(metadata)
metadata.subset <- select(metadata, c(1,10,11,17))
metadata.modified <- metadata %>%
  select(1,10,11,17) %>%
rename(tissue = characteristics_ch1)%>%
rename(metastasis = characteristics_ch1.1)%>%
mutate(tissue = gsub("tissue:", "",tissue))%>%
mutate(metastasis = gsub("metastasis:", "",metastasis))

head(dat)

#reshaping data
dat.long <- dat %>%
  rename(gene = X) %>%
  gather(key = 'samples', value = 'FPKM', -gene)%>%
  head()

#join dataframes = dat.long + metadata.modified
dat.long <- dat.long %>%
  left_join(., metadata.modified, by = c("samples" = "description"))


#explore data
dat.long %>%
  filter(gene == 'BRCA1' | gene == 'BRCA2') %>%
  head()



uniprot_id= c("P17787", "Q9ERK7", "P09484", "A0A3Q1MJN8", "E7F4S7", "H3B5V5", "H9GMJ0", "A4IIS6")

prot_Seq=UniprotR::GetSequences(uniprot_id)

class(prot_Seq)
View(prot_Seq)
library(bio3d)
prot_Seq
prot_Seq$Sequence
prot_Seq[1]

Prot_seq = prot_Seq$Sequence
Prot_seq
pdb <- read.pdb("Prot_seq.pdb")
install.packages("R.utils")
BiocManager::install("Rsubread")


library(ggplot2)
project = read.csv("Book1 (1).csv")
pwd()
setwd()
dr()

project = read.csv("Book (1).csv")
p = ggplot(project, aes(x=Organism, y=Number)) + geom_bar(stat = "identity", fill="blue") + ggtitle("Distribution of Proteins in Bacteria") + theme(legend.position = "none")+coord_flip()

p = p + xlab("Organism") + ylab("Number of Proteins") + geom_text(aes(label = Number), hjust = -0.3, size = 2.5)
p
ggsave("Project.pdf", width = 18, height = 12)

install.packages(compG)

library(pheatmap)
mat=readRDS(expfile)

devtools::install_github("compgenomr/compGenomRData")



library(pheatmap)
expFile=system.file( "extdata","leukemiaExpressionSubset.rds",
                    package="compGenomRData")
mat=readRDS(expFile)
#set the leukemia type annotation for each sample
annotation_col = data.frame(
  LeukemiaType =substr(colnames(mat),1,3))
rownames(annotation_col)=colnames(mat)

pheatmap(mat,show_rownames=FALSE,show_colnames=FALSE,
         annotation_col=annotation_col,
         scale = "none",clustering_method="ward.D2",
         clustering_distance_cols="euclidean")



README.md









