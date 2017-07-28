## Load dataset from file
## Read in data files
data.dir <- "../results/ECOR2HIO_24-96-RNAseq/"
library(magrittr)
hr24 <- readr::read_csv(file = file.path(data.dir,"ECOR2_over_PBS_24hr.csv")) %>% dplyr::rename(SYMBOL = X1)
hr48 <- readr::read_csv(file = file.path(data.dir,"ECOR2_over_PBS_48hr.csv")) %>% dplyr::rename(SYMBOL = X1)
hr96 <- readr::read_csv(file = file.path(data.dir,"ECOR2_over_PBS_96hr.csv")) %>% dplyr::rename(SYMBOL = X1)

## GO GSEA
library(clusterProfiler)
library(org.Hs.eg.db)
## 24hr
hr24 <- bitr(hr24[hr24$baseMean >=1,]$SYMBOL,
               fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = "org.Hs.eg.db") %>%
    dplyr::left_join(.,hr24[hr24$baseMean >=1,], by = 'SYMBOL')

hr24 <- hr24[order(-hr24$log2FoldChange),]

up.list<- sort(hr24$log2FoldChange, decreasing = TRUE)
names(up.list) <- hr24$ENTREZID

hr24.go.gsea <- gseGO(geneList     = up.list,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "BP",
                 nPerm        = 1000,
                 minGSSize    = 50,
                 maxGSSize    = 500,
                 pvalueCutoff = 1,
                 verbose      = TRUE)

write.csv(hr24.go.gsea@result, file = file.path(data.dir,"hr24.GO-GSEA.csv"), row.names = FALSE)

## 48 hr
hr48 <- bitr(hr48[hr48$baseMean >=1,]$SYMBOL,
               fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = "org.Hs.eg.db") %>%
    dplyr::left_join(.,hr48[hr48$baseMean >=1,], by = 'SYMBOL')

hr48 <- hr48[order(-hr48$log2FoldChange),]

up.list<- sort(hr48$log2FoldChange, decreasing = TRUE)
names(up.list) <- hr48$ENTREZID

hr48.go.gsea <- gseGO(geneList     = up.list,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "BP",
                 nPerm        = 1000,
                 minGSSize    = 50,
                 maxGSSize    = 500,
                 pvalueCutoff = 1,
                 verbose      = TRUE)

#save(hr48.go.gsea, file = file.path(data.dir,"hr48.GO-GSEA.RData")) ## used in Figure 5A
write.csv(hr48.go.gsea@result, file = file.path(data.dir,"hr48.GO-GSEA.csv"), row.names = FALSE)

## 96 hr
hr96 <- bitr(hr96[hr96$baseMean >=1,]$SYMBOL,
               fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = "org.Hs.eg.db") %>%
    dplyr::left_join(.,hr96[hr96$baseMean >=1,], by = 'SYMBOL')

hr96 <- hr96[order(-hr96$log2FoldChange),]

up.list<- sort(hr96$log2FoldChange, decreasing = TRUE)
names(up.list) <- hr96$ENTREZID

hr96.go.gsea <- gseGO(geneList     = up.list,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "BP",
                 nPerm        = 1000,
                 minGSSize    = 50,
                 maxGSSize    = 500,
                 pvalueCutoff = 1,
                 verbose      = TRUE)

write.csv(hr96.go.gsea@result, file = file.path(data.dir,"hr96.GO-GSEA.csv"), row.names = FALSE)

## Load dataset from file
## Read in data files
data.dir <- "../results/ECOR2HIO_24-96-RNAseq/"
library(magrittr)
hr24 <- readr::read_csv(file = file.path(data.dir,"ECOR2_over_PBS_24hr.csv")) %>% dplyr::rename(SYMBOL = X1)
hr48 <- readr::read_csv(file = file.path(data.dir,"ECOR2_over_PBS_48hr.csv")) %>% dplyr::rename(SYMBOL = X1)
hr96 <- readr::read_csv(file = file.path(data.dir,"ECOR2_over_PBS_96hr.csv")) %>% dplyr::rename(SYMBOL = X1)

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
## 24hr
hr24 <- bitr(hr24[hr24$baseMean >=1,]$SYMBOL,
               fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = "org.Hs.eg.db") %>%
    dplyr::left_join(.,hr24[hr24$baseMean >=1,], by = 'SYMBOL')

hr24 <- hr24[order(-hr24$log2FoldChange),]

up.list<- sort(hr24$log2FoldChange, decreasing = TRUE)
names(up.list) <- hr24$ENTREZID

hr24.reactome.gsea <- gsePathway(geneList     = up.list,
                 nPerm        = 1000,
                 minGSSize    = 10,
                 pvalueCutoff = 1,
                 verbose      = TRUE)

write.csv(hr24.reactome.gsea@result, file = file.path(data.dir,"hr24.REACTOME-GSEA.csv"), row.names = FALSE)

## 48 hr
hr48 <- bitr(hr48[hr48$baseMean >=1,]$SYMBOL,
               fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = "org.Hs.eg.db") %>%
    dplyr::left_join(.,hr48[hr48$baseMean >=1,], by = 'SYMBOL')

hr48 <- hr48[order(-hr48$log2FoldChange),]

up.list<- sort(hr48$log2FoldChange, decreasing = TRUE)
names(up.list) <- hr48$ENTREZID

hr48.reactome.gsea <- gsePathway(geneList     = up.list,
                 nPerm        = 1000,
                 minGSSize    = 10,
                 pvalueCutoff = 1,
                 verbose      = TRUE)

write.csv(hr48.reactome.gsea@result, file = file.path(data.dir,"hr48.REACTOME-GSEA.csv"), row.names = FALSE)

## 96 hr
hr96 <- bitr(hr96[hr96$baseMean >=1,]$SYMBOL,
               fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = "org.Hs.eg.db") %>%
    dplyr::left_join(.,hr96[hr96$baseMean >=1,], by = 'SYMBOL')

hr96 <- hr96[order(-hr96$log2FoldChange),]

up.list<- sort(hr96$log2FoldChange, decreasing = TRUE)
names(up.list) <- hr96$ENTREZID

hr96.reactome.gsea <- gsePathway(geneList     = up.list,
                 nPerm        = 1000,
                 minGSSize    = 10,
                 pvalueCutoff = 1,
                 verbose      = TRUE)

write.csv(hr96.reactome.gsea@result, file = file.path(data.dir,"hr96.REACTOME-GSEA.csv"), row.names = FALSE)
