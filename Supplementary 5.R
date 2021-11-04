library(magrittr)
library(clusterProfiler)
library(tidyr)
library(dplyr)
library(DOSE)
library(pathview)
library(ggplot2)
library(ggridges)
library(UpSetR)
library(enrichplot)

d <- read.csv("DEgenesEnrichment.csv")
e <- read.csv("DEgenesEnrichmentKEGGID.csv")

geneList <- d[,2]

names(geneList) <- as.character(d[,1])
geneList <- sort(geneList, decreasing = TRUE)
head(geneList)

gene.ls <- names(geneList)[abs(geneList) > 2]
head(gene.ls)


gmt.file <- read.gmt("GO_gmt.gmt")
gmt.file <- gmt.file %>% tidyr::separate(term, c("K.id","Name"), "%")
K.id2gene <- gmt.file %>% dplyr::select(K.id, gene)
K.id2name <- gmt.file %>% dplyr::select(K.id, Name)

ewp <- enricher(gene.ls, TERM2GENE = K.id2gene, TERM2NAME = K.id2name)
head(ewp)

barplot(ewp, showCategory = 20)
dotplot(ewp, showCategory = 20)

y = GSEA(geneList, TERM2GENE=K.id2gene, TERM2NAME=K.id2name, eps = 0, verbose=FALSE)
head((y))
nrow(y)
write.csv(y, "GOGSEAanalysis.csv")
gseaplot(y, "GO:0006334")

y1 <- read.csv("GOGSEAanalysis.csv")

## convert gene ID to Symbol

p1 <- cnetplot(y, foldChange=geneList)
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(y, categorySize="pvalue", foldChange=geneList)
p3 <- cnetplot(y, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))

cnetplot(y, foldChange=geneList)
p1 <- cnetplot(y, node_label="category") 
p2 <- cnetplot(y, node_label="gene") 
p3 <- cnetplot(y, node_label="all") 
p4 <- cnetplot(y, node_label="none")
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])

p1 <- heatplot(y)
p2 <- heatplot(y, foldChange=geneList)
cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])

p1 <- emapplot(y)
p2 <- emapplot(y, pie_scale=1.5)
p3 <- emapplot(y,layout="kk")
p4 <- emapplot(y, pie_scale=1.5,layout="kk") 
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])

ridgeplot(y,
          showCategory = 40,
          label_format = 50,
          core_enrichment = FALSE,
          fill = "p.adjust") + 
  labs(x = "Distribution of logFC")


upsetplot(y, 30, geom = "text")



summary(as.data.frame(y))
head(as.data.frame(y))

gseadist(y, IDs = "GO:0004568", type = "density")

nrow(y1$x)

