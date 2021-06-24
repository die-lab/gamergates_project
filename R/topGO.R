library(topGO)

setwd("~/Desktop/special-eureka/pannzer") #set the correct directory, where geneID2GO and genes_of_interest are.

geneID2GO <- readMappings(file="geneID2GO")
geneNames <- names(geneID2GO)

#the file "gene_of_interest" comes from the file w_q_de_80, from the differential_expression.R
int_genes= read.table("genes_of_interest", sep = "\t")
gene_int_list <- as.vector(int_genes$V1)

geneList <- factor(as.integer(geneNames %in% gene_int_list))
names(geneList) <- geneNames

#I have changed the ontology flag down here for performing test whith a different target.
GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
#ontology: character string specifying the ontology of interest (BP, MF or CC)
#annFUN.gene2GO this function is used when the annotations are provided as a gene-to-GOs mapping

resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
#runTest is used to apply the specified test statistic and method to the data
#classic algorithm = each GO category is tested independently. 
#elim algorithm = removes the genes mapped to significant GO terms from more general (higher level) GO terms.
#weigth algorithm = genes annotated to a GO term receive weights based on the scores of neighboring GO terms.

resultFis
allRes <- GenTable(GOdata, classicFisher = resultFis,ranksOf = "classicFisher", topNodes = 10)

write.table(allRes,file="topGO_BP.txt", quote=FALSE, row.names=FALSE, sep = "\t")




