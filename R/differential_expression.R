#build
##at first set the right working directory (where mycounts, making of mapping stats file, are)
library(NOISeq)

matrix_tab <- read.table("mycounts", sep ="\t",header = TRUE, row.names = 1)
View(matrix_tab)
factors <- data.frame(condizioni= c("w1","w2","w3","q1","q2","q3"))
#mylength<- read.table("mylength",sep="\t")

mydata <- readData(data=matrix_tab, factors = factors)
mydata
head(assayData(mydata)$exprs)


#quality
mysaturation = dat(mydata, k = 0, ndepth = 7, type = "saturation")
explo.plot(mysaturation, toplot = 1, samples = 1:6, yleftlim = NULL, yrightlim = NULL)

mycountsbio = dat(mydata, factor = NULL, type = "countsbio")
explo.plot(mycountsbio, toplot = 1, samples = NULL, plottype = "barplot")


#filter loci with low counts
myfilt10 = filtered.data(matrix_tab, factor = factors$condizioni, norm = FALSE, depth = NULL, method = 1, cv.cutoff = 100, cpm = 10, p.adj = "fdr")
mydata <- readData(data=myfilt10, factors = factors) 
mycountsbio = dat(mydata, factor = NULL, type = "countsbio")
explo.plot(mycountsbio, toplot = 1, samples = NULL, plottype = "barplot")

TMM10 = tmm(assayData(mydata)$exprs, long = 1000, lc = 0) #no length correction is applied
write.table(TMM10, file="tmm10", append = FALSE, eol="\n", quote = FALSE)


#Differential expression between two conditions
two_factors <- data.frame(condizioni =c("w","w","w","q","q","q"))
noiseq<- readData(data=TMM10, factors = two_factors)
mynoiseqbio=noiseqbio(noiseq, k=0.1, norm="n", filter=0, factor="condizioni")
mynoiseqbio_deg= degenes(mynoiseqbio, q = 0.95, M = NULL)
#q=1-FDR
#M if = "up" --> show up-regulated in condition 1; if = "down" --> show down-regulated in condition 1, if = NULL --> show all differentially expressed features
write.table(mynoiseqbio_deg, file="w_q_de", append = FALSE, eol="\n", quote = FALSE)

#plot the average expression value and highlight the feature differentially expressed
DE.plot(mynoiseqbio, q = 0.95, graphic = "expr", log.scale = TRUE)
dev.copy2pdf(file= "w_q_expr.pdf")

#plot the log-FC , M=log2FC, D= |exprCond1 - exprCond2|
DE.plot(mynoiseqbio, q = 0.95, graphic = "MD")
dev.copy2pdf(file= "w_q_DM.pdf")

#Differential expression, less significative, for GO_enrichment
two_factors <- data.frame(condizioni =c("w","w","w","q","q","q"))
noiseq<- readData(data=TMM10, factors = two_factors)
mynoiseqbio=noiseqbio(noiseq, k=0.1, norm="n", filter=0, factor="condizioni")

mynoiseqbio_deg= degenes(mynoiseqbio, q = 0.80, M = NULL)

write.table(mynoiseqbio_deg, file="w_q_de80", append = FALSE, eol="\n", quote = FALSE)

DE.plot(mynoiseqbio, q = 0.80, graphic = "expr", log.scale = TRUE)
dev.copy2pdf(file= "w_q_expr80.pdf")

DE.plot(mynoiseqbio, q = 0.80, graphic = "MD")
dev.copy2pdf(file= "w_q_DM80.pdf")

