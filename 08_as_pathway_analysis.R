# pathway analysis script
# by J.

####~load libraries~~~~~~~~~~####
library(rstudioapi)
library(topGO)
library(genefilter)
library(Rgraphviz)
library(readr)
library(svglite)
library(ggplot2)

####~housekeeping~~~~~~~~~~~~####
rm(list=ls()) #clear the environment
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #set wd to "Scripts" folder
output = "../09_as_pathway_analysis"
dir.create(output) 
setwd(output) #set wd to output folder
GO_file="../02_reference_data/gene2GO.map"

####~load data~~~~~~~~~~~~~~~####
DRIMSeq.dtu = read.csv("../08_as/DRIMSeq_dtu.csv")
salmon = read.csv("../08_as/salmon_gene_transcript_counts.csv")
gene_list = unique(DRIMSeq.dtu$geneID)
gene_universe = unique(salmon$gene_id)

####~analysis~~~~~~~~~~~~~~~~####

###~~make geneList~~~~~~~~~~~####
geneList = factor(as.integer(gene_universe %in% gene_list)) #gets factor vector
names(geneList) = gene_universe

###~~read gene2GO~~~~~~~~~~~~####
gene2GO = readMappings(GO_file)

###~~make topGOdata~~~~~~~~~~####
GOdata = new("topGOdata", 
             ontology = "BP", 
             description = paste("GO term analysis (BP) for AS results"),
             allGenes = geneList,
             annot = annFUN.gene2GO,
             gene2GO = gene2GO)

####~enrichment testing~~~~~~####
test.stat = new("classicCount", testStatistic = GOFisherTest, name = "Fisher test") #define statistical test
resultFisher = getSigGroups(GOdata, test.stat) #get test result
resultFisher.stats = geneData(resultFisher) #get stats for the test - sig GO terms is no. 2 here

####~export results~~~~~~~~~~####

###~~p-value histogram~~~~~~~####
pvalFisher = score(resultFisher) #get P values for each GO term
pvh = hist(pvalFisher, 50, xlab = "p-values")
svglite("p_value_histogram.svg")
pvh
dev.off()

###~~all significant nodes~~~####
GO_sig = GenTable(GOdata, 
                  classic = resultFisher, 
                  orderBy = "classic", 
                  ranksOf = "classic", 
                  topNodes = resultFisher.stats[4]) 
write_csv(GO_sig, "topGO_sig.csv")

###~~nodes graph~~~~~~~~~~~~~####
printGraph(GOdata, 
           resultFisher, 
           firstSigNodes = 5, 
           fn.prefix = "_tGO_", 
           useInfo = "all", 
           pdfSW = TRUE)

###~~print file~~~~~~~~~~~~####
sink("details.txt")
print("Details of topGOdata object")
GOdata #print summary of GOdata object
print("Results summary")
resultFisher #print summary of results
sink()
