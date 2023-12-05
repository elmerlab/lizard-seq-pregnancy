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
output = "../07_dge_pathway_analysis"
dir.create(output) 
setwd(output) #set wd to output folder
GO_file="../02_reference_data/gene2GO.map"
analyses = c("PREGvNON","PREGvPRE","PREGvPOST","POSTvPRE")
count = 00

###~~logfile~~~~~~~~~~~~~~~~~####
log_file=file(paste("06_pathway_analysis_",Sys.Date(),".log",sep=""))
sink(log_file,append=TRUE,type="output")
sink(log_file,append=TRUE,type="message")

####~start loop~~~~~~~~~~~~~~####
for (analysis in analyses) {

  ####~make directory~~~~~~~~~~####
  count = count+1
  if (count < 10) {
    anal_dir=paste("0",count,"_",analysis,sep = "") #if the analysis' index is less than 10, add a trailing 0
  } else {
    anal_dir=paste(count,"_",analysis,sep = "") 
  } 
  dir.create(anal_dir)

  ####~load data~~~~~~~~~~~~~~~####
  
  ###~~gene lists~~~~~~~~~~~~~~####
  gene_universe = read_lines(paste("../06_dge_analysis/",anal_dir,"/gene_universe.txt",sep = ""))
  sig_genes = read_lines(paste("../06_dge_analysis/",anal_dir,"/genes_sig.txt",sep = ""))
  sig_up = read_lines(paste("../06_dge_analysis/",anal_dir,"/genes_sig_up.txt",sep = ""))
  sig_down = read_lines(paste("../06_dge_analysis/",anal_dir,"/genes_sig_down.txt",sep = ""))
  
  ###~~make gene2GO~~~~~~~~~~~~####
  gene2GO = readMappings(GO_file)
  
  ####~start 2nd loop~~~~~~~~~~####
  
  setwd(anal_dir)
  gene_lists=c("sig_genes","sig_up","sig_down")
  for (gene_list in gene_lists) {
  
    ###~~make geneList~~~~~~~~~~~####
    geneList = factor(as.integer(gene_universe %in% eval(as.name(gene_list)))) #gets factor vector
    names(geneList) = gene_universe
  
    ###~~make topGOdata~~~~~~~~~~####
    GOdata = new("topGOdata", 
                 ontology = "BP", 
                 description = paste("GO term analysis (BP)", analysis, sep = ""),
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
    svglite(paste(gene_list,"_p_value_histogram.svg",sep = ""))
    pvh
    dev.off()

    ###~~all significant nodes~~~####
    GO_sig = GenTable(GOdata, 
                      classic = resultFisher, 
                      orderBy = "classic", 
                      ranksOf = "classic", 
                      topNodes = resultFisher.stats[4]) 
    write_csv(GO_sig, paste(gene_list,"_topGO_sig.csv"))

    ###~~nodes graph~~~~~~~~~~~~~####
    printGraph(GOdata, 
               resultFisher, 
               firstSigNodes = 5, 
               fn.prefix = paste(gene_list,"_tGO_",sep=""), 
               useInfo = "all", 
               pdfSW = TRUE)
    
    ###~~print file~~~~~~~~~~~~####
    sink(paste(gene_list,"_details.txt",sep=""))
    print(paste("Details of topGOdata object for:",analysis,gene_list))
    GOdata #print summary of GOdata object
    print(paste("Results summary for:",analysis,gene_list))
    resultFisher #print summary of results
    sink()
    sink(log_file,append=TRUE,type="output")
    sink(log_file,append=TRUE,type="message")
    
  ####~end of 2nd loop~~~~~~~~~####
  }

####~end of loop~~~~~~~~~~~~~####
setwd("..")
}

####~fin~~~~~~~~~~~~~~~~~~~~~####
closeAllConnections()
