# Get functional gene annotation using eggnogg protein mappings
# R script
# by J.

####~load libraries~~~~~~~~~~####
library(phylotools)
library(stringr)
library(readr)
library(rstudioapi)

####~housekeeping~~~~~~~~~~~~####
rm(list=ls()) #clear the environment
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../02_reference_data/")

###~~logfile~~~~~~~~~~~~~~~~~####
log_file=file(paste("02_get_functional_annotation_",Sys.Date(),".log",sep=""))
sink(log_file,append=TRUE,type="output")
sink(log_file,append=TRUE,type="message")
Sys.time()

####~get eggnogg annotation~~####
eggnog = read.delim("../02_reference_data/eggnog.tsv", header = FALSE)
colnames(eggnog) = eggnog[5,] #get correct column headers
eggnog = eggnog[-c(1:5),] #delete empty rows and header row
colnames(eggnog)[1] = "query" #remove hashtag from 1st column name
 
###~~make prot2gene~~~~~~~~~~####
gene2prot = read.delim("../02_reference_data/gene2prot.txt", header = FALSE)
names(gene2prot) = c("SYMBOL","REFSEQ")
prot2gene = gene2prot[,c("REFSEQ","SYMBOL")]

###~~add symbol to eggnogg~~~####
eggnog_symbols = merge(eggnog, prot2gene, by=1, all.x=TRUE)

####~get longest prot only~~~####
all_prot = read.fasta(file = "../02_reference_data/protein.faa")
all_prot$seq.name = gsub(" .*","",all_prot$seq.name) #remove everything but prot ID
all_prot$length = str_count(all_prot$seq.text) #make new column with str length
prot_lengths = all_prot[,c(1,3)] #make new df with lengths and IDs only
prot_lengths = merge(prot_lengths, prot2gene, by=1, all.x=TRUE)
prot_lengths = prot_lengths[order(prot_lengths$SYMBOL, -abs(prot_lengths$length)),] #order df by SYMBOL and prot lengths
prot_longbois = prot_lengths[ !duplicated(prot_lengths$SYMBOL),] #get new df with only longest prots

####~get nogs for longbois~~~####
eggnog_longest = subset(eggnog_symbols, query %in% prot_longbois$seq.name)
row.names(eggnog_longest) = eggnog_longest$SYMBOL #sets gene symbols as row names

####~map file for topGO~~~~~~####
gene2GO = eggnog_longest[,c("SYMBOL","GOs")] #make gene2GO annotation list
write_tsv(gene2GO, file = "../02_reference_data/gene2GO.map", col_names = FALSE)

####~write out annotation~~~~####
eggnog_longest = eggnog_longest[,-22] #removes spare gene symbol column
write.csv(eggnog_longest, file = "../02_reference_data/functional_annotation.csv")

####~fin~~~~~~~~~~~~~~~~~~~~~####
closeAllConnections()