# Alternative splicing analysis script 
# by ~J.
# for analysis of AS in pregnant Z. vivipara oviduct (WV)

####~load libraries~~~~~~~~~~####
library(tximport)
library(readr)
library(GenomicFeatures)
library(stageR)
library(DEXSeq)
library(edgeR)

####~housekeeping~~~~~~~~~~~~####
rm(list=ls()) #clear the environment
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #set wd to Scripts folder

###~~output directory~~~~~~~~####
output = "../08_as/" #specify where the output should go
dir.create(output) #create directory for output
setwd(output) #set the new output directory as the working directory

###~~specify data~~~~~~~~~~~~####
refdata = "../02_reference_data/" #specify where the reference data is kept
saldata = "../03_salmon/" #specify where the salmon quant files are
test_condition = "pregnant" #specify which samples to be compared to all others in sample sheet provided

###~~logfile~~~~~~~~~~~~~~~~~####
log_file=file(paste("07_as_analysis_",Sys.Date(),".log",sep=""))
sink(log_file,append=TRUE,type="output")
sink(log_file,append=TRUE,type="message")

####~load data~~~~~~~~~~~~~~~####

###~~sample sheet~~~~~~~~~~~~####
ss = read.csv(paste(refdata,"sample_sheet.csv",sep = ""), row.names = 3)
ss$sample_id = row.names(ss)
ss_test = subset(ss, Condition == test_condition)
ss_con = subset(ss, Condition != test_condition)
n.small = if (length(row.names(ss_test)) > length(row.names(ss_con))) {length(row.names(ss_con))} else {length(row.names(ss_test))}

##~~~remake ss~~~~~~~~~~~~~~~####
ss_test$sample_group = 2
ss_con$sample_group = 1
ss = rbind(ss_test,ss_con)
ss$sample_group = as.factor(ss$sample_group)
n = length(row.names(ss))

###~~quant files~~~~~~~~~~~~~####
salmonQuantFiles = file.path("../03_salmon",paste(ss$Batch,ss$Barcode,sep = "_"),"quant.sf") #makes a list of filepaths to the quant data
names(salmonQuantFiles) = row.names(ss) #associate filepaths with sampleIDs from ss 
txi = tximport(salmonQuantFiles, type = "salmon", txOut = TRUE, countsFromAbundance = "scaledTPM") #import salmon quant files for DTU
cts = txi$counts
cts = cts[rowSums(cts) > 0,] #get rid of rows with 0 counts

###~~annotation info~~~~~~~~~####
txdb = makeTxDbFromGFF(paste(refdata,"annotation.gff",sep = ""))
txdf = select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID")
tab = table(txdf$GENEID)
txdf$ntx = tab[match(txdf$GENEID, names(tab))]
txdf = txdf[match(rownames(cts),txdf$TXNAME),] #retain only transcripts actually present in the data

###~~combine everything~~~~~~####
counts = data.frame(gene_id=txdf$GENEID,
                    feature_id=txdf$TXNAME,
                    cts)

###~~write out~~~~~~~~~~~~~~~####
write.csv(counts, "salmon_gene_transcript_counts.csv")

####~DRIMSeq~~~~~~~~~~~~~~~~~####
d = dmDSdata(counts=counts,samples=ss)
d
###~~filtering~~~~~~~~~~~~~~~####
d = dmFilter(d,
             min_samps_feature_expr=n.small, min_feature_expr=10, #filter transcripts with at least 10 counts in at least n.small samples
#             min_samps_feature_prop=n.small, min_feature_prop=0.1, #filter genes by proportion of DTU - leave out for now
             min_samps_gene_expr=n, min_gene_expr=10) #filter genes with at least 10 counts in every sample

###~~design~~~~~~~~~~~~~~~~~~####
design_full = model.matrix(~sample_group, data=DRIMSeq::samples(d))

###~~analysis~~~~~~~~~~~~~~~~####
set.seed(1)
system.time({
  d = dmPrecision(d, design=design_full)
  d = dmFit(d, design=design_full)
  d = dmTest(d, coef="sample_group2")
})
res = DRIMSeq::results(d)
res.txp = DRIMSeq::results(d, level="feature")

####~stageR~~~~~~~~~~~~~~~~~~####

###~~screening df~~~~~~~~~~~~####
#detect genes with evidence of DTU
pScreen = res$pvalue
strp = function(x) substr(x,1,15)
names(pScreen) = strp(res$gene_id)

###~~confirmation df~~~~~~~~~####
#detect transcripts for those genes that are DU'd
#construct 1 column matrix of the transcript p values
pConfirmation = matrix(res.txp$pvalue, ncol = 1)
rownames(pConfirmation) = strp(res.txp$feature_id)

###~~make tx2gene~~~~~~~~~~~~####
tx2gene = res.txp[,c("feature_id","gene_id")]
for (i in 1:2) tx2gene[,i] = strp(tx2gene[,i])

###~~analysis~~~~~~~~~~~~~~~~####
stageRObj = stageRTx(pScreen=pScreen, 
                     pConfirmation=pConfirmation,
                     pScreenAdjusted=FALSE,
                     tx2gene=tx2gene)
stageRObj = stageWiseAdjustment(stageRObj, 
                                method = "dtu", 
                                alpha = 0.05)#similar to adj-P, must be specified here and cannot be changed later
suppressWarnings({
  drim.padj = getAdjustedPValues(stageRObj, 
                                 order = FALSE, 
                                 onlySignificantGenes = TRUE)
})

###~~output~~~~~~~~~~~~~~~~~~####
write.csv(drim.padj, "DRIMSeq_dtu.csv")

####~DEXSeq~~~~~~~~~~~~~~~~~~####
#using d object filtered by DRIMSeq
sample.data = DRIMSeq::samples(d)
count.data = round(as.matrix(counts(d)[,-c(1:2)]))
dxd = DEXSeqDataSet(countData = count.data,
                    sampleData = sample.data,
                    design = ~sample + exon + sample_group:exon,
                    featureID = counts(d)$feature_id,
                    groupID = counts(d)$gene_id)
system.time({
  dxd = estimateSizeFactors(dxd)
  dxd = estimateDispersions(dxd, quiet = TRUE)
  dxd = testForDEU(dxd, reducedModel = ~sample + exon)
})

###~~results~~~~~~~~~~~~~~~~~####
dxr = DEXSeqResults(dxd, independentFiltering = FALSE)
qval = perGeneQValue(dxr)
dxr.g = data.frame(gene=names(qval),qval)
columns = c("featureID", "groupID", "pvalue")
dxr = as.data.frame((dxr[,columns]))

####~stageR from DEXSeq~~~~~~####

###~~confirmation df~~~~~~~~~####
pConfirmation = matrix(dxr$pvalue, ncol=1)
dimnames(pConfirmation) = list(strp(dxr$featureID),"transcript")

###~~screening df~~~~~~~~~~~~####
pScreen = qval
names(pScreen) = strp(names(pScreen))

###~~~DTU analysis~~~~~~~~~~~####
stageRObj = stageRTx(pScreen = pScreen,
                     pConfirmation = pConfirmation,
                     pScreenAdjusted = TRUE,
                     tx2gene = tx2gene)
stageRObj = stageWiseAdjustment(stageRObj, method = "dtu", alpha=0.05)
suppressWarnings({
  dex.padj = getAdjustedPValues(stageRObj, order=FALSE, onlySignificantGenes = TRUE)
})

###~~write out results~~~~~~~####
write.csv(dex.padj, "DEXSeq_dtu.csv")

####~fin~~~~~~~~~~~~~~~~~~~~~####
closeAllConnections()
