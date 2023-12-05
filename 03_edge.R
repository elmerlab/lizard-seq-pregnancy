# EdgeR analysis script 
# by ~J.
# for analysis of DGE in pregnant Z. vivipara oviduct (WV)

####~load libraries~~~~~~~~~~####
library(GenomicFeatures)
library(tximport)
library(edgeR)
library(svglite)

####~housekeeping~~~~~~~~~~~~####
rm(list=ls()) #clear the environment
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #set wd to Scripts folder

###~~output directory~~~~~~~~####
output = "../04_edger/" #specify where the output should go
dir.create(output) #create directory for output
setwd(output) #set the new output directory as the working directory

###~~specify data~~~~~~~~~~~~####
refdata = "../02_reference_data/" #specify where the reference data is kept
saldata = "../03_salmon/" #specify where the salmon quant files are

###~~logfile~~~~~~~~~~~~~~~~~####
log_file=file(paste("01_edgeR_",Sys.Date(),".log",sep=""))
sink(log_file,append=TRUE,type="output")
sink(log_file,append=TRUE,type="message")

####~load data~~~~~~~~~~~~~~~####

###~~sample sheet~~~~~~~~~~~~####
ss = read.csv(paste(refdata,"sample_sheet.csv",sep = ""), row.names = 3)

###~~quant files~~~~~~~~~~~~~####

##~~~make tx2gene~~~~~~~~~~~~####
txdb = makeTxDbFromGFF(file= paste(refdata,"annotation.gff",sep = ""), format=c("gff"))
k = keys(txdb,keytype = "TXNAME")
tx2gene = AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

##~~~load quant files~~~~~~~~####
salmonQuantFiles = file.path("../03_salmon",paste(ss$Batch,ss$Barcode,sep = "_"),"quant.sf") #makes a list of filepaths to the quant data
names(salmonQuantFiles) = row.names(ss) #associate filepaths with sampleIDs from ss 
txi = tximport(salmonQuantFiles, type = "salmon", tx2gene = tx2gene) #import salmon quant files for DGE
cts = txi$counts #get gene counts for DGE
write.csv(cts, "salmon_counts.csv") #save gene counts

####~prepare DGEList~~~~~~~~~####
y = DGEList(cts, group = ss$Condition)
y = normLibSizes(y) #normalise DGEList for library size
design = model.matrix(~ group, data = y$samples) #design for filtering
keep = filterByExpr(y, design) 
y = y[keep, ] #should keep only genes with ~10+ reads in at least one group
y = estimateDisp(y, design) #estimate dispersion

###~~get CPM~~~~~~~~~~~~~~~~~####
cpms = edgeR::cpm(y, offset = y$offset, log = FALSE)
write.csv(cpms, "cpm.csv")

###~~get log2 CPM~~~~~~~~~~~~####
logcpm = cpm(y, log = TRUE)
write.csv(logcpm, "log2cpm.csv")

###~~remove batch effect~~~~~####
logcpm.nb = removeBatchEffect(logcpm, batch = ss$Batch)
write.csv(logcpm.nb, "log2cpm_NOBATCH.csv")

####~some basic plots~~~~~~~~####
svglite("mds.svg", width = 4, height = 4)
plotMDS(y) #visualise variation between samples
dev.off()

svglite("bcv.svg", width = 4, height = 4)
plotBCV(y) #visualise dispersion estimates
dev.off()

####~GLM analysis of DGE~~~~~####

###~~fit GLM~~~~~~~~~~~~~~~~~####
design = model.matrix(~ 0 + group, data = y$samples)
colnames(design) = levels(factor(make.names(ss$Condition)))
fit = glmQLFit(y, design)

###~~compare groups~~~~~~~~~~####
my.contrasts = makeContrasts(
  PREGvPRE = pregnant-pre.rep, #compare pregnant to pre rep
  PREGvPOST = pregnant-post.rep, #compare pregnant to post rep
  POSTvPRE = post.rep-pre.rep, #compare post and pre rep
  PREGvNON = pregnant-(pre.rep+post.rep)/2, #compare preg to the average of pre and post rep
  levels = design)
qlf.PREGvPRE = glmQLFTest(fit, contrast=my.contrasts[,"PREGvPRE"])
qlf.PREGvPOST = glmQLFTest(fit, contrast=my.contrasts[,"PREGvPOST"])
qlf.POSTvPRE = glmQLFTest(fit, contrast=my.contrasts[,"POSTvPRE"])
qlf.PREGvNON = glmQLFTest(fit, contrast=my.contrasts[,"PREGvNON"])

###~~get DGE results~~~~~~~~~####
res.PREGvPRE = topTags(qlf.PREGvPRE, n=nrow(y), sort.by = "PValue")
res.PREGvPOST = topTags(qlf.PREGvPOST, n=nrow(y), sort.by = "PValue")
res.POSTvPRE = topTags(qlf.POSTvPRE, n=nrow(y), sort.by = "PValue")
res.PREGvNON = topTags(qlf.PREGvNON, n=nrow(y), sort.by = "PValue")

##~~~write out DGE tables~~~~####
write.csv(res.PREGvPRE, "results_PREGvPRE.csv")
write.csv(res.PREGvPOST, "results_PREGvPOST.csv")
write.csv(res.POSTvPRE, "results_POSTvPRE.csv")
write.csv(res.PREGvNON, "results_PREGvNON.csv")

####~fin~~~~~~~~~~~~~~~~~~~~~####
closeAllConnections()
