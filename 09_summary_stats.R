library(ggVennDiagram)
library(ggvenn)
library(VennDiagram)

####~housekeeping~~~~~~~~~~~~####
rm(list=ls()) #clear the environment
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #set wd to Scripts folder

###~~output directory~~~~~~~~####
output = "../10_summary/" #specify where the output should go
dir.create(output) #create directory for output
setwd(output) #set the new output directory as the working directory

####~load data~~~~~~~~~~~~~~~####
dge_PREGvNON = read.csv("../06_dge_analysis/01_PREGvNON/sig.csv")
dge_PREGvPRE = read.csv("../06_dge_analysis/02_PREGvPRE/sig.csv")
dge_PREGvPOST = read.csv("../06_dge_analysis/03_PREGvPOST/sig.csv")
dge_POSTvPRE = read.csv("../06_dge_analysis/04_POSTvPRE/sig.csv")

fdr_PREGvNON = read.csv("../06_dge_analysis/01_PREGvNON/sig_fdr.csv")
fdr_PREGvPRE = read.csv("../06_dge_analysis/02_PREGvPRE/sig_fdr.csv")
fdr_PREGvPOST = read.csv("../06_dge_analysis/03_PREGvPOST/sig_fdr.csv")
fdr_POSTvPRE = read.csv("../06_dge_analysis/04_POSTvPRE/sig_fdr.csv")

up_PREGvNON = read.csv("../06_dge_analysis/01_PREGvNON/sig_up.csv")
up_PREGvPRE = read.csv("../06_dge_analysis/02_PREGvPRE/sig_up.csv")
up_PREGvPOST = read.csv("../06_dge_analysis/03_PREGvPOST/sig_up.csv")
up_POSTvPRE = read.csv("../06_dge_analysis/04_POSTvPRE/sig_up.csv")

down_PREGvNON = read.csv("../06_dge_analysis/01_PREGvNON/sig_down.csv")
down_PREGvPRE = read.csv("../06_dge_analysis/02_PREGvPRE/sig_down.csv")
down_PREGvPOST = read.csv("../06_dge_analysis/03_PREGvPOST/sig_down.csv")
down_POSTvPRE = read.csv("../06_dge_analysis/04_POSTvPRE/sig_down.csv")

DRIMSeq_dtu = read.csv("../08_as/DRIMSeq_dtu.csv")
DEXSeq_dtu = read.csv("../08_as/DEXSeq_dtu.csv")

####~comparisons~~~~~~~~~~~~~####

###~~all vs. all~~~~~~~~~~~~~####
x = list("PREGvPRE" = dge_PREGvPRE$SYMBOL,
         "PREGvPOST" = dge_PREGvPOST$SYMBOL,
         "POSTvPRE" = dge_POSTvPRE$SYMBOL)
venn.diagram(x = x, filename = "all_sig.png")
allvsall = ggVennDiagram(x)
ggsave("all_vs_all.svg", plot = allvsall)

###~~fdr corrected~~~~~~~~~~~####
x = list("PREGvPRE" = fdr_PREGvPRE$SYMBOL,
         "PREGvPOST" = fdr_PREGvPOST$SYMBOL,
         "POSTvPRE" = fdr_POSTvPRE$SYMBOL)
allvsall_fdr = ggVennDiagram(x)
ggsave("all_vs_all_fdr.svg", plot = allvsall_fdr)

###~~as vs dge~~~~~~~~~~~~~~~####
x = list("DEGs" = fdr_PREGvNON$SYMBOL,
         "DRIMSeq" = unique(DRIMSeq_dtu$geneID),
         "DEXSeq" = unique(DEXSeq_dtu$geneID))
asvsdge = ggVennDiagram(x)
ggsave("as_vs_dge.svg", plot = asvsdge)

###compare AS methodologies~~####
shared_as = unique(DRIMSeq_dtu$geneID[DRIMSeq_dtu$geneID %in% DEXSeq_dtu$geneID])
unique_DRIMSeq = unique(DRIMSeq_dtu$geneID[!(DRIMSeq_dtu$geneID %in% shared_as)])
unique_DEXSeq = unique(DEXSeq_dtu$geneID[!(DEXSeq_dtu$geneID %in% shared_as)])

###~~bar charts~~~~~~~~~~~~~~####
