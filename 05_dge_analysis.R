#DGE analysis script
# by ~J.
#intended for use with output from edgeR2.R script

####~load libraries~~~~~~~~~~####
library(rstudioapi)
library(ggplot2)
library(ggrepel)
library(amap)
library(reshape2)

####~housekeeping~~~~~~~~~~~~####
rm(list=ls()) #clear the environment
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #set wd
count = 00

###~~output dir~~~~~~~~~~~~~~####
output = "../06_dge_analysis/"
dir.create(output) #make folder for output
setwd(output)

###~~specify data~~~~~~~~~~~~####
ss_file="../02_reference_data/sample_sheet.csv" #sample sheet
fa_file="../02_reference_data/functional_annotation.csv" #load eggnog mappings
em_file="../04_edger/cpm.csv" #expression matrix with CPM
de_file_start="../04_edger/results_" #differential expression output from edgeR
analyses = c("PREGvNON","PREGvPRE","PREGvPOST","POSTvPRE")

###~~logfile~~~~~~~~~~~~~~~~~####
log_file=file(paste("05_dge_analysis_",Sys.Date(),".log",sep=""))
sink(log_file,append=TRUE,type="output")
sink(log_file,append=TRUE,type="message")
Sys.time()

####~big loop~~~~~~~~~~~~~~~~####
for (analysis in analyses) {

####~analysis start~~~~~~~~~~####
count = count+1
if (count < 10) {
  anal_dir=paste("0",count,"_",analysis,sep = "") #if the analysis' index is less than 10, add a trailing 0
} else {
  anal_dir=paste(count,"_",analysis,sep = "") 
} 
dir.create(anal_dir)

####~load data~~~~~~~~~~~~~~~####
ss=read.csv(ss_file,row.names=3) #loads sample sheet
fa=read.csv(fa_file,row.names=1) #loads functional annotation
em=read.csv(em_file,row.names=1) #loads expression matrix
de=read.csv(paste(de_file_start,analysis,".csv",sep = ""),row.names=1) #loads differential expression
setwd(anal_dir) #output will now go to the directory for the current analysis

####~parse data~~~~~~~~~~~~~~####

###~~expression matrix~~~~~~~####
em = em[,row.names(ss)] #select and reorder columns in em based on row names in ss
em_scaled=data.frame(t(scale(data.frame(t(em)))))
em_scaled=na.omit(em_scaled)

##~~~write out em~~~~~~~~~~~~###
write.csv(em, file = "em.csv")
write.csv(em_scaled, file = "em_scaled.csv")

###~~master~~~~~~~~~~~~~~~~~~####
master=merge(em,de,by.x=0,by.y=0) #combine DGE results with CPM to make master
row.names(master)=master[,"Row.names"]
names(master)[1]="SYMBOL"
master$mean=rowMeans(master[,2:(nrow(ss)+1)])
master$mlog10p=-log10(master$PValue)
master$sig=as.factor(master$PValue<0.05&abs(master$logFC)>1.0)
master$sigFDR=as.factor(master$FDR<0.1&abs(master$logFC)>1.0)
master=merge(master,fa,by.x=0,by.y=0) #add eggnog functional annotation to master
row.names(master)=master[,"Row.names"]
master=master[,-1]

##~~~sig genes~~~~~~~~~~~~~~~####
master_sig=subset(master,sig==TRUE)
sig_genes=master_sig$SYMBOL
em_sig=em[sig_genes,]
em_scaled_sig=em_scaled[sig_genes,]

##~~~sig genes (FDR)~~~~~~~~~####
master_fdr=subset(master,sigFDR==TRUE)
fdr_genes=master_fdr$SYMBOL
em_fdr=em[fdr_genes,]
em_scaled_fdr=em_scaled[fdr_genes,]

##~~~sig up and sig down~~~~~####
master_sig_up=subset(master_sig,logFC>0)
master_sig_down=subset(master_sig,logFC<0)
master_non_sig=subset(master,sig==FALSE)

##~~~remake master~~~~~~~~~~~####
master_non_sig$direction="ns"
master_sig_up$direction="up"
master_sig_down$direction="down"
master_sig=rbind(master_sig_up,master_sig_down)
master=rbind(master_sig,master_non_sig)
master$direction=factor(master$direction,levels=c("up","down","ns"))

##~~~write out masters~~~~~~~####
write.csv(master, file = "master.csv")
write.csv(master_sig, file = "sig.csv")
write.csv(master_fdr, file = "sig_fdr.csv")
write.csv(master_sig_up, file = "sig_up.csv")
write.csv(master_sig_down, file = "sig_down.csv")

###~~top up and down genes~~~####
order_of_p=order(master_sig_up[,"PValue"],decreasing=FALSE)
master_sig_up=master_sig_up[order_of_p,]
top5_sig_up=master_sig_up[1:5,]
top10_sig_up=master_sig_up[1:10,]
top20_sig_up=master_sig_up[1:20,]
order_of_p=order(master_sig_down[,"PValue"],decreasing=FALSE)
master_sig_down=master_sig_down[order_of_p,]
top5_sig_down=master_sig_down[1:5,]
top10_sig_down=master_sig_down[1:10,]
top20_sig_down=master_sig_down[1:20,]

###~~re-sort master~~~~~~~~~~####
master$direction=factor(master$direction,levels=c("up","down","ns"))
master$sig=factor(master$sig,levels=c("TRUE","FALSE"))
order_of_p=order(master[,"PValue"],decreasing=FALSE)
master=master[order_of_p,]

###~~gene lists~~~~~~~~~~~~~~####
all_genes = row.names(master)
write(all_genes, "gene_universe.txt")
write(sig_genes, "genes_sig.txt")
genes_non_sig = row.names(master_non_sig)
write(genes_non_sig, "genes_non_sig.txt")
genes_sig_up = row.names(master_sig_up)
write(genes_sig_up, "genes_sig_up.txt")
genes_sig_down = row.names(master_sig_down)
write(genes_sig_down, "genes_sig_down.txt")
genes_sig_fdr = row.names(master_fdr)
write(genes_sig_fdr, "genes_sig_fdr.txt")

####~theme~~~~~~~~~~~~~~~~~~~####

js_theme=theme(
  plot.title=element_text(size=14),
  axis.text.x=element_text(size=10),
  axis.text.y=element_text(size=10),
  axis.title.x=element_text(size=18),
  axis.title.y=element_text(size=18)
)

####~make plots~~~~~~~~~~~~~~####

###~~MA~~~~~~~~~~~~~~~~~~~~~~####
ma_plot=ggplot(master,aes(x=log10(mean),y=logFC,colour=direction))+
  geom_point(size=0.9)+
  labs(title="MA plot",x="P",y="Log2FC")+
  theme_bw()+
  geom_vline(xintercept=2,linetype="dashed",colour="grey",linewidth=0.5)+
  geom_hline(yintercept=0,linetype="dashed",colour="grey",linewidth=0.5)+
  scale_colour_manual(values=c("red","blue","black"),labels=c("Up","Down","Non-sig"),name="")+
  geom_label_repel(data=top5_sig_up, aes(label=SYMBOL),show.legend=FALSE)+
  geom_label_repel(data=top5_sig_down, aes(label=SYMBOL),show.legend=FALSE)
ggsave("ma.svg",plot = ma_plot)

###~~volcano~~~~~~~~~~~~~~~~~####
volcano_plot=ggplot(master,aes(x=logFC,y=mlog10p,colour=direction))+
  geom_point()+
  labs(title="Volcano plot",x="Log Fold Change",y="Log P Value")+
  theme_bw()+
  geom_vline(xintercept=-1,linetype="dashed",colour="grey",size=0.5)+
  geom_vline(xintercept=1,linetype="dashed",colour="grey",size=0.5)+
  geom_hline(yintercept=-log10(0.05),linetype="dashed",colour="grey")+
  scale_colour_manual(values=c("red","blue","black"),labels=c("Up","Down","Non-sig"),name="")+
  geom_label_repel(data=top5_sig_up, aes(label=SYMBOL),show.legend=FALSE)+
  geom_label_repel(data=top5_sig_down, aes(label=SYMBOL),show.legend=FALSE)
ggsave("volcano.svg",plot = volcano_plot)


###~~heatmap~~~~~~~~~~~~~~~~~####

##~~~allsig~~~~~~~~~~~~~~~~~~####

#~~~~make matrix~~~~~~~~~~~~~####
hm.matrix=as.matrix(em_scaled_sig) #convert df to matrix to keep gene names
y.dist=amap::Dist(hm.matrix,method="spearman")
y.cluster=hclust(y.dist,method="average")
y.dd=as.dendrogram(y.cluster)
y.dd.reorder=reorder(y.dd,0,FUN="average")
y.order=order.dendrogram(y.dd.reorder)
hm.matrix_clustered=hm.matrix[y.order,]
hm.matrix_clustered=melt(hm.matrix_clustered) #melt the matrix
names(hm.matrix_clustered)=c("gene","sample","expression")

#~~~~make the heatmap~~~~~~~~####
colours=c("blue","grey","red")
colorRampPalette(colours)(200)
hm=ggplot(hm.matrix_clustered,aes(x=sample,y=gene,fill=expression))+
  geom_tile()+
  scale_fill_gradientn(colours=colorRampPalette(colours)(200))+
  ylab("")+
  xlab("")+
  theme(axis.text.y=element_blank(),axis.ticks=element_blank(),
          legend.title=element_blank(),legend.spacing = unit(0.25,"cm"),
          axis.text.x=element_text(angle=45,hjust=1))
ggsave("hm.svg",plot = hm)

###~~boxplots~~~~~~~~~~~~~~~~####

##~~~top 10~~~~~~~~~~~~~~~~~~####
top10=master[1:10,]
candidate_genes=as.vector(row.names(top10)) #get top 10 genes as a vector

#~~~~make gene table~~~~~~~~~####
gene_data=data.frame(t(em_scaled[candidate_genes,]))
gene_data$sample_group=ss$Condition
gene_data.m=melt(gene_data,id.vars = "sample_group")
gene_data.m$sample_group=factor(gene_data$sample_group,levels=c("pre rep","pregnant","post rep")) #reorder

#~~~~make boxplot~~~~~~~~~~~~####
boxplot_candidate_genes=ggplot(gene_data.m,aes(x=variable,y=value,fill=sample_group))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+
  theme(axis.text.x=element_text(angle=45,hjust=1))#must be placed after all other theme, rotates x axis text
ggsave("boxplot_top10.svg", boxplot_candidate_genes)

#~~~~make faceted boxplot~~~~####
faceted_boxplot_top10=ggplot(gene_data.m,aes(y=value,fill=sample_group))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+
  facet_wrap(~variable,ncol=5)
ggsave("faceted_boxplot_top10.svg", faceted_boxplot_top10)

##~~~top 20~~~~~~~~~~~~~~~~~~####
top20=master[1:20,]
candidate_genes=as.vector(row.names(top20)) #get top 20 genes as a vector

#~~~~make gene table~~~~~~~~~####
gene_data=data.frame(t(em_scaled[candidate_genes,]))
gene_data$sample_group=ss$Condition
gene_data.m=melt(gene_data,id.vars = "sample_group")
gene_data.m$sample_group=factor(gene_data$sample_group,levels=c("pre rep","pregnant","post rep")) #reorder

#~~~~make boxplot~~~~~~~~~~~~####
boxplot_candidate_genes_20=ggplot(gene_data.m,aes(x=variable,y=value,fill=sample_group))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+
  theme(axis.text.x=element_text(angle=45,hjust=1))#must be placed after all other theme, rotates x axis text
ggsave("boxplot_top20.svg",plot = boxplot_candidate_genes_20)

#~~~~make faceted boxplot~~~~####
faceted_boxplot_top20=ggplot(gene_data.m,aes(y=value,fill=sample_group))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+
  facet_wrap(~variable,ncol=5)
ggsave("faceted_boxplot_top20.svg", plot = faceted_boxplot_top20)

##~~~top 5 up~~~~~~~~~~~~~~~~####
candidate_genes=as.vector(row.names(top5_sig_up)) #get top genes as a vector

#~~~~make gene table~~~~~~~~~####
gene_data=data.frame(t(em_scaled[candidate_genes,]))
gene_data$sample_group=ss$Condition
gene_data.m=melt(gene_data,id.vars = "sample_group")
gene_data.m$sample_group=factor(gene_data$sample_group,levels=c("pre rep","pregnant","post rep")) #reorder

#~~~~make boxplot~~~~~~~~~~~~####
boxplot=ggplot(gene_data.m,aes(x=variable,y=value,fill=sample_group))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+
  theme(axis.text.x=element_text(angle=45,hjust=1))#must be placed after all other theme, rotates x axis text
ggsave("5up.svg",plot = boxplot)

#~~~~make faceted boxplot~~~~####
boxplot=ggplot(gene_data.m,aes(y=value,fill=sample_group))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+
  facet_wrap(~variable,ncol=5)
ggsave("faceted_5up.svg", plot = boxplot)

##~~~top 10 up~~~~~~~~~~~~~~~####
candidate_genes=as.vector(row.names(top10_sig_up)) #get top genes as a vector

#~~~~make gene table~~~~~~~~~####
gene_data=data.frame(t(em_scaled[candidate_genes,]))
gene_data$sample_group=ss$Condition
gene_data.m=melt(gene_data,id.vars = "sample_group")
gene_data.m$sample_group=factor(gene_data$sample_group,levels=c("pre rep","pregnant","post rep")) #reorder

#~~~~make boxplot~~~~~~~~~~~~####
boxplot=ggplot(gene_data.m,aes(x=variable,y=value,fill=sample_group))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+
  theme(axis.text.x=element_text(angle=45,hjust=1))#must be placed after all other theme, rotates x axis text
ggsave("10up.svg",plot = boxplot)

#~~~~make faceted boxplot~~~~####
boxplot=ggplot(gene_data.m,aes(y=value,fill=sample_group))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+
  facet_wrap(~variable,ncol=5)
ggsave("faceted_10up.svg", plot = boxplot)

##~~~top 20 up~~~~~~~~~~~~~~~####
candidate_genes=as.vector(row.names(top20_sig_up)) #get top genes as a vector

#~~~~make gene table~~~~~~~~~####
gene_data=data.frame(t(em_scaled[candidate_genes,]))
gene_data$sample_group=ss$Condition
gene_data.m=melt(gene_data,id.vars = "sample_group")
gene_data.m$sample_group=factor(gene_data$sample_group,levels=c("pre rep","pregnant","post rep")) #reorder

#~~~~make boxplot~~~~~~~~~~~~####
boxplot=ggplot(gene_data.m,aes(x=variable,y=value,fill=sample_group))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+
  theme(axis.text.x=element_text(angle=45,hjust=1))#must be placed after all other theme, rotates x axis text
ggsave("20up.svg",plot = boxplot)

#~~~~make faceted boxplot~~~~####
boxplot=ggplot(gene_data.m,aes(y=value,fill=sample_group))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+
  facet_wrap(~variable,ncol=5)
ggsave("faceted_20up.svg", plot = boxplot)

##~~~top 5 down~~~~~~~~~~~~~~####
candidate_genes=as.vector(row.names(top5_sig_down)) #get top genes as a vector

#~~~~make gene table~~~~~~~~~####
gene_data=data.frame(t(em_scaled[candidate_genes,]))
gene_data$sample_group=ss$Condition
gene_data.m=melt(gene_data,id.vars = "sample_group")
gene_data.m$sample_group=factor(gene_data$sample_group,levels=c("pre rep","pregnant","post rep")) #reorder

#~~~~make boxplot~~~~~~~~~~~~####
boxplot=ggplot(gene_data.m,aes(x=variable,y=value,fill=sample_group))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+
  theme(axis.text.x=element_text(angle=45,hjust=1))#must be placed after all other theme, rotates x axis text
ggsave("5down.svg",plot = boxplot)

#~~~~make faceted boxplot~~~~####
boxplot=ggplot(gene_data.m,aes(y=value,fill=sample_group))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+
  facet_wrap(~variable,ncol=5)
ggsave("faceted_5down.svg", plot = boxplot)

##~~~top 10 down~~~~~~~~~~~~~####
candidate_genes=as.vector(row.names(top10_sig_down)) #get top genes as a vector

#~~~~make gene table~~~~~~~~~####
gene_data=data.frame(t(em_scaled[candidate_genes,]))
gene_data$sample_group=ss$Condition
gene_data.m=melt(gene_data,id.vars = "sample_group")
gene_data.m$sample_group=factor(gene_data$sample_group,levels=c("pre rep","pregnant","post rep")) #reorder

#~~~~make boxplot~~~~~~~~~~~~####
boxplot=ggplot(gene_data.m,aes(x=variable,y=value,fill=sample_group))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+
  theme(axis.text.x=element_text(angle=45,hjust=1))#must be placed after all other theme, rotates x axis text
ggsave("10down.svg",plot = boxplot)

#~~~~make faceted boxplot~~~~####
boxplot=ggplot(gene_data.m,aes(y=value,fill=sample_group))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+
  facet_wrap(~variable,ncol=5)
ggsave("faceted_10down.svg", plot = boxplot)

##~~~top 20 down~~~~~~~~~~~~~####
candidate_genes=as.vector(row.names(top20_sig_down)) #get top genes as a vector

#~~~~make gene table~~~~~~~~~####
gene_data=data.frame(t(em_scaled[candidate_genes,]))
gene_data$sample_group=ss$Condition
gene_data.m=melt(gene_data,id.vars = "sample_group")
gene_data.m$sample_group=factor(gene_data$sample_group,levels=c("pre rep","pregnant","post rep")) #reorder

#~~~~make boxplot~~~~~~~~~~~~####
boxplot=ggplot(gene_data.m,aes(x=variable,y=value,fill=sample_group))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+
  theme(axis.text.x=element_text(angle=45,hjust=1))#must be placed after all other theme, rotates x axis text
ggsave("20down.svg",plot = boxplot)

#~~~~make faceted boxplot~~~~####
boxplot=ggplot(gene_data.m,aes(y=value,fill=sample_group))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+
  facet_wrap(~variable,ncol=5)
ggsave("faceted_20down.svg", plot = boxplot)

####~close loop~~~~~~~~~~~~~~####
setwd("..")
}

####~end of script~~~~~~~~~~~####
closeAllConnections()