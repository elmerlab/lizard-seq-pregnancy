#DGE analysis script
# by ~J.
#intended for use with output from edgeR2.R script

####~load libraries~~~~~~~~~~####
library(rstudioapi)
library(reshape2)
library(ggplot2)
library(ggrepel)

####~housekeeping~~~~~~~~~~~~####

rm(list=ls()) #clear the environment
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #set wd
dir.create("../05_density_pca") #make folder for output
setwd("../05_density_pca") #set output folder as wd

###~~specify data~~~~~~~~~~~~####

ss_file=("../02_reference_data/sample_sheet.csv") #sample sheet
em_file=("../04_edger/cpm.csv") #expression matrix with CPM
log2cpm_file=("../04_edger/log2cpm.csv")

###~~logfile~~~~~~~~~~~~~~~~~####

log_file=file(paste("02_density_and_pca_",Sys.Date(),".log",sep=""))
sink(log_file,append=TRUE,type="output")
sink(log_file,append=TRUE,type="message")
Sys.time()

####~load data~~~~~~~~~~~~~~~####

em=read.csv(em_file,row.names=1) #loads expression matrix
ss=read.csv(ss_file,row.names=3) #loads sample sheet
log2cpm=read.csv(log2cpm_file,row.names = 1)

####~parse data~~~~~~~~~~~~~~####

#select and reorder columns in em based on row names in ss
em = em[,row.names(ss)]

#make a log10 expression matrix
log10cpm = log10(em)

#make a scaled expression matrix
em_scaled=data.frame(t(scale(data.frame(t(em)))))
em_scaled=na.omit(em_scaled)

#melt the em matrices
em.m=melt(em)
em_scaled.m=melt(em_scaled)
log2cpm.m=melt(log2cpm)
log10cpm.m=melt(log10cpm)

###~~save data files~~~~~~~~~####
write.csv(ss,file="ss.csv")
write.csv(em,file="em.csv")
write.csv(em_scaled,file="em_scaled.csv")

####~theme~~~~~~~~~~~~~~~~~~~####

js_theme=theme(
  plot.title=element_text(size=14),
  axis.text.x=element_text(size=10),
  axis.text.y=element_text(size=10),
  axis.title.x=element_text(size=18),
  axis.title.y=element_text(size=18)
)

####~make plots~~~~~~~~~~~~~~####

###~~density~~~~~~~~~~~~~~~~~####

##~~~faceted~~~~~~~~~~~~~~~~~####
density_plot=ggplot(em.m,aes(x=log10(value),colour=variable))+
  geom_density(alpha=0.75)+
  facet_wrap(~variable,ncol=6)+ #make ncol a variable determined by sample no for generalisability
  js_theme+
  theme(strip.background=element_rect(fill="transparent",linewidth=0),
        legend.position="none")+
  labs(x="Log10(CPM)",y="Density")

#~~~~save plot~~~~~~~~~~~~~~~####
svg("density_faceted.svg")
print(density_plot)
dev.off()

##~~~combined~~~~~~~~~~~~~~~~####
density_plot2=ggplot(em.m,aes(x=log10(value),colour=variable))+
  geom_density(alpha=0.75)+
  js_theme+
  theme(strip.background=element_rect(fill="transparent",linewidth=0),
        legend.position="right")+
  labs(x="Log10(CPM)",y="Density")

#~~~~save plot~~~~~~~~~~~~~~~####
svg("density_overlapping.svg")
print(density_plot2)
dev.off()

###~~boxplots~~~~~~~~~~~~~~~~####

##~~~scaled CPM~~~~~~~~~~~~~~####
boxplot=ggplot(em_scaled.m,aes(y=value,fill=variable))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+
  labs(title = "Expression (All Samples)", y = "Scaled CPM")
ggsave("density_boxplot_scaled_cpm.svg", plot = boxplot)

##~~~log2 CPM~~~~~~~~~~~~~~~~####
boxplot=ggplot(log2cpm.m,aes(y=value,fill=variable))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+
  labs(title = "Expression (All Samples)", y = "Log2 CPM")
ggsave("density_boxplot_log2cpm.svg", plot = boxplot)

##~~~log10 CPM~~~~~~~~~~~~~~~####
boxplot=ggplot(log10cpm.m,aes(y=value,fill=variable))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+
  labs(title = "Expression (All Samples)", y = "Log10 CPM")
ggsave("density_boxplot_log10cpm.svg", plot = boxplot, width = 10, height = 8)

###~~PCAs~~~~~~~~~~~~~~~~~~~~####
em_matrix=t(as.matrix(sapply(em_scaled,as.numeric)))
pca=prcomp(em_matrix)
pca_coord=data.frame(pca$x)

##~~~PCs 1 & 2~~~~~~~~~~~~~~~####

#~~~~make axis labels~~~~~~~~####
vars=apply(pca$x,2,var)
prop_x=round(vars["PC1"]/sum(vars),4)*100
prop_y=round(vars["PC2"]/sum(vars),4)*100
x_axis_label=paste("PC1","(",prop_x,"%)",sep="")
y_axis_label=paste("PC2","(",prop_y,"%)",sep="")

#~~~~make plot~~~~~~~~~~~~~~~####
pca_plot_1_2=ggplot(pca_coord,aes(x=PC1,y=PC2,colour=ss$Condition))+
  geom_point()+
  geom_text_repel(aes(label=row.names(ss)),show.legend=FALSE)+
  labs(title="PCA",x=x_axis_label,y=y_axis_label)
ggsave(file="pca_1_2.svg", plot = pca_plot_1_2, width = 10, height = 8)

##~~~PCs 3 & 4~~~~~~~~~~~~~~~####

#~~~~make axis labels~~~~~~~~####
vars=apply(pca$x,2,var)
prop_x=round(vars["PC3"]/sum(vars),4)*100
prop_y=round(vars["PC4"]/sum(vars),4)*100
x_axis_label=paste("PC3","(",prop_x,"%)",sep="")
y_axis_label=paste("PC4","(",prop_y,"%)",sep="")

#~~~~make plot~~~~~~~~~~~~~~~####
pca_plot_3_4=ggplot(pca_coord,aes(x=PC3,y=PC4,colour=ss$Condition))+
  geom_point()+
  geom_text_repel(aes(label=row.names(ss)),show.legend=FALSE)+
  labs(title="PCA",x=x_axis_label,y=y_axis_label)
ggsave(file="pca_3_4.svg", plot = pca_plot_3_4, width = 10, height = 8)

##~~~PCs 5 & 6~~~~~~~~~~~~~~~####

#~~~~make axis labels~~~~~~~~####
vars=apply(pca$x,2,var)
prop_x=round(vars["PC5"]/sum(vars),4)*100
prop_y=round(vars["PC6"]/sum(vars),4)*100
x_axis_label=paste("PC5","(",prop_x,"%)",sep="")
y_axis_label=paste("PC6","(",prop_y,"%)",sep="")

#~~~~make plot~~~~~~~~~~~~~~~####
pca_plot_5_6=ggplot(pca_coord,aes(x=PC5,y=PC6,colour=ss$Condition))+
  geom_point()+
  geom_text_repel(aes(label=row.names(ss)),show.legend=FALSE)+
  labs(title="PCA",x=x_axis_label,y=y_axis_label)
ggsave(file="pca_5_6.svg", plot = pca_plot_5_6, width = 10, height = 8)


####~end of script~~~~~~~~~~~####
closeAllConnections()
