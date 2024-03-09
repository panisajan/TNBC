
name = int$name
name2 = gsub('\\_.*','', int$name)
df = norm_counts[name,]
df_plot = stack(as.data.frame(df))
df_plot$class = rep(ss_tb_ca$group,each=length(name))
df_plot$name = rep(name,15)
df_plot$name = gsub('\\_.*','', df_plot$name)
df_plot$name = factor(df_plot$name,level=name2)
A_bca = ggplot(df_plot, aes(x=name, y=values, col=class,shape=class)) +  
  theme_light() +
  theme(axis.text=element_text(size=10),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  theme(axis.text.x = element_text(angle = 70, vjust = 0.5, hjust=0.5))+
  geom_point(alpha=0.5,size=3)+theme(axis.title.x=element_blank())


interest = int$name
vsd <- vst(dds_All, blind=FALSE)
pcaData = plotPCA(vsd[interest], intgroup=c("group"), returnData=TRUE) 
percentVar <-  round(100*attr(pcaData, "percentVar"))
pca_plot_bca = ggplot(pcaData, aes(PC1, PC2, color=group, shape=group)) +
  geom_point(size=4,alpha=0.6) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme_light()

df_plot_bca = norm_counts_r[interest,]
name2 = gsub('\\_.*','', interest)
rownames(df_plot_bca) = name2

corMatrix <- cor(df_plot_bca,use="c")
sample = dplyr::select(ss_tb_ca,group)
annoCol<-list(group=c(high="#76dddf", low="#f8766d"))

heat_bca = pheatmap(corMatrix,annotation_col=sample,annotation_colors = annoCol,legend =TRUE,annotation_legend = FALSE,
                    annotation_names_col = FALSE,fontsize_row = 7,fontsize_col = 7,breaks=seq(-1, 1, length.out=101))

hs = ggarrange(heat_bca$gtable,pca_plot_bca,ncol = 2,heights = 0.5,#widths = c(0.7,0.7),
               labels = c("a",'b')) 