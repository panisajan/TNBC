# setwd("C:/Users/Asus/Desktop/AomamxAsus/breast cancer")

## load packages

library(limma)
library(caTools)
library(caret)
library(pheatmap)
library(praznik)
library(DESeq2)
library(dplyr)
library(magrittr)
library(ggvenn)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(igraph)
library(ggpubr)
library(diffuStats)
library(readxl)

## open file

tb_ca_raw_data = read.csv('ca_raw_data2.csv',header = T,sep = ",",row.names = 1)
ss_tb_ca <- read.csv('TempO-Seq_clinical data_mod.csv',header = T,sep = ",",row.names = 2)[,1:2]
colnames(ss_tb_ca) = c('no','group')

#### delete dup genes #### 

### remove strings from gene in a new column
vsd.2 <- as.data.frame(tb_ca_raw_data)
head(vsd.2[1:5])
### add a column of genes with numeric strings removed
df.exp.1 <- vsd.2 %>% mutate(gene.symbol = gsub('\\_.*','', rownames(vsd.2))) %>%
  dplyr::select('gene.symbol', everything())
### any duplicates
sum(duplicated(df.exp.1$gene.symbol)) # 2834 genes with duplicates
dups.1 <- df.exp.1[ c(which( df.exp.1$gene.symbol %in% df.exp.1[ duplicated(df.exp.1$gene.symbol) ,][['gene.symbol']] )) ,]
dim(dups.1)
dups.1 <- dups.1 %>% arrange(gene.symbol) #order

uniq.1 <- df.exp.1[ - which(df.exp.1$gene.symbol %in% dups.1$gene.symbol) ,]
sum(duplicated(uniq.1$gene.symbol))

### for dups, select ones with highest maximum mean
dups.2 <- dups.1 %>% mutate(mean = apply(dups.1[,-1], 1, mean)) %>% dplyr::select('gene.symbol', 'mean', everything())
### sort them based on mean column (highest to lowest) and remove duplicate genes
dups.2 <- dups.2[ order(dups.2$mean, decreasing = T) ,]
uniq.2 <- dups.2[ - which( duplicated(dups.2$gene.symbol) ) ,] # 2198 unique genes retained
sum(duplicated(uniq.2$gene.symbol))
### remove mean column
uniq.2 <- uniq.2[,-2]

### combine two unique into one dataframe
df.exp.2 <- rbind(uniq.1, uniq.2)
dim(df.exp.2) # 19703genes 15samples

#delete dup col
df.exp.2 <- df.exp.2[,-1]

# change file into matrix
All_CA <- as.matrix(df.exp.2)

ss_tb_ca$group = factor(ss_tb_ca$group)
ss_tb_ca$group <- relevel(ss_tb_ca$group, ref = "low")
condition_All_CA <- ss_tb_ca$group

# making sure the row names in colData matches to column names in counts_data
all(colnames(All_CA) %in% rownames(ss_tb_ca))
# are they in the same order?
all(colnames(All_CA) == rownames(ss_tb_ca))
All_CA = All_CA[,rownames(ss_tb_ca)]
all(colnames(All_CA) == rownames(ss_tb_ca))


dds_All <- DESeqDataSetFromMatrix(countData=All_CA, colData=ss_tb_ca, design=~group)
# run DE
dds_All_wald <-DESeq(dds_All) ## for wald test

norm_counts_r = counts(dds_All_wald, normalized = TRUE)
norm_counts = log2(counts(dds_All_wald, normalized = TRUE)+1) ## log(x+1)

res <- results(dds_All_wald)
dat_res <- data.frame(res)
dat_res$name = rownames(dat_res)

## padj < 0.05
res_fill <- subset(dat_res,dat_res$padj < 0.05)
res_fill = res_fill[order(res_fill$pvalue),]
# write.csv(res_fill,'toptable_adj005(15).csv')
summary(abs(res_fill$log2FoldChange))

## volcano

library(ggrepel)
dat_res$diffexpressed <- "-"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
dat_res$diffexpressed[dat_res$log2FoldChange > 1 & dat_res$padj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
dat_res$diffexpressed[dat_res$log2FoldChange < -1 & dat_res$padj < 0.05] <- "DOWN"
dat_res$delabel <- NA
dat_res$delabel[dat_res$diffexpressed != "-"] <- sub("_.*", "",rownames(dat_res)[dat_res$diffexpressed != "NO"])
# dat_res = na.omit(dat_res)
dat_res$diffexpressed <- factor(dat_res$diffexpressed, levels = c('UP','DOWN','-'))
min_lfc = min(dat_res$log2FoldChange,na.rm = TRUE)
max_lfc = max(dat_res$log2FoldChange,na.rm = TRUE)
# cols <- c("UP" = "red", "DOWN" = "blue", "NO" = "grey")
B = ggplot(data=dat_res, aes(x=log2FoldChange, y=-log10(pvalue), fill=diffexpressed, label=delabel)) + 
  geom_point(aes(color = diffexpressed), alpha = 0.25,size=3)+ xlim(-10, max_lfc)+
  # geom_label_repel(aes(label = delabel),
  #                  box.padding   = 0.35, 
  #                  point.padding = 0.5,
  #                  segment.color = 'grey50')+
  scale_color_manual(values=c("red", "blue", "black")) +
  #geom_text_repel() +
  geom_text(aes(color=diffexpressed),size=2.6+0.3,hjust= 0,vjust=-0.25) + 
  #scale_colour_discrete(l = 40)+
  theme_minimal(base_size = 22)+
  theme(legend.position="right")+ labs(x = "logFC")
#B
ggsave("volcano15.png", B, width = 16, height = 8, dpi = 400, bg = 'white')

## compute MI score

df <- as.data.frame(cbind(class = ss_tb_ca$group,t(norm_counts)))
x = as.data.frame(sapply(df[,-1], as.numeric))
y = as.factor(df[,1])
test = data.frame(miScores(x,y))
test$name = rownames(test)

test = test[order(test$miScores.x..y.,decreasing = TRUE),]
high_mi = test[test$miScores.x..y. > 0.50,] ### The highest score is 0.5004024
mi_0 = test[test$miScores.x..y. == 0,]
## Venn diagram of padj<0.05 vs Top MI
s1 <- list('MI' = high_mi$name,
           'DESeq2' = rownames(res_fill))
ggvenn(s1, show_elements = F, label_sep = "\n",fill_color = c('thistle1','lightblue1'),text_size = 5)
#ggsave('venn_MI_Wald.png',dpi=300,width = 6,height = 6,bg='white')

## COnsider top MI 56 genes
MI_dat = dat_res[high_mi$name,]
MI_dat = MI_dat[order(MI_dat$padj),]
# write.csv(MI_dat,'MI_toptable.csv')
summary(abs(MI_dat$log2FoldChange))

name = MI_dat$name
name = res_fill$name[1:10]
name2 = gsub('\\_.*','', res_fill$name[1:10])


MI_dat = dat_res[high_mi$name,]
MI_dat = MI_dat[order(MI_dat$padj),]
name = MI_dat$name
name2 = gsub('\\_.*','', MI_dat$name)


int = res_fill[intersect(high_mi$name,rownames(res_fill)),]
