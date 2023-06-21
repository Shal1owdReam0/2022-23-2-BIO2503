library(DESeq2)
library(dplyr)
#导入counts数据矩阵，以行为基因，列为样本
count <- read.delim("./Martrix.txt",header=TRUE,row.names=1)
## 过滤在所有重复样本中小于1的基因
count <- count[rowMeans(count)>1,]
data <- read.delim("./sample_info.txt",header = TRUE,row.names = 1)
data[,1] <- as.factor(data$condition)
all(rownames(data) %in% colnames(count))
all(rownames(data) == colnames(count))
dds <-  DESeqDataSetFromMatrix(countData = count,colData = data,design = ~ condition)
dim(dds)
dds <- dds[rowSums(counts(dds)) > 1,]
nrow(dds) 
dep <- DESeq(dds)
res <- results(dep)
diff = res
diff <- na.omit(diff) 
dim(diff)
diff
write.csv(diff,"./all_diff.csv")
#进一步筛选差异基因，使用Padj值和log2FC进行筛选
foldChange = 1
padj = 0.05
diffsig <- diff[(diff$pvalue < padj & abs(diff$log2FoldChange) > foldChange),]
dim(diffsig)
diffsig
write.csv(diffsig, "./All_diffsig.csv")
setwd(".")
count <- read.delim("./Martrix.txt",header=T,row.names=1)
DEG_data <- read.csv("./all_diff.csv",header = T,row.names = 1)

#绘制火山图
library(ggpubr) 
DEG_data$logP <- -log10(DEG_data$padj)
#基因分为三类:notsig，up，down
DEG_data$Group <- "notsig"
DEG_data$Group[which((DEG_data$padj< 0.05) & DEG_data$log2FoldChange > 1)] = "up"
DEG_data$Group[which((DEG_data$padj < 0.05) & DEG_data$log2FoldChange < -1)] = "down"
table(DEG_data$Group)

#对差异表达基因调整后的p值进行排序
DEG_data <- DEG_data[order(DEG_data$padj),]

up_label <- head(DEG_data[DEG_data$Group =="up",],5)
down_label <- head(DEG_data[DEG_data$Group == "down",],5)


deg_label_gene <- data.frame(gene = c(rownames(up_label),rownames(down_label)),
                        label = c(rownames(up_label),rownames(down_label)))


DEG_data$gene <- rownames(DEG_data)
DEG_data <- merge(DEG_data,deg_label_gene,by = 'gene',all = T)

ggscatter(DEG_data,x ="log2FoldChange",y = "logP",color = "Group",
            palette = c("green","gray","red"),
            label = DEG_data$label,
            repel = T,
            ylab = "-log10(Padj)",
            size = 2) +
theme(element_line(size = 0),element_rect(size = 1.5))+ 
scale_y_continuous(limits = c(0,35))+
scale_x_continuous(limits = c(-15,15))+
geom_hline(yintercept = 1.3,linetype = "dashed")+
geom_vline(xintercept = c(-2,2),linetype = "dashed")

#绘制热图
library(pheatmap)
#筛选差异基因，使用Padj值和LogFc筛选
diff <- DEG_data
foldChange = 1
padj = 0.05
diffsig <- diff[(diff$pvalue < padj & abs(diff$log2FoldChange) > foldChange),]
dim(diffsig)

diff_expr <- count[match(diffsig$gene,rownames(count)),]

##添加标签信息，一共有6个样本，分别为3个实验组和3个对照组
annotation_col <- data.frame(Group = factor(c(rep("Green",3), rep("Purple",3))))
rownames(annotation_col) <- colnames(diff_expr)

diff_expr <- scale(diff_expr,center = F)

pheatmap(diff_expr,
        annotation_col = annotation_col,
        color = colorRampPalette(c("blue","white","red"))(50),
        cluster_cols = F,
        show_rownames = T,
        show_colnames = T,
        scale = "row", 
        fontsize = 2,
        fontsize_row = 6,
        fontsize_col = 6,
        border = FALSE,
        treeheight_row=100,
        cellwidth = 9, cellheight = 9,
        display_numbers = TRUE,
        number_color = "black",
        number_format = "%.2f",
        legend_breaks = c(-1:1),
        legend = T,
        width = 200,
        height = 500)


com1 <- prcomp(as.matrix(t(diff_expr)), center = FALSE, scale = FALSE)
summary(com1)
df <- as.data.frame(com1$x)
head(df)

sample <- as.data.frame(row.names(df))
names(sample) <- "sample" 
df <- data.frame(df,rownames(df))
names(df)[7] <- "sample"
head(df)

df$type <- ""
df$type <- c(rep("green",3),rep("purple",3))
ggplot(data = df,aes(x=PC1,y=PC2,color= sample))+
        stat_ellipse(aes(fill=type),
                type = "norm", geom ="polygon",alpha=0.2,color=NA)+
                geom_point()+labs(x=xlab,y=ylab,color="")+guides(scale="none")