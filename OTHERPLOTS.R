library(ggplot2)
library(ggrepel)
library(ggpubr)
library(dplyr)


setwd("C:/Users/zjf/Desktop/skeleton muscle")

df <- read.table("HDCKOvsWT_DEGs.csv",sep = ",", header = T)
df$log2FoldChange <- as.numeric(df$log2FoldChange)
df$log2FoldChange <- -df$log2FoldChange


df$group ="none"
df$group[which((df$log2FoldChange > 1) & (df$p.value < 0.05))] = "up" 
df$group[which((df$log2FoldChange < -1) & (df$p.value < 0.05))] = "down"

gene <- c("Hdc", "Cxcr2", "Mmp9", "Vegfb", "Cxcl1", "Cxcl2", "Cxcl3","Cxcl5")
df$label=""
df$label[match(gene,df$gene_symbol)] <- gene 


df$color <- ifelse(df$group == "none" & df$label == "", "color1",   #color1非差异基因
                   ifelse(df$group == "up" & df$label == "", "color2",   #color2上调的差异基因
                          ifelse(df$group == "down" & df$label == "", "color3",  #color3下调的差异基因
                                 ifelse(df$group == "up" & df$label != "", "color4", "color5")))) 




df$logpvalue <- -log10(df$p.value)
df$p.value <- as.numeric(df$p.value)

# 绘图
df <- arrange(df, color)  


pdf("test.pdf")

p <- ggscatter(df,
               x="log2FoldChange",  
               y="logpvalue",  
               color = "color",  
               palette = c("#bcbcbc","#ffab84","#8abddc","#be000e","#0051a6"), 
               label = df$label, 
               font.label = c(15,"plain","black"), 
               repel = T ) + 
  labs(title="HDCKO vs WT",  
       x=expression(paste(Log[2], 'Fold Change')), 
       y=expression(paste(-Log[10], 'P-value')))+  
  theme(plot.title = element_text(hjust = 0.5, vjust = 0.5), 
        text = element_text(size = 15),  
        legend.position="none")  
print(p)

dev.off()


table(df$label)

df$sig <- "none"
df$sig[which((df$log2FoldChange > 3) & (df$p.value < 0.05))] = "Y" 
df$sig[which((df$log2FoldChange < -3) & (df$p.value < 0.05))] = "Y" 

DEGsLID <- df[with(df, df$sig == "Y"), ]


genelist <- bitr(DEGsLID$gene_symbol, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Mm.eg.db')
genelist <- pull(genelist,ENTREZID)               
ekegg <- enrichKEGG(gene = genelist, organism = 'mmu')
KEGG  = ekegg@result
KEGG$Description <- sub("- Mus musculus \\(house mouse\\)", "", KEGG$Description)
#计算Rich Factor（富集因子）：
Enrichment_KEGG2 <- mutate(KEGG,
                           RichFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
#计算Fold Enrichment（富集倍数）：
Enrichment_KEGG2 <- mutate(Enrichment_KEGG2, 
                           FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))
head(Enrichment_KEGG2@result)


p1 <- dotplot(Enrichment_KEGG2,x = "GeneRatio",color = "p.adjust",showCategory = 10) +
  facet_grid(rows = vars(category),scales = 'free_y',space = 'free_y') +
  scale_color_gradientn(colors = c('#BF1E27','#FEB466','#F9FCCB','#6296C5','#38489D'))
