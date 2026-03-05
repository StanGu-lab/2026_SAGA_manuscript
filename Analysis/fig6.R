library(dplyr)
library(ggfortify)
library(DESeq2)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(ggrepel)
library(data.table)
library(tibble)
library(MAGeCKFlute)
library(ggvenn)
library(scales)

#########------------------------------------------------------------------------------------------------#########
#########------------------------------------------Pathway analysis--------------------------------------#########
#########------------------------------------------------------------------------------------------------#########
outdir<-"/Users/dwang20/OneDrive - Inside MD Anderson/20250416_Alok_rna_seq/fig/"
res<-fread("/Users/dwang20/OneDrive - Inside MD Anderson/20250416_Alok_rna_seq/human_Finalcounts/human_all_DEG.csv",header=T)%>% as.data.frame()
rownames(res)<-res$V1


res$gene<-res$V1
up<-res[!is.na(res$padj)&res$padj<0.05&res$log2FoldChange>0,]%>% dplyr::select(gene,log2FoldChange)%>% arrange(desc(log2FoldChange))
down<-res[!is.na(res$padj)&res$padj<0.05&res$log2FoldChange<0,]%>% dplyr::select(gene,log2FoldChange)%>% arrange((log2FoldChange))

#df<-rbind(up,down)%>% arrange(desc(log2FoldChange))
#Pathway+GOBP
ranks<- deframe(down) 
#GOBP#Pathway
enrich <- EnrichAnalyzer(ranks,keytype = "Symbol",type = "H",method = "ORT",organism = "hsa",pvalueCutoff = 0.1,
                         limit = c(10, 500),universe = NULL,filter = FALSE,gmtpath = NULL,verbose = TRUE)

re<-enrich@result
#re$ID<-gsub(":","_",re$ID)
#re$source<-sub("_.*", "", re$ID)
#result<-re
re<-re %>%filter(p.adjust<0.05)
result<-re

'''
result <- re %>%
  group_by(source) %>%
  arrange(p.adjust) %>%
  slice_head(n = 5)%>% ungroup()
'''

#result$GeneRatio<-sapply(strsplit(result$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
#result$Description<-paste0("(",result$source,") ",result$Description)
result<-result[order(result$p.adjust,decreasing = TRUE),]
result$Description<-factor(result$Description,levels=result$Description)

p<-ggplot(data = result, aes(x = -log10(p.adjust), y = Description, 
                             color = NES, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("-log10(p.adjust)")+
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 14),
    #axis.title.y = element_text(size = 14), 
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12))

ggsave(filename = "/Users/dwang20/OneDrive - Inside MD Anderson/20250416_Alok_rna_seq/fig/OVCAR3_down_hallmark_pathway.pdf", plot = p, width = 6.0, height = 4)
 
#########------------------------------------------------------------------------------------------------#########
#########----------------------------correlation analysis with TADA2B-KO rna-seq-------------------------#########
#########------------------------------------------------------------------------------------------------#########
path1<-"/Users/dwang20/OneDrive - Inside MD Anderson/Hok_rna_seq/HGS2_DE/"
Taf5l_KO<-fread(paste0(path1,"Taf5l_KO_vs_control.csv"))%>% as.data.frame()
colnames(Taf5l_KO)[2:7]<-paste0("Taf5l_KO_",colnames(Taf5l_KO)[2:7])

path2<-"/Users/dwang20/OneDrive - Inside MD Anderson/20250416_Alok_rna_seq/mouse_Finalcounts/"
drug_treatment<-fread(paste0(path2,"mouse_all_DEG.csv"))%>% as.data.frame()
colnames(drug_treatment)[2:7]<-paste0("drug_",colnames(drug_treatment)[2:7])

#clrp_nejm <- c("#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF", "#7876B1FF", "#6F99ADFF", "#FFDC91FF", "#EE4C97FF","#b15928","#ffff33")
#"#D16014"

df_merge<-merge(Taf5l_KO,drug_treatment,by="V1")
p1<-ggplot(data = df_merge, mapping = aes(x = Taf5l_KO_stat, y = drug_stat)) +
  geom_point(shape = 21, fill ="#20854EFF", color = "white", size = 2) +
  geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) +
  stat_cor(
    method = "spearman", color = "#A63A50",size=5,
    label.x = -40, label.y = (20) * 0.95,
    cor.coef.name="Spearman: R",
    label.sep = "; p = "
  ) +
  #  stat_cor(
  #    method = "pearson", color = "#4C230A",size=5,
  #    label.x = -50, label.y = (30) * 0.85,
  #    cor.coef.name="Pearson: R",
  #    label.sep = "; p = "
  #  ) +
  labs(x = "Taf5l_KO (stat value)", y = "Drug_Treatment (stat value)")+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

ggsave(filename = "/Users/dwang20/Downloads/correlation_Taf5l_Drug_Treatment.pdf", plot = p1, width = 5.0, height = 4)



