library(data.table)
library("AnnotationDbi")
library("org.Hs.eg.db")
library(ggpubr)
library(dplyr)
library(stringr)
library(MAGeCKFlute)
library(clusterProfiler)
library(ggplot2)
library(data.table)
library(dplyr)
library(ggrepel)
library(tibble)
source("/Users/dwang20/OneDrive - Inside MD Anderson/20260216_SAGA_materials/scripts_for_saga_paper/my_function.R") 

#########------------------------------------------------------------------------------------------------#########
#########---------------------------Crispr screen analysis OVCAR3 (with IFNg)----------------------------#########
#########------------------------------------------------------------------------------------------------#########
path<-"/Users/dwang20/OneDrive - Inside MD Anderson/Hok_crispr_OV3_OV90/"
all_files<-list.files(path,full.names = T)
gdata<-ReadRRA(all_files[1])
gdata$Rank = rank(gdata$Score)
geneList= gdata$Score
names(geneList) = gdata$id
gene_label_up<-c("IFNGR2", "IFNGR1", "RFXAP","B2M","STAT1","TAP1","TAPBP","RFX5","SPPL3","TAP2", "IRF1", "JAK1", "NLRC5", "HLA-B", "NFKB1")
gene_label_down<-c("EED", "STUB1", "TRAF3", "PCGF1", "B3GNT5", "EZH2", "BCORL1", "MEN1", "TMEM127", "SUSD6")


all_label_genes<-c(gene_label_down,gene_label_up)
#c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF")
p<-RankView(geneList, top = 0, bottom = 0, 
            genelist = all_label_genes)
ggsave(p,file="/Users/dwang20/OneDrive - Inside MD Anderson/genes_rank_by_log2FC.pdf",width=6.0,height=4.5)


###############------------------------Distribution of sgRNA--------------------------###############
sdata<-ReadsgRRA(all_files[2])
#sdata$group = "no"
sdata <- sdata %>%
  mutate(color = case_when(
    Gene %in% gene_label_up ~ "#CD534CFF",
    #Gene %in% saga_complex_genes ~ "#EFC000FF",
    Gene %in% gene_label_down ~ "#0073C2FF",
    TRUE ~ "gray80" #keep original value
  ))

colorList= sdata$color
names(colorList) = sdata$Gene
gene_label_up_rank<-sdata[sdata$Gene %in% gene_label_up,]%>% arrange(desc(LFC))
gene_label_down_rank<-sdata[sdata$Gene %in% gene_label_down,]%>% arrange(desc(LFC))
#saga_complex_genes_rank<-sdata[sdata$Gene %in% saga_complex_genes,]%>% arrange(desc(LFC))

all_label_genes_rank<-(c(rev(unique(gene_label_down_rank$Gene)),rev(unique(gene_label_up_rank$Gene))))

# Compute density values for LFC
density_values <- density(sdata$LFC, n = 200)
density_df <- data.frame(x = density_values$x, density = density_values$y)

final_plot <- sgRankView_density(sdata, top = 0, bottom = 0,
                                 gene = all_label_genes_rank) +
  theme(axis.text.y = element_text(
    colour = c(rep("#0073C2FF", length(gene_label_down)), rep("#CD534CFF", length(gene_label_up)))
  ))

ggsave(final_plot,file="/Users/dwang20/OneDrive - Inside MD Anderson/pos_neg_distribution.pdf",width=3,height=5.5)

###Density plot
p2<-ggplot(sdata, aes(x = LFC)) +
  geom_density(alpha = 1, fill="black") + 
  scale_fill_gradient2(low = "white", mid = "black", high = "white", midpoint = 0) +
  xlab("")+
  ylab("")+
  theme_bw() +
  scale_x_continuous(limits = c(-8, 8), breaks = seq(-8, 8, by = 4))+
  #xlim(-8,8)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank()) 

ggsave(p2,file="/Users/dwang20/OneDrive - Inside MD Anderson/density.pdf",width=1.5,height=1.2)


#########------------------------------------------------------------------------------------------------#########
#########---------------------------Saga complex in OVCAR3 (with IFNg) crispr screen---------------------#########
#########------------------------------------------------------------------------------------------------#########
path<-"/Users/dwang20/OneDrive - Inside MD Anderson/Hok_crispr_OV3_OV90/"
all_files<-list.files(path,full.names = T)
top_neg_200<-fread(all_files[1])[1:200,]

#In_vitro_1_2_3_4<-"/Users/dwang20/OneDrive - Inside MD Anderson/2025_0214_CRISPR_Screen/In-vitro_1_2_3_4.gene_summary.txt"
OV3_IFNY<-all_files[1]
saga_complex_genes<-fread("/Users/dwang20/Downloads/saga_complex_genes_mouse_human.csv")
saga_genes<-intersect(top_neg_200$id,saga_complex_genes$V2)

gdata<-ReadRRA(all_files[1])
gdata$Rank = rank(gdata$Score)
geneList= gdata$Score
names(geneList) = gdata$id

all_label_genes<-saga_genes
#c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF")

p<-RankView_saga(geneList, top = 0, bottom = 0, 
            genelist = all_label_genes)
ggsave(p,file="/Users/dwang20/OneDrive - Inside MD Anderson/saga_complex_rank_by_log2FC.pdf",width=6.0,height=4.2)



####Density and distribution of saga sgRNA
sdata<-ReadsgRRA("/Users/dwang20/OneDrive - Inside MD Anderson/Hok_crispr_OV3_OV90/OV3_IFNY_S3_11_S3_13_S3_12_S3_14.sgrna_summary.txt")
#sdata$group = "no"
sdata <- sdata %>%
  mutate(color = case_when(
    #Gene %in% gene_label_up ~ "#CD534CFF",
    Gene %in% saga_genes ~ "#008280FF",
    #Gene %in% gene_label_down ~ "#0073C2FF",
    TRUE ~ "gray80" #keep original value
  ))

colorList= sdata$color
names(colorList) = sdata$Gene
#gene_label_up_rank<-sdata[sdata$Gene %in% gene_label_up,]%>% arrange(desc(LFC))
#gene_label_down_rank<-sdata[sdata$Gene %in% gene_label_down,]%>% arrange(desc(LFC))
saga_complex_genes_rank<-sdata[sdata$Gene %in% saga_genes,]%>% arrange(desc(LFC))

all_label_genes_rank<-(c(rev(unique(saga_complex_genes_rank$Gene))))

# Compute density values for LFC
density_values <- density(sdata$LFC, n = 200)
density_df <- data.frame(x = density_values$x, density = density_values$y)
final_plot <- sgRankView_density(sdata, top = 0, bottom = 0,
                                 gene = all_label_genes_rank) +
  theme(axis.text.y = element_text(
    colour = c(rep("#008280FF", length(all_label_genes)))
  ))
#ggsave(final_plot,file="/Users/dwang20/OneDrive - Inside MD Anderson/20250221_In-vitro_ET_ratio_1_to_2_saga_genes_distribution.pdf",width=3,height=1.8)
ggsave(final_plot,file="/Users/dwang20/OneDrive - Inside MD Anderson/OVCAR3_IFNg_saga_genes_distribution.pdf",width=3,height=3)

p2<-ggplot(sdata, aes(x = LFC)) +
  geom_density(alpha = 1, fill="black") + 
  scale_fill_gradient2(low = "white", mid = "black", high = "white", midpoint = 0) +
  xlab("")+
  ylab("")+
  theme_bw() +
  scale_x_continuous(limits = c(-8, 8), breaks = seq(-8, 8, by = 4))+
  #xlim(-8,8)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank()) 
ggsave(p2,file="/Users/dwang20/OneDrive - Inside MD Anderson/density.pdf",width=2.5,height=1.2)







