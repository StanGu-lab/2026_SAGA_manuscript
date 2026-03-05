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
#########-----------------Check the MHC-I gene expression level across the TCGA tumor--------------------#########
#########------------------------------------------------------------------------------------------------#########

path<-"/Users/dwang20/Downloads/data/TCGA/GDC"
all_path<-list.files(path,full.names = T)

#i<-all_path[1]
all_files<-lapply(all_path,function(i){
  
  print(i)
  expr_meta<-readRDS(i)
  exp<-expr_meta@assayData$exprs
  exp<-as.data.frame(t(exp))
  exp$sample_type<-expr_meta$sample_type
  exp_tumor<-exp[exp$sample_type=="Tumor",]
  exp_tumor<-exp_tumor[,-ncol(exp_tumor)]
  exp_tumor$Avg_HLA_B2M <- rowMeans(exp_tumor[, c("HLA.A","HLA.B","HLA.C","B2M")], na.rm = TRUE)
  
  
  tumor_name<-gsub("_GDC_log2TPMplus1_202306_expr_meta.rds","",basename(i))
  exp_tumor$tumor_type<-tumor_name
  exp_tumor<-dplyr::select(exp_tumor,Avg_HLA_B2M,tumor_type)
  return(exp_tumor)
  
})

all_files_rbind<-do.call(rbind,all_files)

#qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',][c(1,2,3,7,8),]
#col_vector <- unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))

###  =========== ggplot style  =========== 
ggstyle = theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           panel.border = element_rect(colour = "black", fill=NA, linewidth=1.2),
                           axis.text.x=element_text(angle = 70, vjust = 1, hjust = 1, size=12, colour = "black", face = "bold"),
                           axis.text.y=element_text(size=12, colour = "black", face = "bold"),
                           text=element_text(size=12, colour = "black", face = "bold"),
                           legend.text=element_text(),
                           legend.position = "none", 
                           plot.title = element_text(size = 14, hjust=0.5,vjust = 0.5, 
                                                     margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"))

col_vector<-c(rep("#999999",5),"#F0027F",rep("#999999",40))


df_sorted <- all_files_rbind %>%
  group_by(tumor_type) %>%
  mutate(median_HLA_B2M = median(Avg_HLA_B2M)) %>%
  arrange(median_HLA_B2M)

all_files_rbind$tumor_type<-factor(all_files_rbind$tumor_type,levels=unique(df_sorted$tumor_type))

ggboxplot(all_files_rbind, x = "tumor_type", y = "Avg_HLA_B2M",
          color = "tumor_type",
          add="jitter",palette = col_vector,
          add.params = list(size = 0.2)
)+
  theme(plot.title=element_text(hjust=0.5))+
  theme(axis.text.x = element_text(angle = 70, hjust = 1, vjust = 1))+
  ylab('Avg Log2(TPM+1) MHC-I Expression')+
  xlab("Tumor type")+
  #stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "wilcox.test")+
  ggstyle


#########------------------------------------------------------------------------------------------------#########
#########---------------Check the MHC-I gene expression level across the CPTAC tumor samples-------------#########
#########------------------------------------------------------------------------------------------------#########

gene_name<-c("HLA-A","HLA-B","HLA-C","B2M")
#gene_name<-"TADA2B"
convert_to_ensembel_ID<-mapIds(org.Hs.eg.db,
                               keys=gene_name, 
                               column="ENSEMBL",
                               keytype="SYMBOL",
                               multiVals="first")

convert_to_ensembel_ID<-as.data.frame(convert_to_ensembel_ID)
convert_to_ensembel_ID$gene_name<-rownames(convert_to_ensembel_ID)

path<-"/Users/dwang20/OneDrive - Inside MD Anderson/LinkedOmicsKB_CPTAC"
all_files<-list.files(path,full.names = T)
#all_files<-all_files[1:9]
#i<-all_files[1]
all_df<-lapply(all_files,function(i){
  
  sample_name<-basename(i)
  print(sample_name)
  df<-fread(paste0(i,"/",sample_name,"_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt"))%>% as.data.frame()
  df$idx <- sub("\\..*", "", df$idx)
  df_MHCI<-df[df$idx %in% convert_to_ensembel_ID$convert_to_ensembel_ID,]
  rownames(df_MHCI)<-df_MHCI$idx
  df_MHCI<-as.data.frame(t(df_MHCI[,-1]))
  df_MHCI$tumor_type<-sample_name
  return(df_MHCI)
  
})

all_files_rbind <- dplyr::bind_rows(all_df) 
all_files_rbind$mean_value<-rowMeans(all_files_rbind[,1:3],na.rm=TRUE)
#all_files_rbind<-do.call(rbind,all_df)

df_sorted <- all_files_rbind %>%
  group_by(tumor_type) %>%
  mutate(median_mhc= median(mean_value)) %>%
  arrange(median_mhc)

all_files_rbind$tumor_type<-factor(all_files_rbind$tumor_type,levels=unique(df_sorted$tumor_type))

col_vector<-c("#FF7F00",rep("#999999",9))
ggboxplot(all_files_rbind, x = "tumor_type", y = "mean_value",
          color = "tumor_type",
          add="jitter",palette = col_vector,
          add.params = list(size = 0.7)
)+
  theme(plot.title=element_text(hjust=0.5))+
  theme(axis.text.x = element_text(angle = 70, hjust = 1, vjust = 1))+
  ylab('MHC-I Protein abundance')+
  xlab("Tumor type")+
  #stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "wilcox.test")+
  ggstyle


#########------------------------------------------------------------------------------------------------#########
#########---------------------------Crispr screen analysis OVCAR3 (without IFNg)-------------------------#########
#########------------------------------------------------------------------------------------------------#########
path<-"/Users/dwang20/OneDrive - Inside MD Anderson/Hok_crispr_OV3_OV90/"
all_files<-list.files(path,full.names = T)
gdata<-ReadRRA(all_files[3])
gdata$Rank = rank(gdata$Score)
geneList= gdata$Score
names(geneList) = gdata$id
gene_label_up<-c("TAPBP","B2M","TAP1","CALR", "SPPL3", "HLA-B", "IRF2", "TAP2", "RFXAP", "PDIA3", "NFKB1", "RFXANK")
gene_label_down<-c("EED", "BCORL1", "PCGF1", "TRAF3", "SUZ12", "MEN1", "IRF2BP2", "UGCG", "DOT1L", "EZH2", "TMEM127")


all_label_genes<-c(gene_label_down,gene_label_up)
#c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF")
p<-RankView(geneList, top = 0, bottom = 0, 
            genelist = all_label_genes)
ggsave(p,file="/Users/dwang20/OneDrive - Inside MD Anderson/genes_rank_by_log2FC.pdf",width=6.0,height=4.5)


###############------------------------Distribution of sgRNA--------------------------###############
sdata<-ReadsgRRA(all_files[4])
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
#########----------------------Pathway analysis for OVCAR3 (without IFNg) Crispr screen------------------#########
#########------------------------------------------------------------------------------------------------#########
df<-fread("/Users/dwang20/OneDrive - Inside MD Anderson/Hok_crispr_OV3_OV90/OV3_no_IFNY_S3_6_S3_8_S3_7_S3_9.gene_summary.txt")
df<-df%>% dplyr::select(id,`neg|lfc`)
df<-distinct(df, id, .keep_all = TRUE)
df<-df[1:150,]
ranks<- deframe(df) 
enrich<- EnrichAnalyzer(ranks,keytype = "Symbol",type = "GOCC+GOBP",method = "ORT",organism = "hsa",pvalueCutoff = 0.1,
                        limit = c(10, 500),universe = NULL,filter = FALSE,gmtpath = NULL,verbose = TRUE)
re<-enrich@result
re$ID<-gsub(":","_",re$ID)
re$source<-sub("_.*", "", re$ID)
re<-re %>%dplyr::filter(p.adjust<0.1)

result <- re %>%
  arrange(p.adjust) %>%
  slice_head(n = 10)%>% ungroup()

#result$Description<-paste0("(",result$source,") ",result$Description)
#result<-distinct(result, Description, .keep_all = TRUE)
result$Description <- factor(result$Description, levels = result$Description[order(-log10(result$p.adjust))])

ggplot(data = result, aes(x = -log10(p.adjust), y = Description)) + 
  geom_bar(stat="identity",fill="#A3320B",width=0.8) +
  #scale_color_gradient(low = "blue", high = "red") +
  ylab("") + 
  xlab("-log10(p.adjust)")+
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 12),
    #axis.title.y = element_text(size = 14), 
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12))


#########------------------------------------------------------------------------------------------------#########
#########----------------------Saga complex in OVCAR3 (without IFNg) crispr screen ----------------------#########
#########------------------------------------------------------------------------------------------------#########
path<-"/Users/dwang20/OneDrive - Inside MD Anderson/Hok_crispr_OV3_OV90/"
all_files<-list.files(path,full.names = T)
top_neg_200<-fread(all_files[3])[1:200,]

#In_vitro_1_2_3_4<-"/Users/dwang20/OneDrive - Inside MD Anderson/2025_0214_CRISPR_Screen/In-vitro_1_2_3_4.gene_summary.txt"
OV3_no_IFNY<-all_files[3]
saga_complex_genes<-fread("/Users/dwang20/Downloads/saga_complex_genes_mouse_human.csv")
saga_genes<-intersect(top_neg_200$id,saga_complex_genes$V2)
gstable=read.table(all_files[3],header=T)

startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=saga_genes
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
#samplelabel='vitro_3,vitro_4_vs_In-vitro_1,In-vitro_2 neg.'

colors=c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
         "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
         "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec)

pdf(file='/Users/dwang20/OneDrive - Inside MD Anderson/RRA_score.pdf',width=4.5,height=4.5)
plotrankedvalues(pvec,targetgenelist,xlab='Gene Rank',ylab='MAGeCK RRA score',main=paste(""))
dev.off()

####Density and distribution of saga sgRNA
sdata<-ReadsgRRA("/Users/dwang20/OneDrive - Inside MD Anderson/Hok_crispr_OV3_OV90/OV3_no_IFNY_S3_6_S3_8_S3_7_S3_9.sgrna_summary.txt")
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
ggsave(final_plot,file="/Users/dwang20/OneDrive - Inside MD Anderson/OVCAR3_no_IFNg_saga_genes_distribution.pdf",width=3,height=3)

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



#########----------------------------------------------------------------------------------#########
#########----------------------HGS2 E:T=1:2 co-culture crispr screen ----------------------#########
#########----------------------------------------------------------------------------------#########
In_vitro_1_2_3_4<-"/Users/dwang20/OneDrive - Inside MD Anderson/20250214_co_culture_CRISPR_Screen/In-vitro_1_2_3_4.gene_summary.txt"
top_neg<-fread(In_vitro_1_2_3_4)[1:50,]
saga_complex_genes<-fread("/Users/dwang20/Downloads/saga_complex_genes_mouse_human.csv")
saga_genes<-intersect(top_neg$id,saga_complex_genes$V3)
all_label_genes<-saga_genes

sdata<-ReadsgRRA("/Users/dwang20/OneDrive - Inside MD Anderson/20250214_co_culture_CRISPR_Screen/2025_crispr_screen/In-vitro_1_2_3_4.sgrna_summary.txt")
sdata <- sdata %>%
  mutate(color = case_when(
    #Gene %in% gene_label_up ~ "#CD534CFF",
    Gene %in% saga_genes ~ "#008280FF",
    #Gene %in% gene_label_down ~ "#0073C2FF",
    TRUE ~ "gray80" #keep original value
  ))

colorList= sdata$color
names(colorList) = sdata$Gene
saga_complex_genes_rank<-sdata[sdata$Gene %in% saga_genes,]%>% arrange(desc(LFC))
all_label_genes_rank<-(c(rev(unique(saga_complex_genes_rank$Gene))))
density_values <- density(sdata$LFC, n = 200)
density_df <- data.frame(x = density_values$x, density = density_values$y)

final_plot <- sgRankView_density(sdata, top = 0, bottom = 0,
                                 gene = all_label_genes_rank) +
  theme(axis.text.y = element_text(
    colour = c(rep("#008280FF", length(all_label_genes)))
  ))

#print(final_plot)
ggsave(final_plot,file="/Users/dwang20/OneDrive - Inside MD Anderson/20250221_In-vitro_ET_ratio_1_to_2_saga_genes_distribution.pdf",width=3,height=1.8)

##sgRNA Density plot
p2<-ggplot(sdata, aes(x = LFC)) +
  geom_density(alpha = 1, fill="black") + 
  scale_fill_gradient2(low = "white", mid = "black", high = "white", midpoint = 0) +
  xlab("")+
  ylab("")+
  theme_bw() +
  scale_x_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 1))+
  #xlim(-8,8)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank()) 
ggsave(p2,file="/Users/dwang20/OneDrive - Inside MD Anderson/20250221_In-vitro_ET_ratio_1_to_2_saga_genes_density.pdf",width=2.5,height=1.2)

gstable=read.table(In_vitro_1_2_3_4,header=T)
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=saga_genes
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
#samplelabel='vitro_3,vitro_4_vs_In-vitro_1,In-vitro_2 neg.'

colors=c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
         "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
         "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec)

pdf(file='/Users/dwang20/OneDrive - Inside MD Anderson/SAGA_RRA_score.pdf',width=4.5,height=4.5)
plotrankedvalues(pvec,targetgenelist,xlab='Gene Rank',ylab='MAGeCK RRA score',main=paste(""))
dev.off()





