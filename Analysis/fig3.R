library(dplyr)
library(data.table)
library(clusterProfiler)
library(ggplot2)

#########------------------------------------------------------------------------------------------------#########
#########-------------------------------------complex heatmap for RNA-seq data---------------------------#########
#########------------------------------------------------------------------------------------------------#########
OVCAR3_all_tpm<-fread("/Users/dwang20/OneDrive - Inside MD Anderson/Hok_rna_seq/OVCAR3_03_Finalcounts/OVCAR3_all_tpm_matrix.csv")%>%as.data.frame()
rownames(OVCAR3_all_tpm)<-OVCAR3_all_tpm$V1
df_need<-OVCAR3_all_tpm%>% select(control_1,control_2,TADA2B_1,TADA2B_2,TAF5L_1,TAF5L_2,TADA1_1,TADA1_2)
df_need<-log2(df_need+1)


TADA1_KO_vs_control<-fread("/Users/dwang20/OneDrive - Inside MD Anderson/Hok_rna_seq/OVCAR3_DE/TADA1_KO_vs_control.csv")%>%as.data.frame()%>% filter(log2FoldChange>0 & padj<0.05)
TAF5L_KO_vs_control<-fread("/Users/dwang20/OneDrive - Inside MD Anderson/Hok_rna_seq/OVCAR3_DE/TAF5L_KO_vs_control.csv")%>%as.data.frame()%>% filter(log2FoldChange>0 & padj<0.05)
TADA2B_KO_vs_control<-fread("/Users/dwang20/OneDrive - Inside MD Anderson/Hok_rna_seq/OVCAR3_DE/TADA2B_KO_vs_control.csv")%>%as.data.frame()%>% filter(log2FoldChange>0 & padj<0.05)

#all_genes<-unique(c(TADA1_KO_vs_control$V1,TAF5L_KO_vs_control$V1,TADA2B_KO_vs_control$V1))
all_genes<-intersect(intersect(TADA1_KO_vs_control$V1,TAF5L_KO_vs_control$V1),TADA2B_KO_vs_control$V1)

APM_genes<-c("TAP1","TAP2","PSMB8","PSMB9","PSMB10", "NLRC5","TAPBP","ERAP1","ERAP2","B2M","HLA-A","HLA-B","HLA-C")
common_APM_genes<-intersect(all_genes,APM_genes)

HALLMARK_INFLAMMATORY_RESPONSE<-read.gmt("/Users/dwang20/OneDrive - Inside MD Anderson/Hallmark_geneSet/human_hallmark_geneSet/HALLMARK_INFLAMMATORY_RESPONSE.v2025.1.Hs.gmt")
HALLMARK_TNFA_SIGNALING_VIA_NFKB<-read.gmt("/Users/dwang20/OneDrive - Inside MD Anderson/Hallmark_geneSet/human_hallmark_geneSet/HALLMARK_TNFA_SIGNALING_VIA_NFKB.v2025.1.Hs.gmt")
inflammatory_gene_list<-unique(c(HALLMARK_INFLAMMATORY_RESPONSE$gene,HALLMARK_TNFA_SIGNALING_VIA_NFKB$gene))
common_inflammatory_gene<-intersect(all_genes,inflammatory_gene_list)


HALLMARK_INTERFERON_ALPHA<-read.gmt("/Users/dwang20/OneDrive - Inside MD Anderson/Hallmark_geneSet/human_hallmark_geneSet/HALLMARK_INTERFERON_ALPHA_RESPONSE.v2025.1.Hs.gmt")
HALLMARK_INTERFERON_GAMMA<-read.gmt("/Users/dwang20/OneDrive - Inside MD Anderson/Hallmark_geneSet/human_hallmark_geneSet/HALLMARK_INTERFERON_GAMMA_RESPONSE.v2025.1.Hs.gmt")
INTERFERON_gene_list<-c(HALLMARK_INTERFERON_ALPHA$gene,HALLMARK_INTERFERON_GAMMA$gene)
common_interferon_gene<-intersect(all_genes,INTERFERON_gene_list)

needed_genes<-unique(c(common_APM_genes,common_inflammatory_gene,common_interferon_gene))

df<-df_need[rownames(df_need) %in% needed_genes,]
df2 <- df[apply(df, 1, sd, na.rm = TRUE) > 0, , drop = FALSE]
mat_row_z <- t(scale(t(df2)))

mat_row_z <- t(scale(t(df2))) 

##order the genes in rownames
mat_row_z2 <- mat_row_z[needed_genes[needed_genes %in% rownames(mat_row_z)], , drop = FALSE]

i <- which(genes == A)[1]
j <- which(genes == B)[1]
num_between <- abs(i - j) - 1   
span_len    <- abs(i - j) + 1

column_split = rep(1:4, each = 2,labels = rep("", 4))
row_split=rep(1:3,c(length(common_APM_genes),
                    length(which(needed_genes == common_inflammatory_gene[1])[1]:which(needed_genes == common_inflammatory_gene[length(common_inflammatory_gene)])[1]),
                    length(which(needed_genes == common_interferon_gene[1])[1]:which(needed_genes == common_interferon_gene[length(common_interferon_gene)])[1])
))

column_ha = HeatmapAnnotation(
  empty = anno_empty(border = FALSE),
  foo = anno_block(gp = gpar(fill = c("#0073C2FF", "#CD534CFF", "#EFC000FF","#868686FF")), 
                   labels = c("Control","TADA2B-KO","TAF5L-KO","TADA1-KO"),
                   labels_gp = gpar(col = "white", fontsize = 14))
)

gene_list = list(
  text1 = common_APM_genes,
  text2 = needed_genes[which(needed_genes == common_inflammatory_gene[1])[1]:which(needed_genes == common_inflammatory_gene[length(common_inflammatory_gene)])[1]],
  text3 = needed_genes[which(needed_genes == common_interferon_gene[1])[1]:which(needed_genes == common_interferon_gene[length(common_interferon_gene)])[1]]
)
# note how we set the width of this empty annotation
row_ha = rowAnnotation(foo = anno_empty(border = FALSE, 
                                        width = max_text_width(unlist(gene_list)) + unit(4, "mm")))

genes_to_show<-c("HLA-A","HLA-C","HLA-B","B2M","NLRC5","NFKB1","IFIT2","NFKB2","IRF7","TNFSF10","IL6","TNFSF9",
                 "NFKBIA","RELB","STAT1","CXCL8","CXCL1","CXCL2","CXCL10","RELA","STAT1","IFIT3","IFIT1","IFITM2","IFITM3")
right_annotation_ha  = rowAnnotation(foo = anno_mark(at = which(rownames(mat_row_z2) %in% genes_to_show),
                                                     labels = rownames(mat_row_z2)[rownames(mat_row_z2)%in%genes_to_show]))


#clrp_nejm <- c("#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF", "#7876B1FF", "#6F99ADFF", "#FFDC91FF", "#EE4C97FF","#b15928","#ffff33")
#clrp_jco <- c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "#7AA6DCFF", "#003C67FF", "#8F7700FF", "#3B3B3BFF", "#A73030FF", "#4A6990FF")
#clrp_npg <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF")

pdf("/Users/dwang20/Downloads/Heatmap3.pdf",height = 10, width = 7)
Heatmap(mat_row_z2, 
        name = "Z-score", 
        col = circlize::colorRamp2(c(-2, 0, 2), c("#1D3461", "white", "#DF2935")),
        #row_order = order(needed_genes,rownames(mat_row_z)),
        column_split = column_split,
        row_split = row_split,
        top_annotation = column_ha,
        left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c("#b15928","#7876B1FF", "#6F99ADFF")),
                                                         labels = c("APM", "Inflammatory Response Genes", "Interferon-Stimulated Genes"), 
                                                         labels_gp = gpar(col = "white", fontsize = 14))),
        right_annotation = right_annotation_ha,
        na_col = "grey",
        cluster_rows = F,
        cluster_columns = F,
        #top_annotation = column_ha,
        #right_annotation = ha,
        #column_split = Sample_metadata$Category,
        #row_names_max_width = unit(10,"in"),
        row_title_gp = gpar(fontsize = 0),
        column_title_gp = gpar(fontsize = 0),
        #column_names_gp = gpar(fontsize =10),
        #row_names_gp = gpar(fontsize =3),
        show_column_names = F,
        show_row_names = F)
dev.off()


#########------------------------------------------------------------------------------------------------#########
#########-------------------------------saga signature based correlation analysis------------------------#########
#########------------------------------------------------------------------------------------------------#########

all_genes<-fread("/Users/dwang20/OneDrive - Inside MD Anderson/Signature_genes/TADA2B_common_signature_for_mmu_hsa.csv")
colnames(all_genes)[2]<-"gene_name"

path2<-"/Users/dwang20/Downloads/data/TCGA/GDC/"
all_files2<-list.files(path2,full.names = T)
#ff<-all_files2[26] ## OV cancer cohort

all_TCGA<-lapply(all_files2, function(ff){
  
  print(ff)
  cancer_type<-gsub("_GDC_log2TPMplus1_202306_expr_meta.rds","",basename(ff))
  df_file<-readRDS(ff)
  exp<-as.data.frame(t(df_file@assayData$exprs))
  
  common_genes <- intersect(all_genes$gene_name, colnames(exp))
  print(length(common_genes))
  exp_need<-exp[,colnames(exp) %in% common_genes]%>% dplyr::select(sort(colnames(.)))  #322
  all_genes_sub <- all_genes[all_genes$gene_name %in% common_genes, ]%>% arrange(gene_name)
  all_genes_sub<-all_genes_sub[!duplicated(all_genes_sub$gene_name), ]
  weighted_exp <- sweep(exp_need, 2, all_genes_sub$weights_second_method, `*`)
  
  #as.vector(rowSums(weighted_exp)/(length(common_genes)^0.5))
  exp_need$Weighted_Sum <- as.vector(rowSums(weighted_exp)/(length(common_genes)^0.5))
  exp_specific<-exp%>% dplyr::select(HLA.A, HLA.B, HLA.C, B2M)
  #exp_specific<-exp%>% dplyr::select(CD8A, CD8B, GZMA, GZMB, PRF1)
  print(ncol(exp_specific))
  exp_specific$sample_id<-rownames(exp_specific)
  
  exp_need<-exp_need%>% dplyr::select(Weighted_Sum)
  exp_need$sample_id<-rownames(exp_need)
  
  meta<-df_file@phenoData@data
  meta_tumor<-meta[meta$sample_type=="Tumor",]
  
  merge_df<-merge(exp_specific,exp_need,by="sample_id")
  merge_df$Avg_MHC <- rowMeans(merge_df[,2:length(exp_specific)], na.rm = TRUE)
  
  merge_df<-merge_df[merge_df$sample_id %in% rownames(meta_tumor),]
  num<-nrow(merge_df)
  
  corr<-corr.test(merge_df$Avg_MHC,merge_df$Weighted_Sum,method="spearman",adjust="fdr")
  corr_r<-corr$r
  corr_p<-corr$p.adj
  
  res<-c(corr_r,corr_p,cancer_type,num)
  names(res)<-c("corr_r","corr_p","cancer_type","num")
  return(res)
})



##Scatter plot for corr
library(smplot2)
scale_to_minus1_1 <- function(x) {
  xmin <- min(x, na.rm = TRUE)
  xmax <- max(x, na.rm = TRUE)
  
  # avoid divide-by-zero if all values are identical
  if (xmin == xmax) {
    return(rep(0, length(x)))  # or x*0 + 0
  }
  
  x_scaled <- (x - xmin) / (xmax - xmin) * 2 - 1
  return(x_scaled)
}

merge_df$Weighted_Sum_scaled <- scale_to_minus1_1(merge_df$Weighted_Sum)

p<-ggplot(data = merge_df, mapping = aes(x = Weighted_Sum_scaled, y = Avg_MHC)) +
  geom_point(shape = 21, fill = "#EFC000FF", color = "white", size = 2) +
  sm_statCorr(
    color = "black", corr_method = "spearman",
    linetype = "dashed"
  )+
  labs(x='TADA2B-KO signature score',y='MHC-I avg expression \n(log2(TPM+1)')+
  #theme_bw() +
  theme(
    panel.border = element_rect(
      color = "black",
      fill  = NA,
      linewidth = 1.0
    )
  )

ggsave(p,file="/Users/dwang20/Downloads/TADA2B_KO_signature_score_MHC-I_avg_expression.pdf",width=4,height=3)


#CTL signature score
#"#F39B7FB2""#7E6148B2""#A690A4""#a6a57a""#7fb7be""#808080"
#"#CD534CFF""#868686FF""#EFC000FF"
merge_df$Weighted_Sum_scaled <- scale_to_minus1_1(merge_df$Weighted_Sum)
p1<-ggplot(data = merge_df, mapping = aes(x = Weighted_Sum_scaled, y = Avg_MHC)) +
  geom_point(shape = 21, fill = "#CD534CFF", color = "white", size = 2) +
  sm_statCorr(
    color = "black", corr_method = "spearman",
    linetype = "dashed"
  )+
  labs(x='TADA2B-KO signature score',y='CTL signature score')+
  #theme_bw() +
  theme(
    panel.border = element_rect(
      color = "black",
      fill  = NA,
      linewidth = 1.0
    )
  )

ggsave(p1,file="/Users/dwang20/Downloads/TADA2B_KO_signature_score_CTL.pdf",width=3.7,height=3)



#CTL signature score
res_df <- as.data.frame(t(as.data.frame(all_TCGA, check.names = FALSE)))
#res_df <-as.data.frame(res_df,stringsAsFactors=F)
res_df$corr_p.adj <- p.adjust(res_df$corr_p, method = "fdr")
rownames(res_df)<-paste0(res_df$cancer_type,"(",res_df$num,")")
res_df$corr_r<-as.numeric(res_df$corr_r)
res_df$cancer_type<-paste0(res_df$cancer_type," (num=",res_df$num,")")
res_df_order<-res_df[order(-res_df$corr_r), ]

TADA2B_corr<-res_df_order
TADA2B_corr$gene_name<-"TADA2B"

TADA1_corr<-res_df_order
TADA1_corr$gene_name<-"TADA1"

TAF5L_corr<-res_df_order
TAF5L_corr$gene_name<-"TAF5L"

df_all<-rbind(rbind(TADA2B_corr,TADA1_corr),TAF5L_corr)
write.csv(df_all,file="/Users/dwang20/Downloads/df_all_saga_TCGA_corr_with_MHCI.csv")



df_all_heatmap <- df_all %>% 
  mutate(significance = ifelse(corr_p.adj > 0.05| is.na(corr_p.adj), "FDR > 0.05", "FDR <= 0.05")) %>%
  mutate(Spearman_corr = as.numeric(corr_r))%>%
  mutate(cancer_type = gsub("num=","",cancer_type))

#tmp <- df_all_heatmap %>% 
#  filter(cancer_type != "TCGA-OV (429)")

df_all_heatmap$cancer_type<-factor(df_all_heatmap$cancer_type,levels=unique(df_all_heatmap$cancer_type))
df_all_heatmap$gene_name<-factor(df_all_heatmap$gene_name,levels=c("TAF5L","TADA1","TADA2B"))


df_all_heatmap$star <- with(df_all_heatmap, ifelse(
  corr_p.adj < 0.001, "***",
  ifelse(corr_p.adj < 0.01, "**",
         ifelse(corr_p.adj < 0.05, "*", ""))
))

ggplot(df_all_heatmap, aes_string(x = x_axis, y = y_axis)) +
  geom_tile(aes_string(fill = Spearman_corr),
            width = 0.95, height = 0.95) +
  geom_text(
    aes(label = star),
    color = "#474A48",
    size = 4,
    angle = 90,
    vjust=1
  ) +
  scale_fill_gradientn(
    colours = c("navy", "white", "firebrick3"),
    limits = c(-1, 1),
    breaks = c(-1, 0, 1),
    na.value = "white"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_bw() +
  theme(
    axis.text.x = ggtext::element_markdown(
      size = 12, angle = 90, hjust = 1, vjust = 0.5,color="black"
    ),
    axis.text.y = element_text(size = 12,color="black"),
    legend.position = "right"
  ) +
  labs(x = "", y = "")+
  scale_x_discrete(labels = function(labels) {
    sapply(labels, function(label) {
      if (grepl("TCGA-OV", label)) {
        paste0('<span style="color:red;">', label, '</span>')
      } else {
        label
      }
    })
  })









