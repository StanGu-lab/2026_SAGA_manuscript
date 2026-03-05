##
library(Seurat) 
library(Matrix)
library(data.table)
library(progeny)
library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyr)
library(readr)
library(pheatmap)
library(tibble)


#SC564_1	Tada2b WT
#SC564_2	Tada2b KO

#########------------------------------------------------------------------------------------------------#########
#########--------------------------------Standard Single-Cell Processing Workflow------------------------#########
#########------------------------------------------------------------------------------------------------#########

path_control<-"/Users/dwang20/OneDrive - Inside MD Anderson/20250930_Tada2b_single_cell/SC564_1_outs/raw_feature_bc_matrix/"
expression_control_matrix <- Read10X(data.dir = path_control)
mouse_control <-CreateSeuratObject(counts = expression_control_matrix, 
                                   assay="RNA",
                                   project = "mouse_control",min.cells = 3, min.features = 200)

path_KO<-"/Users/dwang20/OneDrive - Inside MD Anderson/20250930_Tada2b_single_cell/SC564_2_outs/raw_feature_bc_matrix/"
expression_KO_matrix <- Read10X(data.dir = path_KO)
mouse_KO <-CreateSeuratObject(counts = expression_KO_matrix, 
                              assay="RNA",
                              project = "mouse_KO",min.cells = 3, min.features = 200)

mouse.combined <- merge(mouse_control, y = mouse_KO, add.cell.ids = c("Control", "Tada2b_KO"), project = "Tada2b_project")
table(mouse.combined$orig.ident)

mouse.combined[["percent.mt"]] <- PercentageFeatureSet(mouse.combined, pattern = "^mt-")
VlnPlot(mouse.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#mouse.combined <- subset(mouse.combined, subset = percent.mt < 25)
mouse.combined <- subset(mouse.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
mouse.combined <- NormalizeData(mouse.combined, normalization.method = "LogNormalize", scale.factor = 10000)
mouse.combined <- FindVariableFeatures(mouse.combined, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(mouse.combined)
mouse.combined <- ScaleData(mouse.combined, features = all.genes)
mouse.combined <- RunPCA(mouse.combined, features = VariableFeatures(object = mouse.combined))
ElbowPlot(mouse.combined)

mouse.combined <- FindNeighbors(mouse.combined, dims = 1:20)
mouse.combined <- FindClusters(mouse.combined, resolution = 0.6)

mouse.combined <- RunUMAP(mouse.combined, dims = 1:20)

mouse.combined <- JoinLayers(mouse.combined)

mouse.combined.markers <- FindAllMarkers(mouse.combined, only.pos = TRUE)
mouse.combined.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5
DoHeatmap(mouse.combined, features = top5$gene) + NoLegend()

#mouse.markers <- FindAllMarkers(mouse.combined, only.pos = TRUE)
#mouse.markers %>%
#  group_by(seurat_clusters) %>%
#  dplyr::filter(avg_log2FC > 0)

DimPlot(mouse.combined, reduction = "umap",label = TRUE, pt.size = 0.5,split.by = "orig.ident")
DimPlot(mouse.combined, reduction = "umap",label = TRUE, pt.size = 0.5,group.by = "orig.ident")

library(genesorteR)
gs = sortGenes(mouse.combined[["RNA"]]$data, Idents(mouse.combined))
gs_marker <- getMarkers(gs, quant = 0.99)
spec_scores <- gs$specScore
genes <- rownames(spec_scores)
top10_genes_df <- data.frame()
for (cluster in colnames(spec_scores)){
  cluster_scores <- data.frame(
    gene = genes,
    specScore = spec_scores[, cluster],
    cluster = cluster
  )
  top10_genes <- cluster_scores %>%
    arrange(desc(specScore)) %>%
    head(10)
  top10_genes_df <- rbind(top10_genes_df, top10_genes)
}

write.csv(top10_genes_df, "/Users/dwang20/OneDrive - Inside MD Anderson/20250930_Tada2b_single_cell/top10_genes_df.csv")


VlnPlot(mouse.combined, features = c("Cd8a", "Cd8b1", "Tcf7","Klrb1c","Ifng","Nkg7","Klrc1","Klrc2","Gzmb"))
VlnPlot(mouse.combined, features = c("Vegfa","Nlrp3","Lyve1","C1qc","Ccl4","Il1rn","Ereg","Il1b","Trem2"))
VlnPlot(mouse.combined, features = c("Top2a","Mki67","Cd3e","Cd4","Cd8a"))
VlnPlot(mouse.combined, features = c("Lpcat2", "Cx3cr1", "Lst1"))

gene.sets <- list()
gene.sets$pDC<-c("Cox6a2","Ly6d","Siglech","Ccr9","Atp1b1","Ly6c2","Ly6c1","Ctsl","Bst2","Ly6a","Cd7","Klk1","Tex2","Sell")
gene.sets$DC1<-c("Xcr1","Rab7b","Pianp","Cxx1a","Cxx1b","Ifi205","Mnda","Naaa","Sept3","Stmn1","Gcsam","Ppt1","Rab32")
gene.sets$DC2<-c("Mmp12","Wfdc17","Clec4n","Il1r2","Emb","Ms4a6d","Dab2","Mgl2","Ccr1","Ctsc","Il1b")
gene.sets$DC3<-c("Fscn1","Ccl5","Ccr7","Cacnb3","Anxa3","Tmem123","Fabp5","Epsti1","Gyg","Samsn1")
gene.sets$Mono1<-c("F13a1","Vcan","Plac8","Ifi27l2a","Rnf213","Ms4a6b","Hpse","Ifi47","Ly6c1","Ifi204")
gene.sets$Mono2<-c("Cd300e","Eno3","Spn","Treml4","Adgre4","Itgal","Ace","Hfe","Ddit4","Smpdl3b")
gene.sets$Mono3<-c("Retnlg","G0s2","Lcn2","Csf3r","Il1b","S100a9","Cxcr2","S100a8","Sike1","Mmp9")
gene.sets$Mac1<-c("Retnla","Ccl17","C1qa","Fcrls","Cxcl16","Mdfic","Ccr5","Ccl12","Mafb","Bcl2l11")
gene.sets$Mac2<-c("Gpnmb","Mfge8","Fabp5","Gm14328","Rilpl2","Spp1","Gm6100","Pdgfa","Lamtor3","Cd63-ps","Gm6977","Ank2")
gene.sets$Mac3<-c("Alox15","Saa3","Prg4","Serpinb2","C4b","F5","Slpi","Ptgis","Icam2","Cd5l","Ednrb","Fcna","Selp","Lrg1","Ecm1")
gene.sets$Mac4<-c("Krt79","Chil3","Plet1","Car4","Cidec","Tcf7l2","Krt19","Hsd3b7","Cpne5","Mfsd7c","Il18")
gene.sets$MonoDC<-c("Cd7","Cd209a","Atp1b1","Clec10a","Cd209d","Klrd1","Flt3","Dpp4","Kmo","Ybx3","Cd209e")
gene.sets$Endothelial<-c("Ptprb","Adgrf5","Egfl7","Esam")
gene.sets$cycling<-c("Top2a","Mki67")
gene.sets$CD8T<-c("Cd3e","Cd3d","Cd3g","Cd8a", "Cd8b1")
gene.sets$CD4T<-c("Cd3e","Cd3d","Cd3g","Cd4")
gene.sets$NK<-c("Nkg7","Prf1","Ncr1","Ctsw","Klrb1c")
gene.sets$Neutrophils<-c("S100a9","S100a8","Lrg1","Csf3r","Vasp")
gene.sets$B<-c("Ms4a1","Cd79a","Fcmr","Cd79b","Cd19")
gene.sets$Mast<-c("Kit","Tpsab1","Cpa3")

gene.sets$macro_TRM<-c("Lyve1","Folr2","Timd4","Gas6","Mertk")
gene.sets$mono_mac<-c("Ccr2","S100a8","S100a9","Il1b","Tnf")
gene.sets$lam_macro<-c("Trem2","Tyrobp","Apoe","Lpl","Ctsb","Ctsd","Lgals3")


mouse.combined_anno <- AddModuleScore(mouse.combined, features = gene.sets,
                                      name = c("pDC", "DC1", "DC2","DC3",
                                               "Mono1","Mono2","Mono3",
                                               "Mac1","Mac2","Mac3","Mac4",
                                               "MonoDC","Endothelial","cycling","CD8T",
                                               "CD4T","NK","Neutrophils","B","Mast",
                                               "macro_TRM","mono_mac","lam_macro"
                                      ))

meta_add_module<-mouse.combined_anno@meta.data[,c(6,8:30)]
module_need<-colnames(mouse.combined_anno@meta.data)[8:30]
result <- meta_add_module %>%
  group_by(seurat_clusters) %>%
  summarise(
    across(all_of(module_need), ~ mean(.x, na.rm = TRUE), .names = "{.col}_mean")
  )

write.csv(result,file="/Users/dwang20/OneDrive - Inside MD Anderson/20250930_Tada2b_single_cell/20251001_module_score_0_6.csv")

##annotate cell types

lab <- c(
  "0"="Macrophage",
  "1"="Macrophage",
  "2"="Macrophage",
  "3"="CD8+ T",
  "4"="Cycling myeloid",
  "5"="B cells",
  "6"="CD4+ T",
  "7"="Endothelial",
  "8"="CD4+ T",
  "9"="Macrophage",
  "10"="Cycling myeloid",
  "11"="NKT",
  "12"="DC",
  "13"="MonoDC",
  "14"="others",
  "15"="Cycling T",
  "16"="Neutrophils"
)

#"#197278","#C44536"
mycol = c("mouse_control"="#62929E", "mouse_KO"="#DC493A") 

mouse.combined@meta.data$annot <- lab[as.character(Idents(mouse.combined))]
DimPlot(mouse.combined, reduction = "umap",label = FALSE, pt.size = 0.1,group.by="orig.ident")+ 
  labs(title = NULL)+
  scale_color_manual(values = c(mycol))

#6E8BA4, #837BA9, #BFD68C, #B59FC7, #794966, #9398C8, #56579C, #8E619D, #F2EE96, #C05949, #C86F85, #D68759, #E7C56D

mycol2 = c("B cells"="#E7C56D","CD4+ T" = "#BFD68C","CD8+ T" = "#C05949","Cycling myeloid"="#B59FC7","Cycling T"="#F2EE96",
           "DC"="#6E8BA4","Endothelial"="#794966","Macrophage"="#837BA9","MonoDC"="#C86F85", "Neutrophils"= "#D5B9B2","NKT"="#B8D3D1", 
           "others"="#D68759")    

DimPlot(mouse.combined, reduction = "umap",label = TRUE, pt.size = 0.1,group.by = "annot",
        split.by ="orig.ident",label.size = 4.5,repel = TRUE)+ 
  labs(title = NULL)+
  scale_color_manual(values = c(mycol2))

DimPlot(mouse.combined, reduction = "umap",label = TRUE, pt.size = 0.1,split.by ="orig.ident",label.size = 4,repel = TRUE)


#########------------------------------------------------------------------------------------------------#########
#########------------------------------------------Single-cell Feature plot------------------------------#########
#########------------------------------------------------------------------------------------------------#########

features_to_plot<- c(
  "Igfbp7", "Col4a1","Col4a2", "Cavin3", "Nfib", "Serpinh1", #endothelial
  "Ms4a1","Cd79a","Fcmr","Cd79b","Cd19", #B
  "S100a9","S100a8","Lrg1","Csf3r", #Neutrophils
  "Cd3e","Cd3d","Cd3g","Cd8a", "Cd8b1",#CD8
  "Cd3e","Cd3d","Cd3g","Cd4",#CD4
  "Nkg7","Prf1","Ncr1","Klrb1c",#NK
  "Adgre1","Csf1r", "Mrc1", "Msr1", "Fcgr1", "C1qa","C1qb", #Macro
  "Fscn1","Ccr7","Cacnb3","Anxa3", ##DC
  "Top2a","Mki67", #cycling
  "Cd7","Cd209a","Clec10a","Cd209d"#mMonoDC
)

mouse.combined$annot<-factor(mouse.combined$annot,levels=c("B cells","CD4+ T","CD8+ T","Cycling myeloid","Cycling T",
                                                           "DC","Endothelial","Macrophage","MonoDC", "Neutrophils","NKT","others"))
Idents(mouse.combined) <- "annot"
mycol2 = c("B cells"="#E7C56D","CD4+ T" = "#BFD68C","CD8+ T" = "#C05949","Cycling myeloid"="#B59FC7","Cycling T"="#F2EE96",
           "DC"="#6E8BA4","Endothelial"="#794966","Macrophage"="#837BA9","MonoDC"="#C86F85", "Neutrophils"= "#D5B9B2","NKT"="#B8D3D1", 
           "others"="#D68759") 

pdf("/Users/dwang20/Downloads/marker_Clustered_DotPlot.pdf", width = 5.5, height = 9.5, useDingbats = FALSE)
scCustomize::Clustered_DotPlot(mouse.combined, features = unique(features_to_plot),
                               plot_km_elbow = FALSE,cluster_feature = FALSE,
                               cluster_ident = FALSE,
                               colors_use_idents=mycol2,
                               row_label_size = 12,column_label_size = 12,
                               legend_label_size = 12,legend_title_size = 12,
                               x_lab_rotate=90)
dev.off()


saveRDS(mouse.combined, file = "/Users/dwang20/OneDrive - Inside MD Anderson/20250930_Tada2b_single_cell/20251001_mouse.combined.rds")
mouse.combined<-readRDS("/Users/dwang20/OneDrive - Inside MD Anderson/20250930_Tada2b_single_cell/20251001_mouse.combined.rds")

#########------------------------------------------------------------------------------------------------#########
#########------------------------------Further analysis cycling T; NKT; CD8;CD4--------------------------#########
#########------------------------------------------------------------------------------------------------#########
Idents(mouse.combined)<-mouse.combined$annot
mouse.combined_T_cells<-subset(mouse.combined, idents = c('CD4+ T','CD8+ T','Cycling T','NKT'))
DefaultAssay(mouse.combined_T_cells) <- "RNA"
all.genes <- rownames(mouse.combined_T_cells)
mouse.combined_T_cells <- ScaleData(mouse.combined_T_cells, features = all.genes)

mouse.combined_T_cells <- RunPCA(mouse.combined_T_cells, features = VariableFeatures(object = mouse.combined_T_cells))
ElbowPlot(mouse.combined_T_cells)

mouse.combined_T_cells <- FindNeighbors(mouse.combined_T_cells, dims = 1:20)
mouse.combined_T_cells <- FindClusters(mouse.combined_T_cells, resolution = 0.8)
mouse.combined_T_cells <- RunUMAP(mouse.combined_T_cells, dims = 1:20)
mouse.combined_T_cells <- RunTSNE(mouse.combined_T_cells, dims = 1:20)

DimPlot(mouse.combined_T_cells, reduction = "umap",label = TRUE, pt.size = 0.5,split.by = "orig.ident")+ NoLegend()

library(genesorteR)
gs = sortGenes(mouse.combined_T_cells[["RNA"]]$data, Idents(mouse.combined_T_cells))
gs_marker <- getMarkers(gs, quant = 0.99)
spec_scores <- gs$specScore
genes <- rownames(spec_scores)
top10_genes_df <- data.frame()
for (cluster in colnames(spec_scores)){
  cluster_scores <- data.frame(
    gene = genes,
    specScore = spec_scores[, cluster],
    cluster = cluster
  )
  top10_genes <- cluster_scores %>%
    arrange(desc(specScore)) %>%
    head(10)
  top10_genes_df <- rbind(top10_genes_df, top10_genes)
}
#VlnPlot(HD_seurat, features = top10_genes_df$gene, stack = T) + theme(legend.position = 'none')
write.csv(top10_genes_df, "/Users/dwang20/OneDrive - Inside MD Anderson/20250930_Tada2b_single_cell/marker_genes_for_T_clusters_resolution0_8.csv")

VlnPlot(mouse.combined_T_cells, features = c("Il23r","Rorc","Il7r","Kit","Cd3d","Cd4","Cd8a","Trdc"))

lab <- c(
  "0"="CD8+ T",
  "1"="CD4+ T",
  "2"="CD8+ T",
  "3"="CD4+ T",
  "4"="CD4+ T",
  "5"="CD8+ T",
  "6"="CD4+ T",
  "7"="NK",
  "8"="CD8+ T",
  "9"="Cycling T",
  "10"="ILCs",
  "11"="others"
)

mouse.combined_T_cells@meta.data$annot_sub <- lab[as.character(Idents(mouse.combined_T_cells))]
meat_mouse.combined_T_cells<-mouse.combined_T_cells@meta.data
meat_mouse.combined_T_cells$sample_id<-rownames(meat_mouse.combined_T_cells)
mouse.combined@meta.data$sample_id<-rownames(mouse.combined@meta.data)

md <- mouse.combined@meta.data %>%
  left_join(meat_mouse.combined_T_cells[, c("sample_id","annot_sub")], by = "sample_id") %>%
  mutate(annot2 = coalesce(annot_sub,annot))%>%
  select(-annot_sub)

mouse.combined<-AddMetaData(mouse.combined, md, col.name = NULL)

mycol2 = c("B cells"="#E7C56D","CD4+ T" = "#BFD68C","CD8+ T" = "#C05949","Cycling myeloid"="#B59FC7","Cycling T"="#F2EE96",
           "DC"="#6E8BA4","Endothelial"="#794966","Macrophage"="#837BA9","MonoDC"="#C86F85", "Neutrophils"= "#D5B9B2",
           "NK"="#B8D3D1","ILCs"="#BFB5AF","others"="#D68759") 
DimPlot(mouse.combined, reduction = "umap",label = TRUE, pt.size = 0.1,group.by = "annot2",
        split.by ="orig.ident",label.size = 4.5,repel = TRUE)+ 
  labs(title = NULL)+
  scale_color_manual(values = c(mycol2))


features_to_plot<- c(
  "Igfbp7", "Col4a1","Col4a2", "Cavin3", "Nfib", "Serpinh1", #endothelial
  "Ms4a1","Cd79a","Fcmr","Cd79b","Cd19", #B
  "S100a9","S100a8","Lrg1","Csf3r", #Neutrophils
  "Cd3e","Cd3d","Cd3g","Cd8a", "Cd8b1",#CD8
  "Cd3e","Cd3d","Cd3g","Cd4",#CD4
  "Ncr1","Klre1","Klrb1c",#NK
  "Adgre1","Csf1r", "Mrc1", "Msr1", "Fcgr1", "C1qa","C1qb", #Macro
  "Fscn1","Ccr7","Cacnb3","Anxa3", ##DC
  "Top2a","Mki67", #cycling
  "Cd209a","Clec10a","Cd209d",#mMonoDC
  "Il23r","Rorc","Il1r1","Il7r" ##ILCs
)

mouse.combined$annot2<-factor(mouse.combined$annot2,levels=c("B cells","CD4+ T","CD8+ T","Cycling myeloid","Cycling T",
                                                             "DC","Endothelial","ILCs","Macrophage","MonoDC", "Neutrophils","NK","others"))
Idents(mouse.combined) <- "annot2"
mycol2 = c("B cells"="#E7C56D","CD4+ T" = "#BFD68C","CD8+ T" = "#C05949","Cycling myeloid"="#B59FC7","Cycling T"="#F2EE96",
           "DC"="#6E8BA4","Endothelial"="#794966","Macrophage"="#837BA9","MonoDC"="#C86F85", "Neutrophils"= "#D5B9B2",
           "NK"="#B8D3D1","ILCs"="#BFB5AF","others"="#D68759") 

pdf("/Users/dwang20/Downloads/marker_Clustered_DotPlot.pdf", width = 5.5, height = 10, useDingbats = FALSE)
scCustomize::Clustered_DotPlot(mouse.combined, features = unique(features_to_plot),
                               plot_km_elbow = FALSE,cluster_feature = FALSE,
                               cluster_ident = FALSE,
                               colors_use_idents=mycol2,
                               row_label_size = 12,column_label_size = 12,
                               legend_label_size = 12,legend_title_size = 12,
                               x_lab_rotate=90)
dev.off()
saveRDS(mouse.combined, file = "/Users/dwang20/OneDrive - Inside MD Anderson/20250930_Tada2b_single_cell/20251002_mouse.combined.rds")

#########------------------------------------------------------------------------------------------------#########
#########-----------------------------------------------Cell populations---------------------------------#########
#########------------------------------------------------------------------------------------------------#########

meta_mouse_combined<-mouse.combined@meta.data[!mouse.combined@meta.data$annot2 %in% c("Endothelial","others"),]
all_immune_cells<-meta_mouse_combined%>% group_by(orig.ident)%>%  summarise(count = n())
number_diff_cells<-meta_mouse_combined%>% group_by(orig.ident,annot2)%>%  summarise(count = n())

df_merge<-merge(number_diff_cells,all_immune_cells,by="orig.ident",all.x=T)
df_merge$pro<-(df_merge$count.x)/(df_merge$count.y)
write.csv(df_merge,file="/Users/dwang20/OneDrive - Inside MD Anderson/20250930_Tada2b_single_cell/cell_proportion.csv")
df_merge<-fread("/Users/dwang20/OneDrive - Inside MD Anderson/20250930_Tada2b_single_cell/cell_proportion.csv")

ggplot(df_merge, aes(x = orig.ident, y = pro)) +
  geom_col(aes(color = annot2, fill = annot2), position = position_stack()) +
  scale_color_manual(values = mycol2,name = "Cell type")+
  scale_fill_manual(values = mycol2,name = "Cell type")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1.2),
        axis.text.x=element_text(size=12, colour = "black"),
        axis.text.y=element_text(size=12, colour = "black"),
        text=element_text(size=12, colour = "black"),
        legend.text=element_text(size=12))+
  scale_x_discrete(labels = c(mouse_control="Control", mouse_KO="Tada2b_KO"))+
  xlab("")+
  ylab("Proportions")


#########------------------------------------------------------------------------------------------------#########
#########----------------------------------CD8+ T signature score analysis-------------------------------#########
#########------------------------------------------------------------------------------------------------#########
mouse.combined<-readRDS("/Users/dwang20/OneDrive - Inside MD Anderson/20250930_Tada2b_single_cell/20251002_mouse.combined.rds")
CD8_T_cytotoxic_state<-c("Ccl5", "Klrg1","Prf1","Fcgr4","Cx3cr1","Gzmk","Eomes","Nkg7")
CD8_T_dysfunction_state<-c("Lag3","Pdcd1","Layn","Havcr2","Ctla4")
CD8_T_IFNg<-c("Ifit3","Ifit2","Stat1","Irf7","Isg15","Ifitm3","Oas2","Jak2","Socs1","Mx1")
#CD8_T_proliferation<-c("Lif","Il2","Cenpv","Nme1","Fabp5","Orc6","G0s2","Gck")
T_naive<-c("Tcf7","Lef1","Ccr7","Sell","Mal")
CD8_T_proliferation <- c("Mki67","Top2a","Birc5","Ccnb1","Ccnb2","Cdc20","Cenpf","Nusap1","Tpx2","Prc1","Kif11","Aurkb","Plk1")
#naive_B <- c("Ms4a1","Cd19","Cd79a","Cd79b","Pax5","Ebf1","Ikzf3",
#                "Ighd","Ighm","Fcer2a","Cr2","Cxcr5","Bach2")
#test<-c("Bcl6", "Aicda", "Mef2b","Cxcr4","Cxcr5", "Irf8", "S1pr2", "Bcl2l11")

mouse.combined <- AddModuleScore(mouse.combined,
                                 features = list(T_naive),
                                 name="T_naive")

data_CD8<-mouse.combined@meta.data[grep("CD8",mouse.combined@meta.data$annot2),]
data_CD8$cell_group<-paste0(data_CD8$orig.ident,"_",data_CD8$annot2)
#data_CD8$subgroup<-factor(data_CD8$subgroup,levels=c("Control","ARID4B_KO"))
ggboxplot(data_CD8, x = "cell_group", y = "CD8_T_proliferation1",
          color = "cell_group",
          add = "jitter", add.params = list(size = 0.5))+
  scale_colour_manual(values = c("#BFD68C","#C05949"))+
  stat_compare_means(aes(group = cell_group),label = "p.signif", method = "wilcox.test",label.x = 1.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        #panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        axis.text.x=element_text(size=12, colour = "black", face = "bold"),
        axis.text.y=element_text(size=12, colour = "black", face = "bold"),
        text=element_text(size=12, colour = "black", face = "bold"),
        legend.text=element_text(size=12))+
  ylab("CD8_T_proliferation score")+
  xlab("")+
  scale_x_discrete(labels = c("mouse_control_CD8+ T"="Control", "mouse_KO_CD8+ T"="Tada2b_KO"))+
  theme(legend.position = "none")


