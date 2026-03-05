library(ggplot2)
library(data.table)

df<-fread("/Users/dwang20/OneDrive - Inside MD Anderson/20251118_Nearby_peaks_around_APM_genes/ATAC_OVCAR3_df_overlapping_sig_peaks_with_MHC_I_genes_cutoff_logfc_plus0.bed_result.csv")
count <- aggregate(GIGGLE_score ~ Factor, df, max)
count <- count[order(count$GIGGLE_score,decreasing = T),]
cr <- 20

tissue <- merge(count[c(1:20),],df,by.x=c("Factor","GIGGLE_score"),by.y=c("Factor","GIGGLE_score"))
#  print(tissue)


data <- df[df$Factor %in% as.character(count[c(1:cr),]$Factor),]
data$Factor<-factor(data$Factor,levels=count$Factor[1:20])
#data$Factor <- Factor(data$Factor, levels = as.character(count[c(1:cr),]$Factor))

#t_count <- as.data.frame(table(data$tissue))
#t_count <- t_count[order(t_count$Freq,decreasing = T),]
#ctr <- nrow(t_count)
#if (ctr > 10) {ctr=10}

#data$tissue <- as.character(data$tissue)

data <- data %>% 
  mutate(Factor_col = if_else(Factor %in% c("IRF1", "RELA"), Factor, "other"))

#data$label_colors <- ifelse(data$Factor %in% c("IRF1", "RELA"), "#D64933", "black")

ggstyle = theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           panel.border = element_rect(colour = "black", fill=NA, linewidth=1.0),
                           axis.text.x=element_text(size=12, colour = "black",angle = 60, hjust = 1),
                           axis.text.y=element_text(size=12, colour = "black"),
                           text=element_text(size=12, colour = "black"),
                           legend.text=element_text(),
                           legend.position = "none", 
                           plot.title = element_text(size = 14, hjust=0.5,vjust = 0.5, 
                                                     margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"))


ggplot(data, aes(x=Factor, y=GIGGLE_score,color=Factor_col)) +
  geom_point(size=2) +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  +
  theme_classic()+
  ggstyle+
  scale_color_manual(
    values = c(
      "other" = "#D7CDCC",   
      "IRF1"   = "#BC2C1A",      
      "RELA"  = "#9C528B"      
    ),
    guide = "none"          
  )+
  xlab("")+
  ylab("GIGGLE Score (Similarity)")