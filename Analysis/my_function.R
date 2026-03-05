RankView <- function(rankdata, genelist=NULL, top=10, bottom=10, cutoff=NULL,
                     main=NULL, filename=NULL, width=5, height=4, ...){
  requireNamespace("ggrepel", quietly=TRUE) || stop("need ggrepel package")
  if(length(cutoff)==0) cutoff = CutoffCalling(rankdata, 2)
  if(length(cutoff)==1) cutoff = sort(c(-cutoff, cutoff))
  data = data.frame(Gene = names(rankdata), diff = rankdata, stringsAsFactors=FALSE)
  data$Rank = rank(data$diff)
  data$group = "no"
  #data$group[data$diff>cutoff[2]] = "up"
  #data$group[data$diff<cutoff[1]] = "down"
  #data <- data %>%
  #  mutate(group = case_when(
  #    Gene %in% gene_label_up ~ "up",
  #    Gene %in% gene_label_down ~ "down",
  #    TRUE ~ group #keep original value
  #  ))
  
  data <- data %>%
    mutate(group = case_when(
      Gene %in% gene_label_up ~ "up",
      Gene %in% gene_label_down ~ "down",
      #Gene %in% saga_complex_genes ~ "saga",
      TRUE ~ group #keep original value
    ))
  #dot_size = ifelse(Gene %in% all_label_genes,0.3,0.3)) 
  
  idx=(data$Rank<=bottom) | (data$Rank>(max(data$Rank)-top)) | (data$Gene %in% genelist)
  #"saga"="#EFC000FF",  
  mycolour = c("up"="#CD534CFF","down"="#0073C2FF")
  p = ggplot(data)
  #p = p + geom_jitter(aes(x=diff,y=Rank,color=group),size=1.2)
  p = p + geom_jitter(
    data = data %>% filter(group == "no"),
    aes(x = diff, y = Rank, color = group), size = 0.5
  )
  
  p = p + geom_jitter(
    data = data %>% filter(group != "no"),
    aes(x = diff, y = Rank, color = group), size = 1.2
  )
  
  if(!all(cutoff==0)) p = p + geom_vline(xintercept = cutoff, linetype = "dotted")
  if(sum(idx)>0)
    p = p + geom_label_repel(aes(x=diff, y=Rank,fill=group,label = Gene),
                             data=data[idx,], fontface = 'bold', color = 'white', size = 4.5,
                             box.padding = unit(0.4, "lines"), segment.color = 'grey50',
                             point.padding = unit(0.3, "lines"), segment.size = 0.3,max.overlaps = Inf)
  p = p + scale_color_manual(values=mycolour)
  p = p + scale_fill_manual(values=mycolour)
  p = p + theme(panel.background = element_rect(fill='white', colour='black'))
  p = p + theme(text = element_text(colour="black",size = 14, family = "Helvetica"),
                plot.title = element_text(hjust = 0.5, size=18),
                axis.text = element_text(colour="gray10"))
  p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank())
  p = p + labs(x="Log2FC",y="Rank",title=main)
  p = p + theme(legend.position="none")#+ylim(-1000,7000)
  
  #if(!is.null(filename)){
  #  ggsave(plot=p, filename=filename, units = "in", width=width, height=height, ...)
  #}
  return(p)
}


sgRankView_density <- function(df, gene = NULL,
                               top = 3, bottom = 3,
                               binwidth = 0.3,
                               interval = 0.1){
  
  if (!all(c("sgrna", "Gene", "LFC") %in% colnames(df)))
    stop("Make sure your data contains columns of 'sgrna', 'Gene', and 'LFC' ...")
  
  df = as.data.frame(df, stringsAsFactors = FALSE)
  df = df[order(df$LFC), ]
  df$Gene = as.character(df$Gene)
  tmp = stats::aggregate(df$LFC, by = list(df$Gene), median)
  colnames(tmp) = c("Gene", "mid")
  tmp = tmp[order(tmp$mid), ]
  
  if(top > 0){
    idx = max((nrow(tmp) - top + 1), 1)
    gene = c(gene, tmp$Gene[idx:nrow(tmp)])
  }
  if(bottom > 0){
    gene = c(gene, tmp$Gene[1:min(bottom, nrow(tmp))])
  }
  gene = unique(gene)
  
  subdf = df[df$Gene %in% gene, ]
  if(nrow(subdf) < 2) return(ggplot())
  subdf$Gene = factor(subdf$Gene, levels = gene)
  subdf = subdf[order(subdf$Gene), ]
  subdf$index = rep(1:length(gene), as.numeric(table(subdf$Gene)[gene]))
  subdf$yend <- (binwidth + interval) * subdf$index - interval
  subdf$y <- (binwidth + interval) * (subdf$index - 1)
  color <- rep("pos", dim(subdf)[1])
  median_value<-subdf%>%group_by(Gene)%>% summarise(value=median(LFC))%>% filter(value<0)
  color[which(subdf[, "Gene"] %in% median_value$Gene)] <- "neg"
  
  #color[which(subdf[, 3] < 0)] <- "neg"
  subdf$color <- color
  subdf = subdf[, c("sgrna", "Gene", "LFC", "y", "yend", "color", "index")]
  
  # **Create Black Boxes for Each Gene**
  gene_boxes <- data.frame(y = subdf$y, yend = subdf$yend, Gene = subdf$Gene)
  
  # Merge Density Background into the Ranking View
  p <- ggplot()
  # **Add Density Background Directly Inside sgRankView**
  p <- p + geom_tile(data = density_df, aes(x = x, y = 0, fill = density), height = Inf, alpha = 0.8)
  # Add Black Boxes
  p <- p + geom_rect(data = gene_boxes, aes(xmin = -8, xmax = 8, ymin = y, ymax = yend), 
                     color = "black", fill = NA, size = 0.5)
  # Add Gene Ranking Segments
  p <- p + geom_segment(aes_string("LFC", "y", xend = "LFC", yend = "yend", color = "color"), data = subdf)
  
  # Use Proper Color and Fill Gradients
  p <- p + scale_color_manual(values = c("pos" = "#e41a1c", "neg" = "#377eb8"))
  p <- p + scale_fill_gradient(low = "white", high = "black")  # Density Gradient
  
  p <- p + scale_y_continuous(breaks = seq(min(subdf$y), max(subdf$y), by = binwidth + interval),
                              labels = gene, expand = c(0, 0))
  p <- p + scale_x_continuous(limits = c(-8, 8), breaks = seq(-4, 4, by = 4))
  
  p <- p + labs(x = "Log2(Fold change)", y = NULL)
  p <- p + theme(
    text = element_text(colour = "black", size = 14, family = "Helvetica"),
    plot.title = element_text(hjust = 0.5, size = 18),
    axis.text = element_text(colour = "gray10"),
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = "none"
  )
  
  return(p)
}


plotrankedvalues<-function(val, tglist, ...){
  
  #plot(val,log='y',ylim=c(max(val),min(val)),type='p',lwd=2, ...)
  plot(val,log='y',ylim=c(max(val),min(val)),type='p',col = 'grey50', pch = 20, cex = 0.6, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=1.5,pch=20)
      #text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 1.5,cex=1,col=colors)
  }
}



RankView_saga <- function(rankdata, genelist=NULL, top=10, bottom=10, cutoff=NULL,
                          main=NULL, filename=NULL, width=5, height=4, ...){
  requireNamespace("ggrepel", quietly=TRUE) || stop("need ggrepel package")
  if(length(cutoff)==0) cutoff = CutoffCalling(rankdata, 2)
  if(length(cutoff)==1) cutoff = sort(c(-cutoff, cutoff))
  data = data.frame(Gene = names(rankdata), diff = rankdata, stringsAsFactors=FALSE)
  data$Rank = rank(data$diff)
  data$group = "no"
  #data$group[data$diff>cutoff[2]] = "up"
  #data$group[data$diff<cutoff[1]] = "down"
  #data <- data %>%
  #  mutate(group = case_when(
  #    Gene %in% gene_label_up ~ "up",
  #    Gene %in% gene_label_down ~ "down",
  #    TRUE ~ group #keep original value
  #  ))
  
  data <- data %>%
    mutate(group = case_when(
      #Gene %in% gene_label_up ~ "up",
      #Gene %in% gene_label_down ~ "down",
      Gene %in% saga_genes ~ "saga",
      TRUE ~ group #keep original value
    ))
  #dot_size = ifelse(Gene %in% all_label_genes,0.3,0.3)) 
  
  idx=(data$Rank<=bottom) | (data$Rank>(max(data$Rank)-top)) | (data$Gene %in% genelist)
  #"saga"="#EFC000FF",  
  #mycolour = c("up"="#CD534CFF","down"="#0073C2FF")
  mycolour = c("saga"="#008280FF")
  p = ggplot(data)
  #p = p + geom_jitter(aes(x=diff,y=Rank,color=group),size=1.2)
  p = p + geom_jitter(
    data = data %>% filter(group == "no"),
    aes(x = diff, y = Rank, color = group), size = 0.5
  )
  
  p = p + geom_jitter(
    data = data %>% filter(group != "no"),
    aes(x = diff, y = Rank, color = group), size = 1.2
  )
  
  if(!all(cutoff==0)) p = p + geom_vline(xintercept = cutoff, linetype = "dotted")
  if(sum(idx)>0)
    p = p + geom_label_repel(aes(x=diff, y=Rank,fill=group,label = Gene),
                             data=data[idx,], fontface = 'bold', color = 'white', size = 4.5,
                             box.padding = unit(0.4, "lines"), segment.color = 'grey50',
                             point.padding = unit(0.3, "lines"), segment.size = 0.3,max.overlaps = Inf)
  p = p + scale_color_manual(values=mycolour)
  p = p + scale_fill_manual(values=mycolour)
  p = p + theme(panel.background = element_rect(fill='white', colour='black'))
  p = p + theme(text = element_text(colour="black",size = 14, family = "Helvetica"),
                plot.title = element_text(hjust = 0.5, size=18),
                axis.text = element_text(colour="gray10"))
  p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank())
  p = p + labs(x="Log2FC",y="Rank",title=main)
  p = p + theme(legend.position="none")#+ylim(-1000,7000)
  
  return(p)
}



