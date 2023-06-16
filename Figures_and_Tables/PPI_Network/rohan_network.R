library(ggplot2) 
library(tidyverse) 
library(ggforce) 
library(gridExtra) 
library(scatterpie) 
library(RColorBrewer) 
library(ggrepel)


plotNetWork <- function(xf.nodes, xf.edges, target_key,alternative_key){ 
  
  xdf.nodes <- read.csv(xf.nodes,sep = '\t',stringsAsFactors = F) 
  
  xdf.nodes$listCount <- rowSums(xdf.nodes[c(target_key,alternative_key)] > 0) 
  
  colnames(xdf.nodes)[1] <- 'key' 
  
  xdf.final <- xdf.nodes[rowSums(xdf.nodes[target_key]) > 0,] 
  
  xdf.other <- xdf.nodes[rowSums(xdf.nodes[target_key]) == 0,] 
  
  xdf.edges <- read.csv(xf.edges,header = T,sep='\t', stringsAsFactors = F) 
  
  xdf.edges$edgeType <- factor(xdf.edges$edgeType, levels = c('C','S','G','M')) 
  
  x.fill.color <- c(brewer.pal(9,'Set1')[1:5],brewer.pal(12,'Set3')) 
  
  
  
  
  
  # split xdf.edges 
  
  xdf.edges.final <- xdf.edges[xdf.edges$node1 %in% xdf.final$key & xdf.edges$node2 %in% xdf.final$key,] 
  
  xdf.edges.other <- xdf.edges[!(xdf.edges$node1 %in% xdf.final$key & xdf.edges$node2 %in% xdf.final$key),] 
  
  x.nnn=2 
  
  p <- ggplot()  + 
    geom_segment(data = xdf.edges.final, aes(x=x,y=y,xend=xend,yend=yend, color=edgeType)) +
    geom_segment(data = xdf.edges.other, aes(x=x,y=y,xend=xend,yend=yend, color=edgeType), size=0.3) +
    scale_color_manual(values=brewer.pal(8,'Set2')) + 
    geom_scatterpie(aes(x=x, y=y, r=sqrt(listCount)*0.04 * x.nnn), data=xdf.nodes, cols=c(target_key,alternative_key), color=NA) + 
    coord_fixed() + 
    theme_void()+ 
    geom_text_repel(min.segment.length = Inf, data = xdf.final, aes(x=x, y=y-0.03 * x.nnn, label=key, min.segment.length = Inf), size=3,color='blue', force_pull = -0.08) + 
    scale_fill_manual(values=x.fill.color) + geom_text(min.segment.length = Inf, data = xdf.other, aes(x=x, y=y-0.03* x.nnn, label=key, min.segment.length = Inf), size=3,color='black') 
  
  return(p) 
  
} 







pdf('/lab01/Projects/Rohan_Projects/CNV_Project/2023/CNV_Pipeline/Results/Plots/Network_PPI/20230605.networkx.pdf',width = 12, height = 12) 

xf.nodes <- '/lab01/Projects/Rohan_Projects/CNV_Project/2023/CNV_Pipeline/Results/Plots/Network_PPI/20230521xiaolong/20230605.networkx.nodes.tsv' 

xf.edges <- '/lab01/Projects/Rohan_Projects/CNV_Project/2023/CNV_Pipeline/Results/Plots/Network_PPI/20230521xiaolong/20230605.networkx.edges.tsv' 

target_key = c('target') 

alternative_key = c('ADHD_evidence', 'other_Neuro', 'SFARI') 

print(plotNetWork(xf.nodes, xf.edges, target_key = target_key,alternative_key=alternative_key)) 



dev.off() 