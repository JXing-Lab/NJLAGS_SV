

library(ggplot2) 
library(tidyverse) 
library(ggforce) 
library(gridExtra) 
library(scatterpie) 
library(RColorBrewer) 
library(ggrepel)
library(GGally)
library(shadowtext)




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
    geom_segment(data = xdf.edges.final, aes(x=x*2,y=y*2,xend=xend*2,yend=yend*2, color=edgeType)) +
    geom_segment(data = xdf.edges.other, aes(x=x*2,y=y*2,xend=xend*2,yend=yend*2, color=edgeType), linewidth=0.3) +
    scale_color_manual(values=brewer.pal(8,'Set2')) + 
    geom_scatterpie(aes(x=x*2, y=y*2, r=sqrt(listCount)*0.06 * x.nnn), data=xdf.nodes, cols=c(target_key,alternative_key), color=NA) + 
    coord_fixed() + 
    theme_void()+ 
    geom_text(data = xdf.final, 
                    #min.segment.length = Inf,
                    aes(x=x*2, y=y*2-0.1 * x.nnn, label=key), 
                    linewidth=3,color='blue') + 
    scale_fill_manual(values=x.fill.color) + 
                    geom_text( 
                    data = xdf.other, 
                    #min.segment.length = Inf,
                    aes(x=x, y=y-0.1* x.nnn, label=key), 
                    linewidth=3,color='black') 
  
  return(p) 
  

} 







pdf('/home/rohan/CNV_Project/2023/Submission_Pipeline/Results/Plots/Network_PPI/20230609_truncated.networkx_lineless_large.pdf',width = 12, height = 12) 

xf.nodes <- '/home/rohan/CNV_Project/2023/Submission_Pipeline/Results/Plots/Network_PPI/20230609_truncated.networkx.nodes.tsv' 

xf.edges <- '/home/rohan/CNV_Project/2023/Submission_Pipeline/Results/Plots/Network_PPI/20230609_truncated.networkx.edges.tsv' 

target_key = c('target') 

alternative_key = c('ADHD_evidence', 'other_Neuro', 'SFARI') 

print(plotNetWork(xf.nodes, xf.edges, target_key = target_key,alternative_key=alternative_key)) 



dev.off() 