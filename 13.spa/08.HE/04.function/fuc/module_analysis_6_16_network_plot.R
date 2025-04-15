library(tidyverse)
library(networkD3)
library(igraph)

#########igraph绘制网络图############
library(igraph)
get_edge.data <- function(g_adj,modules){
  node_data_tmp <- E(g_adj)
  #node_data <- get.edgelist(g_adj)
  node_data_tmp1 <- as_ids(node_data_tmp)
  node_data_tmp2 <- data.table::tstrsplit(node_data_tmp1,'[|]')
  node_data <- data.frame(from=node_data_tmp2[[1]],to=node_data_tmp2[[2]])
  node_data$fre <- get.edge.attribute(g_adj,name = 'weight',index = node_data_tmp)
  return(node_data)
}
modules <- readRDS('D:\\project\\pancancer\\module-analysis\\module_analysis_6_16\\modules_clutser_infomap.rds')
adj <- read.table('D:\\project\\pancancer\\module-analysis\\module_analysis_6_16\\filter.list.yes.gene.stat.filter.matrix')
adj <- as.matrix(adj)
g_adj <- graph_from_adjacency_matrix(adj, diag = FALSE, mode = 'undirected', weighted = TRUE)
edge_data <- get_edge.data(g_adj = g_adj,modules = modules)
node_data <- sapply(names(modules),function(x){
  tmp1 <- data.frame(vertices=modules[[x]],group=x)
  list(tmp1)
})
node_data <- Reduce('rbind',node_data)
graph_plot <- graph_from_data_frame(edge_data,directed = F,vertices = node_data)

####plot####
cols_use <- unique(c('#F68282','#B95FBB','#78BE94','#ff9a36','#2FF18B','#B84D64','#faf4cf','#CCB1F1','#25aff5','#A4DFF2',
                     '#7CA878','#AC8F14','#35A132','#8DD3C7','#CFECBB','#F4F4B9','#AF98B5','#8952A0','#F4867C','#C0979F',
                     '#86B1CD','#8DD3C7','#CFECBB','#F4F4B9','#CFCCCF','#D1A7B9','#F4867C','#C0979F','#86B1CD','#CEB28B',
                     '#EDBC63','#C2D567','#CDD796','#F8CDDE','#E9D3DE','#D5CFD6','#C59CC5','#C09CBF','#C9DAC3','#E1EBA0',
                     '#FFED6F',"#99a9cc","#f89e81","#acd485","#dd9bc5","#f6d573","#84c7b3",'#8DD3C7','#CFECBB','#F4F4B9',
                     '#CFCCCF','#D1A7B9'))
###All module
set.seed(50) #生成随机数，这样图的布局就会可重复，而不是每次生成的时候都变
l<-layout.fruchterman.reingold(graph_plot) #设置图的布局方式为弹簧式发散的布局
#具体修改过程
V(graph_plot)$color <- cols_use[as.factor(V(graph_plot)$group)] #根据类型设置颜色,按照类型分组
node_color <- unique(data.frame(group=V(graph_plot)$group,color=V(graph_plot)$color))
V(graph_plot)$label.color <- 'grey' #设置节点标记的颜色
E(graph_plot)$width <- E(graph_plot)$fre #根据频次列设置边宽度
#E(graph_plot)$label <- E(graph_plot)$fre #根据频次列设置边标签
#E(graph_plot)$arrow.size<-0 #设置箭头大小
#生成图
png('D:\\project\\pancancer\\module-analysis\\module_analysis_6_16\\all_module_network.png',res=300,units='in',w=16,h=16)
plot(graph_plot, layout=l,vertex.size=3,vertex.shape='circle',    #节点不带边框none,,圆形边框circle,方块形rectangle  
     vertex.label=NA, #NULL表示不设置，为默认状态，即默认显示数据中点的名称，可以是中文。如果是NA则表示不显示任何点信息	 
     vertex.label.cex=0,    #节点字体大小  
     vertex.label.dist=0,   #标签和节点位置错开
     edge.arrow.size=0,#连线的箭头的大小,若为0即为无向图，当然有些数据格式不支持有向图  
     edge.width = 0.01, #连接线宽度
     edge.label=NA, #不显示连接线标签，默认为频次
     edge.color="grey")
legend('bottomleft',#位置
       legend = node_color$group,#名称
       pt.cex=2,text.font = 1, #点大小，文字大小 
       col = node_color$color,#颜色
       ncol=2,#两列展示
       pch=rep(19,length(node_color$group)),#点的形状
       cex=1.1#设置图例整体的大小
)
dev.off()

#module split
png('D:\\project\\pancancer\\module-analysis\\module_analysis_6_16\\all_module_network_split_1.png',res=300,units='in',h=28,w=30)
par(mfrow = c(4, 5),mar=c(1,1,2,1)+0.1)
for(module_name in unique(V(graph_plot)$group)[1:20]){
  V(graph_plot)$color <- cols_use[as.factor(V(graph_plot)$group)] #根据类型设置颜色,按照类型分组
  node_color <- data.frame(group=V(graph_plot)$group,color=V(graph_plot)$color)
  node_color$color[!(node_color$group %in% module_name)] <- 'grey'
  V(graph_plot)$color <- node_color$color
  V(graph_plot)$label.color <- 'grey' #设置节点标记的颜色
  E(graph_plot)$width <- E(graph_plot)$fre #根据频次列设置边宽度
  #E(graph_plot)$label <- E(graph_plot)$fre #根据频次列设置边标签
  #E(graph_plot)$arrow.size<-0 #设置箭头大小
  #生成图
  plot(graph_plot, layout=l,vertex.size=3,vertex.shape='circle',    #节点不带边框none,,圆形边框circle,方块形rectangle  
       vertex.frame.color=V(graph_plot)$color, #节点边框颜色
       vertex.label=NA, #NULL表示不设置，为默认状态，即默认显示数据中点的名称，可以是中文。如果是NA则表示不显示任何点信息	 
       vertex.label.cex=0,    #节点字体大小  
       vertex.label.dist=0,   #标签和节点位置错开
       edge.arrow.size=0,#连线的箭头的大小,若为0即为无向图，当然有些数据格式不支持有向图  
       edge.width = 0.01, #连接线宽度
       edge.label=NA, #不显示连接线标签，默认为频次
       edge.color="grey")
  title(main = module_name,font=10)
}
dev.off()

png('D:\\project\\pancancer\\module-analysis\\module_analysis_6_16\\all_module_network_split_2.png',res=300,units='in',h=28,w=30)
par(mfrow = c(4, 5),mar=c(1,1,2,1)+0.1)
for(module_name in unique(V(graph_plot)$group)[21:40]){
  V(graph_plot)$color <- cols_use[as.factor(V(graph_plot)$group)] #根据类型设置颜色,按照类型分组
  node_color <- data.frame(group=V(graph_plot)$group,color=V(graph_plot)$color)
  node_color$color[!(node_color$group %in% module_name)] <- 'grey'
  V(graph_plot)$color <- node_color$color
  V(graph_plot)$label.color <- 'grey' #设置节点标记的颜色
  E(graph_plot)$width <- E(graph_plot)$fre #根据频次列设置边宽度
  #E(graph_plot)$label <- E(graph_plot)$fre #根据频次列设置边标签
  #E(graph_plot)$arrow.size<-0 #设置箭头大小
  #生成图
  plot(graph_plot, layout=l,vertex.size=3,vertex.shape='circle',    #节点不带边框none,,圆形边框circle,方块形rectangle  
       vertex.frame.color=V(graph_plot)$color, #节点边框颜色
       vertex.label=NA, #NULL表示不设置，为默认状态，即默认显示数据中点的名称，可以是中文。如果是NA则表示不显示任何点信息	 
       vertex.label.cex=0,    #节点字体大小  
       vertex.label.dist=0,   #标签和节点位置错开
       edge.arrow.size=0,#连线的箭头的大小,若为0即为无向图，当然有些数据格式不支持有向图  
       edge.width = 0.01, #连接线宽度
       edge.label=NA, #不显示连接线标签，默认为频次
       edge.color="grey")
  title(main = module_name,font=10)
}
dev.off()
rm(list=ls())
gc()

#####filter module######
library(igraph)
get_edge.data <- function(g_adj,modules){
  node_data_tmp <- E(g_adj)
  #node_data <- get.edgelist(g_adj)
  node_data_tmp1 <- as_ids(node_data_tmp)
  node_data_tmp2 <- data.table::tstrsplit(node_data_tmp1,'[|]')
  node_data <- data.frame(from=node_data_tmp2[[1]],to=node_data_tmp2[[2]])
  node_data$fre <- get.edge.attribute(g_adj,name = 'weight',index = node_data_tmp)
  return(node_data)
}
modules <- readRDS('D:\\project\\pancancer\\module-analysis\\module_analysis_6_16\\modules_clutser_infomap.rds')
adj <- read.table('D:\\project\\pancancer\\module-analysis\\module_analysis_6_16\\filter.list.yes.gene.stat.filter.matrix')
adj <- as.matrix(adj)
module_gene_number <- sapply(modules,function(x) length(x))
module_filter <- modules[module_gene_number>9]
module_gene <- Reduce('c',module_filter)
index <- match(module_gene,colnames(adj))
adj <- adj[index,index]
g_adj <- graph_from_adjacency_matrix(adj, diag = FALSE, mode = 'undirected', weighted = TRUE)
edge_data <- get_edge.data(g_adj = g_adj,modules = module_filter)
node_data <- sapply(names(module_filter),function(x){
  tmp1 <- data.frame(vertices=modules[[x]],group=x)
  list(tmp1)
})
node_data <- Reduce('rbind',node_data)
graph_plot <- graph_from_data_frame(edge_data,directed = F,vertices = node_data)

####plot####
cols_use <- unique(c('#F68282','#B95FBB','#ff9a36','#2FF18B','#B84D64','#faf4cf','#CCB1F1','#25aff5','#A4DFF2',
                     '#7CA878','#AC8F14','#35A132','#8DD3C7','#CFECBB','#F4F4B9','#AF98B5','#8952A0','#F4867C','#C0979F',
                     '#86B1CD','#8DD3C7','#CFECBB','#F4F4B9','#CFCCCF','#D1A7B9','#F4867C','#C0979F','#86B1CD','#CEB28B',
                     '#EDBC63','#C2D567','#CDD796','#F8CDDE','#E9D3DE','#D5CFD6','#C59CC5','#C09CBF','#C9DAC3','#E1EBA0',
                     '#FFED6F',"#99a9cc","#f89e81","#acd485","#dd9bc5","#f6d573","#84c7b3",'#8DD3C7','#CFECBB','#F4F4B9',
                     '#CFCCCF','#D1A7B9'))
###All module
set.seed(50) #生成随机数，这样图的布局就会可重复，而不是每次生成的时候都变
l<-layout.fruchterman.reingold(graph_plot) #设置图的布局方式为弹簧式发散的布局
#具体修改过程
V(graph_plot)$color <- cols_use[as.factor(V(graph_plot)$group)] #根据类型设置颜色,按照类型分组
node_color <- unique(data.frame(group=V(graph_plot)$group,color=V(graph_plot)$color))
V(graph_plot)$label.color <- 'grey' #设置节点标记的颜色
E(graph_plot)$width <- E(graph_plot)$fre #根据频次列设置边宽度
#E(graph_plot)$label <- E(graph_plot)$fre #根据频次列设置边标签
#E(graph_plot)$arrow.size<-0 #设置箭头大小
#生成图
png('D:\\project\\pancancer\\module-analysis\\module_analysis_6_16\\filter_module_network.png',res=300,units='in',w=16,h=16)
plot(graph_plot, layout=l,vertex.size=2.5,vertex.shape='circle',    #节点不带边框none,,圆形边框circle,方块形rectangle  
     vertex.label=NA, #NULL表示不设置，为默认状态，即默认显示数据中点的名称，可以是中文。如果是NA则表示不显示任何点信息	 
     vertex.label.cex=0,    #节点字体大小  
     vertex.label.dist=0,   #标签和节点位置错开
     edge.arrow.size=0,#连线的箭头的大小,若为0即为无向图，当然有些数据格式不支持有向图  
     edge.width = 0.01, #连接线宽度
     edge.label=NA, #不显示连接线标签，默认为频次
     edge.color="grey")
legend('bottomleft',#位置
       legend = node_color$group,#名称
       pt.cex=2,text.font = 1, #点大小，文字大小 
       col = node_color$color,#颜色
       ncol=1,#两列展示
       pch=rep(19,length(node_color$group)),#点的形状
       cex=1.2#设置图例整体的大小
)
dev.off()

#module split
png('D:\\project\\pancancer\\module-analysis\\module_analysis_6_16\\filter_module_network_split.png',res=300,units='in',h=28,w=30)
par(mfrow = c(4, 5),mar=c(1,1,2,1)+0.1)
for(module_name in unique(V(graph_plot)$group)){
  V(graph_plot)$color <- cols_use[as.factor(V(graph_plot)$group)] #根据类型设置颜色,按照类型分组
  node_color <- data.frame(group=V(graph_plot)$group,color=V(graph_plot)$color)
  node_color$color[!(node_color$group %in% module_name)] <- 'grey'
  V(graph_plot)$color <- node_color$color
  V(graph_plot)$label.color <- 'grey' #设置节点标记的颜色
  E(graph_plot)$width <- E(graph_plot)$fre #根据频次列设置边宽度
  #E(graph_plot)$label <- E(graph_plot)$fre #根据频次列设置边标签
  #E(graph_plot)$arrow.size<-0 #设置箭头大小
  #生成图
  plot(graph_plot, layout=l,vertex.size=2.5,vertex.shape='circle',    #节点不带边框none,,圆形边框circle,方块形rectangle  
       vertex.frame.color=V(graph_plot)$color, #节点边框颜色
       vertex.label=NA, #NULL表示不设置，为默认状态，即默认显示数据中点的名称，可以是中文。如果是NA则表示不显示任何点信息	 
       vertex.label.cex=0,    #节点字体大小  
       vertex.label.dist=0,   #标签和节点位置错开
       edge.arrow.size=0,#连线的箭头的大小,若为0即为无向图，当然有些数据格式不支持有向图  
       edge.width = 0.01, #连接线宽度
       edge.label=NA, #不显示连接线标签，默认为频次
       edge.color="grey")
  title(main = module_name,font=10)
}
dev.off()


###networkD3绘制网络图
modules <- readRDS('D:\\project\\pancancer\\module-analysis\\module_analysis_6_16\\modules_clutser_infomap.rds')
adj <- read.table('D:\\project\\pancancer\\module-analysis\\module_analysis_6_16\\filter.list.yes.gene.stat.filter.matrix')
adj <- as.matrix(adj)
g_adj <- graph_from_adjacency_matrix(adj, diag = FALSE, mode = 'undirected', weighted = TRUE)
get_node <- function(g_adj,modules){
  node_data_tmp <- E(g_adj)
  #node_data <- get.edgelist(g_adj)
  node_data_tmp1 <- as_ids(node_data_tmp)
  node_data_tmp2 <- data.table::tstrsplit(node_data_tmp1,'[|]')
  node_data <- data.frame(source=node_data_tmp2[[1]],target=node_data_tmp2[[2]])
  node_data$value <- get.edge.attribute(g_adj,name = 'weight',index = node_data_tmp)
  modules_gene <- sapply(names(modules),function(x){
    tmp1 <- data.frame(module_name=x,gene=modules[[x]])
    list(tmp1)
  })
  modules_gene <- Reduce('rbind',modules_gene)
  modules_tmp1 <- cbind(modules_gene$module_name[match(node_data$source,modules_gene$gene)],
                        modules_gene$module_name[match(node_data$target,modules_gene$gene)])
  modules_tmp1 <- sapply(1:nrow(modules_tmp1),function(x){
    if(modules_tmp1[x,1]==modules_tmp1[x,2]){
      out_data <- 'one-module'
    } else {
      out_data <- 'two-module'
    }
    out_data
  })
  node_data$color <- modules_tmp1
  return(node_data)
}

edge_data <- get_node(g_adj = g_adj,modules = data)
node_data <- sapply(names(modules),function(x){
  tmp1 <- data.frame(name=modules[[x]],group=x,size=1)
  list(tmp1)
})
node_data <- Reduce('rbind',node_data)

cols_use <- unique(c('#F68282','#B95FBB','#78BE94','#ff9a36','#2FF18B','#B84D64','#faf4cf','#CCB1F1','#25aff5','#A4DFF2',
              '#7CA878','#AC8F14','#35A132','#8DD3C7','#CFECBB','#F4F4B9','#AF98B5','#8952A0','#F4867C','#C0979F',
              '#86B1CD','#8DD3C7','#CFECBB','#F4F4B9','#CFCCCF','#D1A7B9','#F4867C','#C0979F','#86B1CD','#CEB28B',
              '#EDBC63','#C2D567','#CDD796','#F8CDDE','#E9D3DE','#D5CFD6','#C59CC5','#C09CBF','#C9DAC3','#E1EBA0',
              '#FFED6F',"#99a9cc","#f89e81","#acd485","#dd9bc5","#f6d573","#84c7b3",'#8DD3C7','#CFECBB','#F4F4B9',
              '#CFCCCF','#D1A7B9'))
node_data$group_color <- cols_use[as.factor(node_data$group)]

# 处理数据
# 因为networkD3需要的连线数据，是节点文件里的名称的索引。所以，需要做一个名称到索引的转化
Node2index <- list()
Node2index[node_data$name] <- 0:(length(node_data$name)-1)

edge_data <- edge_data %>%
  mutate(source2 = unlist(Node2index[source])) %>%
  mutate(target2 = unlist(Node2index[target]))
edge_data$edge_color <- c('red','grey')[as.factor(edge_data$color)]

# 定义颜色
group2project <- paste(unique(node_data$group),collapse = '","')
color2project <- paste(unique(node_data$group_color),collapse = '","')
my_color <- paste0('d3.scaleOrdinal().domain(["',group2project,'"]).range(["',color2project,'"])')


# 绘图
forceNetwork(Links = edge_data, 
             Nodes = node_data,
             Source = "source2", 
             Target = "target2",
             Value ="value",
             NodeID = "name",
             Group = "group", 
             opacity= 1,        # 透明度
             Nodesize="size",
             zoom = F,       # 是否可以缩放
             opacityNoHover=1,  # 鼠标没有悬浮在节点上时，文字的透明度(0-1)
             colourScale = JS(my_color),   # 节点颜色，JavaScript
             legend=T, 
             fontSize = 0,
             linkColour= edge_data$edge_color
)




##############test########################
get_edge.data <- function(g_adj,modules){
  node_data_tmp <- E(g_adj)
  #node_data <- get.edgelist(g_adj)
  node_data_tmp1 <- as_ids(node_data_tmp)
  node_data_tmp2 <- data.table::tstrsplit(node_data_tmp1,'[|]')
  node_data <- data.frame(from=node_data_tmp2[[1]],to=node_data_tmp2[[2]])
  node_data$fre <- get.edge.attribute(g_adj,name = 'weight',index = node_data_tmp)
  return(node_data)
}
edge_data <- get_edge.data(g_adj = g_adj,modules = modules)
node_data <- sapply(names(modules),function(x){
  tmp1 <- data.frame(vertices=modules[[x]],group=x)
  list(tmp1)
})
node_data <- Reduce('rbind',node_data)
graph_plot <- graph_from_data_frame(edge_data,directed = F,vertices = node_data)


plot(graph_plot,  
     layout=l,  #layout.fruchterman.reingold表示弹簧式发散的布局，
     #其他还有环形布局layout.circle，分层布局layout.reingold.tilford，中心向外发散layout.reingold.tilford(graph,circular=T) ，核心布局layout_as_star，大型网络可视化layout_with_drl
     vertex.size=8,     #节点大小  
     vertex.shape='circle',    #节点不带边框none,,圆形边框circle,方块形rectangle  
     vertex.color="yellow",#设置颜色，其他如red,blue,cyan,yellow等
     vertex.label=NULL, #NULL表示不设置，为默认状态，即默认显示数据中点的名称，可以是中文。如果是NA则表示不显示任何点信息	 
     vertex.label.cex=0.8,    #节点字体大小  
     vertex.label.color='black',  #节点字体颜色,red  
     vertex.label.dist=0.4,   #标签和节点位置错开
     edge.arrow.size=0.3,#连线的箭头的大小,若为0即为无向图，当然有些数据格式不支持有向图  
     edge.width = 0.5, #连接线宽度
     edge.label=NA, #不显示连接线标签，默认为频次
     edge.color="grey")  #连线颜色 

#生成方式2：
set.seed(50) #生成随机数，这样图的布局就会可重复，而不是每次生成的时候都变
l<-layout.fruchterman.reingold(graph_plot) #设置图的布局方式为弹簧式发散的布局
#具体修改过程
V(graph_plot)$color <- cols_use[as.factor(V(graph_plot)$group)] #根据类型设置颜色,按照类型分组
node_color <- unique(data.frame(group=V(graph_plot)$group,color=V(graph_plot)$color))
V(graph_plot)$label.color <- 'grey' #设置节点标记的颜色
E(graph_plot)$width <- E(graph_plot)$fre #根据频次列设置边宽度
#E(graph_plot)$label <- E(graph_plot)$fre #根据频次列设置边标签
#E(graph_plot)$arrow.size<-0 #设置箭头大小
#生成图
png('D:\\project\\pancancer\\module-analysis\\module_analysis_6_16\\network_test.png',res=300,units='in',w=16,h=16)
plot(graph_plot, layout=l,vertex.size=3,vertex.shape='circle',    #节点不带边框none,,圆形边框circle,方块形rectangle  
     vertex.label=NA, #NULL表示不设置，为默认状态，即默认显示数据中点的名称，可以是中文。如果是NA则表示不显示任何点信息	 
     vertex.label.cex=0,    #节点字体大小  
     vertex.label.dist=0,   #标签和节点位置错开
     edge.arrow.size=0,#连线的箭头的大小,若为0即为无向图，当然有些数据格式不支持有向图  
     edge.width = 0.01, #连接线宽度
     edge.label=NA, #不显示连接线标签，默认为频次
     edge.color="grey")
legend('bottomleft',#位置
       legend = node_color$group,#名称
       pt.cex=2,text.font = 1, #点大小，文字大小 
       col = node_color$color,#颜色
       ncol=2,#两列展示
       pch=rep(19,length(node_color$group)),#点的形状
       cex=1.1#设置图例整体的大小
       )
dev.off()

module_name <- 'module-1'
V(graph_plot)$color <- cols_use[as.factor(V(graph_plot)$group)] #根据类型设置颜色,按照类型分组
node_color <- data.frame(group=V(graph_plot)$group,color=V(graph_plot)$color)
node_color$color[!(node_color$group %in% module_name)] <- 'grey'
V(graph_plot)$color <- node_color$color
V(graph_plot)$label.color <- 'grey' #设置节点标记的颜色
E(graph_plot)$width <- E(graph_plot)$fre #根据频次列设置边宽度
#E(graph_plot)$label <- E(graph_plot)$fre #根据频次列设置边标签
#E(graph_plot)$arrow.size<-0 #设置箭头大小
#生成图
png('D:\\project\\pancancer\\module-analysis\\module_analysis_6_16\\network_test1.png',res=300,units='in',w=16,h=16)
plot(graph_plot, layout=l,vertex.size=3,vertex.shape='circle',    #节点不带边框none,,圆形边框circle,方块形rectangle  
     vertex.label=NA, #NULL表示不设置，为默认状态，即默认显示数据中点的名称，可以是中文。如果是NA则表示不显示任何点信息	 
     vertex.label.cex=0,    #节点字体大小  
     vertex.label.dist=0,   #标签和节点位置错开
     edge.arrow.size=0,#连线的箭头的大小,若为0即为无向图，当然有些数据格式不支持有向图  
     edge.width = 0.01, #连接线宽度
     edge.label=NA, #不显示连接线标签，默认为频次
     edge.color="grey")
dev.off()

png('D:\\project\\pancancer\\module-analysis\\module_analysis_6_16\\network_test1.png',res=300,units='in',h=16,w=30)
par(mfrow = c(2, 5),mar=c(1,1,2,1)+0.1)
for(module_name in unique(V(graph_plot)$group)[1:20]){
  V(graph_plot)$color <- cols_use[as.factor(V(graph_plot)$group)] #根据类型设置颜色,按照类型分组
  node_color <- data.frame(group=V(graph_plot)$group,color=V(graph_plot)$color)
  node_color$color[!(node_color$group %in% module_name)] <- 'grey'
  V(graph_plot)$color <- node_color$color
  V(graph_plot)$label.color <- 'grey' #设置节点标记的颜色
  E(graph_plot)$width <- E(graph_plot)$fre #根据频次列设置边宽度
  #E(graph_plot)$label <- E(graph_plot)$fre #根据频次列设置边标签
  #E(graph_plot)$arrow.size<-0 #设置箭头大小
  #生成图
  plot(graph_plot, layout=l,vertex.size=3,vertex.shape='circle',    #节点不带边框none,,圆形边框circle,方块形rectangle  
       vertex.frame.color=V(graph_plot)$color, #节点边框颜色
       vertex.label=NA, #NULL表示不设置，为默认状态，即默认显示数据中点的名称，可以是中文。如果是NA则表示不显示任何点信息	 
       vertex.label.cex=0,    #节点字体大小  
       vertex.label.dist=0,   #标签和节点位置错开
       edge.arrow.size=0,#连线的箭头的大小,若为0即为无向图，当然有些数据格式不支持有向图  
       edge.width = 0.01, #连接线宽度
       edge.label=NA, #不显示连接线标签，默认为频次
       edge.color="grey")
  title(main = module_name,font=6)
}
dev.off()

module_name <- 'module-1'
V(graph_plot)$color <- cols_use[as.factor(V(graph_plot)$group)] #根据类型设置颜色,按照类型分组
node_color <- data.frame(group=V(graph_plot)$group,color=V(graph_plot)$color)
node_color$color[!(node_color$group %in% module_name)] <- 'grey'
V(graph_plot)$color <- node_color$color
V(graph_plot)$label.color <- 'grey' #设置节点标记的颜色
E(graph_plot)$width <- E(graph_plot)$fre #根据频次列设置边宽度
#E(graph_plot)$label <- E(graph_plot)$fre #根据频次列设置边标签
#E(graph_plot)$arrow.size<-0 #设置箭头大小
#生成图
png('D:\\project\\pancancer\\module-analysis\\module_analysis_6_16\\network_test1.png',res=300,units='in',w=16,h=16)
plot(graph_plot, layout=l,vertex.size=3,vertex.shape='circle',    #节点不带边框none,,圆形边框circle,方块形rectangle  
     vertex.label=NA, #NULL表示不设置，为默认状态，即默认显示数据中点的名称，可以是中文。如果是NA则表示不显示任何点信息	 
     vertex.label.cex=0,    #节点字体大小  
     vertex.label.dist=0,   #标签和节点位置错开
     edge.arrow.size=0,#连线的箭头的大小,若为0即为无向图，当然有些数据格式不支持有向图  
     edge.width = 0.01, #连接线宽度
     edge.label=NA, #不显示连接线标签，默认为频次
     edge.color="grey")
dev.off()
