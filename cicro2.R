############################################
######cicro
#############################################
#rm(list = ls())
library(tidyverse);library(reshape);library(graphics);library(RColorBrewer);library(ggplot2)
#library(xlsx)
# Create dataset
#data <- read.table('clipboard',header = T,sep = '\t',quote = "")
data <- read.table(file = 'inputdata.txt',sep="\t", header=T, check.name=F,comment.char = "",quote="")
# Set a number of 'empty bar'
# empty_bar <- 10
# 
# # Add lines to the initial dataset
# to_add <- matrix(NA, empty_bar, ncol(data))
# colnames(to_add) <- colnames(data)
# data <- rbind(data, to_add)
# data$id <- seq(1, nrow(data))
# Get the name and the y position of each label


label_data <- subset(data, select = c("id","log2FC","to","lgP"))
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# Make the plot
ggplot(data, aes(x=as.factor(id), y=abs(log2FC))) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(stat="identity", fill=alpha("gray", 0.3)) +
  ylim(-15,15) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar(start = 0) +
  geom_text(data=label_data, aes(x=id, y=abs(log2FC), label=to, hjust=hjust), 
            color="black", fontface="bold",alpha=0.6, size=2.5, 
            angle= label_data$angle, inherit.aes = FALSE ) 
ggsave('bar.pdf',width = 10,height = 10)

# ##C,P,d or VIP
# data2 <- subset(data, select = c("name","VIP","to","lgP"))
# data2$VIP <- abs(data2$VIP)
# data2_melt <- melt(data2,id.vars = c('name','to'),measure.vars = c('VIP'))
# data2_melt$name <- factor(data2_melt$name,levels = data2$name)
# #data2_melt$value2 <- ifelse(data2_melt$value>2,)
# ggplot(data2_melt,aes(name,variable,fill=value))+geom_tile()+
#   #scale_fill_distiller(palette = )
#   scale_fill_gradientn(colours = c(#colorRampPalette(colors = c("#bf812d"))(1),
#                                    #colorRampPalette(colors = c('grey','WhiteSmoke'))(40),
#                                    #colorRampPalette(colors = c('#f5f5f5','#f5f5f5'))(19),
#                              colorRampPalette(colors = c('#cddbec','#004eb4'))(81)
#                              #colorRampPalette(colors = c("#1a9850"))(2)
#                              ))+
#   coord_polar(theta = "x")+
#   theme_void()+
#   scale_y_discrete(expand = c(0,10))
# ggsave('C_P_VIP.pdf',width = 8,height = 8)

# class
data3 <- subset(data, select = c("name","to","Super.Class"))
data3$id <- 'Y'
data3$to <- factor(data3$to,levels = data3$to)
colors <- colorRampPalette(colors = brewer.pal(n = 9,'Accent'))(13)
colors <- c("#D89200", "#F28F42", "#8E831E", "#6E9326", "#0080BA", "#2961B3", "#4D2580", "#C42631")
ggplot(data3,aes(to,id,fill=Super.Class))+
  geom_tile(show.legend = T)+
  #theme(axis.text.x = element_text(angle = 90))
  scale_fill_manual(values = colors)+
  coord_polar(theta = "x")+theme_void()+
  scale_y_discrete(expand = c(0,11))
ggsave('super_class.pdf',width = 10,height = 10)

###cicir
library(ggraph)
library(igraph)
library(tidyverse)
library(RColorBrewer) 


# create a data frame giving the hierarchical structure of your individuals
d1=data.frame(from="origin", to=c('C','group1',
                                  'P','group2'))
d2=subset(data, select = c("from","to"))
#d2$to <- factor(d2$to)
edges=rbind(d1, d2)

# create a vertices data.frame. One line per object of our hierarchy
vertices = data.frame(
  name = unique(c(as.character(edges$from), as.character(edges$to)))
) 
vertices$name <- factor(vertices$name,levels = unique(vertices$name))
#abt <- "ss_metabolites_lgpadj" ## define the size of points
abt <- "co_metabolites_lgpadj"## define the size of points
abundant <- subset(data, select = c(abt,"to")) ## define the size of points
abundant$CV <- abundant[,abt] ## define the size of points
vertices <- merge(vertices,abundant,by.x="name",by.y='to',all.x=T)
vertices$CV <- ifelse(is.na(vertices$CV),0,vertices$CV)
table(vertices$CV)
#vertices$value <- 
# Let's add a column with the group of each name. It will be useful later to color points
vertices$group = edges$from[ match( vertices$name, edges$to ) ]


#Let's add information concerning the label we are going to add: angle, horizontal adjustement and potential flip
#calculate the ANGLE of the labels
vertices$id=NA
myleaves=which(is.na( match(vertices$name, edges$from) ))
nleaves=length(myleaves)
vertices$id[ myleaves ] = seq(1:nleaves)
vertices$angle= 90 - 360 * vertices$id / nleaves

# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90
vertices$hjust<-ifelse( vertices$angle < -90, 1, 0)

# flip angle BY to make them readable
vertices$angle<-ifelse(vertices$angle < -90, vertices$angle+180, vertices$angle)

# Create a graph object
#edges$from <- factor(edges$from,levels = c('SBP','group1','DBP','group2'))
#vertices$id <- factor(vertices,levels = unique(vertices$name))
mygraph <- graph_from_data_frame( edges, vertices=vertices )

# Make the plot
pdf(paste(abt,'cicor.pdf',sep = "_"),width = 20,height = 20)
ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
  geom_edge_diagonal(colour="grey") +
  scale_edge_colour_distiller(palette = "RdPu") +
  #geom_node_text(aes(x = x*1.15, y=y*1.15,label=name ,colour=group), size=2.7, alpha=1) +
  geom_node_point(aes(filter = leaf, x = x*1.07, y=y*1.07, colour=group, size=sqrt(CV), alpha=1)) +
  scale_colour_manual(values= brewer.pal(9,"Paired") , 30) +
  scale_size_continuous(range = c(0.2,12)) +
  theme_void() +
  theme(
    legend.position="none",
    plot.margin=unit(c(0,0,0,0),"cm"),
  ) +
  expand_limits(x = c(-1.3, 1.3), y = c(-1.3, 1.3))
dev.off() 


# lgP
data4 <- subset(data, select = c("name","to","lgP"))
data4_melt <- melt(data4,id.vars = c('name','to'))
data4_melt$name <- factor(data4_melt$name,levels = data4$name)
#data4_melt$value2 <- ifelse(data4_melt$value>2,)
ggplot(data4_melt,aes(name,variable,fill=value))+geom_tile()+
  #scale_fill_distiller(palette = )
  scale_fill_gradientn(
    colours = (c(#colorRampPalette(colors = 'white')(1),
    rev(colorRampPalette(colors = c('#0ca678','#e6fcf5'))(40))
    #colorRampPalette(colors = c("#1a9850"))(2)
  )))+
  coord_polar(theta = "x")+theme_void()+
  scale_y_discrete(expand = c(0,10))
ggsave('lgP.pdf',width = 10,height = 10)






# Lactobacillus
data4 <- read.table('clipboard',header = T,sep = '\t',quote = "")
data4_melt <- melt(data4,id.vars = c('name','to'))
data4_melt$name <- factor(data4_melt$name,levels = data4$name)
#data4_melt$value2 <- ifelse(data4_melt$value>2,)
ggplot(data4_melt,aes(name,variable,fill=value))+geom_tile()+
  #scale_fill_distiller(palette = )
  scale_fill_gradientn(colours = (c(#colorRampPalette(colors = c("#bf812d"))(1),
    colorRampPalette(colors = c('#fe621f','#fbe6d1'))(40),
    colorRampPalette(colors = 'white')(3),
    colorRampPalette(colors = c('#e0f3f8','#4575b4'))(40)
    #colorRampPalette(colors = c("#1a9850"))(2)
  )))+
  coord_polar(theta = "x")+theme_void()+
  scale_y_discrete(expand = c(0,10))
ggsave('Lacto.pdf',width = 10,height = 10)




