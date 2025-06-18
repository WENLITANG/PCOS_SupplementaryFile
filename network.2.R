library(phyloseq)
library(ggClusterNet)
library(tidyverse)
library(Biostrings)
library(xlsx)
library(igraph)
library(network)
library(sna)
library(tidyverse)

metadata = read.delim("../../PCOS代谢组batch1/PCOS_batch1.代谢组/PCOS_batch1.meta.txt",row.names = 1)
otutab = read.delim("../../PCOS代谢组batch1/PCOS_batch1.代谢组/neg_pos_KEGG.CAS.HMDB_log10.txt", row.names=1)
otutab$ZC04 <- NULL
# taxonomy = xlsx::read.xlsx("inputdata.xlsx", row.names=1, sheetName = "newclass")
taxonomy = fread("../../PCOS代谢组batch1/PCOS_batch1.代谢组/KEGG.CAS.HMDB_log10_ort3.vip_anno_cv_class.txt")
taxonomy$cSuper.Class = gsub(",", " ", taxonomy$cSuper.Class)
taxonomy <- as.data.frame(taxonomy)
rownames(taxonomy) <- taxonomy$V1
# tree  = read_tree("./otus.tree")
# rep = readDNAStringSet("./otus.fa")

#rownames(taxonomy)[rownames(taxonomy) %in% setdiff(rownames(taxonomy), rownames(otutab))] <- taxonomy[setdiff(rownames(taxonomy), rownames(otutab)), "HMDB"]

ps = phyloseq(sample_data(metadata),
              otu_table(as.matrix(otutab), taxa_are_rows=TRUE),
              tax_table(as.matrix(taxonomy[, c("cSuper.Class", "cClass")]))#,
              # phy_tree(tree),
              # refseq(rep)
)

#################################################
##              model_igraph布局               ##
#################################################

#按model上色和聚类
# result = cor_Big_micro(ps = ps,
#                        N = 0,
#                        r.threshold=0.8,
#                        p.threshold=0.05,
#                        method = "spearman"
# )
# 
# #--提取相关矩阵
# cor = result[[1]]
# dim(cor)
# 
# result2 <- model_igraph(cor = cor,
#                         method = "cluster_fast_greedy",
#                         seed = 12
# )
# node = result2[[1]]
# head(node)
# 
# dat = result2[[2]]
# head(dat)
# tem = data.frame(mod = dat$model,col = dat$color) %>%  
#   dplyr::distinct( mod, .keep_all = TRUE)  
# col = tem$col
# names(col) = tem$mod
# 
# #---node节点注释#-----------
# otu_table = as.data.frame(t(vegan_otu(ps)))
# tax_table = as.data.frame(vegan_tax(ps))
# nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax_table)
# head(nodes)
# #-----计算边#--------
# edge = edgeBuild(cor = cor,node = node)
# colnames(edge)[8] = "cor"
# head(edge)
# 
# tem2 = dat %>% 
#   dplyr::select(OTU,model,color) %>%
#   dplyr::right_join(edge,by =c("OTU" = "OTU_1" ) ) %>%
#   dplyr::rename(OTU_1 = OTU,model1 = model,color1 = color)
# head(tem2)
# 
# tem3 = dat %>% 
#   dplyr::select(OTU,model,color) %>%
#   dplyr::right_join(edge,by =c("OTU" = "OTU_2" ) ) %>%
#   dplyr::rename(OTU_2 = OTU,model2 = model,color2 = color)
# head(tem3)
# 
# tem4 = tem2 %>%inner_join(tem3)
# head(tem4)
# 
# edge2 = tem4 %>% mutate(color = ifelse(model1 == model2,as.character(model1),"across"),
#                         manual = ifelse(model1 == model2,as.character(color1),"#C1C1C1")
# )
# head(edge2)
# col_edge = edge2 %>% dplyr::distinct(color, .keep_all = TRUE)  %>% 
#   select(color,manual)
# col0 = col_edge$manual
# names(col0) = col_edge$color
# 
# library(ggnewscale)
# 
# p1 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = color),
#                               data = edge2, size = 1) +
#   scale_colour_manual(values = col0) 
# 
# # ggsave("./cs1.pdf",p1,width = 16,height = 14)
# p2 = p1 +
#   new_scale_color() +
#   geom_point(aes(X1, X2,color =model), data = dat,size = 4) +
#   scale_colour_manual(values = col) +
#   scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
#   theme(panel.background = element_blank()) +
#   theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
#   theme(legend.background = element_rect(colour = NA)) +
#   theme(panel.background = element_rect(fill = "white",  colour = NA)) +
#   theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
# p2
# ggsave("model_igraph_cs2.pdf",p2,width = 16,height = 14)

########################################################
##  按class布局和上色

# path <- "network.2"
# map = sample_data(ps)
# map$Group = "one"
# sample_data(ps) = map
dir.create(path <- "校正数据后network/cSuper.Class_2group_person_r0.8_p0.05_max30")
result = network.2(ps = ps,
                   N = 0,
                   big = TRUE,
                   maxnode = 30,
                   select_layout = TRUE,
                   layout_net = "model_igraph",
                   r.threshold=0.8,
                   p.threshold=0.05,
                   label = FALSE,
                   fill = "cSuper.Class", #"cClass", #"cSuper.Class",
                   group = "group",
                   path = path,
                   #method = "spearman",
                   zipi = TRUE)

# 多组网络绘制到一个面板
p = result[[1]]
# 全部样本网络参数比对
data = result[[2]]
num= 3
# plotname1 = paste(path,"/network_all.jpg",sep = "")
# ggsave(plotname1, p,width = 16*num,height = 16,dpi = 72)

plotname1 = paste(path,"/network_all.pdf",sep = "")
ggsave(plotname1, p,width = 16*num,height = 16,limitsize = FALSE)

tablename <- paste(path,"/co-occurrence_Grobel_net",".csv",sep = "")
write.csv(data,tablename)

