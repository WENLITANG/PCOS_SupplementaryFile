result3 <- fread("diffsp108_mbFlavonoids34_cor.txt")
filtered <- result3[#result3$spearman.p < 0.05 & 
                      result3$q0.2 %in% c("CX_QJ","CX_QJ_mice"), ] %>% droplevels()
filtered_sp <- distinct(filtered[, c("Var1", "q0.2", "direction")])
dec <- filtered_sp$Var1[filtered_sp$direction %in% "decrease"] %>% as.character()
inc <- filtered_sp$Var1[filtered_sp$direction %in% "increase"] %>% as.character()

cx <- read.table("DESeq_l0.1/transformed_counts.txt", header = T, row.names = 1, sep = "\t", quote = "")
qj <- read.table("../DESeq/transformed_counts.QJ.txt", header = T, row.names = 1, sep = "\t", quote = "")
mice <- read.table("../DESeq/transformed_counts.mice.txt", header = T, row.names = 1, sep = "\t", quote = "")
cx_meta <- read.table("PCOS_CX.meta.txt", header = T, row.names = 1, sep = "\t", comment.char = "")
qj_meta <- read.table("../DESeq/metadata.qiaojie.txt", header = T, row.names = 1, sep = "\t")
mice_meta <- read.table("../DESeq/meta.tsv", header = T, row.names = 1, sep = "\t")
rownames(mice_meta) <- gsub("-", ".", rownames(mice_meta))

dcx <- data.frame(group=cx_meta[colnames(cx), "group"], t(cx[as.character(filtered_sp$Var1), ]), check.names = F)
dqj <- data.frame(group=qj_meta[colnames(qj), "group"], t(qj[as.character(filtered_sp$Var1), ]), check.names = F)
dmice <- data.frame(group=mice_meta[colnames(mice), "group"], t(mice[as.character(filtered_sp$Var1), ]), check.names = F)
apply(dcx[, -1], 2, function(x) tapply(x, dcx$group, mean))

dcx$Group <- ifelse(dcx$group == "C", "Ctrl", "PCOS")
dqj$Group <- ifelse(dqj$group == "C", "Ctrl", "PCOS")
dmice$Group <- ifelse(dmice$group == "ctrl", "Ctrl", "PCOS")
dcx$group <- dqj$group <- dmice$group <- NULL
dcx$dataset <- "CX"
dqj$dataset <- "QJ"
dmice$dataset <- "mice"
dta <- rbind(melt(dcx), melt(dqj), melt(dmice))
dta$subgroup <- paste(dta$dataset, dta$Group, sep = "_")
edta <- distinct(dta[, -4])
library(ggplot2)
library(ggnewscale)
dta$dataset <- factor(dta$dataset, levels = rev(c("CX","QJ","mice")))
fwrite(dta, "CX.QJ.overlap_diffsp_deseqcount_bar.txt", sep = "\t")

cx_qj_mice_info <- fread("mergeDESeq.CXnew_QJ_mice.human_species_q0.2.txt")
padj <- cx_qj_mice_info[cx_qj_mice_info$Annotation %in% as.character(filtered_sp$Var1),
                        c("Annotation", "padj.CX", "padj.QJ", "padj")]
padj <- melt(padj)
padj$dataset <- ifelse(padj$variable == "padj.CX", "CX",
                       ifelse(padj$variable == "padj.QJ", "QJ", "mice"))
padj$subgroup <- paste(padj$dataset, padj$Annotation, sep = "_")
colnames(padj) <- gsub("value", "padj", colnames(padj))
lgfc <- cx_qj_mice_info[cx_qj_mice_info$Annotation %in% as.character(filtered_sp$Var1),
                        c("Annotation", "log2FoldChange.CX", "log2FoldChange.QJ", "log2FoldChange")]
lgfc <- melt(lgfc)
lgfc$dataset <- ifelse(lgfc$variable == "log2FoldChange.CX", "CX",
                       ifelse(lgfc$variable == "log2FoldChange.QJ", "QJ", "mice"))
lgfc$dataset <- factor(lgfc$dataset, levels = rev(c("CX","QJ","mice")))
lgfc$subgroup <- paste(lgfc$dataset, lgfc$Annotation, sep = "_")
lgfc_padj <- merge(lgfc, padj[, c("subgroup", "padj")], by = "subgroup")
lgfc_padj$label <- ifelse(lgfc_padj$padj < 0.001, "***",
                          ifelse(lgfc_padj$padj < 0.01, "**",
                                 ifelse(lgfc_padj$padj < 0.05, "*",
                                        ifelse(lgfc_padj$padj < 0.2, "#", ""))))
lgfc_padj$dataset <- factor(lgfc_padj$dataset, levels = rev(c("CX","QJ","mice")))
fwrite(lgfc_padj, "CX.QJ.overlap_diffsp_direction_plot.txt", sep = "\t")

lgfc_padj <- lgfc_padj[order(lgfc_padj$value), ]
splevel <- lgfc_padj$Annotation %>% unique()
lgfc_padj$Annotation <- factor(lgfc_padj$Annotation, levels = rev(splevel))

plot2 <- ggplot(lgfc_padj[lgfc_padj$Annotation %in% dec, ], 
                aes(x=1, Annotation, group = dataset)) + 
  #geom_tile()+
  geom_col(aes(fill=value), position = position_dodge(),
           color = "black")+
  geom_text(aes(label = label, x = 0.5), size=3, 
            hjust = "middle",vjust = "middle",
            position = position_dodge(0.9))+
  theme_void()+
  labs(fill = "Log2FC")+
  scale_fill_gradient2(low = "#0571B0", mid = "white",
                       high = "#CC0824", 
                       midpoint = 0, limits = c(-4,4))

dta$variable <- factor(dta$variable, levels = rev(splevel))
dta$Group <- factor(dta$Group, levels = c("Ctrl", "PCOS"))
dta$subgroup <- factor(dta$subgroup, levels = paste(rep(c("mice", "QJ","CX"), each = 2),
                                                    c("Ctrl", "PCOS"),sep = "_"))
dta_dec <- dta[dta$variable %in% dec, ]
#dta_dec$dataset <- factor(dta_dec$dataset, levels = rev(c("CX","QJ","mice")))
plot1 <- ggplot(dta_dec)+
  geom_col(aes(y=variable,x=16,group = dataset, fill = dataset),
           position = position_dodge())+
  scale_fill_manual(values = c("CX"="#F4F4F4",
                               "QJ"="#E2E2E2",
                               "mice"="#D1D1D1"))+
  new_scale_fill()+
  geom_bar(aes(y=variable, x=value, 
               group = subgroup, fill = Group),
           stat = "summary", fun="mean", color = "black",
           position = position_dodge(), inherit.aes = F)+
  stat_summary(aes(y=variable, x=value, 
                   group = subgroup),
               fun.data = "mean_se", geom = "errorbar", color = "black",
               width = 0.5, position = position_dodge(.9))+
  scale_fill_manual(values = c("PCOS"="#FF6E6E",
                               "Ctrl"="#6298CE"))+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"))+
  xlab("DESeq counts")+
  ylab(NULL)+
  scale_x_reverse()

library(aplot)
plot1 %>% insert_right(plot2, width = 0.05)
#ggsave("figS5c_decrease_spearman0.05.pdf", width = 8, height = 8)
ggsave("figS5c_decrease.pdf", width = 9.6, height = 12)

plot3 <- ggplot(dta[dta$variable %in% inc, ])+
  geom_col(aes(y=variable,x=20,group = dataset, fill = dataset),
           position = position_dodge())+
  scale_fill_manual(values = c("CX"="#F4F4F4",
                               "QJ"="#E2E2E2",
                               "mice"="#D1D1D1"))+
  new_scale_fill()+
  geom_bar(aes(y=variable, x=value, 
               group = subgroup, fill = Group),
           stat = "summary", fun="mean", color = "black",
           position = position_dodge(), inherit.aes = F)+
  stat_summary(aes(y=variable, x=value, 
                   group = subgroup),
               fun.data = "mean_se", geom = "errorbar", color = "black",
               width = 0.5, position = position_dodge(.9))+
  scale_fill_manual(values = c("PCOS"="#FF6E6E",
                               "Ctrl"="#6298CE"))+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"))+
  xlab("DESeq counts")+
  ylab(NULL)+
  scale_y_discrete(position = "right")
plot4 <- ggplot(lgfc_padj[lgfc_padj$Annotation %in% inc, ], 
                aes(x=1, Annotation, group = dataset)) + 
  #geom_tile()+
  geom_col(aes(fill=value), position = position_dodge(),
           color = "black")+
  geom_text(aes(label = label, x = 0.5), size=3, 
            hjust = "middle",vjust = "middle",
            position = position_dodge(0.9))+
  theme_void()+
  labs(fill = "Log2FC")+
  scale_fill_gradient2(low = "#0571B0", mid = "white",
                       high = "#CC0824", 
                       midpoint = 0, limits = c(-4,4))
plot3 %>% insert_left(plot4, width = 0.05)
# ggsave("figS5c_increase_spearman0.05.pdf", width = 7, height = 6)
ggsave("figS5c_increase.pdf", width = 7, height = 7)

############################## sankeyNetwork
library(webshot);library(sankeyD3)
filtered2 <- subset(result3, q0.2 %in% c("CX_QJ", "CX_QJ_mice") &
                      flavonoids.Class %in% c("3_Chalcones and dihydrochalcones", "4_flavonoid metabolites")) 
bac_info <- distinct(filtered2[, c("Var1", "q0.2", "direction")]) %>% droplevels()
bac_info$Group <- paste0()
met_info <- distinct(filtered2[, c("Var2", "flavonoids.Class")]) %>% droplevels()
met_info <- met_info[order(met_info$flavonoids.Class), ]
met_info$Group <- "met"
mb_level <- as.character(met_info$Var2)
corr <- filtered2[, c("Var1", "Var2", "spearman.r", "spearman.p")]
corr <- subset(corr,Var2 %in% as.character(met_info$Var2) & 
                 Var1 %in% as.character(bac_info$Var1))
corr$value <- ifelse(corr$spearman.r>0, 1, -1)
corr$group <- ifelse(corr$spearman.p < 0.05, ifelse(corr$spearman.r>0, "pos","neg"),
                     ifelse(corr$spearman.r>0, "pos_trend", "neg_trend"))
corr$trend <- ifelse(corr$spearman.r>0, "pos", "neg")
rownames(bac_info) <- bac_info$Var1
bac_info <- bac_info[splevel, ]
rownames(met_info) <- met_info$Var2
met_info <- met_info[mb_level, ]
b1 <- rownames(bac_info)[bac_info$direction %in% "decrease"]
b2 <- rownames(bac_info)[bac_info$direction %in% "increase"]
link1 <- data.frame(source = rep(b1, each = length(mb_level)),
                   target = rep(mb_level, length(b1)))
link1$add <- paste(link1$source, link1$target, sep = "_")
link2 <- data.frame(source = rep(mb_level, length(b2)),
                    target = rep(b2, each = length(mb_level)))
link2$add <- paste(link2$target, link2$source, sep = "_")
corr$add <- paste(corr$Var1, corr$Var2, sep = "_")
link1 <- merge(link1, corr[, c("add", "value", "group", "trend")],
               by = "add", all.x = T, sort = F)
link2 <- merge(link2, corr[, c("add", "value", "group", "trend")], 
               by = "add", all.x = T, sort = F)
link <- rbind(link1, link2)
link$add <- NULL

node <- data.frame(Node=c(bac_info$Var1 %>% as.character(),
                          met_info$Var2 %>% as.character()), 
                   Group=c(bac_info$direction,met_info$Group),
                   #ID=0:(nrow(bac_info)+nrow(met_info)-1), 
                   tax=c(bac_info$q0.2,met_info$flavonoids.Class)#,
                   #width=c(b,a)
)
rownames(node) <- node$Node
node <- node[c(b1, mb_level, b2), ]
node$ID <- 0:(nrow(node)-1)

link <- merge(link, node[, c("Node", "ID")], 
              by.x = "target", by.y = "Node", 
              all.x = T, sort = F)
colnames(link) <- gsub("ID", "IDtarget", colnames(link))
link <- merge(link, node[, c("Node", "ID")], 
              by.x = "source", by.y = "Node", 
              all.x = T, sort = F)
colnames(link) <- gsub("^ID$", "IDsource", colnames(link))

node$Node <- factor(node$Node,levels = node$Node)
node$ID <- factor(node$ID,levels = node$ID)
link$value <- abs(link$value)
my_color <- 'd3.scaleOrdinal() .domain(["neg", "neg_trend", "pos", "pos_trend", "CX_QJ_mice", "CX_QJ","3_Chalcones and dihydrochalcones", "4_flavonoid metabolites"]) .range(["#E61F19", "#EB8788", "#19499C", "#86CCDC","#74C0FC", "#B9DDD7","#BFE8F7", "#81C1E0"])'
p <- sankeyNetwork(Links = link, Nodes = node,LinkGroup = 'group',
                   #linkColor = 
                   NodeID='Node',NodeGroup = 'tax',#sinksRight = T,iterations = 0,
                   #colourScale=JS("d3.scaleOrdinal(d3.schemeCategory10);"),
                   colourScale=my_color,
                   Source = 'IDsource', Target = 'IDtarget', 
                   Value = 'value', nodeCornerRadius = 5,#nodePadding = 3,
                   dragX = F,dragY = T,fontSize = 12,
                   nodeWidth = 15,width = 800,showNodeValues=F,
                   fontFamily = 'Arial',height = 800)
saveNetwork(p,file = "Sankey.html")
#gw <- getwd()
webshot("Sankey.html",file = 'Sankey.pdf',delay = 1)
save(link, node, mb_level, splevel, file = "data_for_sankey.RData")
# library(ggplotify)
# as.ggplot(p)
