library(data.table)
cx_new <- fread("Deseq.group_P_vs_C.txt")
diff_sp <- cx_new[cx_new$padj < 0.2, ]
mg <- fread("transformed_counts.txt")
meta <- fread("PCOS_CX.meta.txt")
mb <- fread("neg_pos_KEGG.CAS.HMDB_log10.txt")

sid <- intersect(colnames(mb)[-1], colnames(mg)[-1])
mb2 <- mb[, sid, with = F] %>% t() %>% as.data.frame()
colnames(mb2) <- mb$V1
mg2 <- mg[, sid, with = F] %>% t() %>% as.data.frame()
colnames(mg2) <- mg$V1
mg2 <- mg2[, diff_sp$V1]
library(psych)
res <- corr.test(mg2, mb2, method = "spearman", adjust = "fdr")
p <- res$p[, apply(res$p, 2,function(x) any(x<0.05))]
p <- p[apply(p, 1,function(x) any(x<0.05)), ]

anno <- fread("KEGG.CAS.HMDB_log10_ort3.vip_anno_cv_class.txt")
anno_flavon <- xlsx::read.xlsx("previous/组学差异标志物.xlsx", sheetIndex = 1)
anno_flavon$ms_mz <- paste(anno_flavon$MS2.name, round(anno_flavon$mz, 3), sep = "_")
anno$ms_mz <- paste(anno$MS2.name, round(anno$mz, 3), sep = "_")
anno2 <- merge(anno[, -c(7,8)], anno_flavon[, c(3, 6,7,13,17:31, 91)],
              by="ms_mz")
anno3 <- anno2[anno2$flavonoids.Class %in% c("1_Flavonoid glycosides",
                                 "2_Flavans", "2_Flavones",
                                 "3_Chalcones and dihydrochalcones",
                                 "4_flavonoid metabolites"), ]
p <- res$p[, anno3$V1] %>% melt()
r <- res$r[, anno3$V1] %>% melt()
identical(p$Var1, r$Var1)
identical(p$Var2, r$Var2)
colnames(r)[3] <- "spearman.r"
result <- data.frame(r, spearman.p=p$value)
result2 <- merge(result, anno3[, c(2,4,19:23, 32,34,36:38)],
                by.x = "Var2", by.y = "V1")

cx_qj_mice <- fread("mergeDESeq.CXnew_QJ_mice.human_species.txt")
cx_qj_mice$Annotation2 <- gsub("'", ".",cx_qj_mice$Annotation)
result3 <- merge(result2, cx_qj_mice, by.x = "Var1",
                 by.y = "Annotation2")

cytoscape1 <- fread("cytoscape_1.txt")
cytoscape23 <- fread("cytoscape_2-3.txt")
cytoscape23_numb <- table(cytoscape23$Annotation, cytoscape23$dataset) %>% 
  as.data.frame() %>% dcast(Var1~Var2)
cytoscape23_numb$q0.2 <- ifelse(cytoscape23_numb$CX==0, "QJ_mice",
                                ifelse(cytoscape23_numb$QJ==0, "CX_mice",
                                       ifelse(cytoscape23_numb$mice==0, "CX_QJ", "CX_QJ_mice")))
table(cytoscape23_numb$q0.2)
cytoscape23_info <- merge(cytoscape23_numb, 
                          distinct(cytoscape23[, c("Annotation", "direction")]),
                          by.x = "Var1", by.y = "Annotation")
sp_info <- rbind(data.frame(Var1=cytoscape1$Annotation, 
                            q0.2=cytoscape1$dataset,
                            direction=cytoscape1$direction),
                 cytoscape23_info[, c("Var1", "q0.2", "direction")])
cx_qj_mice_info <- merge(cx_qj_mice, sp_info,
                         by.x="Annotation", by.y = "Var1",
                         all.x = T)
fwrite(cx_qj_mice_info, "mergeDESeq.CXnew_QJ_mice.human_species_q0.2.txt", sep = "\t")

result3 <- merge(result2, cx_qj_mice_info, by.x = "Var1",
                 by.y = "Annotation2")
fwrite(result3, "diffsp108_mbFlavonoids34_cor.txt", sep = "\t")

