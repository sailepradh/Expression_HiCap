#####################################################################################################################################################
setwd("/Volumes/Work_drive/prj/Expression_HiCap/data/raw_external/")
Genes_target_THP1 <- read.table("/Volumes/Work_drive/prj/Expression_HiCap/data/raw_external/THP1_PD_SP", quote="\"", comment.char="")
Genes_target_THP1$V1 <- as.character(Genes_target_THP1$V1)
Genes_target_TAV <- read.table("/Volumes/Work_drive/prj/Expression_HiCap/data/raw_external/TAV_PD_SP", quote="\"", comment.char="")
Genes_target_TAV$V1 <- as.character(Genes_target_TAV$V1)
Genes_target_SMC <- read.table("/Volumes/Work_drive/prj/Expression_HiCap/data/raw_external/SMC_PD_SP", quote="\"", comment.char="")
Genes_target_SMC$V1 <- as.character(Genes_target_SMC$V1)

Comm_probe_genes = Reduce(intersect,  list(Genes_target_THP1$V1, 
                                           Genes_target_TAV$V1,
                                           Genes_target_SMC$V1))

write.table(unique(Comm_probe_genes),"Common_ProbeGenes_targetted.txt",sep="\t", quote=FALSE,
            row.names=FALSE, col.names = FALSE)

#####################################################################################################################################################
setwd("/Volumes/Work_drive/prj/Expression_HiCap/data/raw_external/expression/")
genes_expressed_THP1 <- read.delim("/Volumes/Work_drive/prj/Expression_HiCap/data/raw_external/expression/genes_expressed_THP1.txt")
colnames (genes_expressed_THP1) <- c("gene","FPKM_THP1_Rep_1","FPKM_THP1_Rep_2")
genes_expressed_THP1$gene <- as.character(genes_expressed_THP1$gene)


genes_expressed_TAV <- read.delim("/Volumes/Work_drive/prj/Expression_HiCap/data/raw_external/expression/genes_expressed_TAV.txt")
colnames (genes_expressed_TAV) <- c("gene","transcipt_ID","FPKM_TAV_Rep_1","FPKM_TAV_Rep_2")
genes_expressed_TAV$gene <- as.character(genes_expressed_TAV$gene)


genes_expressed_SMC <- read.delim("/Volumes/Work_drive/prj/Expression_HiCap/data/raw_external/expression/genes_expressed_SMC.txt")
colnames (genes_expressed_SMC) <- c("gene","transcript_ID","FPKM_SMC_Rep_1","FPKM_SMC_Rep_2")
genes_expressed_SMC$gene <- as.character(genes_expressed_SMC$gene)


Comm_expressed_genes = Reduce(intersect,  list(genes_expressed_THP1$gene, 
                                               genes_expressed_TAV$gene,
                                               genes_expressed_SMC$gene))

write.table(unique(Comm_expressed_genes),"Common_Expressed_gene.txt",sep="\t", quote=FALSE,
            row.names=FALSE, col.names = FALSE)

#####################################################################################################################################################
setwd("/Volumes/Work_drive/prj/Expression_HiCap/data/raw_external/expression/")
genes_Unexpressed_THP1 <- read.delim("/Volumes/Work_drive/prj/Expression_HiCap/data/raw_external/expression/genes_unexpressed_THP1.txt")
colnames (genes_Unexpressed_THP1) <- c("gene","FPKM_THP1_Rep_1","FPKM_THP1_Rep_2")
genes_Unexpressed_THP1$gene <- as.character(genes_Unexpressed_THP1$gene)


genes_Unexpressed_TAV <- read.delim("/Volumes/Work_drive/prj/Expression_HiCap/data/raw_external/expression/genes_unexpressed_TAV.txt")
colnames (genes_Unexpressed_TAV) <- c("gene","transcipt_ID","FPKM_TAV_Rep_1","FPKM_TAV_Rep_2")
genes_Unexpressed_TAV$gene <- as.character(genes_Unexpressed_TAV$gene)


genes_Unexpressed_SMC <- read.delim("/Volumes/Work_drive/prj/Expression_HiCap/data/raw_external/expression/genes_unexpressed_SMC.txt")
colnames (genes_Unexpressed_SMC) <- c("gene","transcript_ID","FPKM_SMC_Rep_1","FPKM_SMC_Rep_2")
genes_Unexpressed_SMC$gene <- as.character(genes_Unexpressed_SMC$gene)


Comm_Unexpressed_genes = Reduce(intersect,  list(genes_Unexpressed_THP1$gene, 
                                               genes_Unexpressed_TAV$gene,
                                               genes_Unexpressed_SMC$gene))

write.table(unique(Comm_Unexpressed_genes),"Common_UnExpressed_gene.txt",sep="\t", quote=FALSE,
            row.names=FALSE, col.names = FALSE)

Dubious_genes_expression = intersect(Comm_expressed_genes,Comm_Unexpressed_genes ) 
#####################################################################################################################################################
setwd("/Volumes/Work_drive/prj/Expression_HiCap/data/raw_external/")

Common_exprss_and_probe <- intersect(Comm_probe_genes,Comm_expressed_genes)
Common_Unexprss_and_probe <- intersect(Comm_probe_genes,Comm_Unexpressed_genes)


write.table(unique(Common_exprss_and_probe),"Common_Expr_probe_gene.txt",sep="\t", quote=FALSE,
            row.names=FALSE, col.names = FALSE)
write.table(unique(Comm_Unexpressed_genes),"Common_UnExpr_probe_gene.txt",sep="\t", quote=FALSE,
            row.names=FALSE, col.names = FALSE)
#####################################################################################################################################################
######################################################################################################################################################

tmp <- Genes_target_THP1[which(Genes_target_THP1$V1 %in% Common_exprss_and_probe),]
SP_FPKM_THP1 <- tmp[order(tmp[,1]),][,c(1,3,4)]

tmp<- Genes_target_TAV[which(Genes_target_TAV$V1 %in% Common_exprss_and_probe),]
SP_FPKM_TAV <- tmp[order(tmp[,1]),][,c(1,3,4)]

tmp <- Genes_target_SMC[which(Genes_target_SMC$V1 %in% Common_exprss_and_probe),]
SP_FPKM_SMC <- tmp[order(tmp[,1]),][,c(1,3,4)]

head (SP_FPKM_SMC)
tail (SP_FPKM_SMC)

head (SP_FPKM_TAV)
tail (SP_FPKM_TAV)

head (SP_FPKM_THP1)
tail (SP_FPKM_THP1)

All_interaction_combine <- cbind(SP_FPKM_THP1,SP_FPKM_TAV[,c(2,3)],SP_FPKM_SMC[,c(2,3)])
colnames (All_interaction_combine) <- c("Genes","SP_RPKM_Rep1_THP1","SP_RPKM_Rep2_THP1",
                                        "SP_RPKM_Rep1_HAEC","SP_RPKM_Rep2_HAEC",
                                        "SP_RPKM_Rep1_HSMC","SP_RPKM_Rep2_HSMC")


row_sub = apply(All_interaction_combine[2:7], 1, function(row) all(row !=0 ))
All_interaction_combine_clean <- All_interaction_combine[row_sub,]

library(factoextra)
log_TF <- log10(All_interaction_combine_clean[,2:7])
head (log_TF)
tr <- t(log_TF)
colnames(tr) <- All_interaction_combine_clean[, 1]

ir.pca <- prcomp(tr  ,
                 center = TRUE,
                 scale. = TRUE)
ir.genes <- row.names(tr)
cal= c("mTHP1_Rep","mTHP1_Rep",
       "HAEC_Rep", "HAEC_Rep",
       "HSMC_Rep","HSMC_Rep")

fviz_eig(ir.pca)
print(ir.pca)
plot(ir.pca, type = "l")
summary(ir.pca)
plot(ir.pca$x[,1:2],xlab = "PC1", ylab = "PC2", col =as.factor(cal), pch =19)
biplot(ir.pca)


#fviz_pca_biplot(ir.pca,
#                col.var = "#2E9FDF", # Variables color
#                col.ind = "#696969", # Individuals color
#                repel = TRUE     # Avoid text overlapping
#)

lab= c("mTHP1","mTHP1","HAEC","HAEC","HSMC","HSMC")

p = fviz_pca_ind(ir.pca, pointsize = 5, habillage=lab, addEllipses=TRUE, ellipse.level=0.95)+
  labs(title ="PCA")+
  xlim(-100,100) + ylim (-100, 100)

p+ scale_color_brewer(palette="Dark2") +
  theme(legend.text=element_text(size=20),
        legend.position = 'right',
        legend.key.size = unit(1.5, 'lines'),
        axis.title.x = element_text(size=30,face="bold", margin = margin(t = 20, r = 0, b =0, l = 0)),
        axis.text.x = element_text(size=20),
        axis.title.y = element_text(size=30,face="bold",margin = margin(t = 0, r = 20, b =0, l = 0)),
        axis.text.y = element_text(size=20),
        plot.title = element_text(colour="grey20",size=30,hjust=0.5),
        plot.margin = unit(c(2, 2, 2, 2), "cm"))+
  ggtitle("\nPCA plot of interaction profiles of expressed genes \n in all three cell types and replicates \n")
