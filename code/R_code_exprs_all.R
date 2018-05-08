### Expresseion PCA

tmp <- genes_expressed_THP1[which(genes_expressed_THP1$gene %in% Common_exprss_and_probe),]
Exprs_FPKM_THP1 <- tmp[order(tmp[,1]),][,c(1,2,3)]

tmp<- genes_expressed_TAV[which(genes_expressed_TAV$gene %in% Common_exprss_and_probe),]
Exprs_FPKM_TAV <- tmp[order(tmp[,1]),][,c(1,3,4)]

tmp <- genes_expressed_SMC[which(genes_expressed_SMC$gene%in% Common_exprss_and_probe),]
Exprs_FPKM_SMC <- tmp[order(tmp[,1]),][,c(1,3,4)]

remove_SMC <- Exprs_FPKM_SMC$gene[duplicated(Exprs_FPKM_SMC$gene)]
remove_TAV <- Exprs_FPKM_TAV$gene[duplicated(Exprs_FPKM_TAV$gene)]
remove_these <- union (remove_SMC,remove_TAV)

tmp <- Exprs_FPKM_THP1[which(!Exprs_FPKM_THP1$gene %in% remove_these),]
Exprs_FPKM_THP1 <- tmp[order(tmp[,1]),][,c(1,2,3)]

tmp <- Exprs_FPKM_TAV[which(!Exprs_FPKM_TAV$gene %in% remove_these),]
Exprs_FPKM_TAV <- tmp[order(tmp[,1]),][,c(1,2,3)]

tmp <- Exprs_FPKM_SMC[which(!Exprs_FPKM_SMC$gene %in% remove_these),]
Exprs_FPKM_SMC <- tmp[order(tmp[,1]),][,c(1,2,3)]


head (Exprs_FPKM_SMC)
tail (Exprs_FPKM_SMC)

head (Exprs_FPKM_TAV)
tail (Exprs_FPKM_TAV)

head (Exprs_FPKM_THP1)
tail (Exprs_FPKM_THP1)

All_interaction_combine <- cbind(Exprs_FPKM_THP1,Exprs_FPKM_TAV[,c(2,3)],Exprs_FPKM_SMC[,c(2,3)])
colnames (All_interaction_combine) <- c("Genes","Exprs_RPKM_Rep1_THP1","Exprs_RPKM_Rep2_THP1",
                                        "Exprs_RPKM_Rep1_HAEC","Exprs_RPKM_Rep2_HAEC",
                                        "Exprs_RPKM_Rep1_HSMC","Exprs_RPKM_Rep2_HSMC")


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
  labs(title ="PCA")

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
  ggtitle("\nPCA plot of  expressed genes \n in all three cell types and replicates \n")


