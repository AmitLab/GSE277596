
mo_mean <- mo_mean[names(head(sort(apply(mo_mean, 1, max), decreasing=T), 6000)), ]
g1 <- intersect(pid_PSA, colnames(mo_mean))
g2 <- intersect(pid_OA, colnames(mo_mean))
x = rowSums(mo_mean[, g1]) / length(g1) * min(length(g1), length(g2))
y = rowSums(mo_mean[, g2]) / length(g2) * min(length(g1), length(g2)) 

outtable1 <- as.data.frame(rownames(mo_mean))
outtable1[, "pV"] <- apply(mo_mean, 1, FUN = function(X){t <- wilcox.test(x= as.numeric(X[g1]), y=as.numeric(X[g2]), paired= F)$p.value})
outtable1[, "Aside"] <- log(x+0.1, 2)
outtable1[, "Bside"] <- log(y+0.1, 2)
outtable1[, "foldp1"] <- outtable1[, "Aside"] - outtable1[, "Bside"]
rownames(outtable1) <- outtable1[, 1]
outtable1 <- outtable1[, -1]
outtable1[, "pVNANA"] <- outtable1[, "pV"]
outtable1[is.na(outtable1[, "pVNANA"]), "pVNANA"] <- 1
outtable1[, "pVCor"] <- p.adjust(outtable1[, "pVNANA"] , "BH")
tmp <- outtable1[, c("Aside", "Bside", "foldp1", "pV", "pVCor")]
colnames(tmp) <- c("PsA_mean", "OA_mean", "fold_change", "pV", "pVCor")
write.table(file="fig3_Mo_PsA_vs_OA.txt", tmp, sep="\t", quote=F)
xmax <- 1.2*max(abs(outtable1$foldp1))
difgenes <- intersect(rownames(outtable1[abs(outtable1$foldp1) > 1, , drop=F]) , rownames(outtable1[outtable1$pVCor < 0.05, , drop=F]))


showgenes <- c("PSME2", "OPTN", "TAP1", "GBP1", "CXCL10", "IFIT3", "NLRC5", "IL18BP", "MX1", "STAT1", "ANKH", "ECM1", "ENHO")
showgenes <- intersect(showgenes, rownames(outtable1))
showgenes <- setdiff(showgenes, bl)
pdf(file = paste(id, "_volcano_Mo_PsA_vs_OA_gl5.pdf", sep=""), width=10, height=10, paper='special', useDingbats = FALSE )
par(mar = c(5, 5, 1, 1))
plot(outtable1$foldp1, -log(outtable1$pVNANA, 10), cex = 0.5, xlim = c(-xmax, xmax), pch = 16, xlab="OA <- fold change -> PsA", ylab="-log 10 of FDR PValue", col="gray", bty="n", cex.lab=2.1, cex.axis=2.1)
points(outtable1[difgenes, "foldp1"], -log(outtable1[difgenes, "pVNANA"], 10), xlim = c(-xmax, xmax), pch=16, cex = 2.1, col = "red")
textplot( outtable1[showgenes, "foldp1"], -log(outtable1[showgenes, "pVNANA"], 10), cex=2.1, words = showgenes, new = F)
abline(v=1, lty="dotted")
abline(v=-1, lty="dotted")
dev.off()
