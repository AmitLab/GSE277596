
mtpath <- 'meta_SPID_PSA_v8_Blood.txt'
excl <- c(333, 334)
mat_id <- 'PsA_v8_Blood'
col_pop_Imm <- c('chartreuse4', 'aquamarine', 'darkorange', 'purple', 'burlywood', 'darkgreen', 'deepskyblue4', 
 'firebrick1', 'chocolate4', 'mediumseagreen', 'magenta', 
 'darkturquoise', 'lightsalmon', 'darkred', 'deeppink', "darkblue",
 'grey11', 'grey55', 'grey99')
order_pop_Imm <- c("Mo_Blood", "Mo_CD16", "Mo", "Mf", "mDC", "cDC1", "cDC2", 
 "NK", "pDC", "B", "Plasma", 
 "Treg_Imm", "T_Effector_GNLY", "T_Effector", "Treg", "T_Naive",
 "Doublets", "KRT19", "UN")
names(col_pop_Imm) <- order_pop_Imm

sc_mat <- scdb_mat(mat_id)
sc_cl<-scdb_mc("test_mc")
sc_2d <- scdb_mc2d("test_mc_2dproj")
cell_stats <- sc_mat@cell_metadata
cell_stats <- merge(cell_stats, sc_cl@mc, by=0)
rownames(cell_stats) <- cell_stats[, "Row.names"]
cell_stats <- cell_stats[, -1]
cell_stats <- cell_stats[! cell_stats[, "y"] %in% excl, ]
cells <- sample(names(sc_cl@mc[!sc_cl@mc %in% excl]))

meta <- read.table(mtpath, header=T, row.names = 1, sep="\t")
meta <- meta[intersect(rownames(meta), rownames(sc_cl@n_bc)),]

cells_Blood <- rownames(cell_stats[cell_stats[, "Diagnosis"] == "Control", ])
cells_PsA <- sample(rownames(cell_stats[cell_stats[, "Diagnosis"] == "PsA", ]), length(cells_Blood))
cells_OA <- sample(rownames(cell_stats[cell_stats[, "Diagnosis"] == "OA", ]), length(cells_Blood))
png(file="colors_2D_1_BSO.png", width=2000, height= 2000)
par(mfrow=c(2, 2))
plot(sc_2d@sc_x[cells_Blood], sc_2d@sc_y[cells_Blood], col= as.character(col_pop_Imm[ann[sc_cl@mc[cells_Blood],"Populations"]]), pch=16, ylab="", xlab=paste(length(cells_Blood), "cells from Blood"), xaxt='n', yaxt='n', cex=1, bty="n")
plot(sc_2d@sc_x[cells_PsA], sc_2d@sc_y[cells_PsA], col= as.character(col_pop_Imm[ann[sc_cl@mc[cells_PsA],"Populations"]]), pch=16, ylab="", xlab=paste(length(cells_PsA), "cells from PsA"), xaxt='n', yaxt='n', cex=1, bty="n")
plot(sc_2d@sc_x[cells_OA], sc_2d@sc_y[cells_OA], col= as.character(col_pop_Imm[ann[sc_cl@mc[cells_OA],"Populations"]]), pch=16, ylab="", xlab=paste(length(cells_OA), "cells from OA"), xaxt='n', yaxt='n', cex=1, bty="n")
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
legend(x = "topleft", legend=names(col_pop_Imm[pop]), col=col_pop_Imm[pop], pch=16, cex=2)
dev.off()


pid2use <- c(pid_Blood, pid_OA, pid_PSA)
meta_CD45 <- meta[meta$Gating == 'CD45' & meta$Tissue %in% c('SF', 'Blood') & meta$PID %in% pid2use, ]
pIDs <- as.vector(unique(meta_CD45$PID))
nbcs <- sc_cl@n_bc
anns <- ann[setdiff(rownames(ann), excl), ]
ctps <- sort(as.vector(unique(anns$Populations)))
ctps <- setdiff(ctps, c('UN'))
res1 <- foreach(i = ctps, .combine = cbind) %do% rowSums(nbcs[, anns[anns$Populations == i, ]$mcID, drop =F])
colnames(res1) <- ctps
head(meta_CD45[meta_CD45[, "PID"] %in% pid_NAIVE, ])

pid2Diagnosis <- unique(meta_CD45[, c("PID", "Diagnosis")]) 
rownames(pid2Diagnosis) <- pid2Diagnosis[, "PID"]

res2 <- foreach(i = pid2use, .combine = rbind) %do% colSums(res1[rownames(meta_CD45[meta_CD45$PID == i, ]), , drop =F])
rownames(res2) <- pid2use
head(res2, n=2)
res2<- res2[order(pid2Diagnosis[rownames(res2), "Diagnosis"]), ]
res2<- res2[rowSums(res2) > 100, ]
res2_pct <-round(t(apply(res2, 1, function(x) x/sum(x))) * 100, 2)
head(res2_pct, 2)
write.table(file="res2_pct.txt", res2_pct, sep="\t", quote=F)

pid_Blood <- intersect(pid_Blood, rownames(res2_pct))
pid_OA <- intersect(pid_OA, rownames(res2_pct))
pid_PSA <- intersect(pid_PSA, rownames(res2_pct))

pdf(file=paste0(outdir, "/", fig, "_box_plot_Diagnosis_in_group.pdf"), width=15, height=9, paper='special', useDingbats = FALSE )
par(mfrow=c(2, 5))
group = "B"
for( group in colnames(res2_pct)){
	names <- c(rep("Blood", times=length(pid_Blood)), 
	rep("OA", times=length(pid_OA)), 
	rep("PsA", times=length(pid_PSA)))

	value <- res2_pct[c(pid_Blood, pid_OA, pid_PSA), group]
	data <- data.frame(names, value)
	data$names <- factor(data$names , levels=c("Blood", "OA", "PsA"))
	#png(file=paste0(outdir, "/", id, "_", gsub(pattern="/", replacement="_", x=group) , "_PsA_SF_vs_Blood_box_dist.png"), width=800, height= 800)
	boxplot(data$value ~ data$names, outline=F, frame=F, ylab="Cells percent", xlab=group, ylim = c(0, 1.3*max(13, res2_pct[, group])), col=c("gray80", "gray60", "gray40", "gray20")) 
	ytop <- 1.2*max(13, res2_pct[, group])
	res1 <- wilcox.test(data[data[, "names"] == "Blood", "value"], data[data[, "names"] == "OA", "value"], paired= F)$p.value
	segments(x0 = 1.25, y0 = ytop, x1 = 1.75, y1 = ytop)
	text(1.5, ytop*1.02, labels = pv2stars(res1))
	res1 <- wilcox.test(data[data[, "names"] == "Blood", "value"], data[data[, "names"] == "PsA", "value"], paired= F)$p.value
	segments(x0 = 1.25, y0 = ytop*0.95, x1 = 2.75, y1 = ytop*0.95)
	text(2, ytop*1.02*0.95, labels = pv2stars(res1))
	res1 <- wilcox.test(data[data[, "names"] == "OA", "value"], data[data[, "names"] == "PsA", "value"], paired= F)$p.value
	segments(x0 = 2.25, y0 = ytop*0.9, x1 = 2.75, y1 = ytop*0.9)
	text(2.5, ytop*1.02*0.9, labels = pv2stars(res1)) 

	levelProportions <- table(data$names)/nrow(data)
	levelProportions
	mycol <- rep("lightblue1", times=4)
	mylevels <- levels(data$names)
	mylevels
	for(i in 1:length(mylevels)){
		thislevel <- mylevels[i]
		thisvalues <- data[data$names==thislevel, "value"]
		myjitter <- jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/2)
		points(myjitter, thisvalues, pch = 21, cex=1.5, col="gray", bg=mycol[i]) 
	}
}
dev.off()

