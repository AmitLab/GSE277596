
pid2use <- c(pid_Blood, pid_OA, pid_NAIVE, pid_NONNAIVE)
meta_CD45 <- meta[meta$Gating == 'CD45' & meta$Tissue %in% c('SF', 'Blood') & meta$PID %in% pid2use, ]
pIDs <- as.vector(unique(meta_CD45$PID))
nbcs <- mc@n_bc
anns <- ann[setdiff(rownames(ann), excl), ]
ctps <- sort(as.vector(unique(anns$Populations)))
ctps <- setdiff(ctps, c('UN'))
res1 <- foreach(i = ctps, .combine = cbind) %do% rowSums(nbcs[, anns[anns$Populations == i, ]$mcID, drop =F])
colnames(res1) <- ctps
head(meta_CD45[meta_CD45[, "PID"] %in% pid_NAIVE, ])

meta_CD45[, "grp"] <- "Control"
meta_CD45[meta_CD45[, "PID"] %in% pid_OA, "grp"] <- "OA"
meta_CD45[meta_CD45[, "PID"] %in% pid_NAIVE, "grp"] <- "Naive"
meta_CD45[meta_CD45[, "PID"] %in% pid_NONNAIVE, "grp"] <- "NonNaive"
meta_CD45$grp <- factor(meta_CD45$grp , levels=c("Control", "OA", "Naive", "Non Naive"))
table(meta_CD45$grp)
pid2grp <- unique(meta_CD45[, c("PID", "grp")]) 
rownames(pid2grp) <- pid2grp[, "PID"]

res2 <- foreach(i = pid2use, .combine = rbind) %do% colSums(res1[rownames(meta_CD45[meta_CD45$PID == i, ]), , drop =F])
rownames(res2) <- pid2use
head(res2, n=2)
res2<- res2[order(pid2grp[rownames(res2), "grp"]), ]

res2<- res2[rowSums(res2) > 100, ]

res2_pct <-round(t(apply(res2, 1, function(x) x/sum(x))) * 100, 2)
head(res2_pct, 2)

col_pop_Imm[colnames(res2_pct)]
space <- c(rep(0.1, times=length(intersect(pid_Blood, rownames(res2_pct)))), 1, 
 rep(0.1, times=length(intersect(pid_OA, rownames(res2_pct)))-1), 1, 
 rep(0.1, times=length(intersect(pid_NAIVE, rownames(res2_pct)))-1), 1, 
 rep(0.1, times=length(intersect(pid_NONNAIVE, rownames(res2_pct)))-1)) #pid_Blood, pid_OA, pid_NAIVE, pid_NONNAIVE
barplot(t(res2_pct), col=col_pop_Imm[colnames(res2_pct)], space=space, las=2)

pdf(file=paste0(outdir, "/", "fig2_low_level_cells_prop.pdf"), width=15, height=7, paper='special', useDingbats = FALSE )
par(mar = c(8, 5, 1, 1))
 barplot(t(res2_pct), col=col_pop_Imm[colnames(res2_pct)], space=space, las=2, ylab='Fraction')
dev.off()