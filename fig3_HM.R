library(reshape2)
library(foreach)
library(cowplot)
library(ComplexHeatmap)
library(circlize)

read_exp <- function(fpath, pop){
    dat <- read.table(fpath, sep = ' ', header = T, row.names = 1)
    dat
}

Mo <- read_exp('/ver_20240215/Mo_mean_lmt.txt', 'Mo')
Mf <- read_exp('/ver_20240215/Mf_mean_lmt.txt', 'Mf')
cDC2 <- read_exp('/ver_20240215/cDC2_mean_lmt.txt', 'cDC2')
cDC1 <- read_exp('/ver_20240215/cDC1_mean_lmt.txt', 'cDC1')

OA <- c('OA_001','OA_002','OA_003','OA_004','OA_006','OA_007','OA_008','OA_009')
psa <- grep("PsA", colnames(Mo), value=T)
psa <- setdiff(psa, c("PsA_003.1","PsA_008.1","PsA_006.1","PsA_012.1","PsA_022.1","PsA_033","PsA_033.2"))

extract_dat <- function(dat, pop, OA, psa){
    dat_OA <- dat[, colnames(dat) %in% OA]
    dat_psa <- dat[, colnames(dat) %in% psa]
    colnames(dat_OA) <- paste0('0_OA|', colnames(dat_OA))
    colnames(dat_psa) <- paste0('1_PsA|', colnames(dat_psa))
    res <- cbind(dat_OA, dat_psa)
    colnames(res) <- paste0(colnames(res), '|', pop)
    res
}

dat_Mo <- extract_dat(Mo, '0_Mo', OA, psa)
dat_Mf <- extract_dat(Mf, '1_Mf', OA, psa)
dat_cDC2 <- extract_dat(cDC2, '2_cDC2', OA, psa)
dat_cDC1 <- extract_dat(cDC1, '3_cDC1', OA, psa)

genes <- Reduce(intersect, list(rownames(dat_Mo), rownames(dat_Mf), rownames(dat_cDC2), rownames(dat_cDC1)))
dat <- cbind(dat_Mo[genes,], dat_Mf[genes,], dat_cDC2[genes,], dat_cDC1[genes,])
head(dat)

geneord <- read.delim(file="geneord.txt", sep="\t", stringsAsFactors=F, header=F)[,]
dat_Now <- dat[intersect(geneord, rownames(dat_Now)), ]
idx <- colnames(dat_Now)
dat_Now <- t(apply(dat_Now, 1, function(x) scale(x, scale = T)))
colnames(dat_Now) <- idx

col_ann <- data.frame(colnames(dat_Now))
col_ann$Pop <- foreach(i = colnames(dat_Now), .combine = c) %do% strsplit(i, split = '\\|')[[1]][3]
col_ann$Grp <- foreach(i = colnames(dat_Now), .combine = c) %do% strsplit(i, split = '\\|')[[1]][1]
col_ann$Grp[col_ann$Grp != '0_OA'] <- '1_PsA'
col_ann$Sec <- paste0(col_ann$Pop, '|', col_ann$Grp)
toSplit <- table(col_ann$Sec)

ann_top = HeatmapAnnotation(Pop = col_ann$Pop, Grp = col_ann$Grp, col= list(Pop = c("0_Mo" = "darkorange", "1_Mf" = "purple", "2_cDC2" = "deepskyblue4", "3_cDC1" = "darkgreen"), Grp = c("0_OA" = "#00A158", "1_PsA" = "#BE202D")))
selected4 <- c("GBP4", "STAT1", "CXCL10", "CXCL9", "SOCS1", "ISG15", "PSME2", "PSMB10", "CCR2", "IL4I1", "GSDMD", "JAK2", "TNFSF13B", "EIF3I", "EIF4EBP1",  "IL21R", "MMP12", "IFTIM3", "IL18BP", "TAP1", "SALAMF7", "CD48", "IDO1", "JAK3", "TNFSF10", "IFITM1", "SDHAF3", "SCL25A11", "GOLOGA5", "SLAMF7", "RAB2B", "TIFAB", "NLRC5", "APOL3", "OPTN", "FOLR2", "ITGB5", "MERTK", "STAB1", "TREM2", "BNIP3")
selected4 <- intersect(selected4, rownames(dat_Now)) 

ann_right = rowAnnotation(foo = anno_mark(at = match(selected4, rownames(dat_Now)), labels = selected4, labels_gp = gpar(fontsize=30)))

breakpoints <- c(-4, 0, 4)  # Replace with your actual scale limits
colors <- c("blue", "white", "red")  # Replace with your desired colors

pdf(file=paste0(outdir, "/", "fig3_all_nonames_all_lmt.pdf"), width=20, height=15, paper='special', useDingbats = FALSE )
Heatmap(dat_Now, 
        col = color_fun,
        top_annotation = ann_top,
        right_annotation = ann_right,
        cluster_columns = T, 
        cluster_column_slices = F, 
		show_row_names=F,
        cluster_rows = F,
        clustering_distance_rows = 'spearman', 
        column_split = rep(names(toSplit), as.vector(toSplit)),
        row_names_gp = gpar(fontsize = 6), 
        column_names_gp = gpar(fontsize = 10),
		column_dend_height = unit(0,"mm"))
dev.off()

