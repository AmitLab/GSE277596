

library(metacell)

cl_compute_e_gc= function(cl_mc, us, norm_by_mc_meansize = T){
 f_g_cov = rowSums(us) > 10
 
 e_gc = t(tgs_matrix_tapply(us[f_g_cov,], cl_mc, function(y) {exp(mean(log(1+y)))-1}))
 rownames(e_gc) = rownames(us)[f_g_cov]
 
 if (norm_by_mc_meansize) {
 mc_meansize = tapply(colSums(us), cl_mc, mean)
 e_gc = t(t(e_gc)/as.vector(mc_meansize))
 }
 return(e_gc)
}

read_large_umis_n <- function(mat_id, bs = 1e4, cells = NULL, norm_to = 2000, do_log = F, genes = NULL) {
	mat = scdb_mat(mat_id)
	if (is.null(cells)) { cells = mat@cells}
	ncells = length(cells)
	umis = NULL
	if (is.null(genes)) { genes = rownames(mat@mat)}
	for (i in seq_len(ncells %/% bs + 1)) {
		from = (i - 1) * bs + 1
		to = min(i * bs, ncells)
		if( to > from ){ 
			print(from)
			print(to)
			mt <- as.matrix(mat@mat[genes, cells[from:to]])
			mtn <- sweep(mt, 2, colSums(mt), "/") * norm_to
			if (do_log == T){
				umis = cbind(umis, log(mtn+1, 2))
			} else {
				umis = cbind(umis, mtn)
			}
			mt <- NULL
			mtn <- NULL
		}
	}
	umis
}
min_cell_count <- 30

sc_mat <- scdb_mat(mat_id)
sc_cl<-scdb_mc("test_mc")
cell_stats <- sc_mat@cell_metadata
cell_stats <- merge(cell_stats, sc_cl@mc, by=0)
rownames(cell_stats) <- cell_stats[, "Row.names"]
cell_stats <- cell_stats[, -1]
cell_stats[, "group"] <- ann[cell_stats[, "y"], "Populations"]

for (g in pop){
	good_PID <- names(which(table(cell_stats[cell_stats[, "group"] %in% g, "PID"]) > min_cell_count))
	cells <- rownames(cell_stats[cell_stats[, "group"] %in% g & cell_stats[, "PID"] %in% good_PID, ])
	umis <- read_large_umis_n(mat_id, cells=cells, genes=rownames(sc_cl@mc_fp)) 
	cl_pid <- cell_stats[cells, "PID"]
	names(cl_pid) <- cells
	g_mean <- cl_compute_e_gc(cl_pid, umis, norm_by_mc_meansize=F)
	write.table(file=gsub(pattern= " ", replacement= "_", x= paste0(g, "_mean_n.txt")), g_mean, quote=F)
}
