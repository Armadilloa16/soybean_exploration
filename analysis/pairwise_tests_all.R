library(parallel)

df_pairs = data.frame()
ind = 0
pb = txtProgressBar(min = 0, max = (length(gene_nms)^2 - length(gene_nms))/2, initial = 0, style = 3)
int_p_val = function(g1, g2){
  f = as.formula(paste('Yield ~', g1, '*', g2))
  m = lm(f, data = pav_table)
  int_nm = paste0(g1, ':', g2)
  if (int_nm %in% row.names(summary(m)$coefficients)) {
    return(summary(m)$coefficients[int_nm, "Pr(>|t|)"])
  } else {
    return(NA)
  }
}
for (i in 1:(length(gene_nms)-1)) {
  gene_nms_tmp = gene_nms[(i+1):length(gene_nms)]
  df_pairs = rbind(
    df_pairs,
    data.frame(gene.1 = factor(gene_nms[i], levels = gene_nms),
               gene.2 = factor(gene_nms_tmp, levels = gene_nms),
               p.val = unlist(mclapply(
                 gene_nms_tmp,
                 function(x){return(int_p_val(gene_nms[i], x))},
                 mc.cores = 7
               ))
    )
  )
  ind = ind + length(gene_nms_tmp)
  setTxtProgressBar(pb, ind)
}

# For purposes of file size optimisation
df_pairs$gene.1 = as.numeric(df_pairs$gene.1)
df_pairs$gene.2 = as.numeric(df_pairs$gene.2)

write.csv(df_pairs, file = file.path('output', 'pairwise_tests_all.csv'), row.names = FALSE)
