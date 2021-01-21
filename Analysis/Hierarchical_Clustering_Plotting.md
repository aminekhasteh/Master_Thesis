Hierarchical Clustering of the PRSs matrices
================

# Dendograms

We take a look at the Dendogram of the phenotypes of the residuals of
the PRSs.

``` r
knitr::opts_chunk$set(echo = TRUE, warning = FALSE)
op <- par(mfrow=c(1,2), mar=c(3,1,1,1))
for(prs_matrix in names(results1)){
                clusters <- results1[[prs_matrix]]$clust.var
                dend <- as.dendrogram(clusters)
                # order it the closest we can to the order of the observations:
                dend <- rotate(dend, 1:length(clusters$labels))
                # Color the branches based on the clusters:
                dend <- color_branches(dend, k = 3) #, groupLabels=iris_species)
                # We hang the dendrogram a bit:
                dend <- hang.dendrogram(dend, hang_height = 0.1)
                # reduce the size of the labels:
                dend <- set(dend, "labels_cex", 0.5)
                # And plot:
                par(mar = c(3, 3, 3, 7))
                if(prs_matrix == 'PCAresults_p_val_1.rds'){
                                jpeg(paste0("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Plots/Dendograms/Dendogram_Pheno_",
                                            gsub("PCAresults_","",gsub(".rds","",prs_matrix)),".jpg"),
                                     width = 1200, height = 1500)
                                plot(dend, cex.main=1.25,
                                     main =paste0('Dendogram of PRSs at phenotypes level',"\n",
                                                  dim(clusters[[1]])[1],' individuals',"\n",
                                                  'at P-value threshold of ',
                                                  unlist(strsplit(prs_matrix, "[_.]"))[4]),
                                     horiz =  TRUE,
                                     nodePar = list(cex = .007)
                                )
                                dev.off()
                                # circlize_dendrogram(dend)
                                
                } else {
                                jpeg(paste0("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Plots/Dendograms/Dendogram_Pheno_",
                                            gsub("PCAresults_","",gsub(".rds","",prs_matrix)),".jpg"),
                                     width = 1200, height = 1500)
                                plot(dend, cex.main=1.25,
                                     main =paste0('Dendogram of PRSs at phenotypes level',"\n",
                                                  dim(clusters[[1]])[1],' individuals',"\n",
                                                  'at P-value threshold of ',
                                                  paste0(unlist(strsplit(prs_matrix, "[_.]"))[4],'.',
                                                         unlist(strsplit(prs_matrix, "[_.]"))[5])),
                                     horiz =  TRUE,
                                     nodePar = list(cex = .007)
                                )
                                dev.off()
                                # circlize_dendrogram(dend)
                }
}
par(op)
```

We take a look at the Dendogram of the phenotypes of the residuals of
PRSs at an individual level.

``` r
knitr::opts_chunk$set(echo = TRUE, warning = FALSE)
op <- par(mfrow=c(1,2), mar=c(3,1,1,1))
for(prs_matrix in names(results1)){
                clusters <- results1[[prs_matrix]]$clust.ind
                dend <- as.dendrogram(clusters)
                # order it the closest we can to the order of the observations:
                dend <- rotate(dend, 1:length(clusters$labels))
                # Color the branches based on the clusters:
                dend <- color_branches(dend, k = 3) #, groupLabels=iris_species)
                # We hang the dendrogram a bit:
                dend <- hang.dendrogram(dend, hang_height = 0.1)
                # reduce the size of the labels:
                dend <- set(dend, "labels_cex", 0.5)
                # And plot:
                par(mar = c(3, 3, 3, 7))
                if(prs_matrix == 'PCAresults_p_val_1.rds'){
                                jpeg(paste0("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Plots/Dendograms/Dendogram_Ind_",
                                            gsub("PCAresults_","",gsub(".rds","",prs_matrix)),".jpg"),
                                     width = 1200, height = 1500)
                                plot(dend, cex.main=1.25,
                                     main =paste0('Dendogram of PRSs at individuals level',"\n",
                                                  dim(clusters[[1]])[1],' individuals',"\n",
                                                  'at P-value threshold of ',
                                                  unlist(strsplit(prs_matrix, "[_.]"))[4]),
                                     horiz =  TRUE,
                                     nodePar = list(cex = .007)
                                )
                                dev.off()
                                # circlize_dendrogram(dend)
                                
                } else {
                                jpeg(paste0("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Plots/Dendograms/Dendogram_Ind_",
                                            gsub("PCAresults_","",gsub(".rds","",prs_matrix)),".jpg"),
                                     width = 1200, height = 1500)
                                plot(dend, cex.main=1.25,
                                     main =paste0('Dendogram of PRSs at individuals level',"\n",
                                                  dim(clusters[[1]])[1],' individuals',"\n",
                                                  'at P-value threshold of ',
                                                  paste0(unlist(strsplit(prs_matrix, "[_.]"))[4],'.',
                                                         unlist(strsplit(prs_matrix, "[_.]"))[5])),
                                     horiz =  TRUE,
                                     nodePar = list(cex = .007)
                                )
                                dev.off()
                                # circlize_dendrogram(dend)
                }
}
par(op)
```
