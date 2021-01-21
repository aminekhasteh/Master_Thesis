Exploring with PRSs and ROS/MAP
================

## Packages and importing datasets:

# Hierarchical Clustering

We take a look at the Dendogram of the phenotypes of the residuals of
the PRSs.

``` r
knitr::opts_chunk$set(echo = TRUE, warning = FALSE)

for(prs_matrix in names(results1)){
  matrix <- as.data.frame(results1[[prs_matrix]]$residuals)
  matrix <- matrix[,-dim(matrix)[2]]
  names(matrix) <- c(paste0('V',c(1:(dim(matrix)[2])))) # Need to come up with a better way of labeling the phenotypes
  clusters <- hclust(dist(t(matrix)))
  if(prs_matrix == 'PCAresults_p_val_1.rds'){
    plot(clusters, col='lightblue', xlab='Phenotypes', 
         main =paste0('Dendogram of PRSs of ',dim(matrix)[2],' phenotypes, at P-value threshold of ',
                                                        paste0(unlist(strsplit(prs_matrix, "[_.]"))[3],' ',
                                                               unlist(strsplit(prs_matrix, "[_.]"))[4])))
    
  } else {
    plot(clusters, col='lightblue', xlab = 'Phenotypes',
         main =paste0('Dendogram of PRSs of ',dim(matrix)[2],' phenotypes, at P-value threshold of ',
                                                        paste0(unlist(strsplit(prs_matrix, "[_.]"))[3],' ',
                                                               unlist(strsplit(prs_matrix, "[_.]"))[4],'.',
                                                               unlist(strsplit(prs_matrix, "[_.]"))[5])))
  }
}
```

<img src="EDA_files/figure-gfm/dd-1.png" style="display: block; margin: auto;" /><img src="EDA_files/figure-gfm/dd-2.png" style="display: block; margin: auto;" /><img src="EDA_files/figure-gfm/dd-3.png" style="display: block; margin: auto;" /><img src="EDA_files/figure-gfm/dd-4.png" style="display: block; margin: auto;" /><img src="EDA_files/figure-gfm/dd-5.png" style="display: block; margin: auto;" /><img src="EDA_files/figure-gfm/dd-6.png" style="display: block; margin: auto;" /><img src="EDA_files/figure-gfm/dd-7.png" style="display: block; margin: auto;" /><img src="EDA_files/figure-gfm/dd-8.png" style="display: block; margin: auto;" /><img src="EDA_files/figure-gfm/dd-9.png" style="display: block; margin: auto;" /><img src="EDA_files/figure-gfm/dd-10.png" style="display: block; margin: auto;" />

We take a look at the Dendogram of the phenotypes of the residuals of
PRSs at an individual level.

``` r
knitr::opts_chunk$set(echo = TRUE, warning = FALSE)

for(prs_matrix in names(results1)){
  matrix <- as.data.frame(results1[[prs_matrix]]$residuals)
  matrix <- matrix[,-dim(matrix)[2]]
  matrix <- t(matrix)
  names(matrix) <- c(paste0('V',c(1:(dim(matrix)[1])))) # Need to come up with a way of categorizing the particpants
  clusters <- hclust(dist(t(matrix)))
  if(prs_matrix == 'PCAresults_p_val_1.rds'){
    plot(clusters, col='lightblue', xlab = 'Individuals',
         main =paste0('Dendogram of PRSs of individuals of ',dim(matrix)[1],' individuals, at P-value threshold of ',
                                                        paste0(unlist(strsplit(prs_matrix, "[_.]"))[3],'.',
                                                               unlist(strsplit(prs_matrix, "[_.]"))[4])))
    
  } else {
    plot(clusters, col='lightblue', xlab = 'Individuals',
         main =paste0('Dendogram of PRSs of individuals of ',dim(matrix)[1],' individuals, at P-value threshold of ',
                                                        paste0(unlist(strsplit(prs_matrix, "[_.]"))[3],' ',
                                                               unlist(strsplit(prs_matrix, "[_.]"))[4],'.',
                                                               unlist(strsplit(prs_matrix, "[_.]"))[5])))
  }
}
```

<img src="EDA_files/figure-gfm/dd1-1.png" style="display: block; margin: auto;" /><img src="EDA_files/figure-gfm/dd1-2.png" style="display: block; margin: auto;" /><img src="EDA_files/figure-gfm/dd1-3.png" style="display: block; margin: auto;" /><img src="EDA_files/figure-gfm/dd1-4.png" style="display: block; margin: auto;" /><img src="EDA_files/figure-gfm/dd1-5.png" style="display: block; margin: auto;" /><img src="EDA_files/figure-gfm/dd1-6.png" style="display: block; margin: auto;" /><img src="EDA_files/figure-gfm/dd1-7.png" style="display: block; margin: auto;" /><img src="EDA_files/figure-gfm/dd1-8.png" style="display: block; margin: auto;" /><img src="EDA_files/figure-gfm/dd1-9.png" style="display: block; margin: auto;" /><img src="EDA_files/figure-gfm/dd1-10.png" style="display: block; margin: auto;" />
