library(tidyverse)
library(RColorBrewer)
library(ggrepel)

#' convert long tibble to pheatmap data input
#'
#' @param tbl.long tibble in long format
#' @param rowVar column name used as row for pheatmap matrix
#' @param colVar column name used as column for pheatmap matrix
#' @param valueVar column name used as the values for pheatmap matrix
#' @param colAnnVars column name(s) used for additonal column annotation
#' @param rowAnnVars column name(s) used for additonal row annotation
#' @return a list with 3 components: mat, rowAnn, colAnn
#' @export
tbl2hmap <- function(tbl.long, rowVar, colVar, valueVar,
                   colAnnVars = NULL, rowAnnVars = NULL) {
  Mat0 <- dplyr::select(tbl.long, one_of(c(rowVar, colVar, valueVar)))%>%
    tidyr::spread_(key = colVar, value=valueVar)
  Mat <- select(Mat0, -one_of(rowVar)) %>% as.matrix()
  rownames(Mat) <- unlist(select(Mat0,one_of(rowVar)))
  
  colAnn <- rowAnn <- NULL
  
  if(!is.null(colAnnVars)) {
    colAnn0 <- dplyr::select(tbl.long,one_of( c(colVar,colAnnVars))) %>%
      unique()
    colAnn <- select(colAnn0, - one_of(colVar)) %>% as.data.frame()
    rownames(colAnn) <- unlist(select( colAnn0,one_of(colVar)))
  }
  
  if(!is.null(rowAnnVars)) {
    rowAnn0 <- dplyr::select(tbl.long,one_of(c( rowVar,rowAnnVars))) %>%
      unique()
    rowAnn <- select(rowAnn0, -one_of( rowVar)) %>% as.data.frame()
    rownames(rowAnn) <- unlist(select(rowAnn0,one_of(rowVar)))
  }
  
  list(mat = Mat, rowAnn =rowAnn, colAnn = colAnn)
}

save_pheatmap_pdf <- function(x, filename, width, height) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# Sets the minimum (0), the maximum (15), and the increasing steps (+1) for the color scale
# Note: if some of your genes are outside of this range, they will appear white on the heatmap
breaksList <- seq(0, 1, by = 0.001)


# -s 20k -l 60k -p 95
HG002_all <- read.table('~/Downloads/Pangenomics/HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.dist.tsv', sep = '\t', header = T)

# -s 100k -l 300k -p 98
HG002_all <- read.table('~/Downloads/Pangenomics/HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.dist.tsv', sep = '\t', header = T)


# Filtering
if (FALSE) {
  HG002_all_filtered <- HG002_all %>%
    filter(!(
      `group.a` %in% c('Peregrine_HiFi_20kb.alt.utgs', 'Peregrine_HiFi_25kb.alt.utgs') | 
        `group.b` %in% c('Peregrine_HiFi_20kb.alt.utgs', 'Peregrine_HiFi_25kb.alt.utgs')
    ))
  imageD12 <- 'HG002.CMS.NoUtgs.D12.png'
  imageD23 <- 'HG002.CMS.NoUtgs.D23.png'
  imageD13 <- 'HG002.CMS.NoUtgs.D13.png'
  
  HG002_all_filtered <- HG002_all %>%
    filter(!(
      `group.a` %in% c('Peregrine_HiFi_20kb.alt.utgs', 'Peregrine_HiFi_25kb.alt.utgs', 'FALCON_Unzip.alt', 'Peregrine_HiFi_20kb.alt', 'Peregrine_HiFi_25kb.alt', 'HiCanu.alt') | 
        `group.b` %in% c('Peregrine_HiFi_20kb.alt.utgs', 'Peregrine_HiFi_25kb.alt.utgs', 'FALCON_Unzip.alt', 'Peregrine_HiFi_20kb.alt', 'Peregrine_HiFi_25kb.alt', 'HiCanu.alt')
    ))
  imageD12 <- 'HG002.CMS.NoUtgs.NoAlt.D12.png'
  imageD23 <- 'HG002.CMS.NoUtgs.NoAlt.D23.png'
  imageD13 <- 'HG002.CMS.NoUtgs.NoAlt.D13.png'
} else {
  HG002_all_filtered <- HG002_all
  imageD12 <- 'HG002.CMS.All.D12.png'
  imageD23 <- 'HG002.CMS.All.D23.png'
  imageD13 <- 'HG002.CMS.All.D13.png'
}

HG002_all_wide <- tbl2hmap(HG002_all_filtered, 'group.a', 'group.b', 'jaccard')

meta_df <- read.table("~/git/HG002_assemblies_assessment/data/HGRC_bakeoff_HG002_assemblies_metadata.tsv", sep = '\t', header = T) %>%
  select(AbbreviatedName, InstitutionOrCons, TopLevel, Haplotype, ContigAlgorithm, Type)
meta_df$AbbreviatedName <- gsub(' ', '_', meta_df$AbbreviatedName)

pheatmap::pheatmap(
  #cor(as.matrix(D), method = "spearman"),
  HG002_all_wide$mat,
  fontsize = 18, fontsize_row = 18, #height = 20,
  #breaks = breaksList,
  cluster_cols = T, clustering_distance_rows = "euclidean",
  cluster_rows = T, clustering_distance_cols = "euclidean",
  #border=TRUE, border_color = "black", gaps_col = 1:6,
  color=colorRampPalette(c("#114477", "#FFFFFF"))(1000),
  #color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
  
  annotation = meta_df %>%
    dplyr::select(
      AbbreviatedName,
      InstitutionOrCons,
      TopLevel,
      Haplotype,
      ContigAlgorithm
      ) %>%
    column_to_rownames(var = "AbbreviatedName"),
  #annotation_colors=list(
  #  Ploidy = c(Hyperdiploid = '#F8766D', `Near-to-diploid`='#00BFC4'),
  #  `Cell line` = c(`17` = '#D89000', `30` = '#DA62DB', `32` = '#55CC7D')
)

# HG002_all45.jaccard.heatmap.11x11.2.pdf
# HG002_all45.s20k.l60k.p95.jaccard.heatmap.clustering.1900x1700.png
# HG002_all45.s20k.l60k.p95.jaccard.heatmap.no_clustering.1900x1700.png
# HG002_all45.s20k.l60k.p95.jaccard.heatmap.NoUtgs.NoAlt.clustering.1900x1700.png
# HG002_all45.s20k.l60k.p95.jaccard.heatmap.NoUtgs.NoAlt.no_clustering.1900x1700.png

# HG002_all45.s100k.l300k.p98.jaccard.heatmap.clustering.1900x1700.png
# HG002_all45.s100k.l300k.p98.jaccard.heatmap.no_clustering.1900x1700.png
# HG002_all45.s100k.l300k.p98.jaccard.heatmap.NoUtgs.NoAlt.clustering.1900x1700.png
# HG002_all45.s100k.l300k.p98.jaccard.heatmap.NoUtgs.NoAlt.no_clustering.1900x1700.png


# Classical (Metric) Multidimensional Scaling

HG002_all_wide_euclidean <- tbl2hmap(HG002_all_filtered, 'group.a', 'group.b', 'euclidean')

#https://stats.stackexchange.com/questions/87681/performing-pca-with-only-a-distance-matrix
#http://www.gastonsanchez.com/visually-enforced/how-to/2013/01/23/MDS-in-R/
#https://stats.stackexchange.com/questions/31908/what-is-percentage-of-variance-in-pca
fit <- cmdscale(as.dist(HG002_all_wide_euclidean$mat), eig=TRUE, k = 3)

#https://stats.stackexchange.com/questions/22019/how-to-calculate-the-r-squared-value-and-assess-the-model-fit-in-multidimensiona/22020
plot(cumsum(fit$eig) / sum(fit$eig), 
     type="h", lwd=5, las=1, 
     xlab="Number of dimensions", 
     ylab=expression(R^2))
plot(fit$eig, 
       type="h", lwd=5, las=1, 
       xlab="Number of dimensions", 
       ylab="Eigenvalues")

fit_df <- as.data.frame(fit$points) %>%
  rownames_to_column(var = "AbbreviatedName")

  
fit_and_meta_df <- merge(fit_df, meta_df)


ggplot(data = fit_and_meta_df, aes(x = V1, y = V2, label = AbbreviatedName, shape = TopLevel, color = Haplotype)) +
  geom_point(size=2) + 
  geom_text_repel(
    size=3,
    show.legend  = FALSE, # to hide the `a` from the legend
    max.overlaps=Inf
    ) +
  xlab(paste0("Dimension 1 (", trunc(fit$eig[1]/sum(fit$eig) * 100 * 10^2)/10^2, "%)"))+
  ylab(paste0("Dimension 2 (", trunc(fit$eig[2]/sum(fit$eig) * 100 * 10^2)/10^2, "%)"))
ggsave(imageD12, width = 21, height = 14,  units = "cm", dpi = 300,  bg = "transparent")

ggplot(data = fit_and_meta_df, aes(x = V2, y = V3, label = AbbreviatedName, shape = TopLevel, color = Haplotype)) +
  geom_point(size=2) + 
  geom_text_repel(
    size=3,
    show.legend  = FALSE, # to hide the `a` from the legend
    max.overlaps=Inf
  ) +
  xlab(paste0("Dimension 2 (", trunc(fit$eig[2]/sum(fit$eig) * 100 * 10^2)/10^2, "%)"))+
  ylab(paste0("Dimension 3 (", trunc(fit$eig[3]/sum(fit$eig) * 100 * 10^2)/10^2, "%)"))
ggsave(imageD23, width = 21, height = 14,  units = "cm", dpi = 300,  bg = "transparent")

ggplot(data = fit_and_meta_df, aes(x = V1, y = V3, label = AbbreviatedName, shape = TopLevel, color = Haplotype)) +
  geom_point(size=2) + 
  geom_text_repel(
    size=3,
    show.legend  = FALSE, # to hide the `a` from the legend
    max.overlaps=Inf
  ) +
  xlab(paste0("Dimension 1 (", trunc(fit$eig[1]/sum(fit$eig) * 100 * 10^2)/10^2, "%)"))+
  ylab(paste0("Dimension 3 (", trunc(fit$eig[3]/sum(fit$eig) * 100 * 10^2)/10^2, "%)"))
ggsave(imageD13, width = 21, height = 14,  units = "cm", dpi = 300,  bg = "transparent")




# All chromosomes
# Classical (Metric) Multidimensional Scaling
meta_df <- read.table("~/git/HG002_assemblies_assessment/data/HGRC_bakeoff_HG002_assemblies_metadata.tsv", sep = '\t', header = T) %>%
  select(AbbreviatedName, InstitutionOrCons, TopLevel, Haplotype, ContigAlgorithm, Type)
meta_df$AbbreviatedName <- gsub(' ', '_', meta_df$AbbreviatedName)

base_dir='~/Downloads/Pangenomics/HG002_bakeoff'
for (N in c(seq(1, 22), 'XY', '1to22', 'All')){
  if (N == 'All') {
    path_chrN = file.path(base_dir, paste0('HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.chr', N, '.dist.tsv'))
  }else{
    path_chrN = file.path(base_dir, paste0('HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.Y.x100.chr', N, '.dist.tsv'))
  }

  print(path_chrN)
  
  HG002_all <- read.table(path_chrN, sep = '\t', header = T)
  
  if (TRUE) {
    if (FALSE) {
      HG002_all_filtered <- HG002_all %>%
        filter(!(
          `group.a` %in% c('Peregrine_HiFi_20kb.alt.utgs', 'Peregrine_HiFi_25kb.alt.utgs') | 
            `group.b` %in% c('Peregrine_HiFi_20kb.alt.utgs', 'Peregrine_HiFi_25kb.alt.utgs')
        ))
      imageD12 = file.path(base_dir, paste0('HG002.LongAlignments.NoUtgs.chr', N, '.CMS.D1vsD2.pdf'))
      imageD23 = file.path(base_dir, paste0('HG002.LongAlignments.NoUtgs.chr', N, '.CMS.D12sD3.pdf'))
      imageD13 = file.path(base_dir, paste0('HG002.LongAlignments.NoUtgs.chr', N, '.CMS.D1vsD3.pdf'))
      path_heatmap = file.path(base_dir, paste0('HG002.LongAlignments.NoUtgs.chr', N, '.JaccardHeatmap.Clustering.pdf'))
      title = paste0('Long alignments - No Utgs - chr', N)
    } else {
      HG002_all_filtered <- HG002_all %>%
        filter(!(
          `group.a` %in% c('Peregrine_HiFi_20kb.alt.utgs', 'Peregrine_HiFi_25kb.alt.utgs', 'FALCON_Unzip.alt', 'Peregrine_HiFi_20kb.alt', 'Peregrine_HiFi_25kb.alt', 'HiCanu.alt') | 
            `group.b` %in% c('Peregrine_HiFi_20kb.alt.utgs', 'Peregrine_HiFi_25kb.alt.utgs', 'FALCON_Unzip.alt', 'Peregrine_HiFi_20kb.alt', 'Peregrine_HiFi_25kb.alt', 'HiCanu.alt')
        ))
      imageD12 = file.path(base_dir, paste0('HG002.LongAlignments.NoUtgs.NoAlt.chr', N, '.CMS.D1vsD2.pdf'))
      imageD23 = file.path(base_dir, paste0('HG002.LongAlignments.NoUtgs.NoAlt.chr', N, '.CMS.D12sD3.pdf'))
      imageD13 = file.path(base_dir, paste0('HG002.LongAlignments.NoUtgs.NoAlt.chr', N, '.CMS.D1vsD3.pdf'))
      path_heatmap = file.path(base_dir, paste0('HG002.LongAlignments.NoUtgs.NoAlt.chr', N, '.JaccardHeatmap.Clustering.pdf'))
      title = paste0('Long alignments - No Utgs and No alt - chr', N)
    }
  } else {
    HG002_all_filtered <- HG002_all
    imageD12 = file.path(base_dir, paste0('HG002.LongAlignments.AllData.chr', N, '.CMS.D1vsD2.pdf'))
    imageD23 = file.path(base_dir, paste0('HG002.LongAlignments.AllData.chr', N, '.CMS.D12sD3.pdf'))
    imageD13 = file.path(base_dir, paste0('HG002.LongAlignments.AllData.chr', N, '.CMS.D1vsD3.pdf'))
    path_heatmap = file.path(base_dir, paste0('HG002.LongAlignments.AllData.chr', N, '.JaccardHeatmap.Clustering.pdf'))
    title = paste0('Long alignments - All data - chr', N)
  }
  
  # Jaccard heatmap
  HG002_all_wide <- tbl2hmap(HG002_all_filtered, 'group.a', 'group.b', 'jaccard')


  ph <- pheatmap::pheatmap(
    HG002_all_wide$mat,
    fontsize = 18, fontsize_row = 18, #height = 20,
    #breaks = breaksList,
    cluster_cols = T, clustering_distance_rows = "euclidean",
    cluster_rows = T, clustering_distance_cols = "euclidean",
    
    #border=TRUE, border_color = "black", gaps_col = 1:6,
    
    color=colorRampPalette(c("#114477", "#FFFFFF"))(1000),
    #color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
    
    main= title,
    
    annotation = meta_df %>%
      dplyr::select(
        AbbreviatedName,
        InstitutionOrCons,
        TopLevel,
        Haplotype,
        ContigAlgorithm
      ) %>%
      column_to_rownames(var = "AbbreviatedName"),
    #annotation_colors=list(
    #  Ploidy = c(Hyperdiploid = '#F8766D', `Near-to-diploid`='#00BFC4'),
    #  `Cell line` = c(`17` = '#D89000', `30` = '#DA62DB', `32` = '#55CC7D')
  )
  save_pheatmap_pdf(ph, path_heatmap, 20, 18)  


}


# Classical (Metric) Multidimensional Scaling
HG002_all_wide_euclidean <- tbl2hmap(HG002_all_filtered, 'group.a', 'group.b', 'euclidean')

fit <- cmdscale(as.dist(HG002_all_wide_euclidean$mat), eig=TRUE, k = 3)

fit_df <- as.data.frame(fit$points) %>%
  rownames_to_column(var = "AbbreviatedName")
fit_and_meta_df <- merge(fit_df, meta_df)

# TODO: TO ADD WITH AND WITHOUT LABELS
ggplot(data = fit_and_meta_df, aes(x = V1, y = V2, label = AbbreviatedName, shape = TopLevel, color = Haplotype)) +
  geom_point(size=2) + 
  geom_text_repel(
    size=3,
    show.legend  = FALSE, # to hide the `a` from the legend
    max.overlaps=Inf
  ) +
  xlab(paste0("Dimension 1 (", trunc(fit$eig[1]/sum(fit$eig) * 100 * 10^2)/10^2, "%)"))+
  ylab(paste0("Dimension 2 (", trunc(fit$eig[2]/sum(fit$eig) * 100 * 10^2)/10^2, "%)")) +
  ggtitle(title) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(imageD12, width = 21, height = 14,  units = "cm", dpi = 300,  bg = "transparent")

ggplot(data = fit_and_meta_df, aes(x = V2, y = V3, label = AbbreviatedName, shape = TopLevel, color = Haplotype)) +
  geom_point(size=2) + 
  geom_text_repel(
    size=3,
    show.legend  = FALSE, # to hide the `a` from the legend
    max.overlaps=Inf
  ) +
  xlab(paste0("Dimension 2 (", trunc(fit$eig[2]/sum(fit$eig) * 100 * 10^2)/10^2, "%)"))+
  ylab(paste0("Dimension 3 (", trunc(fit$eig[3]/sum(fit$eig) * 100 * 10^2)/10^2, "%)")) +
  ggtitle(title) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(imageD23, width = 21, height = 14,  units = "cm", dpi = 300,  bg = "transparent")

ggplot(data = fit_and_meta_df, aes(x = V1, y = V3, label = AbbreviatedName, shape = TopLevel, color = Haplotype)) +
  geom_point(size=2) + 
  geom_text_repel(
    size=3,
    show.legend  = FALSE, # to hide the `a` from the legend
    max.overlaps=Inf
  ) +
  xlab(paste0("Dimension 1 (", trunc(fit$eig[1]/sum(fit$eig) * 100 * 10^2)/10^2, "%)"))+
  ylab(paste0("Dimension 3 (", trunc(fit$eig[3]/sum(fit$eig) * 100 * 10^2)/10^2, "%)")) +
  ggtitle(title) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(imageD13, width = 21, height = 14,  units = "cm", dpi = 300,  bg = "transparent")

library(ggfortify)
autoplot(fit, label = TRUE, label.size = 5)

x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
     main="Metric MDS", type="n") 
text(x, y, labels = row.names(fit$points), cex=.7) 



# NOT SURE THE MATH IS PRESERVED
library(factoextra) # install.packages("factoextra")
res.pca <- prcomp(HG002_all_wide_euclidean$mat, scale = F)
#fviz_eig(res.pca)
fviz_pca_ind(res.pca,
             #col.ind = "cos2", # Color by the quality of representation
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


HG002_Erik <- read.table('~/Downloads/Pangenomics/all.l50k.k256.l256.odgi.asm.dist.tsv', sep = '\t', header = T)
HG002_Erik_wide_euclidean <- tbl2hmap(HG002_Erik, 'group.a', 'group.b', 'euclidean')

HG002_Erik_wide_euclidean_fit <- cmdscale(as.dist(HG002_Erik_wide_euclidean$mat), eig=TRUE, k = 6)

library(ggfortify)
autoplot(HG002_Erik_wide_euclidean_fit, label = TRUE, label.size = 5)

library(factoextra) # install.packages("factoextra")
res.pca <- prcomp(HG002_Erik_wide_euclidean$mat, scale = T)
#fviz_eig(res.pca)
fviz_pca_ind(res.pca,
             #col.ind = "cos2", # Color by the quality of representation
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

