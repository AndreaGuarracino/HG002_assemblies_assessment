library(tidyverse)
library(RColorBrewer)
library(ggrepel)

# convert long tibble to pheatmap data input
# @param tbl.long tibble in long format
# @param rowVar column name used as row for pheatmap matrix
# @param colVar column name used as column for pheatmap matrix
# @param valueVar column name used as the values for pheatmap matrix
# @param colAnnVars column name(s) used for additonal column annotation
# @param rowAnnVars column name(s) used for additonal row annotation
# @return a list with 3 components: mat, rowAnn, colAnn
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

meta_df <- read.table("~/git/HG002_assemblies_assessment/data/HGRC_bakeoff_HG002_assemblies_metadata.tsv", sep = '\t', header = T) %>%
  select(AbbreviatedName, InstitutionOrCons, TopLevel, Haplotype, ContigAlgorithm, Type, ContigAlgorithmBroader)
meta_df$AbbreviatedName <- gsub(' ', '_', meta_df$AbbreviatedName)
meta_df$Haplotype <- factor(
  meta_df$Haplotype,
  levels = c('hap1', 'hap2', 'maternal', 'paternal', 'primary', 'pseudo-hap ', 'alternate', 'unitigs')
)

base_dir <- '~/Downloads/Pangenomics/HG002_bakeoff'

################################################################################
# Jaccard and PCA analysis and plots
#====================================
#dir.create(file.path(base_dir, 'Heatmaps', 'NoUtgs', 'ByChromosome'), recursive = T)
dir.create(file.path(base_dir, 'Heatmaps', 'NoUtgs.NoAlt', 'ByChromosome'), recursive = T)
dir.create(file.path(base_dir, 'Heatmaps', 'NoUtgs.NoAltExceptPeregrine', 'ByChromosome'), recursive = T)
dir.create(file.path(base_dir, 'Heatmaps', 'AllAssemblies', 'ByChromosome'), recursive = T)
#dir.create(file.path(base_dir, 'ClassicalMultidimensionalScaling', 'NoUtgs', 'ByChromosome'), recursive = T)
dir.create(file.path(base_dir, 'ClassicalMultidimensionalScaling', 'NoUtgs.NoAlt', 'ByChromosome'), recursive = T)
dir.create(file.path(base_dir, 'ClassicalMultidimensionalScaling', 'NoUtgs.NoAltExceptPeregrine', 'ByChromosome'), recursive = T)
dir.create(file.path(base_dir, 'ClassicalMultidimensionalScaling', 'AllAssemblies', 'ByChromosome'), recursive = T)

#gg_color_hue <- function(n) {
#  hues = seq(15, 375, length = n + 1)
#  hcl(h = hues, l = 65, c = 100)[1:n]
#}
haplotype_colors <- c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3", "#C77CFF", "#44758B") #gg_color_hue(6) plus 2 other colors

for (N in c('All', 'XY', '1to22', seq(1, 22))){
  # Input matrix of distances
  if (N == 'All') {
    path_chrN <- file.path(base_dir, paste0('HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.chr', N, '.dist.tsv'))
  }else{
    path_chrN <- file.path(base_dir, paste0('HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.Y.x100.chr', N, '.dist.tsv'))
  }

  # Read the matrix
  HG002_all <- read.table(path_chrN, sep = '\t', header = T)
  
  # NOTE: the loop doesn't work for all datasets together. Set one dataset at a time
  ###for (dataset in c('AllAssemblies', 'NoUtgs.NoAlt', 'NoUtgs.NoAltExceptPeregrine')) {
  for (dataset in c('AllAssemblies')) {
  #for (dataset in c('NoUtgs.NoAlt')) {
  #for (dataset in c('NoUtgs.NoAltExceptPeregrine')) {
    if (dataset == 'AllAssemblies') {
      HG002_all_filtered <- HG002_all
      
      title <- paste0('AllAssemblies - chr', N)
      
      if (N %in%  c('All', 'XY', '1to22')){
        heatmaps_dir <- file.path(base_dir, 'Heatmaps', 'AllAssemblies')
      } else {
        heatmaps_dir <- file.path(base_dir, 'Heatmaps', 'AllAssemblies', 'ByChromosome')
      }
      path_heatmap <- file.path(heatmaps_dir, paste0('JaccardHeatmap.LongAlignments.AllAssemblies.chr', N, '.Clustering.pdf'))
      
      if (N %in%  c('All', 'XY', '1to22')){
        CMS_dir <- file.path(base_dir, 'ClassicalMultidimensionalScaling', 'AllAssemblies')
      } else {
        CMS_dir <- file.path(base_dir, 'ClassicalMultidimensionalScaling', 'AllAssemblies', 'ByChromosome')
      }
      prefix <- 'CMS.LongAlignments.AllAssemblies.chr'
    } else if (dataset == 'NoUtgs'){
      HG002_all_filtered <- HG002_all %>%
        filter(!(
          `group.a` %in% c('Peregrine_HiFi_20kb.alt.utgs', 'Peregrine_HiFi_25kb.alt.utgs') | 
            `group.b` %in% c('Peregrine_HiFi_20kb.alt.utgs', 'Peregrine_HiFi_25kb.alt.utgs')
        ))
      
      title <- paste0('No Utgs - chr', N)
      
      if (N %in%  c('All', 'XY', '1to22')){
        heatmaps_dir <- file.path(base_dir, 'Heatmaps', 'NoUtgs')
      } else {
        heatmaps_dir <- file.path(base_dir, 'Heatmaps', 'NoUtgs', 'ByChromosome')
      }
      path_heatmap <- file.path(heatmaps_dir, paste0('JaccardHeatmap.LongAlignments.NoUtgs.NoAlt.chr', N, '.Clustering.pdf'))
      
      
      if (N %in%  c('All', 'XY', '1to22')){
        CMS_dir <- file.path(base_dir, 'ClassicalMultidimensionalScaling', 'NoUtgs')
      } else {
        CMS_dir <- file.path(base_dir, 'ClassicalMultidimensionalScaling', 'NoUtgs', 'ByChromosome')
      }
      prefix <- 'CMS.LongAlignments.NoUtgs.chr'
    } else if (dataset == 'NoUtgs.NoAltExceptPeregrine') {
      HG002_all_filtered <- HG002_all %>%
        filter(!(
          `group.a` %in% c('Peregrine_HiFi_20kb.alt.utgs', 'Peregrine_HiFi_25kb.alt.utgs', 'FALCON_Unzip.alt', 'HiCanu.alt') | 
            `group.b` %in% c('Peregrine_HiFi_20kb.alt.utgs', 'Peregrine_HiFi_25kb.alt.utgs', 'FALCON_Unzip.alt', 'HiCanu.alt')
        ))
      
      title <- paste0('No Utgs and No alt except Peregrine - chr', N)
      
      if (N %in%  c('All', 'XY', '1to22')){
        heatmaps_dir <- file.path(base_dir, 'Heatmaps', 'NoUtgs.NoAltExceptPeregrine')
      } else {
        heatmaps_dir <- file.path(base_dir, 'Heatmaps', 'NoUtgs.NoAltExceptPeregrine', 'ByChromosome')
      }
      path_heatmap <- file.path(heatmaps_dir, paste0('JaccardHeatmap.LongAlignments.NoUtgs.NoAltExceptPeregrine.chr', N, '.Clustering.pdf'))
      
      if (N %in%  c('All', 'XY', '1to22')){
        CMS_dir <- file.path(base_dir, 'ClassicalMultidimensionalScaling', 'NoUtgs.NoAltExceptPeregrine')
      } else {
        CMS_dir <- file.path(base_dir, 'ClassicalMultidimensionalScaling', 'NoUtgs.NoAltExceptPeregrine', 'ByChromosome')
      }
      prefix <- 'CMS.LongAlignments.NoUtgs.NoAltExceptPeregrine.chr'
    } else if (dataset == 'NoUtgs.NoAlt'){
      HG002_all_filtered <- HG002_all %>%
        filter(!(
          `group.a` %in% c('Peregrine_HiFi_20kb.alt.utgs', 'Peregrine_HiFi_25kb.alt.utgs', 'FALCON_Unzip.alt', 'Peregrine_HiFi_20kb.alt', 'Peregrine_HiFi_25kb.alt', 'HiCanu.alt') | 
            `group.b` %in% c('Peregrine_HiFi_20kb.alt.utgs', 'Peregrine_HiFi_25kb.alt.utgs', 'FALCON_Unzip.alt', 'Peregrine_HiFi_20kb.alt', 'Peregrine_HiFi_25kb.alt', 'HiCanu.alt')
        ))
      
      title <- paste0('No Utgs and No alt - chr', N)
      
      if (N %in%  c('All', 'XY', '1to22')){
        heatmaps_dir <- file.path(base_dir, 'Heatmaps', 'NoUtgs.NoAlt')
      } else {
        heatmaps_dir <- file.path(base_dir, 'Heatmaps', 'NoUtgs.NoAlt', 'ByChromosome')
      }
      path_heatmap <- file.path(heatmaps_dir, paste0('JaccardHeatmap.LongAlignments.NoUtgs.NoAlt.chr', N, '.Clustering.pdf'))
      
      if (N %in%  c('All', 'XY', '1to22')){
        CMS_dir <- file.path(base_dir, 'ClassicalMultidimensionalScaling', 'NoUtgs.NoAlt')
      } else {
        CMS_dir <- file.path(base_dir, 'ClassicalMultidimensionalScaling', 'NoUtgs.NoAlt', 'ByChromosome')
      }
      prefix <- 'CMS.LongAlignments.NoUtgs.NoAlt.chr'
    } else {
      print('No dataset')
      next
    }
    
    imageD12 <- file.path(CMS_dir, paste0(prefix, N, '.D1vsD2.noLabels.pdf'))
    imageD23 <- file.path(CMS_dir, paste0(prefix, N, '.D2vsD3.noLabels.pdf'))
    imageD13 <- file.path(CMS_dir, paste0(prefix, N, '.D1vsD3.noLabels.pdf'))
    imageD14 <- file.path(CMS_dir, paste0(prefix, N, '.D1vsD4.noLabels.pdf'))
    imageD15 <- file.path(CMS_dir, paste0(prefix, N, '.D1vsD5.noLabels.pdf'))
    imageD12_withLabels <- file.path(CMS_dir, paste0(prefix, N, '.D1vsD2.withLabels.pdf'))
    imageD23_withLabels <- file.path(CMS_dir, paste0(prefix, N, '.D2vsD3.withLabels.pdf'))
    imageD13_withLabels <- file.path(CMS_dir, paste0(prefix, N, '.D1vsD3.withLabels.pdf'))
    imageD14_withLabels <- file.path(CMS_dir, paste0(prefix, N, '.D1vsD4.withLabels.pdf'))
    imageD15_withLabels <- file.path(CMS_dir, paste0(prefix, N, '.D1vsD5.withLabels.pdf'))

    print(path_chrN)
    print(title)
    
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
          ContigAlgorithmBroader
        ) %>%
        rename(`Institution or Consortium` = InstitutionOrCons) %>%
        rename(`Top Level` = TopLevel) %>%
        rename(`Contig Algorithm` = ContigAlgorithmBroader) %>%
        column_to_rownames(var = "AbbreviatedName"),
      annotation_colors=list(
        `Haplotype` = c(
          `hap1` = haplotype_colors[1],
          `hap2` = haplotype_colors[2],
          `maternal` = haplotype_colors[3],
          `paternal` = haplotype_colors[4],
          `primary` = haplotype_colors[5],
          `pseudo-hap ` = haplotype_colors[6],
          `alternate` = haplotype_colors[7],
          `unitigs` = haplotype_colors[8]
          ),
        `Contig Algorithm` = c(
          Diploid = brewer.pal(n = 8, name = "Dark2")[1],
          Haploid = brewer.pal(n = 8, name = "Dark2")[2],
          `Partial diploid` = brewer.pal(n = 8, name = "Dark2")[3],
          `Reference based` = brewer.pal(n = 8, name = "Dark2")[4],
          `Trio diploid` = brewer.pal(n = 8, name = "Dark2")[5]
        )
        #  `Cell line` = c(`17` = '#D89000', `30` = '#DA62DB', `32` = '#55CC7D')
      )
    )
    save_pheatmap_pdf(ph, path_heatmap, 21, 18)
  }
  
  
  # Classical (Metric) Multidimensional Scaling
  HG002_all_wide_euclidean <- tbl2hmap(HG002_all_filtered, 'group.a', 'group.b', 'euclidean')
  
  fit <- cmdscale(as.dist(HG002_all_wide_euclidean$mat), eig=TRUE, k = 5)
  
  fit_df <- as.data.frame(fit$points) %>%
    rownames_to_column(var = "AbbreviatedName")
  fit_and_meta_df <- merge(fit_df, meta_df) %>%
    rename(`Top Level` = TopLevel) %>%
    rename(`Contig Algorithm` = ContigAlgorithmBroader)

  plotD1D2 <- ggplot(data = fit_and_meta_df, aes(x = V1, y = V2, label = AbbreviatedName, shape = `Top Level`, color = Haplotype)) +
    scale_color_manual(values=haplotype_colors) +
    geom_point(size=2) + 
    xlab(paste0("Dimension 1 (", trunc(fit$eig[1]/sum(fit$eig) * 100 * 10^2)/10^2, "%)"))+
    ylab(paste0("Dimension 2 (", trunc(fit$eig[2]/sum(fit$eig) * 100 * 10^2)/10^2, "%)")) +
    ggtitle(title)
  plotD2D3 <- ggplot(data = fit_and_meta_df, aes(x = V2, y = V3, label = AbbreviatedName, shape = `Top Level`, color = Haplotype)) +
    scale_color_manual(values=haplotype_colors) +
    geom_point(size=2) + 
    xlab(paste0("Dimension 2 (", trunc(fit$eig[2]/sum(fit$eig) * 100 * 10^2)/10^2, "%)"))+
    ylab(paste0("Dimension 3 (", trunc(fit$eig[3]/sum(fit$eig) * 100 * 10^2)/10^2, "%)")) +
    ggtitle(title)
  plotD1D3 <- ggplot(data = fit_and_meta_df, aes(x = V1, y = V3, label = AbbreviatedName, shape = `Top Level`, color = Haplotype)) +
    scale_color_manual(values=haplotype_colors) +
    geom_point(size=2) + 
    xlab(paste0("Dimension 1 (", trunc(fit$eig[1]/sum(fit$eig) * 100 * 10^2)/10^2, "%)"))+
    ylab(paste0("Dimension 3 (", trunc(fit$eig[3]/sum(fit$eig) * 100 * 10^2)/10^2, "%)")) +
    ggtitle(title)
  plotD1D4 <- ggplot(data = fit_and_meta_df, aes(x = V1, y = V4, label = AbbreviatedName, shape = `Top Level`, color = Haplotype)) +
    scale_color_manual(values=haplotype_colors) +
    geom_point(size=2) + 
    xlab(paste0("Dimension 1 (", trunc(fit$eig[1]/sum(fit$eig) * 100 * 10^2)/10^2, "%)"))+
    ylab(paste0("Dimension 4 (", trunc(fit$eig[4]/sum(fit$eig) * 100 * 10^2)/10^2, "%)")) +
    ggtitle(title)
  plotD1D5 <- ggplot(data = fit_and_meta_df, aes(x = V1, y = V5, label = AbbreviatedName, shape = `Top Level`, color = Haplotype)) +
    scale_color_manual(values=haplotype_colors) +
    geom_point(size=2) + 
    xlab(paste0("Dimension 1 (", trunc(fit$eig[1]/sum(fit$eig) * 100 * 10^2)/10^2, "%)"))+
    ylab(paste0("Dimension 5 (", trunc(fit$eig[5]/sum(fit$eig) * 100 * 10^2)/10^2, "%)")) +
    ggtitle(title)
  
  if (N %in%  c('All')) {
    ggsave(plot = plotD1D2, imageD12, width = 21, height = 14,  units = "cm", dpi = 300,  bg = "transparent")
    ggsave(plot = plotD2D3, imageD23, width = 21, height = 14,  units = "cm", dpi = 300,  bg = "transparent")
    ggsave(plot = plotD1D3, imageD13, width = 21, height = 14,  units = "cm", dpi = 300,  bg = "transparent")
    ggsave(plot = plotD1D4, imageD14, width = 21, height = 14,  units = "cm", dpi = 300,  bg = "transparent")
    ggsave(plot = plotD1D5, imageD15, width = 21, height = 14,  units = "cm", dpi = 300,  bg = "transparent")
  }

  plotD1D2 <- plotD1D2 + theme(plot.title = element_text(hjust = 0.5)) + geom_text_repel(
    size=3,
    max.iter=40000,
    max.time=2,
    show.legend  = FALSE, # to hide the `a` from the legend
    max.overlaps=Inf
  )
  plotD2D3 <- plotD2D3 + theme(plot.title = element_text(hjust = 0.5)) + geom_text_repel(
    size=3,
    max.iter=40000,
    max.time=2,
    show.legend  = FALSE, # to hide the `a` from the legend
    max.overlaps=Inf
  )
  plotD1D3 <- plotD1D3 + theme(plot.title = element_text(hjust = 0.5)) + geom_text_repel(
    size=3,
    max.iter=40000,
    max.time=2,
    show.legend  = FALSE, # to hide the `a` from the legend
    max.overlaps=Inf
  )
  plotD1D4 <- plotD1D4 + theme(plot.title = element_text(hjust = 0.5)) + geom_text_repel(
    size=3,
    max.iter=40000,
    max.time=2,
    show.legend  = FALSE, # to hide the `a` from the legend
    max.overlaps=Inf
  )
  plotD1D5 <- plotD1D5 + theme(plot.title = element_text(hjust = 0.5)) + geom_text_repel(
    size=3,
    max.iter=40000,
    max.time=2,
    show.legend  = FALSE, # to hide the `a` from the legend
    max.overlaps=Inf
  )
  ggsave(plot = plotD1D2, imageD12_withLabels, width = 21, height = 14,  units = "cm", dpi = 300,  bg = "transparent")
  ggsave(plot = plotD2D3, imageD23_withLabels, width = 21, height = 14,  units = "cm", dpi = 300,  bg = "transparent")
  ggsave(plot = plotD1D3, imageD13_withLabels, width = 21, height = 14,  units = "cm", dpi = 300,  bg = "transparent")
  ggsave(plot = plotD1D4, imageD14_withLabels, width = 21, height = 14,  units = "cm", dpi = 300,  bg = "transparent")
  ggsave(plot = plotD1D5, imageD15_withLabels, width = 21, height = 14,  units = "cm", dpi = 300,  bg = "transparent")
}
################################################################################

################################################################################
# Systematic analysis of PCA's principal components
#==================================================
dir.create(file.path(base_dir, 'ClassicalMultidimensionalScaling', 'NoUtgs.NoAltExceptPeregrine', 'BestPairsMatPatSeparationInPCA'), recursive = T)

## Write PCA matrixes
for (N in c('All', 'XY', '1to22', seq(1, 22))){
  # Input matrix of distances
  if (N == 'All') {
    path_chrN <- file.path(base_dir, paste0('HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.chr', N, '.dist.tsv'))
  }else{
    path_chrN <- file.path(base_dir, paste0('HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.Y.x100.chr', N, '.dist.tsv'))
  }
  
  # Read the matrix
  HG002_all <- read.table(path_chrN, sep = '\t', header = T)
  
  HG002_all_filtered <- HG002_all %>%
    filter(!(
      `group.a` %in% c('Peregrine_HiFi_20kb.alt.utgs', 'Peregrine_HiFi_25kb.alt.utgs', 'FALCON_Unzip.alt', 'HiCanu.alt') | 
        `group.b` %in% c('Peregrine_HiFi_20kb.alt.utgs', 'Peregrine_HiFi_25kb.alt.utgs', 'FALCON_Unzip.alt', 'HiCanu.alt')
    ))
  prefix <- 'CMS.LongAlignments.NoUtgs.NoAltExceptPeregrine.chr'
  
  # Classical (Metric) Multidimensional Scaling
  HG002_all_wide_euclidean <- tbl2hmap(HG002_all_filtered, 'group.a', 'group.b', 'euclidean')
  
  fit <- cmdscale(as.dist(HG002_all_wide_euclidean$mat), eig=TRUE, k = nrow(HG002_all_wide_euclidean$mat) - 1)
  
  fit_df <- as.data.frame(fit$points) %>%
    rownames_to_column(var = "AbbreviatedName")
  fit_and_meta_df <- merge(fit_df, meta_df %>% select(AbbreviatedName, Haplotype)) 
  fit_and_meta_mat_pat_df <- fit_and_meta_df %>% filter(Haplotype == 'maternal' | Haplotype == 'paternal')
  
  path_PCA_matrix <- file.path(base_dir, paste0(prefix, N, '.PCA.tsv'))
  print(path_PCA_matrix)
  write.table(fit_and_meta_mat_pat_df,path_PCA_matrix, sep = '\t', quote = F, row.names = F)
}

# Go to the `HG002_bakeoff_SVM_on_PCA.ipynb` notebook and produce the 'BestPairsMatPatSeparationInPCA.tsv' file

BestPairsMatPatSeparationInPCA_df <-read.table(
  file.path(base_dir, 'BestPairsMatPatSeparationInPCA.tsv'),
  sep = '\t',
  header = T
)

for(i in 1:nrow(BestPairsMatPatSeparationInPCA_df)) {
  row <- BestPairsMatPatSeparationInPCA_df[i,]
  
  N <- as.character(row$N)
  x <- row$PCx
  y <- row$PCy
  slope <- row$slope
  intercept <- row$intercept
  print(N)

  # Input matrix of distances
  if (N == 'All') {
    path_chrN <- file.path(base_dir, paste0('HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.chr', N, '.dist.tsv'))
  }else{
    path_chrN <- file.path(base_dir, paste0('HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.Y.x100.chr', N, '.dist.tsv'))
  }
  
  # Read the matrix
  HG002_all <- read.table(path_chrN, sep = '\t', header = T)
  
  HG002_all_filtered <- HG002_all %>%
    filter(!(
      `group.a` %in% c('Peregrine_HiFi_20kb.alt.utgs', 'Peregrine_HiFi_25kb.alt.utgs', 'FALCON_Unzip.alt', 'HiCanu.alt') | 
        `group.b` %in% c('Peregrine_HiFi_20kb.alt.utgs', 'Peregrine_HiFi_25kb.alt.utgs', 'FALCON_Unzip.alt', 'HiCanu.alt')
    ))
  title <- paste0('No Utgs and No alt except Peregrine - chr', N)
  CMS_dir <- file.path(base_dir, 'ClassicalMultidimensionalScaling', 'NoUtgs.NoAltExceptPeregrine', 'BestPairsMatPatSeparationInPCA')
  prefix <- 'CMS.LongAlignments.NoUtgs.NoAltExceptPeregrine.chr'
  
  imageDXDY <- file.path(CMS_dir, paste0(prefix, N, '.BestMatPatSeparation.D', x, 'vsD', y, '.noLabels.pdf'))
  imageDXDY_withLabels <- file.path(CMS_dir, paste0(prefix, N, '.BestMatPatSeparation.D', x, 'vsD', y, '.withLabels.pdf'))

  # Classical (Metric) Multidimensional Scaling
  HG002_all_wide_euclidean <- tbl2hmap(HG002_all_filtered, 'group.a', 'group.b', 'euclidean')
  
  fit <- cmdscale(as.dist(HG002_all_wide_euclidean$mat), eig=TRUE, k = nrow(HG002_all_wide_euclidean$mat) - 1)
  
  fit_df <- as.data.frame(fit$points) %>%
    rownames_to_column(var = "AbbreviatedName")
  fit_and_meta_df <- merge(fit_df, meta_df) %>%
    rename(`Top Level` = TopLevel) %>%
    rename(`Contig Algorithm` = ContigAlgorithmBroader)
  
  #n=9
  PCx <- paste0('V', x)
  PCy <- paste0('V', y)
  plotDXDY <- ggplot(data = fit_and_meta_df, aes_string(x = PCx, y = PCy, label = "`AbbreviatedName`", shape = "`Top Level`", color = "`Haplotype`")) +
    scale_color_manual(values=haplotype_colors) +
    geom_point(size=2) + 
    xlab(paste0("Dimension ", x, " (", trunc(fit$eig[x]/sum(fit$eig) * 100 * 10^2)/10^2, "%)")) +
    ylab(paste0("Dimension ", y, " (", trunc(fit$eig[y]/sum(fit$eig) * 100 * 10^2)/10^2, "%)")) +
    ggtitle(title) + 
    geom_abline(slope = slope, intercept = intercept, size = 0.2)
  plotDXDY
  ggsave(plot = plotDXDY, imageDXDY, width = 21, height = 14,  units = "cm", dpi = 300,  bg = "transparent")
  
  plotDXDY <- plotDXDY + theme(plot.title = element_text(hjust = 0.5)) + geom_text_repel(
    size=3,
    max.iter=40000,
    max.time=2,
    show.legend  = FALSE, # to hide the `a` from the legend
    max.overlaps=Inf
  )
  ggsave(plot = plotDXDY, imageDXDY_withLabels, width = 21, height = 14,  units = "cm", dpi = 300,  bg = "transparent")
}
################################################################################


# Old tests
if (FALSE) {
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
}
