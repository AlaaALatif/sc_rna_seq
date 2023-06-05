# Install required packages
#install.packages("Seurat")
#install.packages("SeuratObject)
#install.packages("remotes")
#remotes::install_github("mojaveazure/seurat-disk")
#remotes::install_github("LungCellAtlas/FastCAR")
#install.packages("tidyverse")
#install.packages("harmony")
#ensure you module load gsl/2.6+ before rsc start
#install.packages("gsl")
#BiocManager::install("powerTCR")
#devtools::install_github("ncborcherding/scRepertoire")
#install.packages("SoupX")
#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(ggplot2)
library(tidyverse)
library(Matrix)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(qlcMatrix)
library(FastCAR)
library(scales)
library(harmony)
library(cowplot)
library(purrr)
library(plyr)
library(patchwork)
#library(scRepertoire)
library(SoupX)
library(DoubletFinder)
sample_ids = c(1,2,3,5,6,7,8,9)
sample_names = c("NC1Pre", "NC2Pre", "EC1Pre", "EC2Pre", "NC2Post","NC3Post", "EC1Post", "EC2Post")

# read input GEX data collectively using lists
celldata_path_lst = c('/krummellab/data1/alaa/data/hiv/rutishauser/rutishauser_ATI_gi_RR3692_01/outs/per_sample_outs/rutishauser_ATI_gi_RR3692_01/count/sample_filtered_feature_bc_matrix',
                      '/krummellab/data1/alaa/data/hiv/rutishauser/rutishauser_ATI_gi_RR3692_02/outs/per_sample_outs/rutishauser_ATI_gi_RR3692_02/count/sample_filtered_feature_bc_matrix',
                      '/krummellab/data1/alaa/data/hiv/rutishauser/rutishauser_ATI_gi_RR3692_03/outs/per_sample_outs/rutishauser_ATI_gi_RR3692_03/count/sample_filtered_feature_bc_matrix',
                      '/krummellab/data1/alaa/data/hiv/rutishauser/rutishauser_ATI_gi_RR3692_05/outs/per_sample_outs/rutishauser_ATI_gi_RR3692_05/count/sample_filtered_feature_bc_matrix',
                      '/krummellab/data1/alaa/data/hiv/rutishauser/rutishauser_ATI_gi_RR3692_06/outs/per_sample_outs/rutishauser_ATI_gi_RR3692_06/count/sample_filtered_feature_bc_matrix',
                      '/krummellab/data1/alaa/data/hiv/rutishauser/rutishauser_ATI_gi_RR3692_07/outs/per_sample_outs/rutishauser_ATI_gi_RR3692_07/count/sample_filtered_feature_bc_matrix',
                      '/krummellab/data1/alaa/data/hiv/rutishauser/rutishauser_ATI_gi_RR3692_08/outs/per_sample_outs/rutishauser_ATI_gi_RR3692_08/count/sample_filtered_feature_bc_matrix',
                      '/krummellab/data1/alaa/data/hiv/rutishauser/rutishauser_ATI_gi_RR3692_09/outs/per_sample_outs/rutishauser_ATI_gi_RR3692_09/count/sample_filtered_feature_bc_matrix')
fulldata_path_lst = c('/krummellab/data1/alaa/data/hiv/rutishauser/rutishauser_ATI_gi_RR3692_01/outs/multi/count/raw_feature_bc_matrix',
                      '/krummellab/data1/alaa/data/hiv/rutishauser/rutishauser_ATI_gi_RR3692_02/outs/multi/count/raw_feature_bc_matrix',
                      '/krummellab/data1/alaa/data/hiv/rutishauser/rutishauser_ATI_gi_RR3692_03/outs/multi/count/raw_feature_bc_matrix',
                      '/krummellab/data1/alaa/data/hiv/rutishauser/rutishauser_ATI_gi_RR3692_05/outs/multi/count/raw_feature_bc_matrix',
                      '/krummellab/data1/alaa/data/hiv/rutishauser/rutishauser_ATI_gi_RR3692_06/outs/multi/count/raw_feature_bc_matrix',
                      '/krummellab/data1/alaa/data/hiv/rutishauser/rutishauser_ATI_gi_RR3692_07/outs/multi/count/raw_feature_bc_matrix',
                      '/krummellab/data1/alaa/data/hiv/rutishauser/rutishauser_ATI_gi_RR3692_08/outs/multi/count/raw_feature_bc_matrix',
                      '/krummellab/data1/alaa/data/hiv/rutishauser/rutishauser_ATI_gi_RR3692_09/outs/multi/count/raw_feature_bc_matrix')
cluster_path_lst = c("/krummellab/data1/alaa/data/hiv/rutishauser/rutishauser_ATI_gi_RR3692_01/outs/per_sample_outs/rutishauser_ATI_gi_RR3692_01/count/analysis/clustering/gene_expression_graphclust/clusters.csv",
                     "/krummellab/data1/alaa/data/hiv/rutishauser/rutishauser_ATI_gi_RR3692_02/outs/per_sample_outs/rutishauser_ATI_gi_RR3692_02/count/analysis/clustering/gene_expression_graphclust/clusters.csv",
                     "/krummellab/data1/alaa/data/hiv/rutishauser/rutishauser_ATI_gi_RR3692_03/outs/per_sample_outs/rutishauser_ATI_gi_RR3692_03/count/analysis/clustering/gene_expression_graphclust/clusters.csv",
                     "/krummellab/data1/alaa/data/hiv/rutishauser/rutishauser_ATI_gi_RR3692_05/outs/per_sample_outs/rutishauser_ATI_gi_RR3692_05/count/analysis/clustering/gene_expression_graphclust/clusters.csv",
                     "/krummellab/data1/alaa/data/hiv/rutishauser/rutishauser_ATI_gi_RR3692_06/outs/per_sample_outs/rutishauser_ATI_gi_RR3692_06/count/analysis/clustering/gene_expression_graphclust/clusters.csv",
                     "/krummellab/data1/alaa/data/hiv/rutishauser/rutishauser_ATI_gi_RR3692_07/outs/per_sample_outs/rutishauser_ATI_gi_RR3692_07/count/analysis/clustering/gene_expression_graphclust/clusters.csv",
                     "/krummellab/data1/alaa/data/hiv/rutishauser/rutishauser_ATI_gi_RR3692_08/outs/per_sample_outs/rutishauser_ATI_gi_RR3692_08/count/analysis/clustering/gene_expression_graphclust/clusters.csv",
                     "/krummellab/data1/alaa/data/hiv/rutishauser/rutishauser_ATI_gi_RR3692_09/outs/per_sample_outs/rutishauser_ATI_gi_RR3692_09/count/analysis/clustering/gene_expression_graphclust/clusters.csv")
seurat_objs_lst <- list()
#toc = Seurat::Read10X(celldata_path_lst[[i]])
#tod = Seurat::Read10X(fulldata_path_lst[[i]])
#cluster_info = read.csv(cluster_path_lst[[i]])
#projection_info = read.csv("/krummellab/data1/alaa/data/hiv/rutishauser/rutishauser_ATI_gi_RR3692_01/outs/per_sample_outs/rutishauser_ATI_gi_RR3692_01/count/analysis/umap/gene_expression_2_components/projection.csv")
#sc = SoupChannel(tod, toc)
#sc = setClusters(sc, setNames(cluster_info$Cluster, rownames(cluster_info)))
#sc = autoEstCont(sc)
#decont_counts = adjustCounts(sc)
#contaminationChanceCutoff = 0.05
#emptyDropletCutoff = 175
#toc = Seurat::Read10X(celldata_path_lst[[i]])
#tod = Seurat::Read10X(fulldata_path_lst[[i]])
#cluster_info = read.csv(cluster_path_lst[[i]])
# generate ambient RNA profile for removal
#sc <- SoupChannel(toc, tod)
#sc = setClusters(sc, setNames(cluster_info$Cluster, rownames(cluster_info)))
for (i in 1:length(celldata_path_lst)) {
  print(i)
  toc = Seurat::Read10X(celldata_path_lst[[i]])
  tod = Seurat::Read10X(fulldata_path_lst[[i]])
  cluster_info = read.csv(cluster_path_lst[[i]])
  #projection_info = read.csv("/krummellab/data1/alaa/data/hiv/rutishauser/rutishauser_ATI_gi_RR3692_01/outs/per_sample_outs/rutishauser_ATI_gi_RR3692_01/count/analysis/umap/gene_expression_2_components/projection.csv")
  sc = SoupChannel(tod, toc)
  sc = setClusters(sc, setNames(cluster_info$Cluster, rownames(cluster_info)))
  sc = autoEstCont(sc)
  #ambient_profile = determine.background.to.remove(fullCellMatrix = fulldata, cellMatrix = celldata, 
  #                                                 emptyDropletCutoff = emptyDropletCutoff, 
  #                                                 contaminationChanceCutoff = contaminationChanceCutoff)
  # remove ambient RNA from cell data
  decont_celldata = adjustCounts(sc)
  #decont_celldata = remove.background(celldata, ambient_profile)
  # create Seurat object using post-QC count matrix
  sample_name = sample_names[[i]]
  s_obj = CreateSeuratObject(counts = decont_celldata, project = sample_name, min.cells = 3, min.features = 200)
  seurat_objs_lst[[sample_name]] = s_obj#c(seurat_objs_lst[[sample_name]], s_obj)
}
print(names(seurat_objs_lst))
print(seurat_objs_lst[[1]])
# read input data individually
#celldata_path1 = '/krummellab/data1/alaa/data/hiv/rutishauser/rutishauser_ATI_gi_RR3692_01/outs/per_sample_outs/rutishauser_ATI_gi_RR3692_01/count/sample_filtered_feature_bc_matrix'
#fulldata_path1 = '/krummellab/data1/alaa/data/hiv/rutishauser/rutishauser_ATI_gi_RR3692_01/outs/multi/count/raw_feature_bc_matrix'
#celldata_path2 = '/krummellab/data1/alaa/data/hiv/rutishauser/rutishauser_ATI_gi_RR3692_03/outs/per_sample_outs/rutishauser_ATI_gi_RR3692_03/count/sample_filtered_feature_bc_matrix'
#fulldata_path2 = '/krummellab/data1/alaa/data/hiv/rutishauser/rutishauser_ATI_gi_RR3692_03/outs/multi/count/raw_feature_bc_matrix'
#list.files(celldata_path1)
#list.files(celldata_path2)
# contaminant (ambient RNA) removal
#cellMatrix1 = read.cell.matrix(celldata_path1)
#fullMatrix1 = read.full.matrix(fulldata_path1)
#cellMatrix2 = read.cell.matrix(celldata_path2)
#fullMatrix2 = read.full.matrix(fulldata_path2)
# generate plots of ambient RNA profile
#ambProfile = describe.ambient.RNA.sequence(fullCellMatrix = fullMatrix, 
#                                           start = 10, 
#                                           stop = 500, 
#                                           by = 10, 
#                                           contaminationChanceCutoff = 0.05)
#plot.ambient.profile(ambProfile)
# specify hyperparameters
#recommend.empty.cutoff(ambProfile)
#contaminationChanceCutoff = 0.05
#emptyDropletCutoff = 175
# generate ambient RNA profile for remval
#ambient_profile1 = determine.background.to.remove(fullCellMatrix = fullMatrix1, cellMatrix = cellMatrix1, emptyDropletCutoff = emptyDropletCutoff, contaminationChanceCutoff = contaminationChanceCutoff)
#ambient_profile2 = determine.background.to.remove(fullCellMatrix = fullMatrix2, cellMatrix = cellMatrix2, emptyDropletCutoff = emptyDropletCutoff, contaminationChanceCutoff = contaminationChanceCutoff)
# remove ambient RNA from cell data
# cellMatrix1 = remove.background(cellMatrix1, ambient_profile1)
# cellMatrix2 = remove.background(cellMatrix2, ambient_profile2)
# create Seurat object for processing and analysis
# huati1 <- Read10X(data.dir = data_dir)
# huati1n2 = CreateSeuratObject(counts = cbind(cellMatrix1, cellMatrix2), project = 'huATI_1_and_2', min.cells = 3, min.features = 200)
# huati = CreateSeuratObject(counts = cell_matrices, project = 'huATI', min.cells = 3, min.features = 200)
#seurat_objs_lst
huati = merge(seurat_objs_lst[[1]], y = seurat_objs_lst[2:length(seurat_objs_lst)], 
              add.cell.ids = names(seurat_objs_lst), project = "Rutishauser_huATI")
huati[["pct.mito"]] = PercentageFeatureSet(huati, pattern="^MT-")
table(huati$orig.ident)
# filter out low-quality cells based on mitochondrial expression + RNA (depth) count + gene (feature) count
s_huati <- subset(huati, pct.mito <= 20 & nCount_RNA >= 1500 & nFeature_RNA >= 800)
table(s_huati$orig.ident)
s_huati.split <- SplitObject(huati, split.by = "orig.ident")
length(s_huati.split)
for (i in 1:length(s_huati.split)) {
  ex = s_huati.split[[i]]
  ex = NormalizeData(ex, normalization.method = "LogNormalize", scale.factor = 10000)
  ex = FindVariableFeatures(ex, selection.method="vst", nfeatures = 2000)
  all.genes <- rownames(ex)
  ex = ScaleData(ex, features=all.genes, verbose = TRUE)
  ex = RunPCA(ex, features = VariableFeatures(object = ex))
  # find significant PCs
  stdv = ex[["pca"]]@stdev
  sum.stdv = sum(stdv)
  pct.stdv = (stdv / sum.stdv) * 100
  cumulative = cumsum(pct.stdv)
  co1 = which(cumulative > 90 & pct.stdv < 5)[1]
  co2 = sort(which((pct.stdv[1:length(pct.stdv) - 1] - pct.stdv[2:length(pct.stdv)]) > 0.1),
             decreasing=T)[1] + 1
  min.pc = min(co1, co2)
  print(min.pc)
  # dimensionality reduction and clustering
  ex = RunUMAP(ex, dims = 1:min.pc)
  ex = FindNeighbors(object = ex, dims = 1:min.pc)
  ex = FindClusters(object = ex, resolution=c(0.5))
  # pK identification -- computing bimodality coefficient distribution
  sweep.list = paramSweep_v3(ex, PCs = 1:min.pc, num.cores=34)
  sweep.stats = summarizeSweep(sweep.list)
  bcmvn = find.pK(sweep.stats)
  # identify optimal pK from distribution
  bcmvn.max = bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pK = bcmvn.max$pK
  optimal.pK = as.numeric(levels(optimal.pK))[optimal.pK]
  print(optimal.pK)
  # homotypic doublet proportion estimate
  cluster_assignments = ex@meta.data$seurat_clusters
  homotypic.prop = modelHomotypic(cluster_assignments)
  # assume 7.5% doublet formation rate
  nExp.poi = round(0.075 * nrow(ex@meta.data))
  nExp.poi.adj = round(nExp.poi * (1 - homotypic.prop))
  print(nExp.poi.adj)
  # remove doublets
  ex = doubletFinder_v3(seu = ex,
                        PCs = 1:min.pc,
                        pK = optimal.pK,
                        nExp = nExp.poi.adj)
  metadata = ex@meta.data
  colnames(metadata)[8] = "doublet_finder"
  ex@meta.data = metadata
  # subset singlets and save
  ex_singlets = subset(ex, doublet_finder == "Singlet")
  s_huati.split[[i]] = ex_singlets
  remove(ex_singlets)
}
s_huati = merge(s_huati.split[[1]], y = s_huati.split[2:length(s_huati.split)], 
              add.cell.ids = names(s_huati.split), project = "Rutishauser_huATI")
s_huati[["pct.mito"]] = PercentageFeatureSet(s_huati, pattern="^MT-")
table(s_huati$orig.ident)
# Run downstream processing without integration
s_huati = NormalizeData(s_huati, normalization.method = "LogNormalize", scale.factor = 10000)
s_huati = FindVariableFeatures(s_huati, selection.method="vst", nfeatures = 2000)
all.genes <- rownames(s_huati)
s_huati = ScaleData(s_huati, features=all.genes, verbose = TRUE)
s_huati = RunPCA(s_huati, features = VariableFeatures(object = s_huati))
# visualize differences between datasets in raw PCs
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = s_huati, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = s_huati, features = "PC_1", group.by = "orig.ident", pt.size = .1)
plot_grid(p1,p2)
# integrate datasets with batch correction, visualize convergence of integration
options(repr.plot.height = 2.5, repr.plot.width = 6)
sc_huati <- s_huati %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)
# visualize differences between datasets post-correction
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = sc_huati, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = sc_huati, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
plot_grid(p1,p2)
# perform downstream analysis using corrected harmony embeddings
use.corrected_pcs = 1:15
# run UMAP dimensionality reduction
sc_huati = sc_huati %>% RunUMAP(reduction = "harmony", dims = use.corrected_pcs)
# find nearest neighbors
sc_huati = FindNeighbors(sc_huati, reduction = 'harmony', dims = use.corrected_pcs)
# assign clusters
sc_huati = FindClusters(sc_huati, resolution = c(0.5))
# visualize UMAP embeddings of each corrected data
options(repr.plot.height = 4, repr.plot.width = 10)
DimPlot(sc_huati, reduction = "umap", group.by = "orig.ident", pt.size = .1, split.by = 'orig.ident')
# visualize UMAP embeddings of integrated data
options(repr.plot.height = 4, repr.plot.width = 6)
DimPlot(sc_huati, reduction = "umap", label = TRUE, pt.size = .1)
output_seurat_rds_path = "/krummellab/data1/alaa/data/hiv/rutishauser/outputs/s_huati_seurat_integrated_v2.rds"
saveRDS(sc_huati, output_seurat_rds_path)
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
sc_huati.markers = FindAllMarkers(sc_huati, only.pos = TRUE, 
                                  min.pct = 0.25, 
                                  logfc.threshold = 0.25)
print(sc_huati.markers %>%
  group_by(cluster) %>%
  slice_max(n = 9, order_by = avg_log2FC), n=162)

SaveH5Seurat(sc_huati, output_seurat_path)

# subset for CD8s only
output_cd8_metapath = "/krummellab/data1/alaa/data/hiv/rutishauser/outputs/cd8_s_huati_subset_seurat_integrated_meta_v2.csv"
rownames(sc_huati[grep("CD3", rownames(sc_huati)),])
cd8_s_huati = subset(sc_huati, seurat_clusters %in% c(3, 5))
rownames(cd8_s_huati)
#cd8_sc_huati = subset(sc_huati, subset = CD3 > 1)
# visualize UMAP embeddings of integrated data -- CD8s ONLY
options(repr.plot.height = 4, repr.plot.width = 6)
DimPlot(cd8_s_huati, reduction = "umap", label = TRUE, pt.size = .1)
cd8_s_huati_meta = cd8_s_huati[[]]
unique(cd8_s_huati_meta$orig.ident)
write.csv(cd8_s_huati_meta, file = output_cd8_metapath)

#huati.split
#huati = merge(seurat_objs_lst[[1]][[1]], y = seurat_objs_lst[[2]], add.cell.ids = c("S1", "S2"), project = "Rutishauser_huATI")
# QC, removal of low-quality/dying cells; characterized by extensive mitochondrial content
# calculate % counts originating from mito-related set of features
#huati1[["pct.mito"]] = PercentageFeatureSet(huati1, pattern="^MT-")
#huati1n2[["pct.mito"]] = PercentageFeatureSet(huati1n2, pattern="^MT-")
# visualize QC metrics independently using Violin plots
#VlnPlot(huati, features = c("nFeature_RNA", "nCount_RNA", "pct.mito"))
# visualize QC metrics co-dependently using scatter plots
#count_v_mito_plot = FeatureScatter(huati, feature1 = "nCount_RNA", feature2 = "pct.mito")
# + theme(legend.position = "none")
#count_v_mito_plot = count_v_mito_plot + labs(title = NULL)
#count_v_genes_plot = FeatureScatter(huati, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# + theme(legend.position = "none")
#count_v_genes_plot = count_v_genes_plot + labs(title = NULL)
#qc_plot2 <- count_v_mito_plot + count_v_genes_plot
#qc_plot2 + plot_annotation(title = 'QC Metrics Correlation', subtitle = 'rutishauser_huATI')
# generate tible of QC metrics
#qc_metrics = as_tibble(huati[[]], rownames = "Cell.Barcode")
# histogram of count depths
#qc_metrics %>%
#  ggplot(aes(nCount_RNA)) + 
#  geom_histogram(binwidth = 50, fill="white", colour="black") +
#  geom_vline(aes(xintercept=mean(nCount_RNA)), color="purple", linetype="dashed") +
#  ggtitle("Distribution of Count Depths")
# zoomed-in histogram of RNA count depths per cell
#qc_metrics %>%
#  ggplot(aes(nCount_RNA)) + 
#  geom_histogram(binwidth = 25, fill="white", colour="black") +
#  geom_vline(aes(xintercept=mean(nCount_RNA)), color="purple", linetype="dashed") +
#  xlim(c(0, 6500)) +
#  ggtitle("Distribution of Depth (molecule) per Cell")
# histogram of gene counts per cell
#qc_metrics %>%
#  ggplot(aes(nFeature_RNA)) + 
#  geom_histogram(binwidth = 5, fill="white", colour="black") +
#  geom_vline(aes(xintercept=mean(nFeature_RNA)), color="black", linetype="dashed") +
#  annotate("text", x=2200, y=60, label="mean") +
#  geom_vline(aes(xintercept=800), color="red", label="cutoff") +
#  annotate("text", x=500, y=60, label="cutoff") +
#  ggtitle("Distribution of Gene (feature) Counts per Cell")
# relationship between QC metrics
#qc_metrics %>%
#  ggplot(aes(y=nFeature_RNA, x=nCount_RNA, colour=pct.mito)) + 
#  geom_point() +
#  scale_x_continuous(breaks = scales::pretty_breaks(n = 12)) +
#  scale_y_continuous(breaks = scales::pretty_breaks(n = 12)) +
#  scale_color_gradientn(colours=c('darkblue', 'darkgreen', 'yellow'), limits=c(0,30), oob=scales::squish) +
#  ggtitle("Relationship between depth (mRNA) count, gene count, and fractional mitochondrial expression levels") +
#  theme(plot.title = element_text(size = 8, face = "bold"))
# Add clonetype information
add_clonotype <- function(tcr_prefix, seurat_obj, type="t"){
  tcr <- read.csv(paste(tcr_prefix,"filtered_contig_annotations.csv", sep="/"))
  
  # Remove the -1 at the end of each barcode.
  # Subsets so only the first line of each barcode is kept,
  # as each entry for given barcode will have same clonotype.
  tcr <- tcr[!duplicated(tcr$barcode), ]
  
  # Only keep the barcode and clonotype columns. 
  # We'll get additional clonotype info from the clonotype table.
  tcr <- tcr[,c("barcode", "raw_clonotype_id")]
  names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"
  
  # Clonotype-centric info.
  clono <- read.csv(paste(tcr_prefix,"clonotypes.csv", sep="/"))
  
  # Slap the AA sequences onto our original table by clonotype_id.
  tcr <- merge(tcr, clono[, c("clonotype_id", "cdr3s_aa")])
  names(tcr)[names(tcr) == "cdr3s_aa"] <- "cdr3s_aa"
  
  # Reorder so barcodes are first column and set them as rownames.
  tcr <- tcr[, c(2,1,3)]
  rownames(tcr) <- tcr[,1]
  tcr[,1] <- NULL
  colnames(tcr) <- paste(type, colnames(tcr), sep="_")
  # Add to the Seurat object's metadata.
  clono_seurat <- AddMetaData(object=seurat_obj, metadata=tcr)
  return(clono_seurat)
}
#tcr_data_dir <- '/krummellab/data1/alaa/data/hiv/rutishauser/rutishauser_ATI_gi_RR3692_01/outs/per_sample_outs/rutishauser_ATI_gi_RR3692_01/vdj_t'
#s_huati1 <- add_clonotype(tcr_data_dir, s_huati1, "t")
# T cell numbers identified
#table(!is.na(s_huati1$t_clonotype_id))
# subset for T cells only
#s_huati1 <- subset(s_huati1, cells=colnames(s_huati1)[(!is.na(s_huati1$t_clonotype_id))])
# visualize RNA molecule counts
#RidgePlot(s_huati1, features = "nCount_RNA") + scale_x_continuous(limits = c(0, 9000))
# visualize RNA features (genes)
#RidgePlot(s_huati1, features = "nFeature_RNA") + scale_x_continuous(limits = c(0, 5000))
# visualize mitochondrial expression
#RidgePlot(s_huati1, features = "pct.mito") + scale_x_continuous(limits = c(0, 50))
# ???
#VlnPlot(s_huati1, features = c("nFeature_RNA", "nCount_RNA","pct.mito"), ncol = 1, pt.size = 0.3)
# scatter between RNA counts and mitochondrial percentage expression
#FeatureScatter(s_huati1, feature1="nCount_RNA", feature2="pct.mito")
# scatter between RNA counts and feature counts
#FeatureScatter(s_huati1, feature1="nCount_RNA", feature2="nFeature_RNA", pt.size=0.5)
# QC filters on pct.mito, RNA count, and Feature count
#s_huati1 <- subset(huati1, pct.mito <= 20 & nCount_RNA >= 1500 & nFeature_RNA >= 800)
#s_huati1n2 <- subset(huati1n2, pct.mito <= 20 & nCount_RNA >= 1500 & nFeature_RNA >= 800)
#s_huati <- subset(huati, pct.mito <= 20 & nCount_RNA >= 1500 & nFeature_RNA >= 800)
# expression normalization
#s_huati1 <- NormalizeData(s_huati1, normalization.method = "LogNormalize", scale.factor = 10000)
#s_huati1n2 <- NormalizeData(s_huati1n2, normalization.method = "LogNormalize", scale.factor = 10000)
s_huati <- NormalizeData(s_huati, normalization.method = "LogNormalize", scale.factor = 10000)
s_huati <- FindVariableFeatures(s_huati, selection.method="vst", nfeatures = 2000)
all.genes <- rownames(s_huati)
s_huati <- ScaleData(s_huati, features=all.genes, verbose = TRUE)
s_huati <- RunPCA(s_huati, features = VariableFeatures(object = s_huati))
# feature selection
#s_huati1 <- FindVariableFeatures(s_huati1, selection.method="vst", nfeatures = 2000)
#s_huati1n2 <- FindVariableFeatures(s_huati1n2, selection.method="vst", nfeatures = 2000)
# rescale data based on selected features
#all.genes <- rownames(s_huati1)
#s_huati1 <- ScaleData(s_huati1, features=all.genes)
#all.genes <- rownames(s_huati1n2)
#s_huati1n2 <- ScaleData(s_huati1n2, features=all.genes, verbose = TRUE)
# Identify principle components
#s_huati1 <- RunPCA(s_huati1, features = VariableFeatures(object = s_huati1))
#s_huati1n2 <- RunPCA(s_huati1n2, features = VariableFeatures(object = s_huati1n2))
# generate variable distinguishing batches
#s_huati1n2@meta.data$sample_id <- c(rep("S1", 10424), rep("S2", 6534))



# use principle components to find nearest neighbors (TODO: confirm with elbow/jackstraw to identify top PCs)
use.pcs = 1:15
s_huati1 <- FindNeighbors(s_huati1, dims = use.pcs)
# find clusters using Louvain algorithm
s_huati1 <- FindClusters(s_huati1, resolution = c(0.5))
# run UMAP for dimensionality reduction and visualization
s_huati1 <- RunUMAP(s_huati1, dims = use.pcs)
# generate UMAP plot
DimPlot(s_huati1, reduction = "umap", label = TRUE)
# TODO: QC, contaminants removal (SoupX, DecontX, FastCAR, CellBender)
# TODO: QC, doublet removal (Doubletfinder)
# TODO: norm, expression normalization (SCTransform, shifted log by Ahlmann-Eltze)
# TODO: agg, data integration (Harmony -> 4-CCA -> LIGER)
# TODO: select, feature selection