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
# remove ambient RNA using SoupX
for (i in 1:length(celldata_path_lst)) {
  print(i)
  toc = Seurat::Read10X(celldata_path_lst[[i]])
  tod = Seurat::Read10X(fulldata_path_lst[[i]])
  cluster_info = read.csv(cluster_path_lst[[i]])
  # instantiate soup channel and populate with cluster information
  sc = SoupChannel(tod, toc)
  sc = setClusters(sc, setNames(cluster_info$Cluster, rownames(cluster_info)))
  sc = autoEstCont(sc)
  # remove ambient RNA from cell data
  decont_celldata = adjustCounts(sc)
  # create Seurat object using post-QC count matrix
  sample_name = sample_names[[i]]
  s_obj = CreateSeuratObject(counts = decont_celldata, project = sample_name, min.cells = 3, min.features = 200)
  seurat_objs_lst[[sample_name]] = s_obj#c(seurat_objs_lst[[sample_name]], s_obj)
}
# merge seurat objects for cell-quality filtering
huati = merge(seurat_objs_lst[[1]], y = seurat_objs_lst[2:length(seurat_objs_lst)], 
              add.cell.ids = names(seurat_objs_lst), project = "Rutishauser_huATI")
huati[["pct.mito"]] = PercentageFeatureSet(huati, pattern="^MT-")
table(huati$orig.ident)
# filter out low-quality cells based on mitochondrial expression + RNA (depth) count + gene (feature) count
s_huati <- subset(huati, pct.mito <= 20 & nCount_RNA >= 1500 & nFeature_RNA >= 800)
table(s_huati$orig.ident)
s_huati.split <- SplitObject(huati, split.by = "orig.ident")
length(s_huati.split)
# remove doublets
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
# merge sample data
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
#p1 <- DimPlot(object = s_huati, reduction = "pca", pt.size = .1, group.by = "orig.ident")
#p2 <- VlnPlot(object = s_huati, features = "PC_1", group.by = "orig.ident", pt.size = .1)
#plot_grid(p1,p2)
# integrate datasets with batch correction, visualize convergence of integration
options(repr.plot.height = 2.5, repr.plot.width = 6)
sc_huati <- s_huati %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)
# visualize differences between datasets post-correction
options(repr.plot.height = 5, repr.plot.width = 12)
#p1 <- DimPlot(object = sc_huati, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
#p2 <- VlnPlot(object = sc_huati, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
#plot_grid(p1,p2)
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
#DimPlot(sc_huati, reduction = "umap", group.by = "orig.ident", pt.size = .1, split.by = 'orig.ident')
# visualize UMAP embeddings of integrated data
options(repr.plot.height = 4, repr.plot.width = 6)
#DimPlot(sc_huati, reduction = "umap", label = TRUE, pt.size = .1)
# save integrated, processed object
output_seurat_rds_path = "/krummellab/data1/alaa/data/hiv/rutishauser/outputs/s_huati_seurat_integrated_v2.rds"
sc_huati = readRDS(output_seurat_rds_path)
table(sc_huati$orig.ident)
table(sc_huati$seurat_clusters)
#saveRDS(sc_huati, output_seurat_rds_path)
# save metadata of processed GEX data
output_gex_metapath = "/krummellab/data1/alaa/data/hiv/rutishauser/outputs/s_huati_seurat_integrated_meta_v3.csv"
s_huati_meta = sc_huati[[]]
write.csv(s_huati_meta, file = output_gex_metapath)
# subset on T-cells
output_t_metapath = "/krummellab/data1/alaa/data/hiv/rutishauser/outputs/tcells_s_huati_seurat_integrated_meta_v3.csv"
rownames(sc_huati[grep("CD3", rownames(sc_huati)),])
t_s_huati = subset(sc_huati, CD8A > 0 | CD8B > 0 | CD4 > 0 | CD3E > 0 | CD3D > 0 | CD3G > 0 | CD3E > 0)
t_s_huati_meta = t_s_huati[[]]
write.csv(t_s_huati_meta, file = output_t_metapath)
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
sc_huati.markers = FindAllMarkers(sc_huati, only.pos = TRUE, 
                                  min.pct = 0.25, 
                                  logfc.threshold = 0.25)
print(sc_huati.markers %>%
        group_by(cluster) %>%
        slice_max(n = 9, order_by = avg_log2FC), n=162)
# subset for CD8s only
output_cd8_metapath = "/krummellab/data1/alaa/data/hiv/rutishauser/outputs/cd8_s_huati_seurat_integrated_meta_v3.csv"
#rownames(sc_huati[grep("CD3", rownames(sc_huati)),])
cd8_s_huati = subset(sc_huati, CD8A > 0 | CD8B > 0)
table(cd8_s_huati$orig.ident)
table(cd8_s_huati$seurat_clusters)
rownames(cd8_s_huati)
rownames(cd8_s_huati[grep("CD4", rownames(cd8_s_huati)),])
# save integrated, processed data consisting of CD8 subset by strict gene expression
output_cd8_seurat_rds_path = "/krummellab/data1/alaa/data/hiv/rutishauser/outputs/cd8_s_huati_seurat_integrated_v3.rds"
#saveRDS(cd8_s_huati, output_cd8_seurat_rds_path)
#cd8_sc_huati = subset(sc_huati, subset = CD3 > 1)
# visualize UMAP embeddings of integrated data -- CD8s ONLY
options(repr.plot.height = 4, repr.plot.width = 6)
DimPlot(cd8_s_huati, reduction = "umap", label = TRUE, pt.size = .1)
cd8_s_huati_meta = cd8_s_huati[[]]
unique(cd8_s_huati_meta$orig.ident)
cd8_s_huati_meta
write.csv(cd8_s_huati_meta, file = output_cd8_metapath)
