# Load required library
library(jsonlite)
library(ggplot2)
library(tidyverse)
library(Matrix)
library(Seurat)
library(SeuratObject)
library(scales)
library(harmony)
library(cowplot)
library(purrr)
library(plyr)
library(patchwork)
library(SoupX)
library(DoubletFinder)



# Read configuration
config <- fromJSON("../cfg/quality_control.json")
project_name = config$project_name
celldata_path = config$input_filtered_count_filepath
fulldata_path = config$input_raw_count_filepath
cluster_path = config$input_analysis_cluster_filepath
feature_selection_threshold = config$feature_selection_threshold
normalization_scale_factor = config$normalization_scale_factor
feature_qc_threshold = config$feature_qc_threshold
rna_qc_threshold = config$rna_qc_threshold
pct_mitochondrial_expression_qc_threshold = config$pct_mitochondrial_expression_qc_threshold
num_cores = config$num_cores
cluster_resolution = config$cluster_resolution
output_qc_plot_filepath = config$output_qc_plot_filepath
output_seurat_object_filepath = config$output_seurat_object_filepath
# Print data filepaths
print(celldata_path)
print(fulldata_path)
print(cluster_path)
# load count data
countdata = Seurat::Read10X(data.dir = celldata_path)
# create Seurat object from count data
celldata = CreateSeuratObject(counts = countdata, project = 'huATI_19')
# display initial number of cells
table(celldata$orig.ident)
# filter out ambient RNA using SoupX
fulldata = Seurat::Read10X(fulldata_path)
cluster_info = read.csv(cluster_path)
# instantiate soup channel and populate with cluster information
sc = SoupChannel(tod = fulldata, toc = countdata)
sc = setClusters(sc, setNames(cluster_info$Cluster, rownames(cluster_info)))
sc = autoEstCont(sc = sc)
# remove ambient RNA from cell data
decont_celldata = adjustCounts(sc = sc)
# create Seurat object from count data
celldata = CreateSeuratObject(counts = decont_celldata, project = project_name)
# display total number of cells after ambient RNA removal
table(celldata$orig.ident)
# get fractional mitochondria expression levels
celldata[["pct.mito"]] = PercentageFeatureSet(celldata, pattern="^MT-")
# normalize counts
ex = NormalizeData(celldata, normalization.method = "LogNormalize", scale.factor = normalization_scale_factor)
# feature selection
ex = FindVariableFeatures(ex, selection.method="vst", nfeatures = feature_selection_threshold)
# fetch gene names
all.genes <- rownames(ex)
# scale expression levels
ex = ScaleData(ex, features=all.genes, verbose = TRUE)
# dimensionality reduction using PCA
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
ex = FindClusters(object = ex, resolution=c(cluster_resolution))
# pK identification -- computing bimodality coefficient distribution
sweep.list = paramSweep_v3(ex, PCs = 1:min.pc, num.cores=num_cores)
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
celldata = subset(ex, doublet_finder == "Singlet")
# display total number of cells after doublet removal
table(celldata$orig.ident)
# retrieve table of QC metrics
celldata_qc = as_tibble(celldata[[]], rownames = "Cell.Barcode")
head(celldata_qc)
# relationship between QC metrics
qc_plot <- celldata_qc %>%
  ggplot(aes(y=nFeature_RNA, x=nCount_RNA, colour=pct.mito)) + 
  geom_point() +
  geom_hline(aes(yintercept=feature_qc_threshold), linetype='dashed', color='red') +
  annotate("text", x=85000, y=feature_qc_threshold+100, 
           label=feature_qc_threshold, size=4, color="red") +
  geom_vline(aes(xintercept=rna_qc_threshold), linetype='dashed', color='red') +
  annotate("text", x=rna_qc_threshold+5500, y=8000, label=rna_qc_threshold, size=4, color="red") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 9)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 12)) +
  scale_color_gradientn(colours=c('darkblue', 'darkgreen', 
                                  'yellow', 'orange', 'red'), 
                        limits=c(0,pct_mitochondrial_expression_qc_threshold), 
                        oob=scales::squish) +
  ggtitle("Relationship between depth (mRNA) count, gene count, and fractional mitochondrial expression levels") +
  theme(plot.title = element_text(size = 9, face = "bold"))
# filter out low-quality cells based on mitochondrial expression + RNA (depth) count + gene (feature) count
celldata <- subset(celldata, pct.mito <= pct_mito | nCount_RNA >= molecule_threshold | nFeature_RNA >= feature_threshold)
table(celldata$orig.ident)
saveRDS(celldata, file = output_seurat_object_filepath)
ggsave(filename = output_qc_plot_filepath, 
       qc_plot, dpi=350)