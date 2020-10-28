###########SEURAT analysis
library(Seurat)
###Reading single cell data matrix in format H5 
library(hdf5r)

C51_HD <- Read10X_h5('GSM4475048_C51_filtered_feature_bc_matrix.h5')
C52_HD <- Read10X_h5('GSM4475049_C52_filtered_feature_bc_matrix.h5')
C100_HD <- Read10X_h5('GSM4475050_C100_filtered_feature_bc_matrix.h5')
C148_HD <- Read10X_h5('GSM4475051_C148_filtered_feature_bc_matrix.h5')
C149_HD <- Read10X_h5('GSM4475052_C149_filtered_feature_bc_matrix.h5')
C152_HD <- Read10X_h5('GSM4475053_C152_filtered_feature_bc_matrix.h5')
C141_mildCOVID <- Read10X_h5('GSM4339769_C141_filtered_feature_bc_matrix.h5')
C142_mildCOVID <- Read10X_h5('GSM4339770_C142_filtered_feature_bc_matrix.h5')
C144_mildCOVID <- Read10X_h5('GSM4339772_C144_filtered_feature_bc_matrix.h5')
C143_severeCOVID <- Read10X_h5('GSM4339771_C143_filtered_feature_bc_matrix.h5')
C145_severeCOVID <- Read10X_h5('GSM4339773_C145_filtered_feature_bc_matrix.h5')
C146_severeCOVID <- Read10X_h5('GSM4339774_C146_filtered_feature_bc_matrix.h5')

#creating individual Seurat objects
control_51 <- CreateSeuratObject(counts = C51_HD, project = "control", min.cells = 3, min.features = 200)
control_52 <- CreateSeuratObject(counts = C52_HD, project = "control", min.cells = 3, min.features = 200)
control_100 <- CreateSeuratObject(counts = C100_HD, project = "control", min.cells = 3, min.features = 200)
control_148 <- CreateSeuratObject(counts = C148_HD, project = "control", min.cells = 3, min.features = 200)
control_149 <- CreateSeuratObject(counts = C149_HD, project = "control", min.cells = 3, min.features = 200)
control_152 <- CreateSeuratObject(counts = C152_HD, project = "control", min.cells = 3, min.features = 200)
mildcovid_141 <- CreateSeuratObject(counts = C141_mildCOVID, project = "mc", min.cells = 3, min.features = 200)
mildcovid_142 <- CreateSeuratObject(counts = C142_mildCOVID, project = "mc", min.cells = 3, min.features = 200)
mildcovid_144 <- CreateSeuratObject(counts = C144_mildCOVID, project = "mc", min.cells = 3, min.features = 200)
severecovid_143 <- CreateSeuratObject(counts = C143_severeCOVID, project = "sc", min.cells = 3, min.features = 200)
severecovid_145 <- CreateSeuratObject(counts = C145_severeCOVID, project = "sc", min.cells = 3, min.features = 200)
severecovid_146 <- CreateSeuratObject(counts = C146_severeCOVID, project = "sc", min.cells = 3, min.features = 200)

#merging the Seurat objects
all <- merge(control_51, y = c(control_52, control_100, control_148, control_149, control_152 , mildcovid_141, mildcovid_142, mildcovid_144, severecovid_143, severecovid_145, severecovid_146 ), add.cell.ids = c("ct51", "ct52","ct100","ct148","ct149","ct152","mc141", "mc142", "mc144", "sc143","sc145","sc146"), project = "covid")

########### description of complete object with pooled experiments
all
An object of class Seurat 
23742 features across 90696 samples within 1 assay 
Active assay: RNA (23742 features, 0 variable features)
###########
#description of the object by their cell origin
list <- SplitObject(all, split.by = "orig.ident")
############## description of combined object
list
$control
An object of class Seurat 
23742 features across 39900 samples within 1 assay 
Active assay: RNA (23742 features, 0 variable features)
$mc
An object of class Seurat 
23742 features across 9710 samples within 1 assay 
Active assay: RNA (23742 features, 0 variable features)
$sc
An object of class Seurat 
23742 features across 41086 samples within 1 assay 
Active assay: RNA (23742 features, 0 variable features)
##############
#normalization
list <- lapply(X = list, FUN = function(x) {
x <- NormalizeData(x)
x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
#find common anchors 
anchors <- FindIntegrationAnchors(object.list = list, dims = 1:20)
combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(combined) <- "integrated"
#scale data
combined <- ScaleData(combined, verbose = FALSE)
#dimensionnal reductions
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:20)
#quick graphical representation
DimPlot(combined, reduction = "umap")   
DimPlot(combined, reduction = "umap",group.by ="orig.ident")
#saving graph under linux
pdf(file = "dimplotumap.pdf",  width = 4, height = 4)
DimPlot(combined, reduction = "umap",cols=c("green","blue","red")) 
dev.off()
#save the Seurat objet in *.rda format
save(combined,file="all.rda")
# anchor description
anchors
> anchors
An AnchorSet object containing 38738 anchors between 3 Seurat objects 
# script for graph of combined group UMAP by cell origin
jpeg(file = "dimplotumap.jpeg",  width = 440, height = 440)
DimPlot(combined, reduction = "umap",cols=c("green","blue","red")) 
dev.off()
# script for Dotplot of different markers with their cell origin 
jpeg(file = "DOTPLOTMARKERS.jpeg",  width = 700, height = 300)
DotPlot(combined, features = rev(c("KRT8","CD3E","CD8A","IL7R","CSF1R","CSF2RA","CSF3R","CD14","FCGR3A","NKG7","MS4A1","HBB","PF4","CD34","PTPRC")),cols=c("green","blue","red"), dot.scale = 8, split.by = "orig.ident",assay = "RNA") +
RotatedAxis()
dev.off()
# script for FeaturePlot on UMAP reduction with split by cell origin
jpeg(file = "fcgr3a.jpg",  width = 900, height = 300)
FeaturePlot(combined, features = c("FCGR3A"),min.cutoff = "q9",cols=c("#CCFFFF","darkblue"),split.by= "orig.ident")
dev.off()
# script drawing violinplot on count slot and split by cell origin
jpeg(file = "PPARGVLN.jpg", width = 300, height = 300)
VlnPlot(combined, features = c("PPARG"), slot = "counts", log = TRUE,split.by= "orig.ident",pt.size=0,cols=c("green","blue","red"))
dev.off()
# script drawing scatter biplot with split color on cell origin
jpeg(file = "biplotPPARG&NR3C1.jpg",  width = 250, height = 250)
FeatureScatter(combined, feature1 = "rna_PPARG", feature2 = "rna_NR3C1", slot="counts",
cols=c("green","blue","red"))
dev.off()
# script for UMAP with different color according origin of cells
jpeg(file = "splitumap.jpg",  width = 900, height = 300)
DimPlot(combined, reduction = "umap",group.by ="orig.ident",split.by = "orig.ident",cols=c("green","blue","red"))
dev.off()
#### CD14+/CD16+ subseting in Seurat
mohd<-WhichCells(combined, idents = "control", expression = CD14 > 4 & FCGR3A > 4 & PPARG  >= 2, slot ="counts")
length(mohd)
539
mosc<-WhichCells(combined, idents = "sc", expression = CD14 > 4 & FCGR3A > 4 & PPARG < 2, slot ="counts")
length(mosc)
1134
mixmono<-c(mohd,mosc)
length(mixmono)
[1] 1693
monocytes <- SubsetData(object = combined, cells = mixmono)
###### single cell trajectory figures 5 and 6
### monocle 
#creating monocle object
mat<-as.matrix(expression)
colonnes<-as.data.frame(colnames(mat))
colonnes$group<-meta$group
colnames(colonnes)<-c("id","group")
row.names(colonnes)<-colonnes$id
pd <- new("AnnotatedDataFrame", data = colonnes)
genes<-as.data.frame(row.names(evt))
colnames(genes)<-"gene_short_name"
row.names(genes)<-genes$gene_short_name
fd <- new("AnnotatedDataFrame", data = genes)
cds <- newCellDataSet(mat, phenoData = pd,featureData = fd)
#defined cell hierarchy
 cth <- newCellTypeHierarchy()
PPARG_id <- row.names(subset(fData(cds), gene_short_name == "PPARG"))
> cth <- addCellType(cth, "PPARG_HIGH", classify_func = function(x) { x[PPARG_id,] >= 8})
> cth <- addCellType(cth, "PPARG_MEDIUM", classify_func = function(x) { x[PPARG_id,] < 8 & x[PPARG_id,] > 2} )
> cth <- addCellType(cth, "PPARG_LOW", classify_func = function(x) { x[PPARG_id,] <= 2 })
> cds <- classifyCells(cds, cth)
> table(pData(cds)$CellType)
meta<-as.data.frame(colnames(table))
meta$group<-Idents(monocytes)
colnames(meta)<-c("identifier","group")
row.names(meta)<-meta$identifier
#data transformation
cds@expressionFamily = negbinomial.size()
cds@lowerDetectionLimit = 0
my_feat <- fData(cds)
my_feat$id<-my_feat$gene_short_name
head(my_feat)
#dispersion table and variance 
disp_table <- dispersionTable(cds)
table(disp_table$mean_expression>=0.1)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
plot_ordering_genes(cds)
plot_pc_variance_explained(cds, return_all = FALSE)
cds <- reduceDimension(cds, max_components = 2, num_dim = 30, reduction_method = 'tSNE', verbose = TRUE)
#clustering cells
cds <- clusterCells(cds)
table(pData(cds)$CellType)
#ggplot graph cell type
pie <- ggplot(pData(cds),
aes(x = factor(1), fill = factor(CellType))) + geom_bar(width = 1)
pie + coord_polar(theta = "y") +
theme(axis.title.x = element_blank(), axis.title.y = element_blank())
#cluster cells with hiearchy
cds <- clusterCells(cds,cth) 
plot_cell_clusters(cds, 1, 2, color = "Cluster") +  facet_wrap(~CellType)
plot_cell_clusters(cds, 1, 2, color = "CellType")
save(cds,file="monohiearchy.rda")
#differential expressed genes against cell type
diff_test_res <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = "~CellType")
expressed_genes <-  row.names(subset(fData(cds),num_cells_expressed >= 10))
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
my_ordering_genes <- row.names(diff_test_res)[order(diff_test_res$qval)][1:1000]
cds2 <- setOrderingFilter(cds, ordering_genes = my_ordering_genes)
#dimensional reduction and gene selection
cds2 <- reduceDimension(cds, method = 'DDRTree')
gene_to_cluster <- row.names(diff_test_res)[order(diff_test_res$qval)][1:50] 
conca<-c("PPARG",gene_to_cluster)
cds2 <- orderCells(cds2)
#build graphs on trajectory 
plot_cell_clusters(cds2, color_by = 'as.factor(CellType)',markers="FABP4")
 plot_cell_trajectory(cds2, color_by = "CellType",)
plot_cell_trajectory(cds2, color_by = "CellType",markers="HLA-DQA2",markers_linear = TRUE,show_branch_points=FALSE)
plot_cell_trajectory(cds2, color_by = "Pseudotime")
plot_cell_trajectory(cds2, color_by = "CellType") +  facet_wrap(~group)
my_pseudotime_cluster <- plot_pseudotime_heatmap(cds2[conca,],cores = 8,
show_rownames = TRUE,return_heatmap = TRUE,cluster_rows = TRUE)
plot_genes_in_pseudotime(cds2[c("CD9","LY6D","LAMA3","ANTXR1","DPT","CTGF","PLPP3","JAM2","SERPINE1"),],
cell_size = 2, color_by = "group",ncol = 2)
#save data
write.table(pData(cds),file="phenotype.txt")
write.table(diff_test_res,file="difexpCD9trajectory.txt")
save(cds,file="monoclecds.rda")
save(cds2,file="monoclecdsfiltre.rda")
### end of the code
