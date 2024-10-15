#R4.0
###读入fcs文件
library(cytofWorkflow)
files=list.files("./CyTOF","*fcs")
samp <- read.flowSet(files = files,path = "./CyTOF/")
read.csv("./descrip.csv")->des        
colnames(samp)[1:42]=as.character(des$desc[1:42])         
# construct SingleCellExperiment
library(CATALYST)
panel=read_xls("./CyTOF/panel.xls")
read_xls("./CyTOF/md.xls")->md
sce <- prepData(samp, panel, md, features = panel$fcs_colname)

#####质控
pdf('./CyTOF/all_markers_density.pdf',width=15,height=9)
p <- plotExprs(sce, color_by = "patient_id")
p$facet$params$ncol <- 6
print(p)
dev.off()

##细胞数量
n_cells(sce) 
pdf("./CyTOF/cellnumber.pdf")
plotCounts(sce,color_by = "sample_id")
dev.off()
##MDS
CATALYST::plotMDS(sce, color_by = "patient_id")
###表达差异热图
plotExprHeatmap(sce)




####按20簇聚类
pro='basic_cluster_k40'
set.seed(1234)
sce <- cluster(sce, features = "type",
               xdim = 10, ydim = 10, maxK = 40, seed = 1234)

pdf(paste0(pro,'_cluster_plotExprHeatmap_row_clust_F.pdf'))
p<-plotExprHeatmap(sce, features = "type", fun="median",
                by = "cluster_id", k = "meta40", 
                row_clust = F,
                bars = TRUE, perc = TRUE)
dev.off()

#pdf(paste0(pro,'_plotClusterExprs.pdf'))
#plotClusterExprs(sce, k = "meta15", features = "type")
#dev.off()

pdf(paste0(pro,'_cluster_plotMultiHeatmap.pdf'),width=15,height=8)
plotMultiHeatmap(sce, 
                 hm1 = "type",   k = "meta40", 
                 row_anno = FALSE, bars = TRUE, perc = TRUE)
dev.off()

###
set.seed(1234)
sce <- runDR(sce, "TSNE", cells = 1e5, features = "type")
sce <- runDR(sce, "UMAP", cells = 1e5, features = "type")
plotDR(sce, "UMAP", color_by = "FOXP3")
pdf('./TSNE_K20.pdf',useDingbats=F)
plotDR(sce, "TSNE", color_by = "meta30") 
dev.off()
pdf('./UMAP_K20.pdf',useDingbats=F)
plotDR(sce, "UMAP", color_by = "meta40") 
dev.off()

plotAbundances(sce, k = "meta20", by = "sample_id")

plotDR(sce, "TSNE", color_by = "meta30", facet_by = "patient_id")

plotDR(sce, "UMAP", color_by = "meta40")

plotMedExprs(sce, k = "meta30", 
    facet_by = "cluster_id", shape_by = "patient_id")

pdf("./CyTOF/marker.pdf",useDingbats=F)
for(i in rownames(sce)){
    p<-plotDR(sce, "UMAP", color_by = i ,a_pal = rev(hcl.colors(10, "Spectral")))
    print(p)
}
dev.off()



  