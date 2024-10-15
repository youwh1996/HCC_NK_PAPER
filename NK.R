####R-3.6.3
####Seurat 3.1.5
#####monocle 2.14.0
###root账户下跑
cwl-runner  --outdir /home/cns3/rawdata/N1/  --parallel  /home/cns3/BD/rhapsody_wta_1.8.cwl /home/cns3/BD/template_wta_1.8.yml
####
library(Seurat)

####预处理

read.csv("/home/youwh/data/sc_HCC_NK/N1/N1_RSEC_MolsPerCell.csv",skip=6,header=T,row.names=1,check.names=F)->N1
paste("N1_",rownames(N1),sep="")->rownames(N1)
CreateSeuratObject(counts=t(N1),project="N1")->N1
N1[["percent.mt"]] <- PercentageFeatureSet(N1, pattern ="^MT-")
plot1 <- FeatureScatter(N1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(N1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
N1 <- subset(N1, subset = nFeature_RNA > 200 & nFeature_RNA < 2000  & percent.mt<50 )



read.csv("/home/youwh/data/sc_HCC_NK/N2/N2_RSEC_MolsPerCell.csv",skip=6,header=T,row.names=1,check.names=F)->N2
paste("N2_",rownames(N2),sep="")->rownames(N2)
CreateSeuratObject(counts=t(N2),project="N2")->N2
N2[["percent.mt"]] <- PercentageFeatureSet(N2, pattern ="^MT-")
plot1 <- FeatureScatter(N2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(N2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
N2 <- subset(N2, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mt<50  )


read.csv("/home/youwh/data/sc_HCC_NK/N3/N3_RSEC_MolsPerCell.csv",skip=6,header=T,row.names=1,check.names=F)->N3
paste("N3_",rownames(N3),sep="")->rownames(N3)
CreateSeuratObject(counts=t(N3),project="N3")->N3
N3[["percent.mt"]] <- PercentageFeatureSet(N3, pattern ="^MT-")
plot1 <- FeatureScatter(N3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(N3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
N3 <- subset(N3, subset = nFeature_RNA > 200 & nFeature_RNA < 2000  & percent.mt<50 )

read.csv("/home/youwh/data/sc_HCC_NK/N4/N4_RSEC_MolsPerCell.csv",skip=6,header=T,row.names=1,check.names=F)->N4
paste("N4_",rownames(N4),sep="")->rownames(N4)
CreateSeuratObject(counts=t(N4),project="N4")->N4
N4[["percent.mt"]] <- PercentageFeatureSet(N4, pattern ="^MT-")
plot1 <- FeatureScatter(N4, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(N4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
N4 <- subset(N4, subset = nFeature_RNA > 200 & nFeature_RNA < 2000  & percent.mt<50 )

read.csv("/home/youwh/data/sc_HCC_NK/N5/N5_RSEC_MolsPerCell.csv",skip=6,header=T,row.names=1,check.names=F)->N5
paste("N5_",rownames(N5),sep="")->rownames(N5)
CreateSeuratObject(counts=t(N5),project="N5")->N5
N5[["percent.mt"]] <- PercentageFeatureSet(N5, pattern ="^MT-")
plot1 <- FeatureScatter(N5, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(N5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
N5 <- subset(N5, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mt<50 )

read.csv("/home/youwh/data/sc_HCC_NK/T1/T1_RSEC_MolsPerCell.csv",skip=6,header=T,row.names=1,check.names=F)->T1
paste("T1_",rownames(T1),sep="")->rownames(T1)
CreateSeuratObject(counts=t(T1),project="T1")->T1
T1[["percent.mt"]] <- PercentageFeatureSet(T1, pattern ="^MT-")
plot1 <- FeatureScatter(T1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(T1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
T1 <- subset(T1, subset = nFeature_RNA > 200 & nFeature_RNA < 2000   & percent.mt<50)

read.csv("/home/youwh/data/sc_HCC_NK/T2/T2_RSEC_MolsPerCell.csv",skip=6,header=T,row.names=1,check.names=F)->T2
paste("T2_",rownames(T2),sep="")->rownames(T2)
CreateSeuratObject(counts=t(T2),project="T2")->T2
T2[["percent.mt"]] <- PercentageFeatureSet(T2, pattern ="^MT-")
plot1 <- FeatureScatter(T2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(T2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
T2 <- subset(T2, subset = nFeature_RNA > 200 & nFeature_RNA < 2000   & percent.mt<50  )


read.csv("//home/youwh/data/sc_HCC_NK/T3/T3_RSEC_MolsPerCell.csv",skip=6,header=T,row.names=1,check.names=F)->T3
paste("T3_",rownames(T3),sep="")->rownames(T3)
CreateSeuratObject(counts=t(T3),project="T3")->T3
T3[["percent.mt"]] <- PercentageFeatureSet(T3, pattern ="^MT-")
plot1 <- FeatureScatter(T3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(T3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
T3 <- subset(T3, subset = nFeature_RNA > 200 & nFeature_RNA < 2000   & percent.mt<50 )

read.csv("/home/youwh/data/sc_HCC_NK/T4/T4_T190117_RSEC_MolsPerCell.csv",skip=6,header=T,row.names=1,check.names=F)->T4
paste("T4_",rownames(T4),sep="")->rownames(T4)
CreateSeuratObject(counts=t(T4),project="T4")->T4
T4[["percent.mt"]] <- PercentageFeatureSet(T4, pattern ="^MT-")
plot1 <- FeatureScatter(T4, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(T4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
T4 <- subset(T4, subset = nFeature_RNA > 200 & nFeature_RNA < 2000  &  percent.mt<50 )

read.csv("/home/youwh/data/sc_HCC_NK/T5/T5_combined_RSEC_MolsPerCell.csv",skip=6,header=T,row.names=1,check.names=F)->T5
paste("T5_",rownames(T5),sep="")->rownames(T5)
CreateSeuratObject(counts=t(T5),project="T5")->T5
T5[["percent.mt"]] <- PercentageFeatureSet(T5, pattern ="^MT-")
plot1 <- FeatureScatter(T5, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(T5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
T5 <- subset(T5, subset = nFeature_RNA > 200 & nFeature_RNA < 2000   & percent.mt<50 )
#######

merge(N1,list(N2,N3,N4,N5,T1,T2,T3,T4,T5))->all
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 4000)
all <- ScaleData(all, features = rownames(all))
all <- RunPCA(all, features = VariableFeatures(object = all))
all <- FindNeighbors(all, dims = 1:50)
all <- FindClusters(all, resolution = 0.2)
all <- RunUMAP(all, dims = 1:50)
DimPlot(all,reduction="umap",label=T)
FeaturePlot(all,features=c("CD3D",'IL7R','NCAM1'),reduction='umap',order=T)
###CD3-CD56+/CD7+ IL7R-
X11()
FeaturePlot(all,features=c("CD7","NCAM1","CD3D","IL7R"),reduction="umap")
subset(all,idents=c("0","1","2","3","4","5","7","9","11","12","14","21"))->NK
######harmony去除批次
library(harmony)
NK@meta.data$patient=factor(NK@meta.data$orig.ident)
levels(NK@meta.data$patient)=c("P1","P2","P3","P4","P5","P1","P2","P3","P4","P5")
NK@meta.data$origin=factor(NK@meta.data$orig.ident)
levels(NK@meta.data$origin)=c(rep("N",5),rep("T",5))
NK <- NormalizeData(NK, normalization.method = "LogNormalize", scale.factor = 10000)
NK <- FindVariableFeatures(NK, selection.method = "vst", nfeatures = 3000)
NK <- ScaleData(NK, features = rownames(NK))
NK <- RunPCA(NK, features = VariableFeatures(object = NK))
RunHarmony(NK,"patient")->NK
NK <- FindNeighbors(NK, dims = 1:30,reduction = "harmony")
NK <- FindClusters(NK, resolution = 0.3)
NK <- RunUMAP(NK, dims = 1:30,reduction = "harmony")
NK <- RunTSNE(NK, dims = 1:30,reduction = "harmony")
subset(NK,idents=c(as.character(0:8)))->NK
saveRDS(NK,file="./rawdata_NK/NK.RDS")


pdf("./rawdata_NK/image/umap_all.pdf",useDingbats=F)
DimPlot(NK,reduction="umap",label=T,pt.size=0.5,shape.by=NULL)+scale_color_manual(values=c("#F0CB5C", "#CF5531", "#C8A374", "#F4A9CE", "#FC7970", "#E1EEB9", "#57B66B", "#D0349A", "#F4AC3E"))
dev.off()
pdf("./rawdata_NK/image/umap_origin.pdf",useDingbats=F)
DimPlot(NK,reduction="umap",label=T,pt.size=0.5,group.by="origin",shape.by=NULL)+scale_color_manual(values=c("#3EC8B5","#FE8B82"))
dev.off()
pdf('./rawdata_NK/image/umap_patient.pdf',useDingbats=F)
DimPlot(NK,reduction="umap",label=T,pt.size=0.5,group.by="patient",shape.by=NULL)+scale_color_manual(values=c("#25CDFE", "#46B6D1",'#FFAEC2',"#726FAF", "#E874C2"))
dev.off()
#cell cycle
CellCycleScoring(NK,g2m.features=cc.genes$g2m.genes,s.features=cc.genes$s.genes)->NK
pdf("./rawdata_NK/image/cellcycle.pdf",useDingbats=F)
DimPlot(NK,reduction="umap",label=T,group.by="Phase",pt.size=0.5)
dev.off()
###样本比例
data.frame(sample=NK@meta.data$orig.ident,cluster=NK@meta.data$seurat_clusters)->per
per$cluster<-as.character(per$cluster)
per$sample<-as.character(per$sample)
data.frame(table(per))->per
per$percent<-c()
for(i in 1:nrow(per)){
    per$percent[i]=per$Freq[i]/sum(per$Freq[which(per$sample==per$sample[i])])
}

pdf("./rawdata_NK/image/percentage_NK.pdf",useDingbats=F)
ggplot(per,aes(x=sample,y=percent,fill=cluster))+
geom_col(position="stack")+
scale_fill_d3()+
coord_flip()+
theme_bw()+
theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())
dev.off()
####癌组织和癌旁 桑基图
library(ggalluvial)
library(reshape)
data.frame(sample=NK@meta.data$origin,cluster=NK@meta.data$seurat_clusters)->per
per$cluster<-as.character(per$cluster)
per$sample<-as.character(per$sample)
data.frame(table(per))->per
per$percent<-c()
for(i in 1:nrow(per)){
    per$percent[i]=per$Freq[i]/sum(per$Freq[which(per$sample==per$sample[i])])
}

per[which(per$sample=="N5" | per$sample=="T5"),]->P5
per$order=(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8))
pdf("./rawdata_NK/image/sankey.pdf",useDingbats=F)
ggplot(per,aes(x = sample, stratum =cluster, 
        alluvium = order,y = percent,
        fill = cluster, label = cluster)) +
  scale_x_discrete(expand = c(.1, .1)) +
  scale_fill_manual(values=c("#F0CB5C", "#CF5531", "#C8A374", "#F4A9CE", "#FC7970", "#E1EEB9", "#57B66B", "#D0349A", "#F4AC3E"))+
  geom_flow() +
  geom_stratum(alpha = .5) +
  geom_text(stat = "stratum", size = 2.5) +
  theme(panel.border=element_rect(color='black', fill=NA), 
       panel.grid.major =element_blank(), 
       panel.grid.minor = element_blank(),
       panel.background = element_blank(), 
       axis.line = element_line(colour = "black"),
       axis.text=element_text(colour="black")) 
dev.off()
####不同亚群比
per$group<-c(rep("Normal",5),rep("Tumor",5)) 
per$patient<-c("P1","P2","P3","P4","P5","P1","P2","P3","P4","P5")
#C0
per[which(per$cluster=="0"),]->C0
ggplot(C0,aes(x=group,y=percent))+
geom_bar(position ="stack")+
#geom_point(aes(color=patient))+scale_color_manual(values=c("#FDB462", "#F0027F",'turquoise2',"#80B1D3", "#00AF99"))+
theme_bw()

########
FindAllMarkers(NK,only.pos=T,min.pct=0.25,logfc.threshold=0.25,test.use="MAST")->deg
library(dplyr)
data.frame(deg %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC))->top10
##heatmap
cts <- GetAssayData(NK, slot = "data")

apply(cts[,rownames(NK@meta.data)[which(NK@meta.data$seurat_clusters=="0")]],1,mean)->C0
apply(cts[,rownames(NK@meta.data)[which(NK@meta.data$seurat_clusters=="1")]],1,mean)->C1
apply(cts[,rownames(NK@meta.data)[which(NK@meta.data$seurat_clusters=="2")]],1,mean)->C2
apply(cts[,rownames(NK@meta.data)[which(NK@meta.data$seurat_clusters=="3")]],1,mean)->C3
apply(cts[,rownames(NK@meta.data)[which(NK@meta.data$seurat_clusters=="4")]],1,mean)->C4
apply(cts[,rownames(NK@meta.data)[which(NK@meta.data$seurat_clusters=="5")]],1,mean)->C5
apply(cts[,rownames(NK@meta.data)[which(NK@meta.data$seurat_clusters=="6")]],1,mean)->C6
apply(cts[,rownames(NK@meta.data)[which(NK@meta.data$seurat_clusters=="7")]],1,mean)->C7
apply(cts[,rownames(NK@meta.data)[which(NK@meta.data$seurat_clusters=="8")]],1,mean)->C8
cbind(C0,C1,C2,C3,C4,C5,C6,C7,C8)->AC
AC[top10$gene,]->AC
pdf("./rawdata_NK/image/top10heatmap.pdf",useDingbats=F)
pheatmap(AC,cluster_rows = F,cluster_cols = F,scale="row",fontsize_row=6,breaks=unique(c(seq(-1,1, length=100))),color=colorRampPalette(c("blue", "white", "#EE7600"))(100),border_color="NA")
dev.off()
 




######
pdf("./rawdata_NK/image/CD7_umap.pdf",useDingbats=F)
FeaturePlot(NK,features="CD7",reduction="umap",pt.size=0.5,order=T)+ scale_color_viridis_c()
dev.off()

pdf("./rawdata_NK/image/NCAM1_umap.pdf",useDingbats=F)
FeaturePlot(NK,features="NCAM1",reduction="umap",pt.size=0.5,order=T)+ scale_color_viridis_c()
dev.off()

pdf("./rawdata_NK/image/FGFBP2_umap.pdf",useDingbats=F)
FeaturePlot(NK,features="FGFBP2",reduction="umap",pt.size=0.5,order=T)+ scale_color_viridis_c()
dev.off()
pdf("./rawdata_NK/image/CXCR6_umap.pdf",useDingbats=F)
FeaturePlot(NK,features="CXCR6",reduction="umap",pt.size=0.5,order=T)+ scale_color_viridis_c()
dev.off()
pdf("./rawdata_NK/image/EOMES_umap.pdf",useDingbats=F)
FeaturePlot(NK,features="EOMES",reduction="umap",pt.size=0.5,order=T)+ scale_color_viridis_c()
dev.off()
pdf("./rawdata_NK/image/CX3CR1_umap.pdf",useDingbats=F)
FeaturePlot(NK,features="CX3CR1",reduction="umap",pt.size=0.5,order=T)+ scale_color_viridis_c()
dev.off()
pdf("./rawdata_NK/image/FCGR3A_umap.pdf",useDingbats=F)
FeaturePlot(NK,features="FCGR3A",reduction="umap",pt.size=0.5,order=T)+ scale_color_viridis_c()
dev.off()
pdf("./rawdata_NK/image/MKI67_umap.pdf",useDingbats=F)
FeaturePlot(NK,features="MKI67",reduction="umap",pt.size=0.5,order=T)+ scale_color_viridis_c()
dev.off()
pdf("./rawdata_NK/image/IL7R_umap.pdf",useDingbats=F)
FeaturePlot(NK,features="IL7R",reduction="umap",pt.size=0.5,order=T)+ scale_color_viridis_c()
dev.off()
pdf("./rawdata_NK/image/ITGA1_umap.pdf",useDingbats=F)
FeaturePlot(NK,features="ITGA1",reduction="umap",pt.size=0.5,order=T)+ scale_color_viridis_c()
dev.off()
pdf("./rawdata_NK/image/KLRC1_umap.pdf",useDingbats=F)
FeaturePlot(NK,features="KLRC1",reduction="umap",pt.size=0.5,order=T)+ scale_color_viridis_c()
dev.off()
pdf("./rawdata_NK/image/CD69_umap.pdf",useDingbats=F)
FeaturePlot(NK,features="CD69",reduction="umap",pt.size=0.5,order=T)+ scale_color_viridis_c()
dev.off()
pdf("./rawdata_NK/image/TBX21_umap.pdf",useDingbats=F)
FeaturePlot(NK,features="TBX21",reduction="umap",pt.size=0.5,order=T)+ scale_color_viridis_c()
dev.off()
######免疫检查点：CTLA4、PDCD1、TIGIT、HAVCR2、LAG3
pdf("./rawdata_NK/image/CTLA4_umap.pdf",useDingbats=F)
FeaturePlot(NK,features="CTLA4",reduction="umap",pt.size=0.5,order=T)+ scale_color_viridis_c()
dev.off()
pdf("./rawdata_NK/image/PDCD1_umap.pdf",useDingbats=F)
FeaturePlot(NK,features="PDCD1",reduction="umap",pt.size=0.5,order=T)+ scale_color_viridis_c()
dev.off()
pdf("./rawdata_NK/image/TIGIT_umap.pdf",useDingbats=F)
FeaturePlot(NK,features="TIGIT",reduction="umap",pt.size=0.5,order=T)+ scale_color_viridis_c()
dev.off()
pdf("./rawdata_NK/image/LAG3_umap.pdf",useDingbats=F)
FeaturePlot(NK,features="LAG3",reduction="umap",pt.size=0.5,order=T)+ scale_color_viridis_c()
dev.off()
pdf("./rawdata_NK/image/HAVCR2_umap.pdf",useDingbats=F)
FeaturePlot(NK,features="HAVCR2",reduction="umap",pt.size=0.5,order=T)+ scale_color_viridis_c()
dev.off()
##Cytotoxicity gene expression signature (GZMA, GZMB, GZMH, GZMK, GZMM, PRF1, GNLY, and NKG7) 
list("Cytotoxicity"=c("GZMA","GZMB","GZMH","GZMK","PRF1","GNLY","NKG7"))->Cytotoxicity
AddModuleScore(NK,features=Cytotoxicity,name='Cytotoxicity')->NK
pdf("./rawdata_NK/image/Cytotoxicity.pdf",useDingbats=F)
FeaturePlot(NK,features="Cytotoxicity1",reduction="umap",pt.size=0.5,order=T)+ scale_color_viridis_c()
dev.off()

list("Exausted"=c("PDCD1",'CTLA4','HAVCR2','LAG3','TIGIT','KLRG1'))->Exau
AddModuleScore(NK,features=Exau,name='Exau')->NK
FeaturePlot(NK,features="Exau1",reduction="umap",pt.size=0.5,order=T)+ scale_color_viridis_c()


####HCC3 signature
c("IL1R1","IL23R","LTBR","IL2RB","RORC","TOX","TOX2","IKZF2","ID2","AHR","TNFSF13B","TNFSF13","CCL20","IL22","TNFSF11","NCR2","NCR1","KIT",
"TNFSF4","CD2","HLA-DRB1","JAG2","ADAM10","DLL1","RBPJ","TCF7","TYROBP","FGR","PLCG2","SYK","VAV3","LAT2","LYN","LCK","ZAP70","PLCG1","VAV1","VAV2",
"FCER1G","KLRF2","PECAM1","AMICA1","PRAM1","EREG","SKAP2","NRP1","SH2D1B","XCL1","XCL2","KLRC1","KLRD1","LIF","FES","CLNK","SNCA","PTPN6","SIGLEC7","PTGDR",
"TLE1","TLE3","TCF4","NLRP2","IL26","RUNX3","TNFRSF18","NLRP7","TNFRSF25","EPCAM")->HCC3
list("HCC3"=HCC3)->HCC3
AddModuleScore(NK,features=HCC3,name="HCC3")->NK
FeaturePlot(NK,features="HCC31",reduction="umap",pt.size=0.5,order=T)+ scale_color_viridis_c()





####
FeaturePlot(NK,features=c("GZMB","PRF1","SELL","IL7R","PCNA","LST1","FOS","JUN","AREG","ISG15"),reduction="umap",order=T)
FeaturePlot(Normal,features=c("GZMB","PRF1","SELL","IL7R","PCNA","LST1","FOS","JUN","AREG","ISG15"),reduction="umap",order=T)
FeaturePlot(Tumor,features=c("GZMB","PRF1","SELL","IL7R","PCNA","LST1","FOS","JUN","AREG","ISG15"),reduction="umap",order=T)

#####NK activating receptors
FeaturePlot(NK,features=c("KLRK1","HCST","KLRF1","NCR1","NCR3","CD226","FASL","KLRC2"),reduction="umap",order=T)
FeaturePlot(Tumor,features=c("KLRK1","HCST","KLRF1","NCR1","NCR3","CD226","FASL","KLRC2"),reduction="umap",order=T)
FeaturePlot(Normal,features=c("KLRK1","HCST","KLRF1","NCR1","NCR3","CD226","FASL","KLRC2"),reduction="umap",order=T)

####NK inhibit receptors
FeaturePlot(NK,features=c("CRTAM","KLRC1","KLRD1","KLRB1","CD160","TIGIT","KIR2DL3","KIR2DL4"),reduction="umap",order=T)
FeaturePlot(Tumor,features=c("CRTAM","KLRC1","KLRD1","KLRB1","CD160","TIGIT","KIR2DL3","KIR2DL4"),reduction="umap",order=T)
FeaturePlot(Normal,features=c("CRTAM","KLRC1","KLRD1","KLRB1","CD160","TIGIT","KIR2DL3","KIR2DL4"),reduction="umap",order=T)


#####chemokine genes
FeaturePlot(NK,features=c("XCL1","XCL2","CCL3","CCL4","CC4L2","CCL5"),reduction="umap",order=T)


#########monocle 
FindAllMarkers(cNK,only.pos=T,logfc.threshold=0.25,min.pct=0.25)->deg

library(monocle)
as(as.matrix(lrNK@assays$RNA@counts),'sparseMatrix')->dat1
pd <- new('AnnotatedDataFrame', data = lrNK@meta.data)
fData <- data.frame(gene_short_name = row.names(dat1), row.names = row.names(dat1))
fd <- new('AnnotatedDataFrame', data = fData)
#fData(dat)$use_for_ordering=fData(dat)$num_cells_expressed> (0.05 * ncol(dat))

  #Construct monocle cds
  dat<- newCellDataSet(dat1,
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())

  rm(dat1)
  dat<-estimateSizeFactors(dat)
  dat<-estimateDispersions(dat)
  dat<-detectGenes(dat,min_expr=1)
 expressed_genes = row.names(subset(fData(dat),num_cells_expressed>= 30))
clustering_DEG_genes<-differentialGeneTest(dat[expressed_genes,],fullModelFormulaStr="~seurat_clusters",cores=10)
ordering_genes<-row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:500]


dat<-reduceDimension(dat,max_compnents=2,norm_method="log",num_dim=30,reduction_method="tSNE",verbose=T)
dat<-clusterCells(dat,verbose=T)
#data.frame(deg %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC))->top50
#ordering_genes<-unique(deg$gene)
dat<-setOrderingFilter(dat,ordering_genes=ordering_genes)
dat<-reduceDimension(dat,method="DDRTree")
dat<-orderCells(dat)
saveRDS(dat,file="./monocle_dat_1000.RDS")
pdf("./youwh/trajec1.pdf")
plot_cell_trajectory(dat,color_by="seurat_clusters",cell_link_size=1,cell_size=1)+scale_color_manual(values=c("#E41A1C", "#006837", "#FF7F00"))
dev.off()

pdf("./youwh/trajectime.pdf")
plot_cell_trajectory(dat,color_by="Pseudotime",cell_link_size=1,cell_size=1)
dev.off()

######difussion map
library("destiny")
library(ggplot2)
as.matrix(lrNK@assays$RNA@data)->exp
anno=lrNK@meta.data
pDat=data.frame(cell=colnames(lrNK),cell_type=lrNK@meta.data$seurat_clusters,row.names=colnames(lrNK))
pDat$cell_type<-as.character(pDat$cell_type)

eset <- Biobase::ExpressionSet(exp, phenoData=Biobase::AnnotatedDataFrame(pDat))
dmap <- DiffusionMap(eset,verbose=T,n_pcs=50,k=3)
plot.DiffusionMap(dmap,col_by="cell_type")
plot.DiffusionMap(dmap,dims=c(1,2),col_by="cell_type")+scale_fill_manual(values=c("#279E68", "#D62728", "#AA40FC"))


saveRDS(dmap,file="./ICC/dmap_ICC.RDS")
set.seed(4)
dpt <- DPT(dmap, tips= sample(ncol(eset), 2L))
plot.DPT(dpt)
qplot(DC1,DC2, data=dmap, colour=cell_type)+scale_color_manual(values=c("#F0CB5C", "#FC7970", "#E1EEB9", "#57B66B", "#D0349A", "#F4AC3E"))+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
plot(dpt,col_by="GATA3")

######






#####RNA 速率
samtools view -HS ./rawdata_NK/T5/T5_combined_final.BAM >head.txt
samtools view ./rawdata_NK/T5/T5_combined_final.BAM| grep "MA:Z:*" | sed "s/MA:Z:/UB:Z:/" > ./rawdata_NK/T5/temp.sam
cat head.txt ./rawdata_NK/T5/temp.sam | samtools view -Sb > ./rawdata_NK/T5/T5_changed.bam
rm ./rawdata_NK/T5/temp.sam
velocyto run -m ./ref/hg38_rmsk.gtf ./rawdata_NK/T5/T5_changed.bam ./ref/cellranger/refdata-gex-GRCh38-2020-A/genes/genes.gtf

import loompy
files=['/home/cns3/rawdata_NK/N1/velocyto/N1.loom','/home/cns3/rawdata_NK/N2/velocyto/N2.loom','/home/cns3/rawdata_NK/N3/velocyto/N3.loom',
'/home/cns3/rawdata_NK/N4/velocyto/N4.loom','/home/cns3/rawdata_NK/N5/velocyto/N5.loom','/home/cns3/rawdata_NK/T1/velocyto/T1.loom',
'/home/cns3/rawdata_NK/T2/velocyto/T2.loom','/home/cns3/rawdata_NK/T3/velocyto/T3.loom','/home/cns3/rawdata_NK/T4/velocyto/T4.loom','/home/cns3/rawdata_NK/T5/velocyto/T5.loom']
loompy.combine(files,"merged.loom",key="Accession")

###
colnames(NK)->sample  
gsub("_",":",sample)->sample
sample[which(substr(sample,1,2)=="N1")]=paste("N1_changed_OMI3C",substr(sample[which(substr(sample,1,2)=="N1")],3,nchar(sample[which(substr(sample,1,2)=="N1")])),sep="")
sample[which(substr(sample,1,2)=="N2")]=paste("N2_changed_Y1FMY",substr(sample[which(substr(sample,1,2)=="N2")],3,nchar(sample[which(substr(sample,1,2)=="N2")])),sep="")
sample[which(substr(sample,1,2)=="N3")]=paste("N3_changed_2TCS0",substr(sample[which(substr(sample,1,2)=="N3")],3,nchar(sample[which(substr(sample,1,2)=="N3")])),sep="")
sample[which(substr(sample,1,2)=="N4")]=paste("N4_changed_3YYE3",substr(sample[which(substr(sample,1,2)=="N4")],3,nchar(sample[which(substr(sample,1,2)=="N4")])),sep="")
sample[which(substr(sample,1,2)=="N5")]=paste("N5_changed_4Q609",substr(sample[which(substr(sample,1,2)=="N5")],3,nchar(sample[which(substr(sample,1,2)=="N5")])),sep="")
sample[which(substr(sample,1,2)=="T1")]=paste("T1_changed_H8D0T",substr(sample[which(substr(sample,1,2)=="T1")],3,nchar(sample[which(substr(sample,1,2)=="T1")])),sep="")
sample[which(substr(sample,1,2)=="T2")]=paste("T2_changed_LF1I6",substr(sample[which(substr(sample,1,2)=="T2")],3,nchar(sample[which(substr(sample,1,2)=="T2")])),sep="")
sample[which(substr(sample,1,2)=="T3")]=paste("T3_changed_IYZ0G",substr(sample[which(substr(sample,1,2)=="T3")],3,nchar(sample[which(substr(sample,1,2)=="T3")])),sep="")
sample[which(substr(sample,1,2)=="T4")]=paste("T4_changed_5N0I8",substr(sample[which(substr(sample,1,2)=="T4")],3,nchar(sample[which(substr(sample,1,2)=="T4")])),sep="")
sample[which(substr(sample,1,2)=="T5")]=paste("T5_changed_7S5PW",substr(sample[which(substr(sample,1,2)=="T5")],3,nchar(sample[which(substr(sample,1,2)=="T5")])),sep="")
write.csv(sample,file="./rawdata_NK/velo/cellID_obs.csv")


 cNK@reductions$umap@cell.embeddings->umap
 rownames(umap)=sample
write.csv(umap,file="./rawdata_NK/velo/cNKcell_umap.csv")
cell_cluster=NK@meta.data$seurat_clusters
data.frame(cell_cluster)->cell_cluster
rownames(cell_cluster)=sample
write.csv(cell_cluster,file='./rawdata_NK/velo/cellcluster.csv')

##
import anndata
import scvelo as scv
import pandas as pd
import velocyto as vcy
import numpy as np
import matplotlib as plt
import rpy2
sample_obs=pd.read_csv("./rawdata_NK/velo/cellID_obs.csv")
merged= anndata.read_loom('./rawdata_NK/merged.loom')

umap=pd.read_csv("./rawdata_NK/velo/cNKcell_umap.csv")
cellcluster=pd.read_csv("./rawdata_NK/velo/cNKcellcluster.csv")
umap=umap.rename(columns={'Unnamed: 0':'Cell ID'})
merged=merged[np.isin(merged.obs.index,umap['Cell ID'])]

umap=umap.iloc[:,1:] 
merged.obsm['X_umap']=umap.values
cell_cluster=pd.read_csv("./rawdata_NK/velo/cellcluster.csv")
cellcluster=cellcluster.iloc[:,2]
merged.obs['celltype']=cellcluster.values
###运行RNA速率
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.settings.set_figure_params('scVelo',format='pdf')  # for beautified visualization
scv.pp.filter_and_normalize(merged, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(merged, n_pcs=30, n_neighbors=40)
scv.tl.velocity(merged, mode = "stochastic")
scv.tl.velocity_graph(merged)
scv.pl.velocity_embedding(merged, basis='umap')

scv.pl.velocity_embedding(merged, arrow_length=3, arrow_size=5, dpi=120,basis='umap')

scv.pl.velocity_embedding_stream(merged, basis='umap', title='', smooth=.8, min_mass=0.5,color=['celltype'])


ident_colours=["#377EB8", "#984EA3", "#FFD92F", "#F781BF", "#A65628"]
scv.pl.velocity_embedding_grid(merged, basis='umap', title='',arrow_length=1.5,color=['celltype'],palette = ident_colours,arrow_size=9)

ident_colours=["#377EB8", "#E41A1C", "#006837", "#FF7F00", "#984EA3", "#FFD92F", "#46A040", "#F781BF", "#A65628"]
ident_colours=["#E41A1C", "#006837", "#FF7F00"]
scv.pl.velocity_embedding_grid(merged, basis='umap', title='',arrow_length=3,arrow_size=5,color=['celltype'],palette=ident_colours)

###lrNK
subset(NK,idents=c("1","2","3"))->lrNK
FindAllMarkers(lrNK,logfc.threshold=0.25,min.pct=0.25,only.pos=T)->deg
deg[which(deg$p_val_adj<0.05 & deg$avg_logFC>0.5),]->deg
DoHeatmap(lrNK,group.by="ident",features=unique(deg$gene))

c("FGFBP2","CX3CR1","FCGR3A","HLA-DRA","CD74","CD160","MT-ND1","HLA-DQA2","IFNG","IGHG3","MKI67","CCL4L2","HSPH1")->marker
pdf("./rawdata_NK/feature.pdf",useDingbats=F)
for(i in marker){
    p<-FeaturePlot(NK,features=i,reduction="umap",pt.size=0.5,order=T)+ scale_color_viridis_c()
    print(p)
}
dev.off()
###

NK@meta.data$type=factor(as.character(NK@meta.data$seurat_clusters))
levels(NK@meta.data$type)=c("cNK","lrNK","lrNK","lrNK","cNK","cNK","pNK","cNK","cNK")
NK@active.ident=NK@meta.data$type
names(NK@active.ident)<-colnames(NK)
FindMarkers(NK,ident.1="lrNK",ident.2="cNK",min.pct=0.25,logfc.threshold=0.25)->deg
##火山图
DoHeatmap(NK,cells=colnames(subset(NK,idents=c("cNK","lrNK"))),group.by="ident",features=deg$gene)



##CytoTRACE
library(CytoTRACE)
CytoTRACE(as.matrix(NK@assays$RNA@counts),ncores = 8)->tracescore
tracescore$CytoTRACE->CytoTRACE
data.frame(CytoTRACE=as.numeric(CytoTRACE),cluster=NK@meta.data$seurat_clusters)->dat
ggplot(dat,aes(x=cluster,y=CytoTRACE,fill=cluster))+geom_boxplot()

##lrNK
subset(NK,idents=c("1","2","3"))->lrNK
lrNK <- NormalizeData(lrNK, normalization.method = "LogNormalize", scale.factor = 10000)
lrNK <- FindVariableFeatures(lrNK, selection.method = "vst", nfeatures = 3000)
lrNK <- ScaleData(lrNK, features = rownames(lrNK))
lrNK <- RunPCA(lrNK, features = VariableFeatures(object = lrNK))
RunHarmony(lrNK,"patient")->lrNK
lrNK <- FindNeighbors(lrNK, dims = 1:30,reduction = "harmony")
lrNK <- FindClusters(lrNK, resolution = 0.3)
lrNK <- RunUMAP(lrNK, dims = 1:30,reduction = "harmony")
DimPlot(lrNK,reduction="umap",label=T)

CytoTRACE(as.matrix(lrNK@assays$RNA@counts),ncores=9)->lrscore
data.frame(CytoTRACE=lrscore$CytoTRACE,cluster= lrNK$seurat_clusters)->lr_CytoTRACE
ggplot(lr_CytoTRACE,aes(x=cluster,y=CytoTRACE))+
geom_boxplot(aes(fill=cluster))+
scale_fill_manual(values=c("#E41A1C", "#006837", "#FF7F00"))+
theme_bw()
###cNK
subset(NK,idents=c("0","4","5","7","8"))->cNK
cNK <- NormalizeData(cNK, normalization.method = "LogNormalize", scale.factor = 10000)
cNK <- FindVariableFeatures(cNK, selection.method = "vst", nfeatures = 3000)
cNK <- ScaleData(cNK, features = rownames(cNK))
cNK <- RunPCA(cNK, features = VariableFeatures(object = cNK))
RunHarmony(cNK,"patient")->cNK
cNK <- FindNeighbors(cNK, dims = 1:30,reduction = "harmony")
cNK <- FindClusters(cNK, resolution = 0.2)
cNK <- RunUMAP(cNK, dims = 1:30,reduction = "harmony")
DimPlot(cNK,reduction="umap",label=T)
CytoTRACE(as.matrix(cNK@assays$RNA@counts),ncores=10)->cNKscore
data.frame(CytoTRACE=as.numeric(cNKscore$CytoTRACE),cluster=cNK@meta.data$seurat_clusters)->dat
ggplot(dat,aes(x=cluster,y=CytoTRACE,fill=cluster))+geom_boxplot()+
scale_fill_manual(values=c("#377EB8", "#984EA3", "#FFD92F", "#F781BF", "#A65628"))+
theme_bw()
###
immune<-c("RGS1","IL2RB","KLRC1","PDCD4","CCL3","EOMES","CD74","CD69","TIGIT","TXNIP","RUNX3","NKG7","ITGB2","TBX21","GNLY","PRF1","GZMB","GZMH","FCGR3A","CX3CR1","FGFBP2")
NK@meta.data$celltype1=factor(as.character(NK@meta.data$celltype1))
levels(NK@meta.data$celltype1)=c("cNK","lrNK","lrNK","lrNK","cNK","cNK","cycling","cNK","cNK")
subset(NK,celltype1=="cNK" | celltype1=="lrNK")->NK2
pdf("./rawdata_NK/image/immun.pdf")
for(i in immune){
  p<-VlnPlot(NK2,features=i,group.by="celltype1",pt.size=0)
  print(p)
}
dev.off()
####

deg[which(deg$p_val_adj<0.05 & deg$avg_logFC>0.35),]->deg1
signature<-list()
for(i in as.character(0:8)){
  deg1[which(deg1$cluster==i),]$gene->signature[[i]]
}

 gsva(as.matrix(LIHC_fpkm),gset.idx.list=signature,method="ssgsea",ssgsea.norm=T,parallel.sz=10)->ssgsea_score    

 read.table("./TCGAdata/TCGA-LIHC.survival.tsv.gz",sep="\t",header=T,row.names=1)->survival  
survival$C1=as.numeric(ssgsea_score["1",rownames(survival)])
for(i in 1:nrow(survival)){
  if(survival$C1[i]<= quantile(survival$C1,0.4)){
    survival$group[i]="low"
  }else{
    survival$group[i]="high"
  }
}

fit<-survfit(Surv(OS.time/30,OS)~group,data=survival)
ggsurvplot(fit,pval=T)


for(i in 1:nrow(OS)){
  if(OS$C6[i]<= quantile(OS$C6,0.6)){
    OS$group6[i]="low"
  }else{
    OS$group6[i]="high"
  }
}

####cibersortX
read.csv("./rawdata_NK/CIBERSORTx_Job2_Results.csv",header=T,check.names=F,row.names=1)->cibersortX
read.table("./TCGAdata/TCGA-LIHC.survival.tsv.gz",sep="\t",row.names=1,header=T)->OS
OS[which(substr(rownames(OS),14,15)=="01"),]->OS
rownames(OS)<-substr(rownames(OS),1,15)
OS[rownames(cibersortX),]->OS
OS$C4=cibersortX$'4'
for(i in 1:nrow(OS)){
  if(OS$C4[i]<=quantile(OS$C4,0.5)){
    OS$group[i]="low"
  }else{
    OS$group[i]="high"
  }
}
fit<-survfit(Surv(OS.time/30,OS)~group,data=OS)
ggsurvplot(fit,pval=T)
######ssgsea
FindAllMarkers(NK,only.pos=T,logfc.threshold=0.25,min.pct=0.25,test.use="MAST")
deg[which(deg$p_val_adj<0.05),]->deg1

signature<-list()
for(i in as.character(0:8)){
  deg1[which(deg1$cluster==i),]$gene->signature[[i]]
}
 gsva(as.matrix(LIHC_fpkm),gset.idx.list=signature,method="ssgsea",ssgsea.norm=T,parallel.sz=10)->ssgsea_score    
read.table("./TCGAdata/TCGA-LIHC.survival.tsv.gz",sep="\t",row.names=1,header=T)->OS
OS[which(substr(rownames(OS),14,15)=="01"),]->OS
rownames(OS)<-substr(rownames(OS),1,15)
ssgsea_score[,intersect(rownames(OS),colnames(ssgsea_score))]->ssgsea_score
OS[intersect(rownames(OS),colnames(ssgsea_score)),]->OS
OS$C1=ssgsea_score["8",]

for(i in 1:nrow(OS)){
  if(OS$C2[i]<= quantile(OS$C2,0.58)){
    OS$group2[i]="low"
  }else{
    OS$group2[i]="high"
  }
}
fit<-survfit(Surv(OS.time/30,OS)~group,data=OS)
ggsurvplot(fit,pval=T,color=c("#2E9FDF",'#E7B800'))

####hallmark
readRDS("./rawdata_NK/NK_score.RDS")->hallmark
library(limma)
annotation_col<-data.frame(Type=NK$seurat_clusters,row.names=colnames(NK))
#把0-clu的定义一组，不是0-clu的定义为另外的组
tvalue<-matrix(ncol=9,nrow=50)->tvalue

colnames(tvalue)<-c("0","1","2","3","4","5","6","7","8")
for(i in colnames(tvalue)){
  some_cluster<-annotation_col
  some_cluster$Type<-as.character(some_cluster$Type)
  some_cluster$Type[which(some_cluster$Type!= i)]<-"control"
some_cluster$Type[which(some_cluster$Type== i)]<-"case"
grouP<-as.factor(some_cluster$Type)
desigN <- model.matrix(~ grouP + 0)
rownames(desigN)<-colnames(hallmark)
comparE <- makeContrasts(grouPcase-grouPcontrol,levels=desigN)
fiT <- lmFit(hallmark, desigN)
fiT2 <- contrasts.fit(fiT, comparE)
fiT3 <- eBayes(fiT2)
Diff<-topTable(fiT3,p.value=1,num=Inf)
tvalue[,i]<-Diff[rownames(hallmark),"t"]
}
rownames(tvalue)=rownames(hallmark)
rownames(tvalue)<-substr(rownames(tvalue),10,nchar(rownames(tvalue)))
gsub("_", " ", rownames(tvalue))->rownames(tvalue)
library(pheatmap)
pheatmap(tvalue,cluster_cols=F,cluster_rows=T,show_colnames=T,show_rownames=T,cellwidth=12,cellheight=8,border_color="white",fontsize_row=8,breaks=unique(c(seq(-100,100,length=500))),color=colorRampPalette(c("navy","white","red"))(500))

###lrNK hallmark
subset(NK,idents=c("1","2","3"))->lrNK
hallmark[,colnames(lrNK)]->hallmark
annotation_col=data.frame(Type=lrNK$origin,row.names=colnames(lrNK))
some_cluster<-annotation_col
  some_cluster$Type<-as.character(some_cluster$Type)
  some_cluster$Type[which(some_cluster$Type!= "T")]<-"control"
some_cluster$Type[which(some_cluster$Type== "T")]<-"case"
grouP<-as.factor(some_cluster$Type)
desigN <- model.matrix(~ grouP + 0)
rownames(desigN)<-colnames(hallmark)
comparE <- makeContrasts(grouPcase-grouPcontrol,levels=desigN)
fiT <- lmFit(hallmark, desigN)
fiT2 <- contrasts.fit(fiT, comparE)
fiT3 <- eBayes(fiT2)
Diff<-topTable(fiT3,p.value=1,num=Inf)



df <- data.frame(ID = rownames(hallmark), score = Diff[rownames(hallmark),"t"])

# 按照score的值分组
cutoff <- 10   # 这里的值可以改的，很灵活的
df$group <- cut(df$score, breaks = c(-Inf, -cutoff, cutoff, Inf),labels = c(1,2,3))
# 按照score排序
sortdf <- df[order(df$score),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
head(sortdf)
ggplot(sortdf, aes(ID, score, fill = group)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  scale_fill_manual(values = c('palegreen3', 'snow3', 'dodgerblue4'), guide = FALSE) + 
  
  # 画2条虚线
  geom_hline(yintercept = c(-cutoff,cutoff), 
             color="white",
             linetype = 2, #画虚线
             size = 0.3) + #线的粗细
  
  # 写label
  geom_text(data = subset(df, score < 0),
            aes(x=ID, y= 0, label= paste0(" ", ID), color = group),#bar跟坐标轴间留出间隙
            size = 3, #字的大小
            hjust = "outward" ) +  #字的对齐方式
  geom_text(data = subset(df, score > 0),
            aes(x=ID, y= -0.1, label=ID, color = group),
            size = 3, hjust = "inward") +  
  scale_colour_manual(values = c("black","snow3","black"), guide = FALSE) +
  
  xlab("") +ylab("t value of GSVA score, tumor \n versus non-malignant")+
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_rect(size = 0.6)) + #边框粗细
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) #去除y轴

####scrat
signature<-deg[which(deg$p_val_adj<0.05 & deg$cluster=="0"),]$gene

pdf("./survival1.pdf")
for(j in 20:84){
  signature1=signature[1:j]
  score=gsva(as.matrix(LIHC_fpkm),gset.idx.list=list(signature1),method="ssgsea",ssgsea.norm=T)
  survival$score<-as.numeric(score)
  for(i in 1:nrow(survival)){
  if(survival$score[i]<=quantile(survival$score,0.5)){
    survival$group[i]="low"
  }else{
    survival$group[i]="high"
  }
}
fit<-survfit(Surv(OS.time,OS)~group,data=survival)
p<-ggsurvplot(fit,pval=T)+ggtitle(j)
print(p)
}
dev.off()



for(i in 1:nrow(survival)){
  if(survival$score[i]<=quantile(survival$score,0.5)){
    survival$group[i]="low"
  }else{
    survival$group[i]="high"
  }
}
fit<-survfit(Surv(OS.time,OS)~group,data=survival)
ggsurvplot(fit,pval=T)





 cellphonedb method statistical_analysis meta.txt  count.txt --verbose  --iterations=10 --threads=10 --pvalue 0.05  --output-path=./interaction --counts-data gene_name
mypvals <- read.delim("./rawdata_NK/cellphonedb/onlyTumor1/pvalues.txt", check.names = FALSE)
mymeans <- read.delim("./rawdata_NK/cellphonedb/onlyTumor1/means.txt", check.names = FALSE)
mypvals[!duplicated(mypvals$interacting_pair),]->mypvals
mymeans[!duplicated(mymeans$interacting_pair),]->mymeans
rownames(mypvals)<-mypvals$interacting_pair
pval.f<-mypvals[,12:ncol(mypvals)]

RR=c('Macrophage|lrNK','lrNK|Macrophage','cNK|Macrophage','Macrophage|cNK')

RR=c('Macrophage|lrNK.1','Macrophage|lrNK.2','lrNK.1|Macrophage','lrNK.2|Macrophage')
RR=c('CLEC9A+ DC|lrNK.1','CLEC9A+ DC|lrNK.2','CD1C+ DC|lrNK.1','CD1C+ DC|lrNK.2','LAMP3+ DC|lrNK.1','LAMP3+ DC|lrNK.2','Macrophage|lrNK.1','Macrophage|lrNK.2','CLEC9A+ DC|cNK','CD1C+ DC|cNK','LAMP3+ DC|cNK','Macrophage|cNK')
pval.f[,RR]->pval.f
index<-c()
for(i in 1:nrow(pval.f)){
	if(length(which(pval.f[i,]<0.05))>0){
		index[i]=i
	}else{
		index[i]=NA
	}
}
as.numeric(as.character(na.omit(index)))->index
pval.f[index,]->pval.f
pval.f[PR,]->pval.f
rownames(mymeans)=mymeans$interacting_pair
means.f<-mymeans[,12:ncol(mymeans)]
means.f<-means.f[rownames(pval.f),RR]

df_names = expand.grid(rownames(pval.f), colnames(pval.f))
pval.f1 = unlist(pval.f)
pval.f1[pval.f1==0] = 0.00009
plot.data = cbind(df_names,pval.f1)
pr = unlist(as.data.frame(means.f))
pr[pr==0] = 1
plot.data = cbind(plot.data,log2(pr))
colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')


plot.data[which(plot.data$pair=='CXCR3_CXCL9' | plot.data$pair=='KIR2DL3_CXCL9' | plot.data$pair=='TNFSF14_LTBR' | plot.data$pair=='CXCR3_CCL19') ,]->plot.data1

pdf("./rawdata_NK/image/inetr_cDC1.pdf",useDingbats=F,height=13,width=10)
ggplot(plot.data,aes(x=clusters,y=pair)) +
  geom_point(aes(size=-log10(pvalue),color=mean)) +
   scale_size_continuous(range = c(1,4))+
   scale_color_distiller(palette = "RdBu")+
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=14, colour = "black"),
        axis.text.x = element_text(angle = 90,colour="black", hjust = 1),
        axis.text.y = element_text(size=8, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
dev.off()


read.table("./count_network.txt",sep='\t',header=T,check.names=F)->mynet
library(psych)
library(qgraph)
library(igraph)
library(tidyverse)
net<- graph_from_data_frame(mynet)
plot(net)

allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB",
            "#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3",
            "#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
            "#40E0D0","#5F9EA0","#FF1493",
            "#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
            "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347",
            "#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

karate_groups <- cluster_optimal(net)
coords <- layout_in_circle(net, order =
                             order(membership(karate_groups)))  # 设置网络布局

E(net)$width  <- E(net)$count/10  # 边点权重（粗细）
plot(net, edge.arrow.size=.1, 
     edge.curved=0,
     vertex.color=allcolour,
     vertex.frame.color="#555555",
     vertex.label.color="black",
     layout = coords,
     vertex.label.cex=.7) 


net2=net
length(unique(mynet$SOURCE)) # 查看需要绘制多少张图，以方便布局
pdf("./inter.pdf")
par(mfrow=c(2,7), mar=c(.3,.3,.3,.3))

for (i in 1: length(unique(mynet$SOURCE)) ){
  net1<-net2
  E(net1)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$color <- allcolour[i]

  plot(net1, edge.arrow.size=.1, 
       edge.curved=0.4,
       vertex.color=allcolour,
       vertex.frame.color="#555555",
       vertex.label.color="black",
       layout = coords,
       vertex.label.cex=1) 

}





qusage::read.gmt("./NKcyto.gmt")->NKcyto
data.frame(ont="Notch signaling",gene=hallmark$HALLMARK_NOTCH_SIGNALING)->Notch
geneList<-data.frame(Gene=rownames(deg),FC=as.numeric(deg$avg_logFC))
glist <- geneList[,2]
names(glist) <- as.character(geneList[,1])
 glist <- sort(glist,decreasing = T)
gsea<-GSEA(glist, TERM2GENE=Notch, verbose=T, pvalueCutoff = 1)##
enrichplot::gseaplot2(gsea,1)



exau=c("PDCD1", "LAYN", "LAG3", "HAVCR2", "CD244","CD160")

list("Dysfunctional"=c("CTLA4","TIGIT","HLA-DPA1","HLA-DQB1",
"HLA-DRB1",
"HLA-DRB5",
"LAG3",
"PDCD1",
"CXCR6",
"PIK3R6",
"GOLIM4",
"GZMB",
"HLA-DMA",
"PDE4D",
"CTSA",
"CCL3",
"HAVCR2",
"PTMS",
"CCL5",
"VCAM1",
"BHLHE40",
"CD74"))->Dysfunctional


#####
plot1 <- FeatureScatter(NK, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(NK, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
NK <- subset(NK, subset = nFeature_RNA > 200 & nFeature_RNA < 2000   & percent.mt<15 )



#conda activate pyscenic
library(reticulate)
#use_python("/bioapps/rhel7/python-3.6.6/bin/python3.6",required=T)

use_python("./miniconda3/envs/pyscenic/bin/python3.7",required=T)
py_config()
sc <- import('scanpy', convert = FALSE)
scvi <- import('scvi', convert = FALSE)
scvi$settings$progress_bar_style = 'tqdm'
library(Seurat)
readRDS("./youwh/NK.RDS")->NK
adata <- sc$AnnData(
  X   = Matrix::t(Matrix::as.matrix(GetAssayData(NK,slot='counts'))), #scVI requires raw counts
  obs = NK[[]],
  var = GetAssay(NK)[[]]
)
scvi$data$setup_anndata(adata,batch_key="orig.ident")
model = scvi$model$SCVI(adata)
model$train()
latent = model$get_latent_representation()

latent <- as.matrix(latent)
rownames(latent) = colnames(NK)
NK[['scvi']] <- CreateDimReducObject(embeddings = latent, key = "scvi_", assay = DefaultAssay(NK))

NK <- FindNeighbors(NK, dims = 1:10, reduction = 'scvi')
NK <- FindClusters(NK, resolution =0.4)
NK <- RunUMAP(NK, dims = 1:10, reduction = 'scvi')
saveRDS(NK,file="./rawdata_NK/NK_scvi.RDS")
pdf("./youwh/NK.pdf",useDingbats=F)
DimPlot(NK, reduction = "umap",label=T)
dev.off()
pdf("./youwh/features.pdf",useDingbats=F)
c("CXCR6","MKI67","CX3CR1","FCGR3A","ITGA1")->feature
for(i in feature){
  p<-FeaturePlot(NK,features=i,reduction="umap")
  print(p)
}
dev.off()


###NicheNet
library(nichenetr)
library(tidyverse)
ligand_target_matrix = readRDS("./ref/nichenet/ligand_target_matrix.rds")
c("NOTCH1",'NOTCH2','NOTCH3','NOTCH4','JAG1','JAG2','DLL3','DLL4','WNT4')->geneset_oi
background_expressed_genes=

ILC <- NormalizeData(ILC, normalization.method = "LogNormalize", scale.factor = 10000)
ILC <- FindVariableFeatures(ILC, selection.method = "vst", nfeatures = 4000)
ILC <- ScaleData(ILC, features = rownames(ILC))
ILC <- RunPCA(ILC, features = VariableFeatures(object = ILC))
ILC <- FindNeighbors(ILC, dims = 1:50,reduction = "pca")
ILC <- FindClusters(ILC, resolution = 0.3)
ILC <- RunUMAP(ILC, dims = 1:50,reduction = "pca")
DimPlot(ILC,reduction='umap',label=T)
FindAllMarkers(ILCN,only.pos=T,logfc.threshold=0.25,min.pct=0.25)->deg
FeaturePlot(HCC,features=c('CD3D','CD4','CD8A','NCR1'),order=T,reduction='umap')


require(dplyr)
# arrange data
coords_tbl = bind_cols(cellName = rownames(coords), as.data.frame(coords))

join_vec = setNames(colnames(labelData)[1], nm = colnames(coords_tbl)[1])
cellinfo_tbl = left_join(coords_tbl, labelData, by = join_vec)

density_obj = getDensity3D(cellinfo_tbl$x, cellinfo_tbl$y, cellinfo_tbl$z)
cellinfo_tbl = cellinfo_tbl %>% mutate(density = density_obj)

p_3Ddensity = plot3D(csomap, color_by = "density", title = "3D density")

plot3d(csomap$x, csomap$y, csomap$z,type = "p", col=csomap$label, 
      add = FALSE,  xlim = NULL, ylim = NULL, zlim = NULL, 
      forceClipregion = FALSE)



for(i in colnames(per)){
  for(j in rownames(per)){
    per[j,i]=nrow(a[which(a$mouse==j & a$mc==i),])/nrow(a[which(a$mouse==j),])
  }
}

NK=data.frame(per=per[,1])

###
library(nichenetr)
library(Seurat) # please update to Seurat V4
ligand_target_matrix=readRDS("./ref/nichenet/ligand_target_matrix.rds")
lr_network = readRDS('./ref/nichenet/lr_network.rds')
weighted_networks = readRDS("./ref/nichenet/weighted_networks.rds")

HCC@active.ident=HCC$type

nichenet_output = nichenetr::nichenet_seuratobj_aggregate(
  seurat_obj = HCC, 
  receiver = "cNK", 
  condition_colname = "Tissue", condition_oi = "Tumor", condition_reference = "Normal", 
  sender = c("LAMP3+ DC",'B cell','CD1C+ DC','CLEC9A+ DC'), 
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks, organism = "human")


##
deg=rbind(deg[which(deg$cluster=='2' & deg$avg_logFC>0.4),],deg[which(deg$cluster=='3' & deg$avg_logFC>0.3),],deg[which(deg$cluster=='4'&deg$avg_logFC>0.35),])
cluster_markers<-list()
for(i in 1:7){
  cluster_i <- deg %>% 
    filter(cluster == i-1)
    cluster_markers[[paste("C",i,sep='')]] <- as.numeric(bitr(cluster_i$gene, fromType="SYMBOL", toType= "ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID)
}
res<-compareCluster(geneCluster = cluster_markers, fun = "enrichGO",ont='BP',OrgDb='org.Hs.eg.db')
res<-compareCluster(geneCluster = diff, fun = "enrichKEGG",organism = "hsa")
res<-compareCluster(geneCluster = diff, fun = "enrichGO",ont='BP',OrgDb='org.Hs.eg.db')

dotplot(res)+scale_color_distiller(palette='Reds')+
  theme(axis.text.x = element_text(angle =45, hjust = 1,size = 8),
        axis.text.y = element_text(size = 10),
        strip.text.x = element_text(size = 14),
        axis.title =element_blank(),legend.position = 'none')+
        scale_y_discrete(labels=function(x) str_wrap(x, width=20))+
        scale_x_discrete(labels=function(x) str_wrap(x, width=8))+
        coord_flip()
###
c("GADD45A",'NCOR2','SPEN','GSK3B','RPS6KA6','NOTCH1','ATXN1','MYC','MAML1','MED27','ELAVL1','ATXN1L','KDM1A','HK2','ACVR1B','NOTCH2','ENO3','TACSTD2','CYB5A','PPP1R12A')->rbpjass
data.frame(ont="RBPJ-associated genes",gene=rbpjass)->rbpj
geneList<-data.frame(Gene=rownames(deg),FC=as.numeric(deg$avg_logFC))
glist <- geneList[,2]
names(glist) <- as.character(geneList[,1])
 glist <- sort(glist,decreasing = T)
gsea<-GSEA(glist, TERM2GENE=Notch, verbose=T, pvalueCutoff = 1)##
enrichplot::gseaplot2(gsea,1)

###
library(metacell)
setwd("./rawdata_NK/GSE158485_RAW/")
if(!dir.exists("liverdb")) dir.create("liverdb/")
scdb_init("liverdb/", force_reinit=T)
if(!dir.exists("figs")) dir.create("figs/")
scfigs_init("figs/")
scdb_add_mat("liver",mat)
mcell_plot_umis_per_cell("liver")
nms = c(rownames(mat@mat), rownames(mat@ignore_gmat)) 
nms <- unique(c(rownames(mat@mat), rownames(mat@ignore_gmat)))
pre_nr_term <- c("^RPS","^RPL","^MT-","^MTMR","^MTND","^MTRN","^MTCO","^MRPL","^MRPS","^HBA","^HBB","^MTATP")
pre_nr_genes <- foreach(i=pre_nr_term, .combine = c) %do% grep(i, nms, v=T)
pre_ex_genes <- c("MALAT1", "XIST", "XIST_intron")
pre_bad_genes <- unique(c(pre_nr_genes, pre_ex_genes))
pre_bad_genes
##
genes_RBC <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
cells_RBC <- names(which((apply(mat@mat[intersect(genes_RBC,rownames(mat@mat)),],2,sum))>=1))
genes_Mt<-foreach(i='^MT-', .combine = c) %do% grep(i, nms, v=T)
cells_mt= names(which((Matrix::colSums(mat@mat[genes_Mt,])/Matrix::colSums(mat@mat))>0.2))

mcell_mat_ignore_genes(new_mat_id='liver', mat_id='liver', pre_bad_genes, reverse=F)

cell=c(names(which(apply(mat@mat,2,sum)>8000)),names(which(apply(mat@mat,2,sum)<200)))
unique(c(cells_RBC,cell,cells_mt))->rmcells
mcell_mat_ignore_cells(new_mat_id='liver', mat_id='liver', ig_cells = rmcells, reverse = F)

genes_anchors = c('FOS','FOSB','NFKBIA','NFKBIZ','JUN','ZFP36','ISG15','HMGB2','STMN1','TOP2A','MKI67','MX1','RSAD2')
tab_fn = "./lateral_gmods.txt"
mcell_mat_rpt_cor_anchors(mat_id='liver', gene_anchors = genes_anchors, cor_thresh = 0.1,
                          gene_anti = c(), tab_fn = tab_fn, sz_cor_thresh = 0.2)
gcor_mat = read.table('./lateral_gmods.txt', header=T)
foc_genes = apply(gcor_mat[, intersect(colnames(gcor_mat),genes_anchors)], 1, which.max)
#
mat = scdb_mat("liver")
mcell_add_gene_stat(gstat_id="liver", mat_id="liver", force=T)
mcell_gset_filter_varmean(gset_id="liver_feats", gstat_id="liver", T_vm=0.08, force_new=T)
mcell_gset_filter_cov(gset_id = "liver_feats", gstat_id="liver", T_tot=100, T_top3=2)
mcell_plot_gstats(gstat_id="liver", gset_id="liver_feats")
gset <- scdb_gset("liver_feats")
pst_genes <- names(gset@gene_set)
pst_nr_term <- c("^AC[0-9]+\\.", "^AP[0-9]+\\.", "^AL[0-9]+\\.", "-AS[0-9]*$", "^MIR[0-9]", "^LINC[0-9]", "^SNOR", "^sno", "^SNR",
                 "^TMSB", "^HIST", "^HSP", "^IFI", "^HLA-", "^ATP", "-", ":", "\\.", '^KIAA',
                 "^IGJ", "^IGH", "^IGK", "^IGL", "^DNAJ", "^GZM", "^CCL", "^XCL", '^FTH', '^FTL', '^LGALS')
pst_nr_genes <- foreach(i=pst_nr_term, .combine = c) %do% grep(i, pst_genes, v=T)
pst_ex_genes <- c()
pst_bad_genes <- unique(c(pst_nr_genes, pst_ex_genes, names(foc_genes)))
pst_add_genes <- c()
final_genes <- unique(setdiff(pst_genes, pst_bad_genes), pst_add_genes)
final_genes

gset = gset_new_gset(sets = gset@gene_set[final_genes], desc = "final genes")
scdb_add_gset("liver_feats", gset)


mcell_add_cgraph_from_mat_bknn(mat_id="liver",
                gset_id = "liver_feats",
                graph_id="liver_graph",
                K=50,
                dsamp=T)###一般K为总细胞数的平方根的数目，若细胞数很少，K=20-40

mcell_coclust_from_graph_resamp(
                coc_id="liver_coc1000",
                graph_id="liver_graph",
                min_mc_size=40,
                p_resamp=0.75, n_resamp=1000)

mcell_mc_from_coclust_balanced(
                coc_id="liver_coc1000",
                mat_id= "liver",
                mc_id= "liver_mc",
                K=40, min_mc_size=40, alpha=2)

mcell_plot_outlier_heatmap(mc_id="liver_mc", mat_id = "liver", T_lfc=3)
mcell_mc_split_filt(new_mc_id="liver_mc",
            mc_id="liver_mc",
            mat_id="liver",
            T_lfc=3, plot_mats=F)

mc = scdb_mc("liver_mc")
mc@colors <- colorRampPalette(c("darkgray", "burlywood1", "chocolate4","orange", "red", "purple", "blue","darkgoldenrod3", "cyan"))(ncol(mc@mc_fp))
scdb_add_mc("liver_mc",mc)
mcell_gset_from_mc_markers(gset_id="liver_markers", mc_id="liver_mc")
mcell_mc_plot_marks(mc_id="liver_mc", gset_id="liver_markers", mat_id="liver")
##Projecting metacells and cells in 2D
mcell_mc2d_force_knn(mc2d_id="liver_2dproj",mc_id="liver_mc", graph_id="liver_graph")
tgconfig::set_param("mcell_mc2d_height",1000, "metacell")
tgconfig::set_param("mcell_mc2d_width",1000, "metacell")
mcell_mc2d_plot(mc2d_id="liver_2dproj")
mcell_mc2d_plot_gene(mc2d_id="liver_2dproj",gene="ZNF683")
mcell_mc2d_plot_gene(mc2d_id="liver_2dproj",gene="NCAM1")
mcell_mc2d_plot_gene(mc2d_id="liver_2dproj",gene="CXCR6")
mcell_mc2d_plot_gene(mc2d_id="liver_2dproj",gene="EOMES",max_lfp=2,min_lfp=-2)
mcell_mc2d_plot_gene(mc2d_id="liver_2dproj",gene="ITGA1")

mc_hc <-mcell_mc_hclust_confu(mc_id="liver_mc",graph_id="liver_graph")
mc_sup <- mcell_mc_hierarchy(mc_id="liver_mc",mc_hc=mc_hc, T_gap=0.04)

mcell_mc_plot_hierarchy(mc_id="liver_mc",
                        graph_id="liver_graph",
                        mc_order=mc_hc$order,
                        sup_mc = mc_sup,
                        width=2800,
                        height=2000,
                        min_nmc=2,
                         show_mc_ids=TRUE)


lfp <- log2(mc@mc_fp)
plt = function(gene1, gene2, lfp, colors)
{
    plot(lfp[gene1, ], lfp[gene2, ], pch=21, cex=3, bg=colors, xlab=gene1, ylab=gene2)
    text(lfp[gene1, ], lfp[gene2, ], colnames(lfp))

}
plt(gene1 = 'CD3G', gene2 = 'NCAM1', lfp = lfp, colors = mc@colors)
##
mc@annots=c(rep("ILC1",7),rep("ILC3",3),'ILC2b',rep("ILC2a",4))
mc@annots<-factor(mc@annots)
mc@colors=mc@annots
library(chameleon)
levels(mc@colors)=distinct_colors(4,minimal_saturation = 15,minimal_lightness = 15,maximal_lightness = 100)$name
mc@colors<-as.character(mc@colors)
mc@annots<-as.character(mc@annots)

scdb_add_mc("mc_annot",mc)
scdb_mc2d("liver_2dproj")->mc2d
mcell_mc2d_plot_by_factor(
  mc2d_id='liver_2dproj',
  mat_id='liver',
  meta_field='sub',
  meta_data_vals = NULL,
  single_plot = F,
  filter_values = NULL,
  filter_name = NULL,
  ncols = NULL,
  neto_points = F,
  colors = mc@colors
)



ggplot(per,aes(x=group,y=per,fill=mc.mc))+
geom_col(position="stack")+
scale_fill_d3()+
theme_bw()+
theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())


ggplot(score,aes(x=RBPJ,y=score)) + 
  stat_density2d(aes(alpha = 1,fill = ..density..),geom = "raster", contour = FALSE) +
  scale_fill_gradient (low = "#FFFFFF", high = "#377EB8") +
  ylab('ssgsea score') + xlab('log2(RBPJ FPKM)') +
  stat_smooth(method="lm",se=T) + 
  stat_cor(method = "pearson",size=2,label.x = 0.1,label.y =0.1) + 
  geom_point(colour='#377EB8',size=0.5) + theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="none",axis.line = element_line(colour = "black"), axis.text=element_text(colour="black"))  
#####

dt = readRDS(paste0("D:\\Dropbox\\Scdata\\result\\pca\\pca_patienttissue\\all\\dat_hvg_m2.rds")) # read in file
dt$response[dt$resi_tumor>0.1]<-'NR'
dt$response[is.na(dt$response)==T]<-'R'
min = min(dt$PC1,dt$PC2)
max = max(dt$PC1,dt$PC2)
ggplot(dt, aes(x=PC1,y=PC2,color =tissue)) +
  geom_point(size=3) +
  geom_point(shape = 1,size = 3,colour = "black")+
  theme_classic() +
  scale_color_jco()+
  xlim(min,max) + ylim(min,max)+
  theme(axis.text.x = element_text(angle =45, hjust = 1,size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        strip.text.x = element_text(size = 14,face = "bold"),
        axis.title =element_blank(),legend.position = 'none')
ggsave("D:\\Dropbox\\Scdata\\result\\ms\\Fig1\\pca.tissue.pdf",height =3,width =3.5)




mat <- NormalizeData(mat, normalization.method = "LogNormalize", scale.factor = 10000)
mat <- FindVariableFeatures(mat, selection.method = "vst", nfeatures = 3000)
mat <- ScaleData(mat, features = rownames(mat))
mat <- RunPCA(mat, features = VariableFeatures(object = mat))
mat <- FindNeighbors(mat, dims = 1:40,reduction = "pca")
mat <- FindClusters(mat, resolution = 0.5)
mat <- RunUMAP(mat, dims = 1:40,reduction = "pca")
DimPlot(mat,reduction='umap',label=T,group.by='type')


##

suppression<-c('TIGIT','PDCD1','SIGLEC7','TGFBR3','TGFBR1','TGFBR2','CD96','LAIR1','KLRC1','LILRB1','KIR2DL1','KIR2DL5','KIR2DL2','KIR3DL1','KIR3DL2','KLRB1')





mat <- NormalizeData(mat, normalization.method = "LogNormalize", scale.factor = 10000)
mat <- FindVariableFeatures(mat, selection.method = "vst", nfeatures = 4000)
mat <- ScaleData(mat, features = rownames(mat))
mat <- RunPCA(mat, features = VariableFeatures(object = mat))
mat <- FindNeighbors(mat, dims = 1:20,reduction = "pca")
mat <- FindClusters(mat, resolution = 0.3)
mat <- RunUMAP(mat, dims = 1:20,reduction = "pca")
DimPlot(mat,reduction='umap',label=T)


GC[["percent.mt"]] <- PercentageFeatureSet(GC, pattern ="^MT-")
plot1 <- FeatureScatter(GC, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(GC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
GC[["percent.HB"]] <- PercentageFeatureSet(GC, features=c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ"))

lrNK <- NormalizeData(lrNK, normalization.method = "LogNormalize", scale.factor = 10000)
lrNK <- FindVariableFeatures(lrNK, selection.method = "vst", nfeatures = 4000)
lrNK <- ScaleData(lrNK, features = rownames(lrNK))
lrNK <- RunPCA(lrNK, features = VariableFeatures(object = lrNK))
lrNK <- FindNeighbors(lrNK, dims = 1:10,reduction = "scvi")
lrNK <- FindClusters(lrNK, resolution = 0.3)
lrNK <- RunUMAP(lrNK, dims = 1:10,reduction = "scvi")
DimPlot(lrNK,reduction='umap',label=T,group.by='subtype')

lrNK@active.ident=factor(as.character(lrNK$subtype))
names(lrNK@active.ident)=colnames(lrNK)

ident.colors <- (scales::hue_pal())(n = length(x = levels(x = lrNK)))
names(x = ident.colors) <- levels(x = lrNK)
cell.colors <- ident.colors[Idents(object = lrNK)]
names(x = cell.colors) <- colnames(x = lrNK)
show.velocity.on.embedding.cor(emb = Embeddings(object = lrNK, reduction = "umap"), vel = Tool(object = lrNK, 
    slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
    cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
    do.par = FALSE, cell.border.alpha = 0.1)


plot1 <- FeatureScatter(HCC, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(HCC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
HCC <- subset(HCC, subset = nFeature_RNA > 200 & nFeature_RNA < 7500  & percent.mt<50 & percent.rb<10)


library(circlize)
pdf("./rawdata_NK/circlize_lrNK.2.pdf",useDingbats=F)
chordDiagram(network, grid.col =distinct_colors(16,minimal_saturation = 5,minimal_lightness = 14,maximal_lightness = 100)$name,annotationTrackHeight = c(0.05, 0.05),annotationTrack = c("name", "grid"),direction.type = "arrows", link.arr.length = 0.2,link.visible=network[[1]]=="lrNK")
circos.clear()
dev.off()

####
library(iTALK)
library(Seurat)
library(Matrix)
library(dplyr)
library(SeuratDisk)
LoadH5Seurat("./data/HCC/HCC.h5seurat")->HCC
#early<-subset(HCC,stage=='early')
#advanced<-subset(HCC,stage=='advanced')
colors <-chameleon::distinct_colors(20,minimal_saturation = 12,minimal_lightness = 14,maximal_lightness = 100)$name
iTalk_data <- as.data.frame(t(HCC@assays$RNA@counts))
iTalk_data$cell_type <- as.character(HCC@meta.data$celltype)
iTalk_data$compare_group <- as.character(HCC@meta.data$stage)
cell_types <- as.character(unique(iTalk_data$cell_type))
cell_col <- structure(colors[1:length(cell_types)], names=cell_types)
rm(HCC)
highly_exprs_genes <- rawParse(iTalk_data, top_genes=50, stats="mean")
# 通讯类型
iTalk_res <- NULL
comm_list<-c("checkpoint")
for(comm_type in comm_list){
  res_cat <- FindLR(highly_exprs_genes, datatype='mean count', comm_type=comm_type)
  iTalk_res <- rbind(iTalk_res, res_cat)
}

iTalk_res <- iTalk_res[order(iTalk_res$cell_from_mean_exprs*iTalk_res$cell_to_mean_exprs,decreasing=T),]
#iTalk_res1 <- iTalk_res[which(iTalk_res$cell_to=='lrNK'),]
NetView(iTalk_res,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
iTalk_res1 <- iTalk_res[order(iTalk_res$cell_from_mean_exprs*iTalk_res$cell_to_mean_exprs,decreasing=T),][1:30,]
LRPlot(iTalk_res1,datatype='mean count',cell_col=cell_col,link.arr.lwd=iTalk_res$cell_from_mean_exprs[1:30],link.arr.width=iTalk_res$cell_to_mean_exprs[1:30])



##
deg_lrNK<-DEG(iTalk_data %>% filter(cell_type=='lrNK'),method='MAST',contrast=c('early', 'advanced'))
deg_Macrophage<-DEG(iTalk_data %>% filter(cell_type=='Macrophage'),method='MAST',contrast=c('early', 'advanced'))
deg_mDC<-DEG(iTalk_data %>% filter(cell_type == 'mDC'),method='MAST',contrast=c('early', 'advanced'))
deg_pDC<-DEG(iTalk_data %>% filter(cell_type == 'pDC'),method='MAST',contrast=c('early', 'advanced'))
deg_cDC1<-DEG(iTalk_data %>% filter(cell_type == 'cDC1'),method='MAST',contrast=c('early', 'advanced'))
deg_cDC2<-DEG(iTalk_data %>% filter(cell_type == 'cDC2'),method='MAST',contrast=c('early', 'advanced'))

#comm_list<-c('checkpoint')
res<- NULL
for(cell_type in c('deg_Macrophage','deg_mDC','deg_pDC','deg_cDC1','deg_cDC2')){
  res_cat<-FindLR(deg_lrNK, get(cell_type), datatype='DEG',comm_type='cytokine')
  #res_cat<-FindLR(deg_lrNK,  datatype='DEG',comm_type=comm_type)
  res<-rbind(res,res_cat)
}
# FindLR DEG类型的数据，可以输入一个基因集合，结果为相应基因内的配体-受体列表
# 如果有超过20组配体-受体结果，取前20进行展示
res=res[order(res$cell_from_logFC*res$cell_to_logFC,decreasing=T),]
res[which(res$cell_from=='lrNK'),]->res1
res1=res[order(res1$cell_from_logFC*res1$cell_to_logFC,decreasing=T),]

res1[1:40,]->res1
LRPlot(res1, datatype='DEG', link.arr.lwd=res1$cell_from_mean_exprs,
       cell_col=cell_col, link.arr.width=res1$cell_to_mean_exprs)
####
mypvals <- read.delim("./early/pvalues.txt", check.names = FALSE)
mymeans <- read.delim("./early/means.txt", check.names = FALSE)

mypvals <- read.delim("./all/pvalues.txt", check.names = FALSE)
mymeans <- read.delim("./all/means.txt", check.names = FALSE)

mypvals[!duplicated(mypvals$interacting_pair),]->mypvals
mymeans[!duplicated(mymeans$interacting_pair),]->mymeans
rownames(mypvals)<-mypvals$interacting_pair
pval.f<-mypvals[,12:ncol(mypvals)]
RR=c('early_lrNK|early_Macrophage',
'advanced_lrNK|advanced_Macrophage')

RR=c('Macrophage|lrNK.1','Macrophage|lrNK.2','lrNK.1|Macrophage','lrNK.2|Macrophage')
RR=c('CLEC9A+ DC|lrNK.1','CLEC9A+ DC|lrNK.2','CD1C+ DC|lrNK.1','CD1C+ DC|lrNK.2','LAMP3+ DC|lrNK.1','LAMP3+ DC|lrNK.2','Macrophage|lrNK.1','Macrophage|lrNK.2','CLEC9A+ DC|cNK','CD1C+ DC|cNK','LAMP3+ DC|cNK','Macrophage|cNK')
pval.f[,RR]->pval.f
index<-c()
for(i in 1:nrow(pval.f)){
	if(length(which(pval.f[i,]<0.05))>0){
		index[i]=i
	}else{
		index[i]=NA
	}
}
as.numeric(as.character(na.omit(index)))->index
pval.f[index,]->pval.f
#pval.f[PR,]->pval.f
rownames(mymeans)=mymeans$interacting_pair
means.f<-mymeans[,12:ncol(mymeans)]
means.f<-means.f[rownames(pval.f),RR]

df_names = expand.grid(rownames(pval.f), colnames(pval.f))
pval.f1 = unlist(pval.f)
pval.f1[pval.f1==0] = 0.00009
plot.data = cbind(df_names,pval.f1)
pr = unlist(as.data.frame(means.f))
pr[pr==0] = 1
plot.data = cbind(plot.data,log2(pr))
colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')


plot.data[which(plot.data$pair=='CXCR3_CXCL9' | plot.data$pair=='KIR2DL3_CXCL9' | plot.data$pair=='TNFSF14_LTBR' | plot.data$pair=='CXCR3_CCL19') ,]->plot.data1

pdf("./rawdata_NK/image/inetr_cDC1.pdf",useDingbats=F,height=13,width=10)
ggplot(plot.data,aes(x=clusters,y=pair)) +
  geom_point(aes(size=-log10(pvalue),color=mean)) +
   #scale_size_continuous(range = c(2,4))+
   scale_color_distiller(palette = "RdBu")+
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=14, colour = "black"),
        axis.text.x = element_text(angle = 90,colour="black", hjust = 1),
        axis.text.y = element_text(size=8, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
dev.off()
chordDiagram(early, grid.col =chameleon::distinct_colors(16,minimal_saturation = 5,minimal_lightness = 14,maximal_lightness = 100)$name,annotationTrackHeight = c(0.05, 0.05),annotationTrack = c("name", "grid"),direction.type = "arrows", link.arr.length = 0.2,link.visible=early[[1]]=="early_lrNK")
#####





library(grid)
library(pbapply)
library(circlize)
library(ggsci)
library(plyr)
library(ggplot2)
library(ggrepel)
library(dendextend)
library(ComplexHeatmap)

regulons <- read.csv("./rawdata_NK/scenic/regulons.csv",header=T, sep = ",")
regulon.names <- sapply(strsplit(as.character(regulons$TF), split = "\\)"), function(x) x[1])
regulon.names<-gsub("\\'","",regulon.names)

regulon.sizes <- sapply(strsplit(as.character(regulons$Enrichment.TargetGenes), split = ","), function(x) length(x))
regulon.names2 <- regulon.names[regulon.sizes>=30]
auc<-read.table("./rawdata_NK/scenic/auc_mtx.csv",header=T,row.names=1,sep=",")
colnames(auc)<-gsub("\\.","",colnames(auc))
auc2 <- auc[, which(colnames(auc) %in% regulon.names2)]
pccMat <- cor(auc2)

CSI <- function(r1, r2) {
  delta <- pccMat[r1,r2]
  r.others <- setdiff(colnames(pccMat), c(r1,r2))
  N <- sum(pccMat[r1, r.others] < delta) + sum(pccMat[r2, r.others] < delta)
  M <- length(r.others) * 2
  return(N/M)
}

csiMat <- pblapply(rownames(pccMat), function(i) sapply(colnames(pccMat), function(j) CSI(i, j)))
csiMat <- do.call(rbind, csiMat)
rownames(csiMat) <- rownames(pccMat)

csiMat.binary <- matrix(as.numeric(csiMat >= 0.7), nrow = nrow(csiMat))
colnames(csiMat.binary) <- colnames(csiMat)
rownames(csiMat.binary) <- rownames(csiMat)
csiMat.binary[1:10,1:10]
h = 3
row_dend = as.dendrogram(hclust(dist(csiMat), method = "complete"))
clusters <- cutree(row_dend, h = h) # dendextend::cutree()
row_dend = color_branches(row_dend, h = h, col = pal_d3("category20")(20))




####
library(symphony)
reference = symphony::buildReference(
    NK4@assays$RNA@data,                   # reference expression (genes by cells)
    NK4@meta.data,              # reference metadata (cells x attributes)
    vars = c('patient'),         # variable(s) to integrate over
    K = 100,                   # number of Harmony soft clusters
    verbose = TRUE,            # display verbose output
    do_umap = TRUE,            # run UMAP and save UMAP model to file
    do_normalize = FALSE,      # perform log(CP10k) normalization on reference expression
    vargenes_method = 'vst',   # variable gene selection method: 'vst' or 'mvp'
    vargenes_groups = 'louvain', # metadata column specifying groups for variable gene selection within each group
    topn = 5000,               # number of variable genes (per group)
    theta = 2,                 # Harmony parameter(s) for diversity term
    d = 20,                    # number of dimensions for PCA
    save_uwot_path = './ICB_model', # file path to save uwot UMAP model
    additional_genes = NULL    # vector of any additional genes to force include
)

umap_labels = cbind(NK4@meta.data, reference$umap$embedding)
plotBasic(umap_labels, title = 'Reference', color.by = 'louvain')
query = mapQuery(NK@assays$RNA@data,             # query gene expression (genes x cells)
                 NK@meta.data,        # query metadata (cells x attributes)
                 reference,             # Symphony reference object
                 vars = NULL,           # Query batch variables to harmonize over (NULL treats query as one batch)
                 do_normalize = FALSE,  # perform log(CP10k) normalization on query (set to FALSE if already normalized)
                 do_umap = TRUE)   
 query = knnPredict(query, reference, reference$meta_data$louvain, k = 5)


reference$meta_data$cell_type_pred_knn = NA
reference$meta_data$cell_type_pred_knn_prob = NA
reference$meta_data$ref_query = 'reference'
query$meta_data$ref_query = 'query'

# Add the UMAP coordinates to the metadata
meta_data_combined = rbind(query$meta_data[,c('cell_type_pred_knn','cell_type_pred_knn_prob','ref_query')],reference$meta_data[,c('cell_type_pred_knn','cell_type_pred_knn_prob','ref_query')])
umap_combined = rbind(query$umap,reference$umap$embedding)
umap_combined_labels = cbind(meta_data_combined, umap_combined)

plotBasic(umap_combined_labels, title = 'Reference and query cells', 
          color.by = 'cell_type_pred_knn', facet.by = 'ref_query')+
scale_color_manual()


dat <- dat %>% NormalizeData(verbose = T) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000, verbose =  T) %>% 
  ScaleData(verbose = T) %>%
  RunPCA(npcs = 50, verbose =  T)
dat <- dat %>% RunBBKNN(reduction = "pca", 
                                    run_TSNE = T,
                                    batch_key = "orig.ident")
dat <-dat %>% FindNeighbors(reduction = "pca", k.param = 10, dims = 1:40) %>% 
  FindClusters(resolution = 1, algorithm = 1, graph.name="bbknn")%>% 
  identity()


 DimPlot(dat,group.by='cell_type',reduction='tsne',label=F)+
 scale_color_manual(values=c('B cell'='#000000','T_naive'='#628DE1','Endothelial-Lyve1'='#979CE1','AT2'='#437193','IM'='#A9F1F9','neutrophils'='#8BC1CB',
 'NK'='#3780E9','Monocyte'='#254585','CD4T-Tnfrsf4'='#00007D','CD8T-Ccr9'='#D3B188','SPP1+ Macrophage'='#CDA100',
 'NKT'='#E96BAA','Fibroblast'='#B85880','AM'='#723955','cDC1'='#EEF37C','Ciliated'='#AA5828','gdT'='#1A421C','Endothelial-Kdr'='#86B400','cDC2'='#74876B',
 'Club'='#C40002','SMC'='#E28B60','Profilterating T'='#B4D9E5','AT1'='#466EE7','Profilterating SPP1+ Macrophage'='#6521AE','ILC2'='#2C3A74','mDC'='#54835A','Basophil'='#F08B64','Plasma'='#9D54B6',
 'Mesotheliocyte'='#7AC1A3','pDC'='#D93986'))


plotEnrichment(fgsea_sets[1],
               ranks)
library(ggalluvial)
ggplot(tmp, aes(x = Var1, 
               y=per*100,
               fill = Var2,alluvium = Var2)) +
geom_col()+
geom_flow(alpha = 1, width = 0)+
 scale_fill_manual(values=c('B cell'='#000000','T_naive'='#628DE1','Endothelial-Lyve1'='#979CE1','AT2'='#437193','IM'='#A9F1F9','neutrophils'='#8BC1CB',
 'NK'='#3780E9','Monocyte'='#254585','CD4T-Tnfrsf4'='#00007D','CD8T-Ccr9'='#D3B188','SPP1+ Macrophage'='#CDA100',
 'NKT'='#E96BAA','Fibroblast'='#B85880','AM'='#723955','cDC1'='#EEF37C','Ciliated'='#AA5828','gdT'='#1A421C','Endothelial-Kdr'='#86B400','cDC2'='#74876B',
 'Club'='#C40002','SMC'='#E28B60','Profilterating T'='#B4D9E5','AT1'='#466EE7','Profilterating SPP1+ Macrophage'='#6521AE','ILC2'='#2C3A74','mDC'='#54835A','Basophil'='#F08B64','Plasma'='#9D54B6',
 'Mesotheliocyte'='#7AC1A3','pDC'='#D93986'))+
 theme_bw()+
 ylab("Fraction of cells")+
 xlab('Group')

irGSEA.score(object = endo, assay = "RNA", slot = "data",
     custom = T, geneset = sig, method = c("AUCell"),
     kcdf = 'Gaussian')->endo