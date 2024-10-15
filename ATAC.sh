bowtie2 -p 10 -X 2000 --very-sensitive   \
-x /home/youwh/reference/bowtie2_mm10/mm10 \
-1 HFHC_B_SE.paired.R1.fq.gz -2 HFHC_B_SE.paired.R2.fq.gz  | samtools sort -@ 10 -O bam -o -> SE.bam
samtools index SE.bam
#Remove mitochondrial reads
samtools idxstats SE.bam > SE.idxstats
grep "chrM" SE.idxstats
samtools flagstat SE.bam  > SE.flagstat

samtools view -h SE.bam  | grep -v chrM | samtools sort -O bam -o SE.rmChrM.bam -T .
#Remove duplicates & low-quality alignments
java -jar ~/Bioapp/picard.jar MarkDuplicates QUIET=true INPUT=SE.rmChrM.bam OUTPUT=SE.rmChrM.marked.bam METRICS_FILE=SE.dup.metrics REMOVE_DUPLICATES=false CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=.
head -n 8 SE.dup.metrics | cut -f 7,9 | grep -v ^# | tail -n 2
##Remove multi-mapped reads 
samtools view -h -q 30 SE.rmChrM.marked.bam  > SE.rmChrM.marked.rmMulti.bam
samtools sort -O bam -o SE.rmChrM.marked.rmMulti.sort.bam SE.rmChrM.marked.rmMulti.bam
###Remove ENCODE blacklist regions
bedtools intersect -nonamecheck -v -abam SE.rmChrM.marked.rmMulti.sort.bam -b ~/mm10-blacklist.v2.bed >  SE.rmChrM.marked.rmMulti.encode.bam
samtools sort -O bam -o SE.rmChrM.marked.rmMulti.encode_filtered.bam SE.rmChrM.marked.rmMulti.encode.bam
samtools index SE.rmChrM.marked.rmMulti.encode_filtered.bam 
rm  SE.rmChrM.marked.rmMulti.encode.bam
##Shifting reads
alignmentSieve --numberOfProcessors 8 --ATACshift --bam  SE.rmChrM.marked.rmMulti.encode_filtered.bam  -o  SE.rmChrM.marked.rmMulti.encode_filtered.tmp.bam 
samtools sort -@ 8 -O bam -o SE.rmChrM.marked.rmMulti.encode_filtered.shift.bam  SE.rmChrM.marked.rmMulti.encode_filtered.tmp.bam 
samtools index -@ 8 SE.rmChrM.marked.rmMulti.encode_filtered.shift.bam
rm  SE.rmChrM.marked.rmMulti.encode_filtered.tmp.bam 
##macs2 Peak calling 
macs2 callpeak -f BAMPE -g mm --keep-dup all --cutoff-analysis -n SE \
  -t SE.rmChrM.marked.rmMulti.encode_filtered.shift.bam --outdir ./SE
bamCoverage --numberOfProcessors 8 --binSize 10 --normalizeUsing RPGC \
  --effectiveGenomeSize 2652783500 --bam EE.rmChrM.marked.rmMulti.encode_filtered.shift.bam -o EE.shifted.bw
###homer peak calling


##
computeMatrix reference-point --referencePoint center  -S  EE.shifted.bw  -R ./EE/EE_summits.bed \
                          --skipZeros  -b 1000 -a 1000 -bs 1  -o ./matrix_gene.mat.gz -p 10
plotHeatmap -m ./matrix_gene.mat.gz -out ./EE.pdf --heatmapHeight 7 --heatmapWidth 11 --colorList 'white,red' --plotType lines --perGroup





dat<-list()
for (i in 1:length(mysample)){
  # i=1
  dir=paste0('./',mysample[i])
  test=Read10X(data.dir = dir)
  colnames(test)=paste0('GSE176021_',mysample[i],'_',colnames(test))
  if('Pigr'%in%rownames(test)){
    cat('Now is: ', dir,'\n')
      print(dim(test))
    cat('Cell location is: ','\n')
    print(which(test['Pigr',]!=0)) #Pigr出现在第几个细胞
    cat('Expression is: ','\n')
    print(test['Pigr',test['Pigr',]!=0]) # Pigr表达量
    cat('\n\n')
  }
  CreateSeuratObject(test)->test
  test$PID=paste0('GSE176021',"_",mysample[i])
  test$dataset='GSE176021'
  test$Tissue='primary'
  test$Disease='LUAD'
  dat[[mysample[i]]]=test
}

 paste0(unlist(lapply(myfile,function(x){strsplit(x,'_')[[1]][2]})),'_', unlist(lapply(myfile,function(x){strsplit(x,'_')[[1]][3]})),'_',unlist(lapply(myfile,function(x){strsplit(x,'_')[[1]][4]})))


for(i in 1: length(myfile)){
  #创建每个样本自己的文件夹:
  dir.create(paste0('./',mysample[i]))
  #把文件依次移动到新的文件夹中并命名:
  file.copy(paste0('./',myfile[i]),
        paste0('./',mysample[i],'/',fileformat[i]))
  #删除原文件:      
  unlink(paste0('./',myfile[i]))
}
dat<-list()


for(i in 1:length(aa$SID)){
  dat[[aa$SID[i]]]$PID=paste0('GSE183904_',aa$PID[i])
   dat[[aa$SID[i]]]$dataset='GSE183904'
   dat[[aa$SID[i]]]$Tissue=aa$Tissue[i]
   dat[[aa$SID[i]]]$Disease='GC'
}



cluster1 <- ConsensusClusterPlus::ConsensusClusterPlus(
    d = as.matrix(mat), #数据矩阵，列是样本，行是变量
    maxK = 5, #要评估的最大聚类簇数量
    seed = 1234, #指定随机数种子，用于子样本处理等随机的过程
    reps = 200, #抽取的子样本数量
    pItem = 0.8, #抽样样本的比例
    pFeature = 1, #抽样变量的比例
    clusterAlg = 'km', #选择聚类算法
   distance = 'euclidean', #指定聚类时使用的距离或相关性类型
    title = 'test3', plot = 'pdf' #聚类簇评估结果的输出格式
)
 

Tumor <- NormalizeData(Tumor, normalization.method = "LogNormalize", scale.factor = 10000)
Tumor <- FindVariableFeatures(Tumor, selection.method = "vst", nfeatures = 3000)
Tumor <- ScaleData(Tumor, features = rownames(Tumor))
Tumor <- RunPCA(Tumor, features = VariableFeatures(object = Tumor))
Tumor <- FindNeighbors(Tumor, dims = 1:30,reduction = "pca")
Tumor <- RunUMAP(Tumor, dims = 1:30,reduction = "pca")


fviz_pca_ind(pre.pca,
             geom= "point",
             col.ind = as.factor(batchType),
             addEllipses = TRUE,
             legend.title="Group")

fviz_pca_ind(pre.pca,
             geom= "point",
             col.ind = as.factor(batch1),
             addEllipses = TRUE,
             legend.title="Group")


ggplot(data.frame(reducedDim(sce.combined, "UMAP.HARMONY")), 
       aes(x = UMAP1, y = UMAP2, color = factor(sce.combined$spatial.cluster))) +
  geom_point() +
  labs(color = "Sample") +
  theme_bw()



seurat_to_spe <- function(seu, sample_id, img_id) {
    ## Convert to SCE
    sce <- Seurat::as.SingleCellExperiment(seu)
    colnames(sce)=paste0(sample_id,'_',colnames(sce))
    ## Extract spatial coordinates
    spatialCoords <- as.matrix(
        seu@images[[img_id]]@coordinates[, c("col", "row")])
    
    ## Extract and process image data
    img <- SpatialExperiment::SpatialImage(
        x = as.raster(seu@images[[img_id]]@image))
    
    imgData <- DataFrame(
        sample_id = sample_id,
        image_id = img_id,
        data = I(list(img)),
        scaleFactor = seu@images[[img_id]]@scale.factors$lowres)
    
    # Convert to SpatialExperiment
    spe <- SpatialExperiment(
        assays = assays(sce),
        rowData = rowData(sce),
        colData = colData(sce),
        metadata = metadata(sce),
        reducedDims = reducedDims(sce),
        altExps = altExps(sce),
        sample_id = sample_id,
        spatialCoords = spatialCoords,
        imgData = imgData
    )
    # indicate all spots are on the tissue
    spe$in_tissue <- 1
    spe$sample_id <- sample_id
    # Return Spatial Experiment object
    spe
}















sce.combined$row[sce.combined$sample_name == "P11T"] = 
  100 + sce.combined$row[sce.combined$sample_name == "P11T"]
sce.combined$col[sce.combined$sample_name == "P3T"] = 
  150 + sce.combined$col[sce.combined$sample_name == "P3T"]
sce.combined$row[sce.combined$sample_name == "P5T"] = 
  100 + sce.combined$row[sce.combined$sample_name == "P5T"]
sce.combined$col[sce.combined$sample_name == "P7T"] = 
  150 + sce.combined$col[sce.combined$sample_name == "P7T"]
sce.combined$row[sce.combined$sample_name == "P8T"] = 
  100 + sce.combined$row[sce.combined$sample_name == "P8T"]
  sce.combined$col[sce.combined$sample_name == "P9T"] = 
  150 + sce.combined$col[sce.combined$sample_name == "P9T"]
    sce.combined$row[sce.combined$sample_name == "P9T"] = 
  50 + sce.combined$row[sce.combined$sample_name == "P9T"]
clusterPlot(sce.combined, "sample_name", color = NA) + #make sure no overlap between samples
  labs(fill = "Sample", title = "Offset check")

sce.combined = spatialCluster(sce.combined, use.dimred = "HARMONY", q = 10, nrep = 10000) #use HARMONY
