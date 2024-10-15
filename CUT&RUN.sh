##比对
 module load bedtools/2.30.0
 module load bowtie2/2.4.5
 module load samtools/1.15.1
 cat file.txt | while read sample;
 do
cores=15
ref='/public/workspace/tmpuser/reference/bowtie2_hg38/hg38'
spikeInRef="/public/workspace/tmpuser/reference/Ecoli/Ecoil"
chromSize="/public/workspace/tmpuser/reference/hg38.chrom.sizes"
echo $sample
mkdir -p ./${sample}/align_bowtie2
bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${cores} -x ${ref} -1 ./${sample}/${sample}_R1.fq.gz -2 ./${sample}/${sample}_R2.fq.gz -S ./${sample}/align_bowtie2/${sample}_bowtie2.sam &> ./${sample}/align_bowtie2/${sample}_bowtie2.txt
bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${cores} -x ${spikeInRef} -1 ./${sample}/${sample}_R1.fq.gz -2 ./${sample}/${sample}_R2.fq.gz -S ./${sample}/align_bowtie2/${sample}_bowtie2_spikeIn.sam &> ./${sample}/align_bowtie2/${sample}_bowtie2_spikeIn.txt

seqDepthDouble=`samtools view -F 0x04  ./${sample}/align_bowtie2/${sample}_bowtie2_spikeIn.sam | wc -l`
seqDepth=$((seqDepthDouble/2))
echo $seqDepth > ./${sample}/align_bowtie2/${sample}_bowtie2_spikeIn.seqDepth

picardCMD="java -jar /public/workspace/tmpuser/picard.jar"
$picardCMD SortSam I=./${sample}/align_bowtie2/${sample}_bowtie2.sam O=./${sample}/align_bowtie2/${sample}_bowtie2.sorted.sam SORT_ORDER=coordinate
samtools view -bS -F 0x04 ./${sample}/align_bowtie2/${sample}_bowtie2.sorted.sam  > ./${sample}/align_bowtie2/${sample}_bowtie2.sorted.mapped.bam
mkdir -p ./${sample}/bed
bedtools bamtobed -i  ./${sample}/align_bowtie2/${sample}_bowtie2.sorted.mapped.bam -bedpe > ./${sample}/bed/${sample}_bowtie2.bed
awk '$1==$4 && $6-$2 < 1000 {print $0}' ./${sample}/bed/${sample}_bowtie2.bed > ./${sample}/bed/${sample}_bowtie2.clean.bed
cut -f 1,2,6 ./${sample}/bed/${sample}_bowtie2.clean.bed | sort -k1,1 -k2,2n -k3,3n  > ./${sample}/bed/${sample}_bowtie2.clean.fragments.bed
if [[ "$seqDepth" -gt "1" ]]; then
    
    mkdir -p ./${sample}/bedgraph

    scale_factor=`echo "10000 / $seqDepth" | bc -l`
   # echo "Scaling factor for Igg is: $scale_factor!"
    bedtools genomecov -bg -scale $scale_factor -i ./${sample}/bed/${sample}_bowtie2.clean.fragments.bed -g $chromSize >  ./${sample}/bedgraph/${sample}_bowtie2.clean.fragments.normalized.bedgraph
    
fi
samtools index ./${sample}/align_bowtie2/${sample}_bowtie2.sorted.mapped.bam 
mkdir -p bw                                                                                                       
bamCoverage -b ./${sample}/align_bowtie2/${sample}_bowtie2.sorted.mapped.bam --numberOfProcessors 10 -o ./bw/${sample}_raw.bw
done
##########
 cat file.txt | while read sample;
 do
 seqDepth=`cat ./${sample}/align_bowtie2/${sample}_bowtie2_spikeIn.seqDepth | bc -l`
 scale_factor=`echo "10000 / $seqDepth" | bc -l`
bamCoverage -b ./${sample}/align_bowtie2/${sample}_bowtie2.sorted.mapped.bam --numberOfProcessors 10  --scaleFactor  $scale_factor --binSize 10 --normalizeUsing CPM --effectiveGenomeSize 2913022398 -o ./bw/${sample}_normalized_raw.bw
done


######call peak
#sample=c('1344954','1365407')
seacr="/public/workspace/tmpuser/Bioapp/SEACR-master/SEACR_1.3.sh"
cat sample.txt | while read id
do
echo $id
mkdir ./${id}/SEACR -p
bash $seacr ./${id}/${id}RBPJCD49apos/bedgraph/${id}RBPJCD49apos_bowtie2.clean.fragments.normalized.bedgraph \
  ./${id}/${id}iggCD49apos/bedgraph/${id}iggCD49apos_bowtie2.clean.fragments.normalized.bedgraph \
   non stringent ./${id}/SEACR/${id}CD49apos_seacr_control.peaks
bash $seacr ./${id}/${id}RBPJCD49apos/bedgraph/${id}RBPJCD49apos_bowtie2.clean.fragments.normalized.bedgraph 0.01 non stringent \
        ./${id}/SEACR/${id}CD49apos_seacr_top0.01.peaks
bash $seacr ./${id}/${id}RBPJCD49aneg/bedgraph/${id}RBPJCD49aneg_bowtie2.clean.fragments.normalized.bedgraph \
  ./${id}/${id}iggCD49aneg/bedgraph/${id}iggCD49aneg_bowtie2.clean.fragments.normalized.bedgraph \
   non stringent ./${id}/SEACR/${id}CD49aneg_seacr_control.peaks
bash $seacr ./${id}/${id}RBPJCD49aneg/bedgraph/${id}RBPJCD49aneg_bowtie2.clean.fragments.normalized.bedgraph 0.01 non stringent \
       ./${id}/SEACR/${id}CD49aneg_seacr_top0.01.peaks
done
##macs2
#narrow peak calling
cat sample.txt | while read id
do
mkdir ./${id}/MACS2 -p
macs2 callpeak -t ./${id}/${id}RBPJCD49apos/align_bowtie2/${id}RBPJCD49apos_bowtie2.sorted.mapped.bam  \
               -g hs \
               -f BAMPE  \
               -n  ${id}_CD49apos \
               --outdir ./${id}/MACS2 -q 0.01 -B --SPMR --keep-dup all 
##broad peak calling
macs2 callpeak -t ./${id}/${id}RBPJCD49apos/align_bowtie2/${id}RBPJCD49apos_bowtie2.sorted.mapped.bam \
            -g hs \
            -f BAMPE \
            -n ${id}_CD49apos \
            --outdir ./${id}/MACS2 \
            --broad --broad-cutoff 0.1 -B --SPMR --keep-dup all
###Getting broad peak summits
python3 /public/workspace/tmpuser/CUT-RUNTools-2.0-master/install/get_summits_broadPeak.py ./${id}/MACS2/${id}_CD49apos_peaks.broadPeak | ~/bin/sort-bed - > ./${id}/MACS2/${id}_CD49apos_summits_broad.bed
#
macs2 callpeak -t ./${id}/${id}RBPJCD49aneg/align_bowtie2/${id}RBPJCD49aneg_bowtie2.sorted.mapped.bam  \
               -g hs \
               -f BAMPE  \
               -n  ${id}_CD49aneg \
               --outdir ./${id}/MACS2 -q 0.01 -B --SPMR --keep-dup all 
##broad peak calling
macs2 callpeak -t ./${id}/${id}RBPJCD49aneg/align_bowtie2/${id}RBPJCD49aneg_bowtie2.sorted.mapped.bam \
            -g hs \
            -f BAMPE \
            -n ${id}_CD49aneg \
            --outdir ./${id}/MACS2 \
            --broad --broad-cutoff 0.1 -B --SPMR --keep-dup all
###Getting broad peak summits
python3 /public/workspace/tmpuser/CUT-RUNTools-2.0-master/install/get_summits_broadPeak.py ./${id}/MACS2/${id}_CD49aneg_peaks.broadPeak | ~/bin/sort-bed - > ./${id}/MACS2/${id}_CD49aneg_summits_broad.bed
done
###########HOMER Call peak
cat sample.txt | while read id
do
mkdir ./${id}/HOMER -p
mkdir ./${id}/HOMER/Tag -p
mkdir ./${id}/HOMER/Tag/CD49aneg -p
mkdir ./${id}/HOMER/Tag/CD49apos -p
mkdir ./${id}/HOMER/Tag/CD49aneg_igg -p
mkdir ./${id}/HOMER/Tag/CD49apos_igg -p
makeTagDirectory  ./${id}/HOMER/Tag/CD49aneg/  ./${id}/${id}RBPJCD49aneg/align_bowtie2/${id}RBPJCD49aneg_bowtie2.sorted.mapped.bam
makeTagDirectory  ./${id}/HOMER/Tag/CD49apos/  ./${id}/${id}RBPJCD49apos/align_bowtie2/${id}RBPJCD49apos_bowtie2.sorted.mapped.bam
#makeTagDirectory  ./${id}/HOMER/Tag/CD49aneg_igg/  ./${id}/${id}iggCD49aneg/align_bowtie2/${id}iggCD49aneg_bowtie2.sorted.mapped.bam
#makeTagDirectory  ./${id}/HOMER/Tag/CD49apos_igg/  ./${id}/${id}iggCD49apos/align_bowtie2/${id}iggCD49apos_bowtie2.sorted.mapped.bam
findPeaks ./${id}/HOMER/Tag/CD49aneg/ -style factor -o auto  -tagThreshold 5 -LP 0.05 -L 0.00000001 -C 5 
findPeaks ./${id}/HOMER/Tag/CD49apos/ -style factor -o auto  -tagThreshold 5 -LP 0.05 -L 0.00000001 -C 5 
done
####pos to bed
cat sample.txt | while read id
do
pos2bed.pl ./${id}/HOMER/Tag/CD49aneg/peaks.txt > ./${id}/HOMER/Tag/CD49aneg/peaks.bed
pos2bed.pl ./${id}/HOMER/Tag/CD49apos/peaks.txt > ./${id}/HOMER/Tag/CD49apos/peaks.bed
done

mergePeaks ./1162472/HOMER/Tag/CD49aneg/peaks.txt ./1162472/HOMER/Tag/CD49apos/peaks.txt  \
./1344954/HOMER/Tag/CD49aneg/peaks.txt ./1344954/HOMER/Tag/CD49apos/peaks.txt \
./1345016/HOMER/Tag/CD49aneg/peaks.txt ./1345016/HOMER/Tag/CD49apos/peaks.txt \
./1365407/HOMER/Tag/CD49aneg/peaks.txt ./1365407/HOMER/Tag/CD49apos/peaks.txt \
./1383848/HOMER/Tag/CD49aneg/peaks.txt ./1383848/HOMER/Tag/CD49apos/peaks.txt > merge_peak.txt

pos2bed.pl merge_peak.txt > merge_peak.bed
 awk '{print "peak_"NR"\t"$1"\t"$2"\t"$3"\t""."}' merge_peak.bed > merge.bed.saf

~/subread-2.0.6-Linux-x86_64/bin/featureCounts -T 4 -a  merge.bed.saf -F SAF -p -o counts_subread.txt \
./1162472/1162472RBPJCD49apos/align_bowtie2/1162472RBPJCD49apos_bowtie2.sorted.mapped.bam ./1162472/1162472RBPJCD49aneg/align_bowtie2/1162472RBPJCD49aneg_bowtie2.sorted.mapped.bam \
./1365407/1365407RBPJCD49apos/align_bowtie2/1365407RBPJCD49apos_bowtie2.sorted.mapped.bam ./1365407/1365407RBPJCD49aneg/align_bowtie2/1365407RBPJCD49aneg_bowtie2.sorted.mapped.bam \
./1383848/1383848RBPJCD49apos/align_bowtie2/1383848RBPJCD49apos_bowtie2.sorted.mapped.bam ./1383848/1383848RBPJCD49aneg/align_bowtie2/1383848RBPJCD49aneg_bowtie2.sorted.mapped.bam \
./1344954/1344954RBPJCD49apos/align_bowtie2/1344954RBPJCD49apos_bowtie2.sorted.mapped.bam ./1344954/1344954RBPJCD49aneg/align_bowtie2/1344954RBPJCD49aneg_bowtie2.sorted.mapped.bam \
./1345016/1345016RBPJCD49apos/align_bowtie2/1345016RBPJCD49apos_bowtie2.sorted.mapped.bam ./1345016/1345016RBPJCD49aneg/align_bowtie2/1345016RBPJCD49aneg_bowtie2.sorted.mapped.bam 

##### heatmap 
cat sample.txt | while read id
do
#CD49a pos
#computeMatrix reference-point -S ./${id}/bw/${id}RBPJCD49apos_raw.bw -R ./${id}/HOMER/Tag/CD49apos/peaks.bed  \
 #             --skipZeros -o ./${id}/${id}CD49apos_HOMER.mat.gz -p $cores -a 3000 -b 3000 --referencePoint center 
plotHeatmap -m  ./${id}/${id}CD49apos_HOMER.mat.gz -out ./${id}/${id}_CD49apos_HOMER_heatmap.pdf  --sortUsing sum --startLabel "Peak Start" -\
-endLabel "Peak End" --xAxisLabel "" --regionsLabel "Peaks" --samplesLabel "${id}_CD49apos" --colorList 'white,#339933' 
#CD49a neg 
#computeMatrix reference-point -S  ./${id}/bw/${id}RBPJCD49aneg_raw.bw -R ./${id}/HOMER/Tag/CD49aneg/peaks.bed  \
#             --skipZeros -o ./${id}/${id}CD49aneg_HOMER.mat.gz -p $cores -a 3000 -b 3000 --referencePoint center 
plotHeatmap -m  ./${id}/${id}CD49aneg_HOMER.mat.gz -out ./${id}/${id}_CD49aneg_HOMER_heatmap.pdf  --sortUsing sum --startLabel "Peak Start" -\
-endLabel "Peak End" --xAxisLabel "" --regionsLabel "Peaks" --samplesLabel "${id}_CD49aneg"   --colorList 'white,#339933'  
done




###
getDifferentialPeaksReplicates.pl -t ./1365407/HOMER/Tag/CD49apos ./1162472/HOMER/Tag/CD49apos ./1383848/HOMER/Tag/CD49apos -b ./1365407/HOMER/Tag/CD49aneg ./1162472/HOMER/Tag/CD49aneg ./1383848/HOMER/Tag/CD49aneg > posvsneg.txt 
###computeMatrix heatmap
cat sample.txt | while read id
do
cores=15
computeMatrix  reference-point -S  ./${id}/bw/${id}RBPJCD49apos_raw.bw \
                                ./${id}/bw/${id}iggCD49apos_raw.bw \
                                ./${id}/bw/${id}RBPJCD49aneg_raw.bw \
                                ./${id}/bw/${id}iggCD49aneg_raw.bw \
                            -R /public/workspace/tmpuser/Bioapp/hg38_gene/hg38_gene.bed --beforeRegionStartLength 3000 \
                            --afterRegionStartLength 3000 \
                            --binSize 50 \
                            --referencePoint TSS  -o ./${id}/${id}_allrefPoint_matrix_gene.mat.gz -p $cores
plotHeatmap -m ./${id}/${id}_allrefPoint_matrix_gene.mat.gz -out ./${id}/${id}_all_heatmap.pdf --sortUsing sum --missingDataColor 1  --colorList 'white,#339933'  --heatmapHeight 12 
done

####ON peaks heatmap seacr
cat sample.txt | while read id
do
cores=15
echo $id
awk '{split($6, summit, ":"); split(summit[2], region, "-"); print summit[1]"\t"region[1]"\t"region[2]}'  ./${id}/SEACR/${id}CD49apos_seacr_control.peaks.stringent.bed >  ./${id}/SEACR/${id}CD49apos_seacr_control.peaks.summitRegion.bed
awk '{split($6, summit, ":"); split(summit[2], region, "-"); print summit[1]"\t"region[1]"\t"region[2]}'  ./${id}/SEACR/${id}CD49aneg_seacr_control.peaks.stringent.bed >  ./${id}/SEACR/${id}CD49aneg_seacr_control.peaks.summitRegion.bed
computeMatrix reference-point -S  ./${id}/bw/${id}RBPJCD49apos_raw.bw \                  
              -R  ./${id}/SEACR/${id}CD49apos_seacr_control.peaks.summitRegion.bed  \
              --skipZeros -o ./${id}/${id}CD49apos_SEACR.mat.gz -p $cores -a 3000 -b 3000 --referencePoint center 
computeMatrix reference-point -S  ./${id}/bw/${id}RBPJCD49aneg_raw.bw \
              -R  ./${id}/SEACR/${id}CD49aneg_seacr_control.peaks.summitRegion.bed  \
              --skipZeros -o ./${id}/${id}CD49aneg_SEACR.mat.gz -p $cores -a 3000 -b 3000 --referencePoint center
plotHeatmap -m  ./${id}/${id}CD49apos_SEACR.mat.gz -out ./${id}/${id}_CD49apos_SEACR_heatmap.pdf  --sortUsing sum --startLabel "Peak Start" -\
-endLabel "Peak End" --xAxisLabel "" --regionsLabel "Peaks" --samplesLabel "${id}_CD49apos" --colorMap RdYlBu  

plotHeatmap -m  ./${id}/${id}CD49aneg_SEACR.mat.gz -out ./${id}/${id}_CD49aneg_SEACR_heatmap.pdf  --sortUsing sum --startLabel "Peak Start" -\
-endLabel "Peak End" --xAxisLabel "" --regionsLabel "Peaks" --samplesLabel "${id}_CD49aneg" --colorMap RdYlBu  
done
######MACS2
cat sample.txt | while read id
do
cores=15
echo $id
computeMatrix reference-point -S  ./${id}/bw/${id}RBPJCD49apos_raw.bw \
              -R  ./${id}/MACS2/${id}_CD49apos_summits.bed  \
              --skipZeros -o ./${id}/${id}CD49apos_MACS2.mat.gz -p $cores -a 3000 -b 3000 --referencePoint center
computeMatrix reference-point -S  ./${id}/bw/${id}RBPJCD49aneg_raw.bw \
              -R  ./${id}/MACS2/${id}_CD49aneg_summits.bed  \
              --skipZeros -o ./${id}/${id}CD49aneg_MACS2.mat.gz -p $cores -a 3000 -b 3000 --referencePoint center
plotHeatmap -m  ./${id}/${id}CD49apos_MACS2.mat.gz -out ./${id}/${id}_CD49apos_MACS2_heatmap.pdf  --sortUsing sum --startLabel "Peak Start" -\
-endLabel "Peak End" --xAxisLabel "" --regionsLabel "Peaks" --samplesLabel "${id}_CD49apos"

plotHeatmap -m ./${id}/${id}CD49aneg_MACS2.mat.gz -out ./${id}/${id}_CD49aneg_MACS2_heatmap.pdf  --sortUsing sum --startLabel "Peak Start" -\
-endLabel "Peak End" --xAxisLabel "" --regionsLabel "Peaks" --samplesLabel "${id}_CD49aneg"
done
##



#######################R analysis
##Create the peak x sample matrix.
library(dplyr)
library(stringr)
library(ggplot2)
library(viridis)
library(GenomicRanges)
library(chromVAR) ## For FRiP analysis and differential analysis
library(DESeq2) ## For differential analysis section
library(ggpubr) ## For customizing figures
library(corrplot) #
library(ChIPseeker)
mPeak = GRanges()
for(hist in c('1162472','1365407','1383848','1344954')){
    for(rep in c('CD49apos','CD49aneg')){
        peakRes = read.table(paste0('./',hist, "/HOMER/", '/Tag/', rep, "/peaks.bed"), header = FALSE, fill = TRUE)
        mPeak = GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*") %>% append(mPeak, .)
    }
}
masterPeak = reduce(mPeak)
library(DESeq2)
countMat = matrix(NA, length(masterPeak), length(c('1162472','1365407','1383848','1344954'))*length(c('CD49apos','CD49aneg')))
i = 1
for(hist in c('1162472','1365407','1383848','1344954')){
  for(rep in  c('RBPJCD49apos','RBPJCD49aneg')){
    bamFile = paste0( "./", hist, '/',hist,rep,"/align_bowtie2/", hist,rep, "_bowtie2.sorted.mapped.bam")
    fragment_counts <- getCounts(bamFile, masterPeak, paired = TRUE, by_rg = FALSE, format = "bam")
    countMat[, i] = counts(fragment_counts)[,1]
    i = i + 1
  }
}
colnames(countMat) = paste(c('1162472','1162472','1365407','1365407','1383848','1383848','1344954','1344954'),rep(c('RBPJCD49apos','RBPJCD49aneg'),2), sep = "_")
cbind(countMat,data.frame(masterPeak))->countMat1
rownames(countMat1)=paste0('peak_',1:nrow(countMat1))

###peak annotation
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(clusterProfiler)
 readPeakFile("./merge_peak.bed")->peaK
 peakAnno <- annotatePeak(peaK, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
data.frame(peakAnno@anno)->anno
anno$pos=paste0(anno$seqnames,'_',anno$start,'_',anno$end)
countMat1$pos=paste0(countMat1$Chr,'_',countMat1$Start+1,'_',countMat1$End)
countMat1$SYMBOL=anno$SYMBOL[match(countMat1$pos,anno$pos)]
countMat1$annotation=anno$annotation[match(countMat1$pos,anno$pos)]

###
library(DESeq2)
row.names(countMat1)=countMat1$Geneid
dataS = countMat1[,c(7:16)]
dataS=dataS[,c(1,2,5,6,7,8)]
selectR=rowSums(dataS)>5
 dataS[selectR,]->dataS
condition = factor(rep(c('RBPJCD49apos','RBPJCD49aneg'),3))

dds<-DESeqDataSetFromMatrix(countData=dataS,colData=data.frame(condition),design=~condition)
dds<- estimateSizeFactors(dds)
dds <- DESeq(dds)
res= results(dds)
cbind(res,countMat1[,c(2:6,18,19)])->res
res[!is.na(res$pvalue),]->res

res1=res[which(res$annotation %in% c("Promoter (1-2kb)","Promoter (<=1kb)","Promoter (2-3kb)")),]
diff_genes <- rownames(diff)[abs(diff$log2FC) > log2(1.5) & diff$p.value <0.05]
diff$sig <- ""
diff[diff_genes,]$sig <- diff_genes
diff <- diff[order(diff$sig),]

ggplot(diff,aes(x=log2FC,y=-log10(p.value)))+
geom_point(color = ifelse(diff$sig == "", "grey", "red"))+
geom_vline(xintercept = log2(1.5) , linetype="dashed", color="grey77") + 
geom_vline(xintercept = -log2(1.5), linetype="dashed", color="grey77") +geom_hline(yintercept = -log10(0.05), linetype="dashed", color="grey77") +
 ylab('-log10(pvalue)') + xlab('log2 f')+
xlim(-5,5)+
theme_classic()

###
library(edgeR)
dataS = countMat1[,c(7:14)]

dgelist <- DGEList(counts = dataS, group = condition)
keep <- rowSums(cpm(dgelist) > 1 ) >= 2
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
design <- model.matrix(~condition)
dge <- estimateDisp(dgelist_norm, design, robust = TRUE)
fit <- glmFit(dge, design, robust = TRUE)
lrt <- topTags(glmLRT(fit), n = nrow(dgelist$counts))
lrt<-lrt$table

lrt1=lrt[which(lrt$annotation %in% c("Promoter (1-2kb)","Promoter (<=1kb)","Promoter (2-3kb)")),]













library(BioSeqUtils)
library(ggplot2)
file <- list.files(path = "./",pattern = '.bw',full.names = T)
bw <- loadBigWig(file)
gtf <- rtracklayer::import.gff( '~/Bioapp/refdata-gex-GRCh38-2020-A/genes/genes.gtf',format = "gtf") %>% data.frame() 
trackVisProMax(Input_gtf = gtf,Input_bw = bw,Input_gene = c('GZMB','GZMK','GZMA','CXCR3','CCL5','TNFRSF18'))

computeMatrix scale-regions -S  ./bw/${id}RBPJCD49apos_raw.bw \
                                ./bw/${id}iggCD49apos_raw.bw \
                                ./bw/${id}RBPJCD49aneg_raw.bw \
                                ./bw/${id}iggCD49aneg_raw.bw \
                            -R /public/workspace/tmpuser/Bioapp/hg38_gene/hg38_gene.bed \
                            --beforeRegionStartLength 3000 \
                            --afterRegionStartLength 3000 \
                            --regionBodyLength 5000 \
                            --skipZeros  -o ./${id}_reigion_matrix_gene.mat.gz -p 14
plotHeatmap -m ./${id}_reigion_matrix_gene.mat.gz -out ./${id}_all_heatmap.pdf --sortUsing sum --missingDataColor 1 --colorMap RdYlBu   --heatmapHeight 12 


cwl-runner --parallel --outdir  /home/youwh/Project/Project_EE_immune/HCC_BD/SE/ /home/youwh/reference/v2.0/rhapsody_pipeline_2.0.cwl /home/youwh/reference/pipeline.yml
