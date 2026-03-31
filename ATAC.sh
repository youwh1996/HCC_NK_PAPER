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












  labs(fill = "Sample", title = "Offset check")

sce.combined = spatialCluster(sce.combined, use.dimred = "HARMONY", q = 10, nrep = 10000) #use HARMONY
