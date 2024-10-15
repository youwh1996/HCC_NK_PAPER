library(foreach)
library(Matrix)
library(ggrepel)
library(cowplot)
library(doParallel)
library(plyr)
library(MASS)
theme_set(theme_cowplot())
###B cell(including plasma B cell)
##comparing INV vs pre-INV
###
registerDoParallel(cores = 15)
getDoParWorkers()
##
Bcellid=colnames(LUAD)[which(LUAD$type_global %in% c('B cell','Plasma B'))]
Tcellid=colnames(LUAD)[which(LUAD$type_global %in% c('CD4+ T','CD8+ T','Cycling T'))]
###
mat=scdb_mat("submc_Mac")
annot_Mac=mat@cell_metadata[intersect(colnames(mat@mat),names(mc2d@sc_x)),c('Time','Tissue_mod_only','Group')]
annot_Mac$type=mc@mc[rownames(annot_Mac)]
annot_Mac$type=factor(annot_Mac$type)
levels(annot_Mac$type)=annot$type
annot_Mac$Well_ID=rownames(annot_Mac)
<-annot_Mac
cells_Str$pop='Mf'

anno_Tcell=mat@cell_metadata[Tcellid,c('PatientID','stage')]
anno_Tcell$type=LUAD@meta.data[Tcellid,'type']
anno_Tcell$Well_ID=rownames(anno_Tcell)
cells_Str<-anno_Tcell
cells_Str$pop='Tcell'
colnames(cells_Str)[3]='Disease'
colnames(cells_Str)[1]='PID'



anno_Mac=mat@cell_metadata[Macrophage_id,c('ident','tissue','stage')]
anno_Mac$type=ICC@meta.data[Macrophage_id,'type']
anno_Mac$Well_ID=rownames(anno_Mac)
cells_Str<-anno_Mac
cells_Str$pop='Macrophage'
colnames(cells_Str)[3]='Disease'
colnames(cells_Str)[1]='PID'

anno_Bcell=mat@cell_metadata[Bcellid,c('PatientID','stage')]
anno_Bcell$type=LUAD@meta.data[Bcellid,'type']
anno_Bcell$Well_ID=rownames(anno_Bcell)
cells_Str<-anno_Bcell
cells_Str$pop='Bcell'
colnames(cells_Str)[2]='Disease'
colnames(cells_Str)[1]='PID'
#levels(cells_Str$Disease)=c('INV','pre-INV','pre-INV','pre-INV')
mat=scdb_mat('LUAD')
mat_ds <- mat@mat[,cells_Str$Well_ID]

# Calculate diff genes between hip control and individual clones 
n <- 1000 # number of cells for each sampling
ntime <- 150 # number of times for sampling
table(cells_Str$pop)
pop <- 'Bcell'
diseaseA <- 'INV'; cutA <- 15
diseaseB <- 'adjacent'; cutB <- 15
###
pop <- 'Tcell'
diseaseA <- 'CD4+ Tfh'; cutA <- 15
diseaseB <- 'other'; cutB <- 15


diseaseA <- 'SPF'; cutA <- 15
diseaseB <- 'GF'; cutB <- 15
########## Sample hip cells ########## This is the only function you need to change for your purpose, basically is how to sample cells.
get_n_cells_per_pop_per_disease <- function(cells_Str, n, pop, disease, cut = 15, diag = F, mc = F){
  cells_Now <- cells_Str[cells_Str$pop %in% pop & cells_Str$Disease %in% disease,]
  message(paste(pop, disease, nrow(cells_Now)))
  inf <- table(cells_Now$PID)
  inf <- sort(inf[inf > 0], decreasing = T)
  message('All: \n', paste(inf, collapse = '\t'))
  inf <- inf[inf >= cut]
  # message(length(inf))
  
  n_Ind <- min(cut, ceiling(n/length(inf)))
  message(n_Ind)
  cells_pool <- foreach(i = sort(names(inf)), .combine = c) %do% sample(rownames(cells_Now)[cells_Now$PID %in% i], n_Ind)
  cells_pool <- rep(cells_pool, ceiling(n/length(cells_pool)))
  res <- sample(cells_pool, n)
  if(diag){
    message('Sel: \n', paste(inf, collapse = '\t'))
    message(n_Ind)
    message('N:', length(inf))
    return(length(inf))
  }else if(mc){
    cells_Now <- cells_Now[cells_Now$PID %in% names(inf),]
    res <- structure(as.vector(cells_Now$PID), names = rownames(cells_Now))
    return(res)
  }else{
    return(res)
  }
}
##
get_n_cells_per_pop_per_disease(cells_Str, n, pop, diseaseA, cutA, diag =F)
get_n_cells_per_pop_per_disease(cells_Str, n, pop, diseaseB, cutB, diag =F)
########## Calc of diff ##########
diff_calc <- function(cells_Str, mat_ds, pop, n, ntime, diseaseA, diseaseB, cutA, cutB){
  res <- foreach(i = seq(ntime), .combine = rbind) %dopar% {
    grpA <- get_n_cells_per_pop_per_disease(cells_Str, n, pop, diseaseA, cutA)
    grpB <- get_n_cells_per_pop_per_disease(cells_Str, n, pop, diseaseB, cutB)
    cal_exp_diff(grpA, grpB, mat_ds, i)
  }
}
############ Calculate diff genes between two sets of cells ##########
cal_exp_diff <- function(grpA, grpB, mat_ds, ntime, blck = NULL){
  # wilcox.test
  # kruskal.test
  if(is.null(blck)){
   #blck <- c("^MT-", "^MTMR", "^MTND", "^MTRN", "^MTCO", "^MTATP", "^MRPS","^MRPL",
   #         "^AC[0-9]", "^AL[0-9]", "^AP[0-9]", "^Gm[0-9]", "^MIR", "^SNOR", "^ATP", "^FP[0-9]", "^FO[0-9]","^RPS","^RPL","^HBA","^HBB",
    #          "^IGJ", "^IGH", "^IGK", "^IGL", "^DNAJ", "^IGD")
    # blck <- c("^MT-", "^MTMR", "^MTND", "^MTRN", "^MTCO", "^MTATP", "^MRPS","^MRPL",
    #        "^AC[0-9]", "^AL[0-9]", "^AP[0-9]", "^Gm[0-9]", "^MIR", "^SNOR", "^ATP", "^FP[0-9]", "^FO[0-9]","^RPS","^RPL","^HBA","^HBB"
   #        )
 # }
  if(is.null(blck)){
    blck <- c("Rik[0-9]?$", "^AC[0-9]+\\.", "^AL[0-9]+\\.", "^Gm[0-9]+", "^Snor", "-ps[0-9]?$", "-rs[0-9]?$", "\\.", "-",   
                 "^Dnaj", "^Hist", '^mt-',"^Mir", "^Hsp", "^Atp", "^Uqc",
                 "^A[A-Z]","^B[A-Z]","^D[0-9]")
  }
  # "^Igh", "^Igj", "^Igk", "^Igl"
  message(paste("grpA", length(grpA)))
  message(paste("grpB", length(grpB)))
  data <- mat_ds
  gene <- rownames(data)
  gene_bad <- foreach(i = blck, .combine = c) %do% grep(i, gene, v = T)
  #gene_bad<-c("Sparc",'Col1a1','Col1a2')
  gene <- setdiff(gene, gene_bad)
  message(paste("gene", length(gene)))
  A <- data[gene, grpA]
  B <- data[gene, grpB]
  AB <- cbind(A, B)
  pctA <- Matrix::rowSums(A > 0)/ncol(A)
  pctB <- Matrix::rowSums(B > 0)/ncol(B)
  avgA <- Matrix::rowSums(A)/sum(Matrix::rowSums(A)) * 1e3
  avgB <- Matrix::rowSums(B)/sum(Matrix::rowSums(B)) * 1e3
  wilcox.p <- apply(AB, 1, function(x) wilcox.test(x[1:ncol(A)], x[(ncol(A) + 1):ncol(AB)])[[3]])
  wilcox.p[is.na(wilcox.p)] <- 1
  zstat <- abs(qnorm(wilcox.p/2)) * sign(avgB - avgA)
  wilcox.padj <- p.adjust(wilcox.p, 'BH')
  wilcox.padj[is.na(wilcox.padj)] <- 1
  zstatadj <- abs(qnorm(wilcox.padj/2)) * sign(avgB - avgA)
  ori <- as.data.frame(cbind(gene, pctA, pctB, avgA, avgB, zstat, zstatadj, wilcox.p, wilcox.padj))
  ori$ntime <- ntime
  ori <- ori[order(ori$wilcox.padj),]
  message(paste(dim(ori), collapse = '\t'))
  return(ori)
}
#
diff_stat <- function(diff_res, ntime){
  genes_freq <- table(diff_res$gene)
  genes <- names(genes_freq)
  message(paste("final gene counts: ", length(genes)))
  res <- foreach(i = genes, .combine = rbind) %dopar% {
    now <- diff_res[diff_res$gene == i,]
    pctA.mean <- mean(as.numeric(as.vector(now$pctA)))
    pctB.mean <- mean(as.numeric(as.vector(now$pctB)))
    avgA.mean <- mean(as.numeric(as.vector(now$avgA)))
    avgB.mean <- mean(as.numeric(as.vector(now$avgB)))
    zstat.mean <- mean(as.numeric(as.vector(now$zstat)))
    wilcox.p.mean <- pnorm(abs(zstat.mean), lower.tail = F) * 2
    pctA.sd <- sd(as.numeric(as.vector(now$pctA)))
    pctB.sd <- sd(as.numeric(as.vector(now$pctB)))
    avgA.sd <- sd(as.numeric(as.vector(now$avgA)))
    avgB.sd <- sd(as.numeric(as.vector(now$avgB)))
    zstat.sd <- sd(as.numeric(as.vector(now$zstat)))
    res <- data.frame(pctA.mean, pctB.mean, avgA.mean, avgB.mean, zstat.mean, wilcox.p.mean,
                      pctA.sd, pctB.sd, avgA.sd, avgB.sd, zstat.sd)
    rownames(res) <- i
    res
  }
  res$wilcox.padj.mean <- p.adjust(res$wilcox.p.mean, method = 'BH')
  res <- res[order(res$wilcox.p.mean, decreasing = F),]
}


##
ind_diff <- diff_calc(cells_Str, mat_ds, pop, n, ntime, diseaseA, diseaseB, cutA, cutB)
ind_diff_stat <- diff_stat(ind_diff, ntime)

deg<-list()
n <- 200 # number of cells for each sampling
ntime <- 500 # number of times for sampling

pop='T'
for(i in unique(sc_info$Disease)){
  diseaseB=i
  diseaseA=setdiff(unique(sc_info$Disease),diseaseB)
  cutA=15
  cutB=15
  ind_diff <- diff_calc(cells_Str, mat_ds, pop, n, ntime, diseaseA, diseaseB, cutA, cutB)
  deg[[i]]= diff_stat(ind_diff, ntime)
}

for(i in 1:5){
deg[[i]]$f=log2((deg[[i]]$avgB.mean + 0.01)/(deg[[i]]$avgA.mean + 0.01))
deg[[i]]$q=log(-(log(deg[[i]]$wilcox.p.mean, 10)) + 1)
deg[[i]]$gene=rownames(deg[[i]])
deg[[i]]$pop=names(deg[i])
}

DEG<-cbind(deg[[1]][,c('f','q','gene','pop')],deg[[2]][,c('f','q','gene','pop')],deg[[3]][,c('f','q','gene','pop')],
deg[[4]][,c('f','q','gene','pop')],deg[[5]][,c('f','q','gene','pop')])





##
########## plot volcano
f_cut <- 1.5
q_cut <- 0.05
pct_cut <- 0.05###
ind_diff_stat_good <- ind_diff_stat[ind_diff_stat$pctA.mean > pct_cut | ind_diff_stat$pctB.mean > pct_cut,]
f <- log2((ind_diff_stat_good$avgB.mean + 0.01)/(ind_diff_stat_good$avgA.mean + 0.01))
q <- log(-(log(ind_diff_stat_good$wilcox.p.mean, 10)) + 1)
dat <- data.frame(g = rownames(ind_diff_stat_good), f = f, q = q)
up <- dat[dat$f > log2(f_cut) & dat$q  > log(-log10(q_cut) + 1),]
message('up: ', nrow(up))
dn <- dat[dat$f < -log2(f_cut) & dat$q  > log(-log10(q_cut) + 1),]
message('dn: ', nrow(dn))
rownames(dat) <- rownames(ind_diff_stat_good)
diff_genes <- rownames(dat)[abs(dat$f) > log2(f_cut) & dat$q > log(-log10(q_cut) + 1)]
message('diff: ', length(diff_genes))
dat$sig <- ""
dat[diff_genes,]$sig <- diff_genes
dat <- dat[order(dat$sig),]
pdf("./CD8cd160vsothercd8.pdf",useDingbats=F)
ggplot(dat, aes(f, q, label = sig)) +
  geom_point(color = ifelse(dat$sig == "", "grey", "red")) +
  geom_text_repel(data = dat[dat$sig != "",], col="blue",max.overlaps=16) +
  geom_vline(xintercept = log2(f_cut) , linetype="dashed", color="grey77") + 
  geom_vline(xintercept = -log2(f_cut), linetype="dashed", color="grey77") + ylab('log((-log10 q) + 1)') + xlab('log2 f')+
  xlim(-2,2)+theme_classic()
  #ggtitle("DEG of B cell between INV and Normal")
dev.off()

ggplot(dat,aes(x=,y=))+
geom_point()+
 geom_text_repel(data = dat[dat$label != "",], col="blue",max.overlaps=20) +

dat2$label=''
for(i in 1:nrow(dat2)){
  if(abs(dat2$f_cxcl13[i])>=2 | abs(dat2$f_ctl[i])>=2){
    dat2$label[i]=rownames(dat2)[i]
  }else{
    dat2$label[i]=''
  }
}
 ggplot(dat2,aes(x=f_cxcl13,y=f_ctl,label=label))+ geom_point()+ 
 geom_text_repel(size=5)+
 geom_hline(yintercept = c(0), color="grey",linetype = 2,size = 1)+
 geom_vline(xintercept = c(0), color="grey",linetype = 2,size = 1)+
 scale_color_manual(values=c("CD8-CXCL13 UP"="","CD8-PLCG2 UP"='',"CD8-CTL UP"="",'nosig'='lightgrey'))

 diseaseA='Naive-like'
  diseaseB='CD8-Tem'



aa$state=''
aa$label=''
for(i in 1:nrow(aa)){
  if(aa$CXCL13_p[i]>log(-log10(0.05) + 1) &  aa$CXCL13_fc[i]>1.5 & aa$CTL_fc[i]<0){
    aa$label[i]=rownames(aa)[i]
    aa$state[i]='CXCL13_specific'
  }else if(aa$CTL_p[i]>log(-log10(0.05) + 1) &  aa$CTL_fc[i]>1.5 & aa$CXCL13_fc[i]<0){
    aa$label[i]=rownames(aa)[i]
    aa$state[i]='CTL_specific'
  }else if(aa$CTL_p[i]>log(-log10(0.05) + 1) & aa$CXCL13_p[i]>log(-log10(0.05) + 1) & aa$CTL_fc[i]<(-1.5) & aa$CXCL13_fc[i]<(-1.5)){
    aa$label[i]=rownames(aa)[i]
    aa$state[i]='Naive_specific'
  }else if(aa$CTL_fc[i]>1.5 & aa$CXCL13_fc[i]>1.5 & aa$CTL_p[i]>log(-log10(0.05) + 1) & aa$CXCL13_p[i]>log(-log10(0.05) + 1)){
    aa$label[i]=rownames(aa)[i]
    aa$state[i]='CTL_CXCL13'
  }else{
    aa$label[i]=''
    aa$state[i]='ns'
  }
}

p<-ggplot(aa,aes(x=CXCL13_fc,y=CTL_fc,label=label,color=state))+
  geom_point(aes(color=state),size=3)+
  scale_color_manual(values=c("CXCL13_specific"='#1775B6',"CTL_specific"='#FF9998',"Naive_specific"="#C8AFD4",'ns'="lightgrey",'CTL_CXCL13'='#76AD6D'))+
  geom_text_repel(size=4,max.overlaps=25)+
  geom_hline(yintercept = c(0), color="grey",linetype = 2,size = 1)+
  geom_vline(xintercept = c(0), color="grey",linetype = 2,size = 1)+
  theme_bw()+theme(axis.text.x = element_text(angle = 0, hjust = 0, vjust = 0.5),panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))+
    xlab("CD8-CXCL13 over Naive-like")+ylab("CD8-CTL over Naive-like")+NoLegend()


data.frame(PLCG2_fc=log2((PLCG2$avgB.mean + 0.01)/(PLCG2$avgA.mean + 0.01)),PLCG2_p=log(-(log(PLCG2$wilcox.p.mean+0.0000000000001, 10)) + 1),TEM_fc=log2((TEM$avgB.mean + 0.01)/(TEM$avgA.mean + 0.01)),TEM_p=log(-(log(TEM$wilcox.p.mean+0.0000000000001, 10)) + 1))->bb
rownames(bb)=rownames(PLCG2)
bb$state=''
bb$label=''
for(i in 1:nrow(bb)){
  if(bb$PLCG2_p[i]>log(-log10(0.05) + 1) &  bb$PLCG2_fc[i]>1 & bb$TEM_fc[i]<0){
    bb$label[i]=rownames(bb)[i]
    bb$state[i]='PLCG2_specific'
  }else if(bb$TEM_p[i]>log(-log10(0.05) + 1) &  bb$TEM_fc[i]>1 & bb$PLCG2_fc[i]<0){
    bb$label[i]=rownames(bb)[i]
    bb$state[i]='TEM_specific'
  }else if(bb$TEM_p[i]>log(-log10(0.05) + 1) & bb$PLCG2_p[i]>log(-log10(0.05) + 1) & bb$TEM_fc[i]<(-1) & bb$PLCG2_fc[i]<(-1)){
    bb$label[i]=rownames(bb)[i]
    bb$state[i]='Naive_specific'
  }else if(bb$TEM_fc[i]>1 & bb$PLCG2_fc[i]>1 & bb$TEM_p[i]>log(-log10(0.05) + 1) & bb$PLCG2_p[i]>log(-log10(0.05) + 1)){
    bb$label[i]=rownames(bb)[i]
    bb$state[i]='TEM_PLCG2'
  }else{
    bb$label[i]=''
    bb$state[i]='ns'
  }
}
p1<-ggplot(bb,aes(x=PLCG2_fc,y=TEM_fc,label=label,color=state))+
  geom_point(aes(color=state),size=3)+
  scale_color_manual(values=c("PLCG2_specific"='#ED201E',"TEM_specific"='#F7F393',"Naive_specific"="#C8AFD4",'ns'="lightgrey",'PLCG2_TEM'='#76AD6D'))+
  geom_text_repel(size=4,max.overlaps=25)+
  geom_hline(yintercept = c(0), color="grey",linetype = 2,size = 1)+
  geom_vline(xintercept = c(0), color="grey",linetype = 2,size = 1)+
  theme_bw()+theme(axis.text.x = element_text(angle = 0, hjust = 0, vjust = 0.5),panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))+
    xlab("CD8-PLCG2 over Naive-like")+ylab("CD8-TEM over Naive-like")+NoLegend()

#####
diseaseA <- 'SPF'; cutA <- 15
diseaseB <- 'GF'; cutB <- 15
res<-list()
for(i in unique(cells_Str$pop)[c(1:7,9:13)]){
  print(i)
  ind_diff <- diff_calc(cells_Str, mat_ds, i, n, ntime, diseaseA, diseaseB, cutA, cutB)
  ind_diff_stat <- diff_stat(ind_diff, ntime)
  res[[i]]=ind_diff_stat
  print(paste0(i,'finished'))
}

marker_condition<-data.frame()
for(i in unique(cells_Str$pop)[c(1:7,9:13)]){
  tmpdata<-res[[i]]
  tmpdata[which(tmpdata$pctA.mean>0.1 | tmpdata$pctB.mean >0.1 ),]->tmpdata
  tmpdata$f<-log2((tmpdata$avgB.mean + 0.01)/(tmpdata$avgA.mean + 0.01))
  tmpdata$q<-log(-(log(tmpdata$wilcox.p.mean, 10)) + 1)
  tmpdata$condition<-ifelse(tmpdata$f >0,paste0(i,"_GF"),paste0(i,"_SPF"))
  tmpdata$cluster=i
  tmpdata$gene=rownames(tmpdata)
  tmpdata=tmpdata%>%arrange(desc(f))
  rownames(tmpdata)<-NULL
  marker_condition=marker_condition%>%rbind(tmpdata)
}

marker_condition$sig=""
marker_condition$sig[abs(marker_condition$f) > log2(1.5) & marker_condition$q  > (log(-log10(0.05) + 1))] = "sig"
marker_condition$sig2=paste(marker_condition$cluster,marker_condition$sig,sep = "_")
marker_condition$sig2[str_detect(marker_condition$sig2,"_$")]="not_sig"
marker_condition$sig2=str_replace(marker_condition$sig2,"_sig","")

marker_condition$sig2=factor(marker_condition$sig2,levels = c("not",sort(unique(marker_condition$cluster))))
marker_condition$cluster=factor(marker_condition$cluster,levels = sort(unique(marker_condition$cluster)))
marker_condition=marker_condition%>%arrange(cluster,sig2)

###控制范围
marker_condition$f[marker_condition$f > 3]=3
marker_condition$f[marker_condition$f < c(-3)]= -3
##配色
library(RColorBrewer)
library(scales)
color_ct=c(brewer.pal(12, "Set3")[-c(2,3,9,12)],
           brewer.pal(5, "Set1")[2],
           brewer.pal(3, "Dark2")[1])
names(color_ct)=sort(unique(as.character(marker_condition$cluster)))
###plot
marker_condition %>% ggplot(aes(x=cluster,y=f,color=sig2))+geom_jitter(width = 0.25,size=0.5)+
  scale_color_manual(values = c(color_ct,"not"="#dee1e6"))+
  scale_y_continuous("SPF  VS GF, average log2FC",expand = c(0.02,0))+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,size = 14,color = "black"),
    axis.text.y.left = element_text(size = 14,color = "black"),
    axis.title.x.bottom = element_blank(),
    axis.title.y.left = element_text(size = 16)
  )
###enhanced plot
marker_condition2=marker_condition
marker_condition2$padj_log10_neg= marker_condition2$q
marker_condition2$padj_log10_neg=ifelse(marker_condition2$f > 0,
                                        marker_condition2$padj_log10_neg,
                                        -marker_condition2$padj_log10_neg)
                            
plot.list=list()
marker_condition2$cluster<-as.character(marker_condition2$cluster)
for (ci in sort(unique(as.character(marker_condition2$cluster)))) {
  tmpdf=marker_condition2[marker_condition2$cluster == ci,]
  minabs=abs(min(tmpdf$padj_log10_neg))
  maxabs=max(tmpdf$padj_log10_neg)
  thre=0
  if(minabs < maxabs) {
    tmpdf$padj_log10_neg[tmpdf$padj_log10_neg > minabs] = minabs
    thre=minabs
  }
  if(minabs > maxabs) {
    tmpdf$padj_log10_neg[tmpdf$padj_log10_neg < (-maxabs)] = -maxabs
    thre=maxabs
  }
  if(minabs == maxabs & maxabs == Inf) {
    thre = min(
      abs(
        range(
          tmpdf$padj_log10_neg[tmpdf$padj_log10_neg < Inf & tmpdf$padj_log10_neg > -Inf]
        )
      )
    )
    tmpdf$padj_log10_neg[tmpdf$padj_log10_neg < (-thre)] = -thre
    tmpdf$padj_log10_neg[tmpdf$padj_log10_neg > thre] = thre
  }
  plotdata = tmpdf
  tmpdf=tmpdf[tmpdf$sig2 != "not",] #这里我取了logFC最极端的几个gene来标注文本，实际处理中不一定这样做
  tmpdf=tmpdf%>%arrange(desc(f))
  tmpdf.a=head(tmpdf[which(tmpdf$f>0),],5)
  tmpdf.a$d=thre*2*0.05+(-thre)-tmpdf.a$padj_log10_neg
  tmpdf.b=tail(tmpdf[which(tmpdf$f<0),],5)
  tmpdf.b$d=thre*2*0.95-thre  - tmpdf.b$padj_log10_neg
  textdata.down = tmpdf.b
  textdata.up   = tmpdf.a
  
    tmpplot=plotdata%>%ggplot(aes(x=padj_log10_neg,y=f))+
    geom_point(aes(color=sig2),size=1)+
    geom_hline(yintercept = c(-0.25,0.25),linetype="dashed")+
    geom_text_repel(data = textdata.down,
                    mapping = aes(label=gene),
                    nudge_x=textdata.down$d,
                    direction = "y", hjust = 1,segment.size = 0.2)+
    geom_text_repel(data = textdata.up,
                    mapping = aes(label=gene),
                    nudge_x=textdata.up$d,
                    direction = "y", hjust = 0,segment.size = 0.2)+
    labs(title = ci)+
    scale_color_manual(values = c(color_ct,"not"="#dee1e6"))+
    scale_y_continuous("SPF VS GF, average log2FC",expand = c(0.02,0),limits = c(-3,3))+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      
      axis.ticks.x.bottom = element_blank(),
      axis.text.x.bottom = element_blank(),
      axis.title.x.bottom = element_blank(),
      axis.text.y.left = element_text(size = 14,color = "black"),
      axis.title.y.left = element_text(size = 16),
      
      plot.title = element_text(size = 16,hjust = 0.5)
    )
  
  index=which(ci == sort(unique(as.character(marker_condition2$cluster))))
  if (index!=1) {
    tmpplot=tmpplot+theme(
      axis.title.y.left = element_blank(),
      axis.ticks.y.left = element_blank(),
      axis.text.y.left = element_blank()
    )
  }
  if (index == length(sort(unique(as.character(marker_condition2$cluster))))) {
    segment.df=data.frame(x=c(0 - thre / 5,0 + thre / 5),
                          xend=c(-thre,thre),
                          y=c(-3,-3),
                          yend=c(-3,-3))
    tmpplot=tmpplot+geom_segment(data = segment.df,
                                 mapping = aes(x=x,xend=xend,y=y,yend=yend),
                                 arrow = arrow(length=unit(0.3, "cm")))
    
  }
  plot.list[[get("index")]]=tmpplot
}
patchwork::wrap_plots(plot.list,ncol = 10)&theme(plot.margin = unit(c(0,0,0,0),"cm"))



mctoseurat=function(mat_id,mc_id,mc2d_id){
    require(Seurat)
     mat=scdb_mat(mat_id)
     mc=scdb_mc(mc_id)
     mc2d=scdb_mc2d(mc2d_id)
     dimensional=data.frame(metacell_1=mc2d@sc_x,metacell_2=mc2d@sc_y,row.names=names(mc@mc))
     dimensional[!is.na(dimensional$metacell_1),]->dimensional
     Seuratobj<-Seurat::CreateSeuratObject(rbind(mat@mat[,rownames(dimensional)],mat@ignore_gmat[,rownames(dimensional)]),meta.data=mat@cell_metadata[rownames(dimensional),])
     Seuratobj[['metacell']]=CreateDimReducObject(as.matrix(dimensional),assay = DefaultAssay(Seuratobj))
     
     return(Seuratobj)
}
mctoseurat('GC','mc_annot','GC_2dproj')->GC

c('HK1','HK2','HK3','PFKFB1','PFKFB2','PFKFB3','PFKFB4','PFKL','ALDOA','ALDOB','ALDOC','TPI1','GAPDH','PGK1','PGM1','ENO1','ENO2','PKM','LDHA','LDHB','PDK1','PDK2','PDK3','G6PD','PGD','RPIA','TKT','TALDO1')







ggplot(Diff, aes(logFC, -log10(adj.P.Val), label = sig)) +
  geom_point(color = ifelse(Diff$sig == "", "grey", "red")) +
  geom_text_repel(data = Diff[Diff$sig != "",], col="blue",max.overlaps=15) +
  geom_vline(xintercept = log2(1.5) , linetype="dashed", color="grey77") + 
  geom_vline(xintercept = -log2(1.5), linetype="dashed", color="grey77") + ylab('-log10P') + xlab('log2 f')+
    xlim(-3,3) +
    theme_classic()


diff_genes <- rownames(Diff)[abs(Diff$logFC) > 1.5 & Diff$adj.P.Val < 0.05]
Diff$sig <- ""
Diff[diff_genes,]$sig <- diff_genes
Diff <- Diff[order(Diff$sig),]

####



fat=cells_Str
fat$CellID=rownames(fat)
fat->fat.scatter2
fat1 <- droplevels(fat.scatter2)

labels <- as.character(unique(fat1$stage))

CD8IL17cell=fat1[which(fat1$type=='CD8-TRM-IL17A'),'CellID']
CD8CD160cell=fat1[which(fat1$type=='CD8-TRM-CD160'),'CellID']

otherCD8=fat1[which(fat1$type %in% c('CD8-CXCL13','CD8-GZMK','CD8-TRM-IFNG','CD8-TRM-CD160','CD8-CTL','CD8-ISG15','CD8-MT1X')),'CellID']

IACcell=fat1[which(fat1$stage=='INV' & fat1$type=='CD83+ Activated B'),'CellID']
Normalcell=fat1[which(fat1$stage=='adjacent' & fat1$type=='CD83+ Activated B'),'CellID']
length(IACcell)
length(Normalcell)

# downsample UMIs
all_gene_sets = rep(1,length(rownames(mat@mat)))
names(all_gene_sets) = rownames(mat@mat)
scdb_add_gset("all_mat_genes",gset_new_gset(all_gene_sets,"all_mat_genes"))

set.seed(2305)
dus = gset_get_feat_mat(gset_id = "all_mat_genes",mat_id = 'GC',downsamp = T)
as.matrix(dus)->dus
# downsample cells
ds <- min(c(length(intersect(IACcell,colnames(dus))),length(intersect(Normalcell,colnames(dus)))))
g1 <- sample(intersect(IACcell,colnames(dus)),ds)
g3 <- sample(intersect(Normalcell,colnames(dus)),ds)

ds <- min(c(length(intersect(CD8IL17cell,colnames(dus))),length(intersect(CD8CD160cell,colnames(dus)))))
g1 <- sample(intersect(CD8IL17cell,colnames(dus)),ds)
g3 <- sample(intersect(CD8CD160cell,colnames(dus)),ds)



umis = dus[,c(g1,g3)]
umis_n = sweep(umis,2,colSums(umis), "/") * 1000 # normalized


# Fold-change
logFC = as.vector(apply(umis, 1, function(x) log2((mean(x[g1]))/(mean(x[g2])))))
logFCx = as.vector(apply(umis, 1, function(x) log2(mean(x[g1]))))
logFCy = as.vector(apply(umis, 1, function(x) log2(mean(x[g2]))))


# Wilcoxon-Mann-Whitney test
wilc = p.adjust(as.vector(apply(umis_n, 1, function(x) wilcox.test(x[g1], x[g2])$p.value)), method = "fdr")


pb <- data.frame(rownames(umis_n), logFC, logFCx, logFCy, wilc)
colnames(pb) <- c("gene", "fc", "x", "y", "pval")

pb$pval[pb$pval=="NaN"] <- NA
pb <- na.omit(pb)

pb <- pb[-which(abs(pb$fc) == "Inf"),]
rownames(pb) <- pb$gene

threshold <- 0.05  #quantile(abs(pb$pval),1:20/20)[1]
pb[["mark"]] = ifelse(pb$pval < threshold & abs(pb$fc)>log2(1.5), "red", "#BFBEBE")
#show<-c('TIGIT','CTLA4','TOX','TOX2','IL21R','NFKBIZ','TNFRSF18','FOXP3','CD82','GABARAPL1','PLAC8','LGALS1','ITGA6','XBP1','S100A10','LMNA','ANXA1','S1PR1','KLF3')
show1<- ifelse(pb$mark %in% c("red"), as.character(pb$gene), "")
show<-c(unique(show1),'IL17A')
#pb[["show"]] = ifelse(pb$pval < threshold, as.character(pb$gene), " ")
pb[["show"]] = ifelse(pb$gene %in% show, as.character(pb$gene), " ")
#pb[["show"]] = ifelse(pb$mark == "red", as.character(pb$gene), " ")
mark <- as.character((pb[,"mark"]))
names(mark) <- pb$mark
#pb <- pb[-which(pb$x < log2(100/732) & pb$fc > 0 | pb$y < log2(100/732) & pb$fc < 0),]
pb <- pb[-which(pb$pval == 1),]
#pb$pval[pb$pval == 0] <- 6.891173e-281
pb$lp <- -log10(pb$pval)


# plot
theme_publa <- function(base_size = 8, base_family = "sans", legend_position = "right",
                     title_size=8){
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
  theme(legend.position = legend_position, legend.background = element_blank(),
        strip.background = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "#000000"), axis.ticks = element_line(colour = "#000000"), 
        legend.key = element_blank(), 
        axis.text = element_text(size = base_size, face="plain"), plot.title=element_text(face="plain", size = title_size),
        axis.title = element_text(face="plain", size = base_size), legend.text=element_text(size = base_size),
        legend.title=element_text(face="plain", size = base_size), strip.text=element_text(face="plain", size = base_size)
        )
}
require(ggplot2)

p <- ggplot(pb, aes(x=(y), y=(x), fill=mark,label =show)) + 
  labs(title  = NULL, # plot title
       y      = paste0("CD8_IL17A  log2(average UMI count)") , # x-axis label
       x      = paste0("CD8_CD160 log2(average UMI count)")) + # y-axis label
  #  geom_rect(mapping=aes(xmin=1, xmax=5, ymin=0.3619223, ymax=2.20, fill=F), color="#000000", linetype="dashed",size=0.1) +
  geom_point(shape=21, size=3, color="#BFBEBE",stroke=0.01) +
  geom_text_repel(data = pb[pb$show != " ",],col='#3551AA', size=3.5,max.overlaps=15) + 
  #  geom_text_repel(aes(label=show),size=1.5) +
  #  scale_y_continuous(limits = c(-1,2.25), expand = c(0,0)) +
  geom_abline(slope=1, intercept=0.6, colour="grey", na.rm = FALSE, show.legend = NA) + 
  geom_abline(slope=1, intercept=-0.6, colour="grey", na.rm = FALSE, show.legend = NA)+
    #xlim(-8,4) +
    #ylim(-8,4)+
  scale_fill_manual(values=mark, guide=F) + 
  #scale_fontface_continuous(guide=F) + 
  theme_publa(base_size=8)
#


ggplot(data=tmp,aes(x=CN,y=per,fill=Tissue))+
geom_()+
theme_classic()



threshold <- 0.05  #quantile(abs(pb$pval),1:20/20)[1]
pb[["mark"]] = ifelse(pb$p_val < threshold & abs(pb$avg_log2FC)>log2(1.5), "red", "#BFBEBE")
#show<-c('TIGIT','CTLA4','TOX','TOX2','IL21R','NFKBIZ','TNFRSF18','FOXP3','CD82','GABARAPL1','PLAC8','LGALS1','ITGA6','XBP1','S100A10','LMNA','ANXA1','S1PR1','KLF3')
show1<- ifelse(pb$mark %in% c("red"), as.character(pb$gene), "")
show<-c(unique(show1),c('APOE','MRC1'))
#pb[["show"]] = ifelse(pb$pval < threshold, as.character(pb$gene), " ")
pb[["show"]] = ifelse(pb$gene %in% show, as.character(pb$gene), " ")
#pb[["show"]] = ifelse(pb$mark == "red", as.character(pb$gene), " ")
mark <- as.character((pb[,"mark"]))
names(mark) <- pb$mark
#pb <- pb[-which(pb$x < log2(100/732) & pb$fc > 0 | pb$y < log2(100/732) & pb$fc < 0),]
pb <- pb[which(pb$p_val !=1),]
pb$p_val[pb$p_val == 0] <- 6.891173e-281
pb$lp <- -log10(pb$p_val)

 ggplot(pb, aes(x=(y), y=(x), fill=mark,label =show)) + 
  labs(title  = NULL, # plot title
       y      = paste0("S100P+  log2(average UMI count)") , # x-axis label
       x      = paste0("S100P- log2(average UMI count)")) + # y-axis label
  #  geom_rect(mapping=aes(xmin=1, xmax=5, ymin=0.3619223, ymax=2.20, fill=F), color="#000000", linetype="dashed",size=0.1) +
  geom_point(shape=21, size=3, color="#BFBEBE",stroke=0.01) +
  geom_text_repel(data = pb[pb$show != " ",],col='#3551AA', size=3.5,max.overlaps=10000) + 
  #  geom_text_repel(aes(label=show),size=1.5) +
  #  scale_y_continuous(limits = c(-1,2.25), expand = c(0,0)) +
  geom_abline(slope=1, intercept=0.6, colour="grey", na.rm = FALSE, show.legend = NA) + 
  geom_abline(slope=1, intercept=-0.6, colour="grey", na.rm = FALSE, show.legend = NA)+
    xlim(-8,8) +
    ylim(-8,8)+
  scale_fill_manual(values=mark, guide=F) + 
  #scale_fontface_continuous(guide=F) + 
  theme_publa(base_size=8)


deg[["mark"]] = ifelse(deg$p_val < 0.05 & abs(deg$avg_log2FC)>log2(1.5), "red", "#BFBEBE")
show1<- ifelse(deg$mark %in% c("red"), as.character(deg$gene), "")
show<-c(show1,'MRC1','CD163')
show[show!='']->show

mark <- as.character((deg[,"mark"]))
names(mark) <- deg$mark

deg[["show"]] = ifelse(deg$gene %in% show, as.character(deg$gene), " ")
 ggplot(deg, aes(x=APOE, y=MRC1, fill=mark,label =show)) + 
  labs(title  = NULL, # plot title
       y      = paste0("Mac-MRC1 log2(average UMI count)") , # x-axis label
       x      = paste0("Mac-APOE log2(average UMI count)")) + # y-axis label
  #  geom_rect(mapping=aes(xmin=1, xmax=5, ymin=0.3619223, ymax=2.20, fill=F), color="#000000", linetype="dashed",size=0.1) +
  geom_point(shape=21, size=3, color="#BFBEBE",stroke=0.01) +
  geom_text_repel(data = deg[deg$show != " ",],col='#3551AA', size=3.5,max.overlaps=20) + 
  #  geom_text_repel(aes(label=show),size=1.5) +
  #  scale_y_continuous(limits = c(-1,2.25), expand = c(0,0)) +
  geom_abline(slope=1, intercept=0.6, colour="grey", na.rm = FALSE, show.legend = NA) + 
  geom_abline(slope=1, intercept=-0.6, colour="grey", na.rm = FALSE, show.legend = NA)+
    #xlim(-8,4) +
    #ylim(-8,4)+
  scale_fill_manual(values=mark, guide=F) + 
  #scale_fontface_continuous(guide=F) + 
  theme_publa(base_size=8)



scores <-  rawEnrichmentAnalysis(fpkm, 
                                 signatures = xCell.data$signatures,
                                 genes = xCell.data$genes, 
                                 parallel.sz = 4,
                                 parallel.type = "SOCK")
tscores <-  transformScores(scores,
                              fit.vals = xCell.data$spill$fv,
                              scale = T)


                              scores1 <-  rawEnrichmentAnalysis(exp, 
                                 signatures = sig,
                                 genes = xCell.data$genes, 
                                 parallel.sz = 4,
                                 parallel.type = "SOCK")



gl<-list()
for(i in names(dat)){
  gl[[i]]=rownames(dat[[i]])
}
Reduce(gl,intersect)
for(i in names(dat)){
dat[[i]]=subset(dat[[i]],features=gene)
}


geneglist <- list("CD69" = c("CD69"),
                  "ITGAE" = c("ITGAE"),
                  "TIGIT" = c("TIGIT"),
                  "HAVCR2" = c("HAVCR2"),
                  "CD69_TIGIT" = c("CD69",'TIGIT'),
                  "CD69_HAVCR2" = c("CD69",'HAVCR2'),
                  "TIGIT_HAVCR2" = c("TIGIT",'HAVCR2'),
                  'CD69_TIGIT_HAVCR2'=c("CD69",'HAVCR2','TIGIT'),
                  'CD69_TIGIT_HAVCR2_ITGAE'=c("CD69",'HAVCR2','TIGIT','ITGAE')
                  ) 