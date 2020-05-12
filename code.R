setwd("/raid1/Liuw/fly/");
options(stringsAsFactors = FALSE)
library(WGCNA);
enableWGCNAThreads(6)
dat0=read.csv("fly_entrez.csv",header=TRUE)
datSummary=dat0[,1:1];
dim(dat0)
datExpr = t(dat0[,2: ncol(dat0)]);
no.samples = dim(datExpr)[[1]];
dim(datExpr)
library(preprocessCore)
datExpr=t(normalize.quantiles(as.matrix(dat0[,-1])))
dimnames(datExpr)=list(names(dat0)[-1],dat0[,1])
GeneName= dat0$EntrezID
ArrayName= names(data.frame(dat0[,-1]))
powers=c(seq(1,10,by=1),seq(12,14,by=2));
sft=pickSoftThreshold(datExpr, powerVector=powers,networkType = "signed")
#sft=pickSoftThreshold(datExpr, powerVector=powers,networkType = "signed",corFnc = bicor,corOptions = list(maxPOutliers =0.1))
RpowerTable=sft[[2]]
sizeGrWindow(9, 5);
pdf('choosing power.pdf');
par(mfrow = c(1,2));cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.80,col="red");
dev.off()
# Mean connectivity as a function of the soft-thresholding power
sizeGrWindow(9, 5);
pdf('mean connectivity.pdf');
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red");
dev.off()
softPower =12
Connectivity=softConnectivity(datExpr,corFnc = "cor", corOptions = "use ='p'",power=softPower,type="signed")
#Connectivity=softConnectivity(datExpr,corFnc = "bicor", corOptions =list(maxPOutliers =0.1),
                              #power=softPower,type="signed")
pdf("scale-free.pdf");
scaleFreePlot(Connectivity,nBreaks = 10,truncated = FALSE,removeFirst = FALSE, main = "");
dev.off()
adjacency = adjacency(datExpr,corFnc = "cor", corOptions = "use ='p'",
                      type = "signed", power = softPower)
#adjacency = adjacency(datExpr,corFnc = "bicor", corOptions =list(maxPOutliers =0.1),
#                      type = "signed", power = softPower)
TOM = TOMsimilarity(adjacency,TOMType="signed");dissTOM = 1-TOM
#method="complete"  ?
geneTree = hclust(as.dist(dissTOM), method = "average")#高版本已经用hclust
minModuleSize =30;
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 4, 
                            pamRespectsDendro = FALSE,minClusterSize = minModuleSize,
                            cutHeight=0.99);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

MEList = moduleEigengenes(datExpr, colors = dynamicMods)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average");#
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
MEDissThres = 0.2
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr, dynamicMods, cutHeight = MEDissThres, verbose = 3);
mergedColors = merge$colors;
mergedMEs = merge$newMEs;
sizeGrWindow(12, 9)
pdf("DendroAndColors.pdf")
plotDendroAndColors(geneTree, cbind(dynamicMods, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),dendroLabels = FALSE, 
                    hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleColors = mergedColors
colorOrder = c("grey", standardColors(unique(moduleColors)));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average");#
pdf("METree.pdf")
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
dev.off()
MEList = moduleEigengenes(datExpr, colors = dynamicMods)
nSamples=nrow(datExpr)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = cbind.data.frame(datSummary,corPvalueStudent(as.matrix(geneModuleMembership), 
                                                        nSamples));
write.table(data.frame(ArrayName,MEs),"MEs.csv",row.name=F)
kMEdat=data.frame(geneModuleMembership,MMPvalue)
write.table(data.frame(datSummary,kMEdat),"kME-MMPvalue.csv",row.names=FALSE)
k.in=intramodularConnectivity(adjacency(datExpr,corFnc = "cor", corOptions = "use ='p'",type = "signed", power = softPower),moduleColors,scaleByMax = FALSE)
#k.in=intramodularConnectivity(adjacency(datExpr,corFnc = "bicor", corOptions =list(maxPOutliers =0.1),type = "signed", power = softPower),moduleColors,scaleByMax = FALSE)

datout=data.frame(datSummary, colorNEW=moduleColors, k.in)
write.table(datout, file="OutputCancerNetwork.csv", sep=",", row.names=F)
hubs    = chooseTopHubInEachModule(datExpr, moduleColors)
write.csv(data.frame(module=names(hubs),moduleColor=labels2colors(names(hubs)),hub=hubs),
          "num2color.csv",row.names=F)

gene=read.csv("OutputCancerNetwork.csv",header=T)
library(gProfileR)
for (i in unique(gene$colorNEW)){
  genes=subset(gene$datSummary,gene$colorNEW==i)
  go=gprofiler(genes, 
               organism = "dmelanogaster",numeric_ns="ENTREZGENE")
  write.table(go,"module_enrichment.csv",append =T,row.names=rep(i,nrow(go)),sep=",")}

moduleColors=gene$colorNEW #keep gene order same between gene and dat0
modules=unique(moduleColors)
n=length(modules)
pb <- txtProgressBar(min = 0, max = n, style = 3)
for (p in 1:n){
  inModule = is.finite(match(moduleColors, modules[p]));
  dat2=data.frame(t(datExpr))[inModule,] #make sure it is matrix data
  resamples=lapply(1:1000,function(i) a=t(sample(dat2[,1:nSamples],round(nSamples/2),replace=F)))
  K1=sapply(resamples,softConnectivity,power= softPower,type="signed") #,type="signed"?
  K=softConnectivity(t(dat2[,1:nSamples]),power= softPower,type="signed") #,type="signed"
  #outfile=paste(modules[p],"-edit.txt",sep="")
  write.table(data.frame(mean(cor(K,K1)),apply(cor(K,K1),1,sd)), file = "module-stability.csv", row.names = modules[p], append = TRUE, col.names = FALSE, sep = ", ")
  setTxtProgressBar(pb, p)}#所有loop结果写到一个文件
close(pb) 

#automatic finish the Cytoscape mods exporting
probes = dat0[,1]
n=length(unique(moduleColors))
pb <- txtProgressBar(min = 0, max = n, style = 3)
for (p in 1:n){ modules=unique(moduleColors)[p]
inModule = is.finite(match(moduleColors,modules));modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               #nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,threshold = quantile(abs(modTOM),probs=0.8),nodeNames = modProbes ,nodeAttr = moduleColors[inModule]);
setTxtProgressBar(pb, p)}
close(pb)
