### R code from vignette source 'PRECISE.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: PRECISE
###################################################
options(keep.source = TRUE, width = 60)
PenPC <- packageDescription("PRECISE")


###################################################
### code chunk number 2: PRECISE.Rnw:31-32
###################################################
library(PRECISE)


###################################################
### code chunk number 3: PRECISE.Rnw:39-43
###################################################
pw.array = c("Apoptosis","Breast reactive","Cell cycle","Core reactive"
             ,"DNA damage response","EMT","PI3K/AKT","RAS/MAPK","RTK"
             ,"TSC/mTOR","Hormone receptor","Hormone signaling (Breast)")
length(pw.array)


###################################################
### code chunk number 4: PRECISE.Rnw:46-50
###################################################
data(rppadat)
length(rppadat)
RPPAdat = rppadat[[1]] # RPPA data for Apoptosis pathway
dim(RPPAdat)


###################################################
### code chunk number 5: PRECISE.Rnw:57-66
###################################################
p = ncol(RPPAdat)
B=100
alpha = 0.1
alpha.array = seq(0.0001,0.1,length=100)
#pcfit = fitPC(dat = RPPAdat,alpha=alpha,stable=TRUE
#,alpha.array=alpha.array,B=B,labels=as.character(1:p),verbose=T)
data(pcfit)
pcfit =  pcfit.list[[1]] # the fitPC result for Apoptosis pathway
names(pcfit)


###################################################
### code chunk number 6: PRECISE.Rnw:69-75
###################################################
adj.pc = as(pcfit$fit@graph,"matrix")
addr = matrix(as.numeric(pcfit$rownames),ncol=2)
addr.rev = cbind(addr[,2],addr[,1])
adj.pc[addr] = adj.pc[addr] * pcfit$pcmaxSel[,1]
adj.pc[addr.rev] = adj.pc[addr.rev] * pcfit$pcmaxSel[,2]
adj.pc


###################################################
### code chunk number 7: PRECISE.Rnw:78-80
###################################################
data(ppi)
names(ppi)


###################################################
### code chunk number 8: PRECISE.Rnw:83-108
###################################################
pw="Apoptosis"
genelist = strsplit(colnames(RPPAdat),split=", ")
membership = rep(1:p,lapply(genelist,length))
indigenes = unlist(genelist)
stringscore = ppi$ppi[[which(names(ppi$ppi) == pw)]]
# select edges included in the RPPAdat
stringscore = stringscore[rowSums(cbind(as.numeric(stringscore[,1] %in% indigenes)
                    ,as.numeric(stringscore[,2] %in% indigenes))) == 2,]
address = cbind(match(stringscore[,1],indigenes),match(stringscore[,2],indigenes))
adj.string = matrix(0,ncol=length(indigenes),nrow=length(indigenes))
adj.string[address] = as.numeric(stringscore[,3])/1000
adj.string.cl = matrix(0,nrow=max(membership),ncol=max(membership))

id = 1:length(membership)
w.off = rbind(which(upper.tri(adj.string.cl),arr.ind=T)
              ,which(lower.tri(adj.string.cl),arr.ind=T))
v.off = c(which(upper.tri(adj.string.cl),arr.ind=F)
          ,which(lower.tri(adj.string.cl),arr.ind=F))

for (i in 1:nrow(w.off)) {
  addr1 = id[membership == w.off[i,1]]
  addr2 = id[membership == w.off[i,2]]
  adj.string.cl[v.off[i]] = mean(adj.string[addr1,addr2])
}
Gmat = adj.string.cl/2 + adj.pc/2 ### Make averages 


###################################################
### code chunk number 9: PRECISE.Rnw:113-122
###################################################
data(mRNA)
data(miRNA)
data(Methylation)
dim(mRNA$Data)
dim(miRNA$Data)
dim(Methylation$Data)
head(mRNA$Des)
head(miRNA$Des)
head(Methylation$Des)


###################################################
### code chunk number 10: PRECISE.Rnw:128-135
###################################################
data(rppasample) ## sample names for the KIRC RPPAdat 
covsample = colnames(mRNA$Data)
intsample = intersect(covsample,rppasample)
RPPAdat = RPPAdat[match(intsample,rppasample),]
mRNA$Data = mRNA$Data[,match(intsample,covsample)]
miRNA$Data = miRNA$Data[,match(intsample,covsample)] 
Methylation$Data = Methylation$Data[,match(intsample,covsample)]


###################################################
### code chunk number 11: PRECISE.Rnw:138-144
###################################################
data(anno.miRNA) 
miRNAname =  gsub("mir","miR",miRNA$Des)
anno.miRNA = anno.miRNA[anno.miRNA[,1]%in%miRNAname,c(7,1)] 
## reduce the annotation file with microRNAs in the dataset
anno.miRNA[,2] = gsub("miR","mir",anno.miRNA[,2])
data(anno.methyl)


###################################################
### code chunk number 12: PRECISE.Rnw:147-154
###################################################
rownames(RPPAdat) = intsample
if (!is.null(miRNA$Des))miRNA$Des = as.matrix(miRNA$Des)
dat = getregcovDat(Gmat = Gmat,RPPAdat=RPPAdat,mRNA=mRNA,miRNA=miRNA
          ,Methylation=Methylation,anno.miRNA=anno.miRNA,anno.methyl=anno.methyl)
names(dat)
length(dat$ylist)
length(dat$Xlist)


###################################################
### code chunk number 13: computation
###################################################
bmsfit= getBMS(dat,Gmat)


###################################################
### code chunk number 14: PRECISE.Rnw:160-161
###################################################
names(bmsfit)


###################################################
### code chunk number 15: PRECISE.Rnw:164-167
###################################################
nodes = names(dat$ylist)
netfit = getPosteriors(bmsfit$outlist,nodes)
names(netfit)


###################################################
### code chunk number 16: PRECISE.Rnw:170-171
###################################################
netfit$G>0.5 # median probability model


###################################################
### code chunk number 17: PRECISE.Rnw:174-177
###################################################
length(netfit$intGlist)
nodes[1]
netfit$intGlist[[1]] # the integrative network for nodes[1]


###################################################
### code chunk number 18: computation
###################################################
delta = 0.5
psNet = getPRECISE(bmsfit$outlist,bmsfit$pdlist,nodes,delta)


###################################################
### code chunk number 19: PRECISE.Rnw:187-188
###################################################
names(psNet)


###################################################
### code chunk number 20: PRECISE.Rnw:191-193
###################################################
dim(psNet$net)
head(psNet$net)


###################################################
### code chunk number 21: PRECISE.Rnw:198-200
###################################################
dim(psNet$score.mat)
head(psNet$score.mat)


###################################################
### code chunk number 22: PRECISE.Rnw:203-206
###################################################
length(psNet$net.status)
head(psNet$net.status)
head(apply(psNet$score.mat,1,which.max))


