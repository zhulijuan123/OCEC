#################################################################################################################
###Create working directory
setwd("F:/obesity_cancer/")
##################################################################### 
### General
DATA.FOLD <- paste0("./Data/")
NETWORK.FOLD <- paste0(DATA.FOLD, "R_input/")
if(!file.exists(NETWORK.FOLD)){dir.create(NETWORK.FOLD)} 
NETWORK.FOLD1 <- paste0(DATA.FOLD, "R_output/")
if(!file.exists(NETWORK.FOLD1)){dir.create(NETWORK.FOLD1)} 


##################################################################################################################
###the genes of cancer
addg=list()
n=0
dgbp<- file(paste0(NETWORK.FOLD,'different cancer map gene.txt'), "r")
line=readLines(dgbp,n=1)
while( length(line) != 0 ) {
   n=n+1
   addg[[n]]=line  
   line=readLines(dgbp,n=1)

}
close(dgbp)
###
fadge=list()
for(i in 1:length(addg)){
fadge[[i]]=unlist(strsplit(addg[[i]], "\t"))
}
zcg=fadge


##################################################################### 
######the genes of obesity
addg=list()
n=0
dgbp<- file(paste0(NETWORK.FOLD,'different obesity map gene.txt'), "r")
line=readLines(dgbp,n=1)
while( length(line) != 0 ) {
   n=n+1
   addg[[n]]=line  
   line=readLines(dgbp,n=1)

}
close(dgbp)
###
fadge=list()
for(i in 1:length(addg)){
fadge[[i]]=unlist(strsplit(addg[[i]], "\t"))
}
og=fadge[[1]]


##################################################################################
###gene set enrichment analysis
library(HPO.db)
library(org.Hs.eg.db)
library(clusterProfiler)

################################################################################## 
###Convert genesymbol to ID(obesity)
bog <- bitr(og1, fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
oid=bog[,2]
geneList <- oid
###
ego_result <- enrichGO(geneList, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
write.table(ego_result,paste0(NETWORK.FOLD, 'ego_result.txt'), sep = '\t', row.names = FALSE, quote = FALSE)

################################################################################## 
###Convert genesymbol to ID(cancer)
for(s in 1:length(zcg)){
bcg <- bitr(zcg[[s]], fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
cid=bcg[,2]
geneList <- cid
###
ego_result <- enrichGO(geneList, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
write.table(ego_result,paste0(NETWORK.FOLD, 'cego_result',s,'.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
}


################################################################################## 
###selected MN
jbg=list()
n=0
jbbp<- file(paste0(NETWORK.FOLD,'selected subnet index_MN_adg.txt'), "r")
line=readLines(jbbp,n=1)
while( length(line) != 0 ) {
   n=n+1
   jbg[[n]]=line  
   line=readLines(jbbp,n=1)

}
close(jbbp)
###
fdg=list()
for(j in 1:length(jbg)){
fdg[[j]]=unlist(strsplit(jbg[[j]], " "))
}
###
xb=list()
for(i in 1:length(zcg)){
xb[[i]]=intersect(as.numeric(fdg[[1]]), as.numeric(fdg[[i+2]]))
}

################################################################################## 
###BP terms
bpjg=list()
m=0
bp<- file(paste0(NETWORK.FOLD,'hsa.GO.BP.txt'), "r")
line=readLines(bp,n=1)
while(length(line)!= 0 ) {
   m=m+1
   bpjg[[m]]=line  
   line=readLines(bp,n=1)

}
close(bp)
###
fqu=list()
symbollist=list() 
bpg=list()
for(i in 1:length(bpjg)){
fqu[[i]]=unlist(strsplit(bpjg[[i]], "\t"))
symbollist[[i]]=fqu[[i]][1]
bpg[[i]]=fqu[[i]][3:length(fqu[[i]])]
}
symbollist=unlist(symbollist)
###
sbpg=list()
n=0
for(i in 1:length(bpg)){
if(length(bpg[[i]])>=30&length(bpg[[i]])<=500){
n=n+1
sbpg[[n]]=bpg[[i]]
 }
}

###select BPs
swz=list()
n=0
for(i in 1:length(bpg)){
if(length(bpg[[i]])>=30&length(bpg[[i]])<=500){
n=n+1
swz[[n]]=i
 }
}
swz=unlist(swz)
###
aid=list()
gna=list()
for(i in 1:length(symbollist)){
fs=unlist(strsplit(symbollist[i], "_"))
aid[[i]]=fs[1]
gna[[i]]=fs[2]
}
aid=unlist(aid)
tqaid=aid[swz]
gna=unlist(gna)
tqna=gna[swz]

##################################################################################
###
ego_result <- read.delim(paste0(NETWORK.FOLD,'ego_result.txt'), stringsAsFactors = FALSE)
yz=0.005
pzhi=ego_result['p.adjust']
sid=unlist(ego_result['ID'])
wz=which(pzhi<=yz)
zid1=unlist(sid[wz])
zwxb=list()
for(s in 1:length(zcg)){
ego_result <- read.delim(paste0(NETWORK.FOLD,'cego_result',s,'.txt'), stringsAsFactors = FALSE)
pzhi=ego_result['p.adjust']
sid=unlist(ego_result['ID'])
wz=which(pzhi<=yz)
zid=unlist(sid[wz])
xid=intersect(zid,zid1)
wz=match(xid,tqaid)
xwz=na.omit(wz)
zwxb[[s]]=xwz
}

##################################################################################
###Calculate the OCEC_eDMN
ppn<- read.delim(paste0(NETWORK.FOLD,'human network.txt'),header=F,stringsAsFactors=F)[,c(1,3)]
ppn<- ppn[ppn[,1]!=ppn[,2],]
###
nawz1=which(is.na(ppn[,1])==TRUE)
nawz2=which(is.na(ppn[,2])==TRUE)
nawz=union(nawz1,nawz2)
ppn=ppn[-nawz,]
###
uhppg=unique(unlist(ppn))
###
bpng=list()
bps=rep(0,length(sbpg))
for(i in 1:length(sbpg)){
bpng[[i]]=intersect(uhppg,sbpg[[i]])
bps[i]=length(bpng[[i]])
}
###
eMN.FOLD <- paste0(NETWORK.FOLD1, "expand BP network/")
if(!file.exists(eMN.FOLD)){dir.create(eMN.FOLD)} 
##################################################################################
###random walk on subnetworks
ljjz=matrix(0,length(uhppg),length(uhppg))
for(i in 1:length(uhppg)){
wz1=which(ppn[,1]==uhppg[i])
hzg1=ppn[wz1,2]
wz2=which(ppn[,2]==uhppg[i])
hzg2=ppn[wz2,1]
hzg=c(hzg1,hzg2)
hzgwz=match(hzg,uhppg)
ljjz[i,hzgwz]=1
ljjz[hzgwz,i]=1
}
###
bljz=matrix(0,length(uhppg),length(uhppg))
for(i in 1:length(uhppg)){
bljz[,i]=ljjz[,i]/sum(ljjz[,i])
}
###
sp=list()
r=0.7
for(i in 1:length(bpng)){
pwz=match(bpng[[i]],uhppg)
p=matrix(0,length(uhppg),1)
p[pwz]=1/length(bpng[[i]])
pxh1=p
pxh2=(1-r)*bljz%*%pxh1+r*p
while(max(abs(pxh2-pxh1))>10^(-6)){
pxh1=pxh2
pxh2=(1-r)*bljz%*%pxh2+r*p
 }
sp[[i]]=pxh2
}
###
for(s in 1:length(sp)){
write.table(sp[[s]],file= paste0(NETWORK.FOLD1,'expand BP network/expand BP network probability',s,'.txt'), row.names=F, col.names=F, quote=F, sep='\n')
}

##################################################################### 
###expand subnetworks
N=2
kg=list()
for(i in 1:length(sp)){
psp=sort(unlist(sp[[i]]),decreasing=TRUE)
ks=min(N*length(bpng[[i]]),500)
gwz=match(psp[1:ks],unlist(sp[[i]]))
kg[[i]]=uhppg[gwz]
}
###
for(s in 1:length(kg)){
write.table(kg[[s]],file= paste0(NETWORK.FOLD1,'expand BP network/expand BP gene',s,'.txt'), row.names=F, col.names=F, quote=F, sep='\n')
}

##################################################################### 
###random walk on obesity and cancer
kbog=list()
for(i in 1:length(kg)){
kbog[[i]]=intersect(og,kg[[i]])
}

###cancer
kbcg=list()
for(i in 1:length(zcg)){
kbcg[[i]]=list()
for(j in 1:length(kg)){
kbcg[[i]][[j]]=intersect(zcg[[i]],kg[[j]])
 }
}

##################################################################### 
###RWR expand
library(foreach)
library(iterators)
library(parallel)
library(doParallel)
cl <- makeCluster(8)
registerDoParallel(cl)
###
fun <- function(s){
bphz=ppn[(ppn[,1]%in% kg[[s]])&(ppn[,2]%in% kg[[s]]),]
###
bpfile<-paste0(NETWORK.FOLD1,'expand BP network/expand BP subnetwork interaction',s,'.txt')
sink(bpfile,append = T)
for(j in 1:nrow(bphz))
{
cat(bphz[j,1])
cat('\t')
cat(bphz[j,2])
cat('\n') 
}  
sink() 
}
###
bpn<- foreach(s=1:length(kg)) %dopar% fun(s)
stopCluster(cl)
###
bppn=list()
for(s in 1:length(kg)){
bppn[[s]]<- read.delim(paste0(NETWORK.FOLD1,'expand BP network/expand BP subnetwork interaction',s,'.txt'),header=F,stringsAsFactors=F)
}
###
r=0.7
cl <- makeCluster(8)
registerDoParallel(cl)
fun <- function(i){
ubpg=unique(unlist(bppn[[i]]))
ljjz=matrix(0,length(ubpg),length(ubpg))
for(i1 in 1:length(ubpg)){
wz1=which(bppn[[i]][,1]==ubpg[i1])
hzg1=bppn[[i]][wz1,2]
wz2=which(bppn[[i]][,2]==ubpg[i1])
hzg2=bppn[[i]][wz2,1]
hzg=c(hzg1,hzg2)
hzgwz=match(hzg,ubpg)
ljjz[i1,hzgwz]=1
ljjz[hzgwz,i1]=1
}
###
bljz=matrix(0,length(ubpg),length(ubpg))
for(i2 in 1:length(ubpg)){
bljz[,i2]=ljjz[,i2]/sum(ljjz[,i2])
}
###RWR expand on obesity
kpo=list()
tpwz=match(kbog[[i]],ubpg)
tp=matrix(0,length(ubpg),1)
tp[tpwz]=1/length(kbog[[i]])
###
pxh1=tp
pxh2=(1-r)*bljz%*%pxh1+r*tp
while(max(abs(pxh2-pxh1))>10^(-6)){
pxh1=pxh2
pxh2=(1-r)*bljz%*%pxh2+r*tp
}
kpo=pxh2
###
eDN1.FOLD <- paste0(NETWORK.FOLD1, "expand OB1 network/")
if(!file.exists(eDN1.FOLD)){dir.create(eDN1.FOLD)} 
###
write.table(kpo,file= paste0(NETWORK.FOLD1,'expand OB1 network/expand OB1 network based on expand BP subnetwork',i,'.txt'), row.names=F, col.names=F, quote=F, sep='\n')


###RWR expand on cancer
eDD.FOLD <- paste0(NETWORK.FOLD1, "expand cancer network/")
if(!file.exists(eDD.FOLD)){dir.create(eDD.FOLD)} 
###
kpc=list()
for(s in 1:length(zcg)){
pwz=match(kbcg[[s]][[i]],ubpg)
p=matrix(0,length(ubpg),1)
p[pwz]=1/length(kbcg[[s]][[i]])
###
pxh1=p
pxh2=(1-r)*bljz%*%pxh1+r*p
while(max(abs(pxh2-pxh1))>10^(-6)){
pxh1=pxh2
pxh2=(1-r)*bljz%*%pxh2+r*p
  }
kpc[[s]]=pxh2
###
write.table(kpc[[s]],file= paste0(NETWORK.FOLD1,'expand cancer network/expand cancer network based on expand BP subnetwork',i,'_',s,'.txt'), row.names=F, col.names=F, quote=F, sep='\n')
 }
}
###
bpcn<- foreach(i=1:length(bppn)) %dopar% fun(i)
stopCluster(cl)

###obesity
kogv=list()
for(i in 1:length(bppn)){
kogv[[i]]<- read.delim(paste0(NETWORK.FOLD1,'expand OB1 network/expand OB1 network based on expand BP subnetwork',i,'.txt'),header=F,stringsAsFactors=F)
}
###cancer
kcgv=list()
for(s in 1:length(zcg)){
kcgv[[s]]=list()
for(i in 1:length(bppn)){
kcgv[[s]][[i]]<- read.delim(paste0(NETWORK.FOLD1,'expand cancer network/expand cancer network based on expand BP subnetwork',i,'_',s,'.txt'),header=F,stringsAsFactors=F)
 }
}

##################################################################### 
###OCEC_eDMN
cl <- makeCluster(8)
registerDoParallel(cl)
fun <- function(N){
kkog=list()
kkbcg=list()
###
for(s in 1:length(zcg)){
kkog[[s]]=list()
kkbcg[[s]]=list()
###
for(i in 1:length(bppn)){
pcv=unlist(kcgv[[s]][[i]])
ubpg=unique(unlist(bppn[[i]]))
bcg=setdiff(kg[[i]],ubpg)
nubpg=union(ubpg,bcg)
npcv=c(pcv,rep(0,length(bcg)))
xpcv=rep(0,length(npcv))
gtwz=match(kbog[[i]],nubpg)
if(length(gtwz)>0){
xpcv[gtwz]=npcv[gtwz]+sqrt((length(nubpg)-length(kbog[[i]]))/length(kbog[[i]]))
xpcv[-gtwz]=npcv[-gtwz]-sqrt(length(kbog[[i]])/(length(nubpg)-length(kbog[[i]])))
}
###
psp=sort(xpcv,decreasing=TRUE)
ks=min(N*length(kbcg[[s]][[i]]),length(nubpg))
if(ks>0){
gwz=match(psp[1:ks],xpcv)
kkbcg[[s]][[i]]=nubpg[gwz]
} 
###
pov=unlist(kogv[[i]])
ubpg=unique(unlist(bppn[[i]]))
bcg=setdiff(kg[[i]],ubpg)
nubpg=union(ubpg,bcg)
npov=c(pov,rep(0,length(bcg)))
xpov=rep(0,length(npov))
gtwz=match(kbcg[[s]][[i]],nubpg)
if(length(gtwz)>0){
xpov[gtwz]=npov[gtwz]+sqrt((length(nubpg)-length(kbcg[[s]][[i]]))/length(kbcg[[s]][[i]]))
xpov[-gtwz]=npov[-gtwz]-sqrt(length(kbcg[[s]][[i]])/(length(nubpg)-length(kbcg[[s]][[i]])))
 }
###
psp=sort(xpov,decreasing=TRUE)
ks=min(N*length(kbog[[i]]),length(nubpg))
if(ks>0){
gwz=match(psp[1:ks],xpov)
kkog[[s]][[i]]=nubpg[gwz]
}
}
}
###
ltd=matrix(0,length(zcg),length(bppn))
for(s1 in 1:length(zcg)){
for(i1 in 1:length(bppn)){
hz1=bppn[[i1]][(bppn[[i1]][,1]%in% kkog[[s1]][[i1]])&(bppn[[i1]][,2]%in%kkbcg[[s1]][[i1]]),]
hz2=bppn[[i1]][(bppn[[i1]][,2]%in% kkog[[s1]][[i1]])&(bppn[[i1]][,1]%in%kkbcg[[s1]][[i1]]),]
ltd[s1,i1]=nrow(rbind(hz1,hz2))
}
}
write.table(ltd,file= paste0(NETWORK.FOLD1,'OB1_interactions',N,'.txt'), row.names=F, col.names=F, quote=F, sep='\t')
}
###
khz<- foreach(N=1:10) %dopar% fun(N)
stopCluster(cl)



##################################################################################
###
for(bs in 1:10){
###bs=7
ec=read.table(paste0(NETWORK.FOLD1,'OB1_interactions',bs,'.txt'),header=F,sep='\t')
mec=rep(0,length(zcg))
for(i in 1:length(zcg)){
zxb=intersect(xb[[i]],zwxb[[i]])
cec=unlist(ec[i,zxb])
mec[i]=mean(cec[cec>0&cec<3000])
}
}














