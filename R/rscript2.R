source("/home/lluc/R_functions/functions.R")
source("~/R_functions/plfaplotfunction.R")


setwd("/home/lluc/Dropbox/MUN/phd/results/NL-BELT/Batches21-25combined/R")

mols<-  data.frame(read.csv("mols.csv"))
fames<-  read.csv("fames.csv")
samples<-    read.csv("samples.csv")

cond.rel<-samples$trusted_rel==T
cond.abs<-samples$trusted_abs==T
cond.ref <- samples$Region!="Ref"
cond1 <-cond.abs & cond.ref
cond2 <-cond.rel& cond.ref

sum<-rowSums(mols)

sum

rel<-100*mols/rowSums(mols)

rel
ord<-rda(rel[cond.rel,T])
plot(ord)
ord.plot(ord,site.sep1=samples$Horizon[cond.rel], site.sep2=samples$Site[cond.rel], pch=c(21,21,21,22,22,22,24), pt.bg=1:6, spe.label.type="text", spe.labels=fames$FAME, cex.leg=.5)

ord2<-rda(rel[samples$Horizon=="F"&cond.rel,T])

ord.plot(ord2, site.sep1=samples$Site[samples$Horizon=="F"&cond.rel], site.sep=samples$Region[samples$Horizon=="F"&cond.rel], pch=c(21,21,21,22,22,22,24), pt.bg=1:6, spe.label.type="text", spe.labels=fames$FAME, cex.leg=.5)

ord.l<-rda(rel[samples$Horizon=="L"&cond.rel,T])
ord.plot(ord.l, site.sep1=samples$Site[samples$Horizon=="L"&cond.rel], site.sep=samples$Region[samples$Horizon=="L"&cond.rel], pch=c(21,21,21,22,22,22,24), pt.bg=1:6, spe.label.type="text", spe.labels=fames$FAME, cex.leg=.5)

fac<-samples$Site
horlev<-c("L","F","H","B")
hor<-samples$Horizon
cond1<-is.element(hor,horlev)
hor1<-factor(hor[cond1], ordered=T, levels=horlev)
fac1<-factor(fac[cond1], levels=unique(fac[cond1]))

var1<-sum[cond1]

means<-tapply(var1, list(fac1,hor1), mean)
error<-tapply(var1, list(fac1,hor1), stderr) 

hor.plot(sum[cond.abs]/samples$percent_C[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="bottomleft", xlab="nmol PLFA g-1 TOC", ylab="Soil horizon")

hor.plot(mols$X16.0[cond.abs]/samples$percent_C[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="bottomleft", xlab="nmol PLFA g-1 TOC", ylab="Soil horizon")

hor.plot(mols$X16.0[cond.abs]/samples$percent_C[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6),  pt.bg=c(2,2,2,3,3,3), legpl="bottomright", xlab="nmol PLFA g-1 TOC", ylab="Soil horizon", pch=rep(21:23,2), er.type="se")

hor.plot(mols$X2.2.16.1w9_[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6),  pt.bg=c(2,2,2,3,3,3), legpl="bottomright", xlab="nmol PLFA g-1 TOC", ylab="Soil horizon", pch=rep(21:23,2), er.type="se")

hor.plot((mols$X2.2.16.1w9_+mols$X2.6.16.1w5_)[cond.abs]/samples$percent_C[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="bottomleft", xlab="nmol PLFA g-1 TOC", ylab="Soil horizon")

hor.plot(mols$X2.4.16.1w7[cond.abs]/samples$percent_C[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="bottomleft", xlab="nmol PLFA g-1 TOC", ylab="Soil horizon")

hor.plot(sum[cond.abs]*samples$weight[cond.abs]/samples$area[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="bottomleft", xlab="nmol PLFA m-2", ylab="Soil horizon")

hor.plot(samples$percent_C[cond.abs]*samples$weight[cond.abs]/samples$area[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="bottomleft", xlab="nmol PLFA m-2", ylab="Soil horizon")

hor.plot(mols$X3.6.18.1w9c[cond.abs]*samples$weight[cond.abs]/samples$area[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="topright", xlab="nmol PLFA m-2", ylab="Soil horizon")


hor.plot(mols$X3.11.18.2w6.9[cond.abs]*samples$weight[cond.abs]/samples$area[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="topright", xlab="nmol PLFA m-2", ylab="Soil horizon")

hor.plot(mols$X3.20.18.3w3.6.9[cond.abs]*samples$weight[cond.abs]/samples$area[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="bottomright", xlab="nmol PLFA m-2", ylab="Soil horizon")


hor.plot(mols$X3.3.18.1w9t._[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="topright", xlab="nmol 1:8w9c g d.w.", ylab="Soil horizon")

hor.plot(rowSums(mols[cond.abs,fames$group=="fungi"])/samples$percent_C[cond.abs], hor=samples$Horizon[cond.abs],  c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="topright", xlab="nmol PLFA m-2", ylab="Soil horizon")

hor.plot(rowSums(mols[cond.abs,fames$group=="Gpos"])/samples$percent_C[cond.abs], hor=samples$Horizon[cond.abs],  c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="topright", xlab="nmol PLFA m-2", ylab="Soil horizon", er.type="se")

hor.plot(rowSums(mols[cond.abs,fames$group=="Gneg"])/samples$percent_C[cond.abs], hor=samples$Horizon[cond.abs],  c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="topright", xlab="nmol PLFA m-2", ylab="Soil horizon", er.type="se")

pdf("ratios.pdf", height=4)
par(mfrow=c(1,2))

hor.plot(  rowSums(mols[cond.rel,fames$group=="fungi"])/
             (rowSums(mols[cond.rel,fames$group=="Gneg"])+
                rowSums(mols[cond.rel,fames$group=="Gpos"])), 
           hor=samples$Horizon[cond.rel],  c("L","F","H","B"), fac=samples$Site[cond.rel],
           lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="bottomright", 
           xlab="Fungi:Bacteria (mol:mol, SE)", ylab="Soil horizon", er.type="se", lwd=2)

hor.plot(  rowSums(mols[cond.rel,fames$group=="Gpos"])/rowSums(mols[cond.rel,fames$group=="Gneg"]),
           hor=samples$Horizon[cond.rel],  c("L","F","H","B"), fac=samples$Site[cond.rel],
           lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="bottomright", 
           xlab="G+/G- (mol:mol, SE)", ylab="Soil horizon", er.type="se", lwd=2)
dev.off()
  hor.plot(  rowSums(mols[cond.rel,fames$group=="fungi"])/
             (rowSums(mols[cond.rel,fames$group=="Gneg"])+
                rowSums(mols[cond.rel,fames$group=="Gpos"])), 
           hor=samples$Horizon[cond.rel],  c("L","F","H","B"), fac=samples$Site[cond.rel],
           lty=rep(1,6),  pt.bg=c(2,2,2,3,3,3), legpl="bottomright", pch=rep(21:23,2),
           xlab="nmol PLFA m-2", ylab="Soil horizon", er.type="se", lwd=2)

hor.plot(  rowSums(mols[cond.rel,fames$group=="Gpos"])/rowSums(mols[cond.rel,fames$group=="Gneg"]),
           hor=samples$Horizon[cond.rel],  c("L","F","H","B"), fac=samples$Site[cond.rel],
           lty=rep(1,6),  pt.bg=c(2,2,2,3,3,3), legpl="bottomright", pch=rep(21:23,2),
           xlab="G+/G- (mol:mol, SE)", ylab="Soil horizon", er.type="se", lwd=2)

hor.plot(sums[cond.abs], hor=samples$Horizon[cond.abs],  
         c("L","F","H","B"), fac=samples$Site[cond.abs],
           lty=rep(1,6),  pt.bg=c(2,2,2,3,3,3), legpl="bottomright", pch=rep(21:23,2),
           xlab="nmol PLFA m-2", ylab="Soil horizon", er.type="se", lwd=2)

unsaturation16.0<-(rel$X2.1.16.1w_+rel$X2.2.16.1w9_+rel$X2.3.16.1w_+rel$X2.4.16.1w7+rel$X2.6.16.1w5_)/rel$X16.0

  unsaturation18.0<-
  (rel$X3.3.18.1w9t._+rel$X3.6.18.1w9c+rel$X3.7.18.1w7c)/rel$X3.1.18.0

hor.plot(unsaturation16.0[cond.rel],
           hor=samples$Horizon[cond.rel],  c("L","F","H","B"), fac=samples$Site[cond.rel],
           lty=rep(1,6),  pt.bg=c(2,2,2,3,3,3), legpl="bottomright", pch=rep(21:23,2),
           xlab="nmol PLFA m-2", ylab="Soil horizon", er.type="se", lwd=2)

hor.plot(unsaturation18.0[cond.rel],
         hor=samples$Horizon[cond.rel],  c("L","F","H","B"), fac=samples$Site[cond.rel],
         lty=rep(1,6),  pt.bg=c(2,2,2,3,3,3), legpl="bottomright", pch=rep(21:23,2),
         xlab="nmol PLFA m-2", ylab="Soil horizon", er.type="se", lwd=2)

hor.plot(unsaturation18.0[cond.rel],
         hor=samples$Horizon[cond.rel],  c("L","F","H","B"), fac=samples$Site[cond.rel],
         lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="bottomright", 
         xlab="nmol PLFA m-2", ylab="Soil horizon", er.type="se", lwd=2)

hor.plot(unsaturation16.0[cond.rel],
hor=samples$Horizon[cond.rel],  c("L","F","H","B"), fac=samples$Site[cond.rel],
lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="bottomright", 
xlab="nmol PLFA m-2", ylab="Soil horizon", er.type="se", lwd=2)

pdf("physiological_indices.pdf", width=10, height=6)
par(mfrow=c(2,3))
hor.plot((rel$i.15.0/rel$ai.15.0)[cond.rel],
         hor=samples$Horizon[cond.rel],  c("L","F","H","B"), fac=samples$Site[cond.rel],
         lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="bottomleft", 
         xlab="i15:0/ai15:0 (mol:mol, SE)", ylab="Soil horizon", er.type="se", lwd=2)

hor.plot(((rel$X3.14.cy.19.0+rel$X2.15.cy17.0)/(rel$X2.4.16.1w7+rel$X3.6.18.1w9c))[cond.rel],
         hor=samples$Horizon[cond.rel],  c("L","F","H","B"), fac=samples$Site[cond.rel],
         lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="topright", 
         xlab="cyclopropyl:precursor (mol:mol, SE)", ylab="Soil horizon", er.type="se", lwd=2)

hor.plot(unsaturation16.0[cond.rel],
         hor=samples$Horizon[cond.rel],  c("L","F","H","B"), fac=samples$Site[cond.rel],
         lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="bottomleft", 
         xlab="sum(16:1)/16:0 (mol:mol, SE)", ylab="Soil horizon", er.type="se", lwd=2)

hor.plot((rel$X3.6.18.1w9c/(rel$X3.6.18.1w9c+rel$X3.11.18.2w6.9))[cond.rel],
         hor=samples$Horizon[cond.rel],  c("L","F","H","B"), fac=samples$Site[cond.rel],
         lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="bottomright", 
         xlab="18:1w9c/(18:1w9c + 18:2w6,9), mol:mol, SE", ylab="Soil horizon", er.type="se", lwd=2)

hor.plot((
  rel$X3.11.18.2w6.9/(rel$X3.11.18.2w6.9+rel$X3.20.18.3w3.6.9)
          )[cond.rel],
         hor=samples$Horizon[cond.rel],  c("L","F","H","B"), fac=samples$Site[cond.rel],
         lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="bottomright", 
         xlab="18:2w6,9/(18:2w6,9 + 18:3w3,6,9), mol:mol, SE", ylab="Soil horizon", er.type="se", lwd=2)
dev.off()

hor.plot((rel$X3.6.18.1w9c+rel$X3.11.18.2w6.9+rel$X3.20.18.3w3.6.9)[cond.rel],
         hor=samples$Horizon[cond.rel],  c("L","F","H","B"), fac=samples$Site[cond.rel],
         lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="bottomright", 
         xlab="nmol PLFA m-2", ylab="Soil horizon", er.type="se", lwd=2)
hor.plot(rowSums(mols[cond.rel,fames$group=="fungi"])/
             (rowSums(mols[cond.rel,fames$group=="Gneg"])+
                rowSums(mols[cond.rel,fames$group=="Gpos"])), 
           hor=samples$Horizon[cond.rel],  c("L","F","H","B"), fac=samples$Site[cond.rel],
           lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="bottomright", 
           xlab="nmol PLFA m-2", ylab="Soil horizon", er.type="se", lwd=2)

hor.plot(mols$X3.14.cy.19.0[cond.rel]/mols$X3.6.18.1w9c[cond.rel],
         hor=samples$Horizon[cond.rel],  c("L","F","H","B"), fac=samples$Site[cond.rel],
         lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="bottomright", 
         xlab="nmol PLFA m-2", ylab="Soil horizon", er.type="se", lwd=2)

hor.plot(mols$X3.14.cy.19.0[cond.rel]/mols$X3.6.18.1w9c[cond.rel],
         hor=samples$Horizon[cond.rel],  c("L","F","H","B"), fac=samples$Site[cond.rel],
         lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="bottomright", 
         xlab="nmol PLFA m-2", ylab="Soil horizon", er.type="se", lwd=2)

hor.plot(mols$X2.15.cy17.0[cond.rel]/mols$X2.4.16.1w7[cond.rel],
         hor=samples$Horizon[cond.rel],  c("L","F","H","B"), fac=samples$Site[cond.rel],
         lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="bottomright", 
         xlab="nmol PLFA m-2", ylab="Soil horizon", er.type="se", lwd=2)

hor.plot((rel$X3.23.20.0[cond.rel]+rel$X4.6.22.0[cond.rel]+rel$X4.11.24.0[cond.rel]),
         hor=samples$Horizon[cond.rel],  c("L","F","H","B"), fac=samples$Site[cond.rel],
         lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="bottomright", 
         xlab="nmol PLFA m-2", ylab="Soil horizon", er.type="se", lwd=2)

hor.plot(
  (
    rel$X4.2.20.2_[cond.rel]+rel$X4.3.20.3w6.9.15[cond.rel]+rel$X4.5.20.4w6.9.12.15[cond.rel]+rel$X4.7.20.5w3.6.9.12.15[cond.rel]
    ),
         hor=samples$Horizon[cond.rel],  c("L","F","H","B"), fac=samples$Site[cond.rel],
         lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="bottomright", 
         xlab="nmol PLFA m-2", ylab="Soil horizon", er.type="se", lwd=2)

rel

hor.plot( 
  rowSums(mols[cond.rel,fames$group=="fungi"])/
    (rowSums(mols[cond.rel,fames$group=="Gneg"])+
       rowSums(mols[cond.rel,fames$group=="Gpos"])),
  samples$percent_C[cond.rel], hor=samples$Horizon[cond.rel],  hc("L","F","H","B"), 
  fac=samples$Site[cond.rel], lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, 
  legpl="topright", xlab="nmol PLFA m-2", ylab="Soil horizon", er.type="se")


col=c(1,2)

col<-  rep(col, ncol(means)      )
cond2<-is.na(means)==F & is.na(error)==F

nested<-samples$Region

nested1<-nested[cond1]
means -> means.old
error -> error.old
means.old
nest
nest<-1
for (i in 1:nrow(means.old))
{
  fac<-rownames(means.old)[i]
 tmp <- unique(nested1[which(fac1==fac)])

 nest[i]<-as.character(factor(tmp, levels=unique(nested1)))
}

nest<-as.factor(nest)
cond2
nested1 <- factor(nested1)
means<-data.frame(matrix(nrow= length( levels(nest)), ncol=ncol(means.old)))

rownames(means)<-levels(nest)
colnames(means)<-colnames(means.old)
error <- means
er.type.nested="se"

i<-1
for (i in 1:nrow(means))
{
  means[i,T] <- colMeans( means.old[nest==rownames(means)[i],T])
  for (j in 1:ncol(means)) 
  {
    if(er.type.nested=="se") 
    {
      error[i,j]<-stderr(means.old[nest==rownames(means)[i],j])
    }
    if(er.type.nested=="sd") 
    {
      error[i,j]<-sd(means.old[nest==rownames(means)[i],j])
    }
    if(er.type.nested=="ci") 
    {
      error[i,j]<-CI(means.old[nest==rownames(means)[i],j])
    }
  }
}

means
error
pt.bg<-1:2
cex.pt<-1

cond2<-is.na(means)==F & is.na(error)==F

plot(t(means[,T]), rep(ncol(means):1,nrow(means)), yaxt="n", 
     xlim=c(0,1.2*max(means[cond2]+error[cond2])), tck=0.01, type="n")
i<-1
for(i in 1:nrow(means)) {
  tmp.mean<-as.numeric(as.vector(means[i,T]))
  tmp.er<-as.numeric(as.vector(error[i,T]))
  plotCI(tmp.mean, ncol(means):1, uiw=tmp.er, err="x", pch=pch[i], lty=lty, col=col[i], pt.bg=pt.bg[i], add=T, gap=0, cex=cex.pt)
}
lty<-1

sd(x, na.rm=T)
  x[is.na(x)==F]
       )
x
<-error.old[,2]

var(x[is.na(x==F)]
    
stderr(c(2047.646, 1336.905, NA))
       104.71059
       3   1361.735 2555.529)

summary(sum)
hor.plot(sum[cond.rel], hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=rep(1), nested=samples$Region)
hor.plot(sum[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=1, nested=samples$Region[cond.abs], er.type.nested="sd", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="Sum of PLFA (nmol g-1 d.w.)", ylab="Soil horizon")

hor.plot(mols$X3.11.18.2w6.9[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=1, nested=samples$Region[cond.abs], er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="18:2w6,9 (nmol g-1 d.w.)", ylab="Soil horizon")

hor.plot(
  rowSums(mols[cond.rel, fames$group=="fungi"])/
    (rowSums(mols[cond.rel, fames$group=="Gpos"]+rowSums(mols[cond.rel, fames$group=="Gneg"]))), 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="Fungi:Bacteria (SE)", ylab="Soil horizon")

hor.plot(
  mols$i.15.0[cond.rel]/mols$ai.15.0[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="Fungi:Bacteria (SE)", ylab="Soil horizon")

hor.plot(
  mols$X2.15.cy17.0[cond.rel]/mols$X2.4.16.1w7[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="Fungi:Bacteria (SE)", ylab="Soil horizon")

hor.plot(
  mols$X3.14.cy.19.0[cond.rel]/mols$X3.6.18.1w9.t_[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="Fungi:Bacteria (SE)", ylab="Soil horizon")

hor.plot(
  mols$X3.14.cy.19.0[cond.rel]/(mols$X3.7.18.1w9c_[cond.rel]+mols$X3.6.18.1w9.t_[cond.rel]), 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="Fungi:Bacteria (SE)", ylab="Soil horizon")

hor.plot(
  mols$X3.7.18.1w9c_[cond.rel]/mols$X3.6.18.1w9.t_[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="18:1c/t (SE)", ylab="Soil horizon")


hor.plot(
  rel$X4.3.20.3w6.9.15[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="18:1c/t (SE)", ylab="Soil horizon")

hor.plot(
  rel$X3.9.18.2w_.i19.0_[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="18:1c/t (SE)", ylab="Soil horizon")

hor.plot(
  rel$ai.16.0_[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="18:1c/t (SE)", ylab="Soil horizon")

hor.plot(
  rel$X2.13.i18.0.or.10.Me18.0[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="18:1c/t (SE)", ylab="Soil horizon")

hor.plot(
  rel$X2.15.cy17.0[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="18:1c/t (SE)", ylab="Soil horizon")

hor.plot(
  rel$X3.14.cy.19.0[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="18:1c/t (SE)", ylab="Soil horizon")

hor.plot(
  rel$X2.6.16.1w5_[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="18:1c/t (SE)", ylab="Soil horizon")

hor.plot(
  rel$X3.8.18.2w9.12[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="18:1c/t (SE)", ylab="Soil horizon")

pdf("longchain.pdf")
par(mfrow=c(2,3))

hor.plot(
  rel$X3.23.20.0[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="20:0 (mol%, SE)", ylab="Soil horizon",legsize=0.8)

hor.plot(
  rel$X4.6.22.0[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="22:0 (mol%, SE)", ylab="Soil horizon",legsize=0.8)

hor.plot(
  rel$X4.11.24.0[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="24:0 (mol%, SE)", ylab="Soil horizon",legsize=0.8)

hor.plot(
  rel$X4.3.20.3w6.9.15[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="20:3 (mol%, SE)", ylab="Soil horizon",legsize=0.8)

hor.plot(
  rel$X4.5.20.4w6.9.12.15[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="20:4 (mol%, SE)", ylab="Soil horizon",legsize=0.8)

hor.plot(
  rel$X4.7.20.5w3.6.9.12.15[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="20:5 (mol%, SE)", ylab="Soil horizon",legsize=0.8)
dev.off()

pdf("16:1s.pdf")
par(mfrow=c(2,3))

hor.plot(
  rel$X2.1.16.1w_[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomleft", lwd=2, xlab="16:1 (2-1) (mol%, SE)", ylab="Soil horizon",legsize=0.8)

hor.plot(
  rel$X2.2.16.1w9_[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomleft", lwd=2, xlab="16:1w9 (mol%, SE)", ylab="Soil horizon",legsize=0.8)

hor.plot(
  rel$X2.3.16.1w_[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="topright", lwd=2, xlab="16:1 (2-3) (mol%, SE)", ylab="Soil horizon",legsize=0.8)

hor.plot(
  rel$X2.4a.16.1w_[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="topright", lwd=2, xlab="16:1 (2-4a) (mol%, SE)", ylab="Soil horizon",legsize=0.8)

hor.plot(
  rel$X2.4.16.1w7[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomleft", lwd=2, xlab="16:1w7 (mol%, SE)", ylab="Soil horizon",legsize=0.8)

hor.plot(
  rel$X2.6.16.1w5_[cond.rel], 
dev.off()

pdf("isoanteisos.pdf")
par(mfrow=c(2,3))

hor.plot(
  rel$i.15.0[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomleft", lwd=2, xlab="i-15:0 (mol%, SE)", ylab="Soil horizon",legsize=0.8)

hor.plot(
  rel$ai.15.0[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomleft", lwd=2, xlab="ai-15:0 (mol%, SE)", ylab="Soil horizon",legsize=0.8)

hor.plot(
  rel$i.16.0[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomleft", lwd=2, xlab="i-16:0 (2-3) (mol%, SE)", ylab="Soil horizon",legsize=0.8)

hor.plot(
  rel$ai.16.0_[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="topleft", lwd=2, xlab="ai 16:0 (2-4a) (mol%, SE)", ylab="Soil horizon",legsize=0.8)

hor.plot(
  rel$X2.5.i.17.0_[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomleft", lwd=2, xlab="i17:0 (mol%, SE)", ylab="Soil horizon",legsize=0.8)

hor.plot(
  rel$X2.7._ai.17.0_[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomleft", lwd=2, xlab="ai(?)17:0 (mol%, SE)", ylab="Soil horizon",legsize=0.8)
dev.off()

  
pdf("sufa.pdf")
par(mfrow=c(2,3))
  
  hor.plot(
    rel$X14.0[cond.rel], 
    hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="topleft", lwd=2, xlab="14:0 (mol%, SE)", ylab="Soil horizon",legsize=0.8)
  
hor.plot(
    rel$n.15.0[cond.rel], 
    hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="topleft", lwd=2, xlab="15:0 (mol%, SE)", ylab="Soil horizon",legsize=0.8)
  
  hor.plot(
    rel$X16.0[cond.rel], 
    hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="topleft", lwd=2, xlab="16:0 (mol%, SE)", ylab="Soil horizon",legsize=0.8)
  
hor.plot(
    rel$X3.1.18.0[cond.rel], 
    hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="topleft", lwd=2, xlab="18:0 (mol%, SE)", ylab="Soil horizon",legsize=0.8)
  
  hor.plot(
    rel$X3.6.18.1w9.t_[cond.rel], 
    hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="topleft", lwd=2, xlab="18:1w9t (3-6) (mol%, SE)", ylab="Soil horizon",legsize=0.8)
  
  hor.plot(
    rel$X3.7.18.1w9c_[cond.rel], 
    hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="topright", lwd=2, xlab="18:1w9c (3-7) (mol%, SE)", ylab="Soil horizon",legsize=0.8)
  dev.off()
  

sums <- rowSums(mols)
cond <- is.na(mols$X3.11.18.2w6.9)==F
cond
mols$

f2b<-rowSums(mols[T, fames$X0=="fungi"])/
  rowSums(mols[,fames$X0=="Gpos"|fames$X0=="Gneg"])

hor.plot(mols$X3.11.18.2w6.9[cond], samples$Horizon[cond], c("L","F","H","B"), fac=samples$Site[cond], col=c(1,2), nested=samples$Region[cond], lwd=c(2,2), er.type="se", lty=c(1,2))
cond <- is.na(sums)==F
hor.plot(sums[cond1], samples$Horizon[cond1], c("L","F","H","B"), fac=samples$Region[cond1], col=c(1,2), lwd=2, er.type="se")

f2b<-
  rowSums(mols[T, fames$X0=="fungi"])/
  rowSums(mols[,fames$X0=="Gpos"|fames$X0=="Gneg"])
cond <- is.na(f2b)==F
hor.plot(f2b[cond], samples$Horizon[cond], c("L","F","H","B"), fac=samples$Region[cond], col=c(1,2), lwd=2, er.type="se", lty=c(1,1))



fames
mols
hor.plot()


head(mols)

samples$Horizon
<-
  t(tmp[1:19, 9:ncol(tmp)])
mols <- 
  data.frame(t(tmp[22:nrow(tmp), 9:ncol(tmp)]))
summary(mols)
<-
  as.numeric(mols[,3])


lda.all<-lda(rel[cond.rel&cond.ref,T], samples$Region[cond.rel&cond.ref])
lda.l<-lda(rel[cond.rel&samples$Horizon=="L",T], samples$Region[cond.rel&samples$Horizon=="L"])
lda.f<-lda(rel[cond.rel&samples$Horizon=="F",T], samples$Region[cond.rel&samples$Horizon=="F"])
lda.h<-lda(rel[cond.rel&samples$Horizon=="H",T], samples$Region[cond.rel&samples$Horizon=="H"])
lda.b<-lda(rel[cond.rel&samples$Horizon=="B",T], samples$Region[cond.rel&samples$Horizon=="B"])


plot(lda.b)
lda.values<-data.frame( lda.all$scaling, lda.l$scaling, lda.f$scaling, lda.h$scaling)
lda.values
dev.off()