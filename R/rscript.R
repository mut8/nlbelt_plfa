source("/home/lluc/R_functions/functions.R")
source("~/R_functions/plfaplotfunction.R")


setwd("/home/lluc/Dropbox/MUN/phd/results/NL-BELT/Batches21-25combined/R")

mols<-  data.frame(read.csv("mols.csv"))
fames<-read.csv("fames.csv")
samples<-  
  read.csv("samples.csv")

cond.rel<-samples$trusted_rel==T
cond.abs<-samples$trusted_abs==T
cond.ref <- samples$Region!="Ref"
cond1 <-cond.abs & cond.ref

sum<-rowSums(mols)

rel<-mols/rowSums(mols)


ord<-rda(rel[cond.rel,T])
plot(ord)
ord.plot(ord,site.sep1=samples$Horizon[cond.rel], site.sep2=samples$Site[cond.rel], pch=c(21,21,21,22,22,22,24), pt.bg=1:6, spe.label.type="text", spe.labels=fames$FAME, cex.leg=.5)

fac<-samples$Site
horlev<-c("L","F","H","B")
hor<-samples$Horizon
cond1<-is.element(hor,horlev)
hor1<-factor(hor[cond1], ordered=T, levels=horlev)
fac1<-factor(fac[cond1], levels=unique(fac[cond1]))

var1<-sum[cond1]

means<-tapply(var1, list(fac1,hor1), mean)
error<-tapply(var1, list(fac1,hor1), stderr) 

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

pdf("longchain.pdf", height=11)
par(mfrow=c(4,3))

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
  rel$X2.6.16.1w5_[cond.rel], hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomleft", lwd=2, xlab="16:1w5 (mol%, SE)", ylab="Soil horizon",legsize=0.8)
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