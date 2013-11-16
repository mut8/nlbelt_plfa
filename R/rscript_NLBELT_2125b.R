source("/home/lluc/R_functions/functions.R")
source("~/R_functions/plfaplotfunction.R")


setwd("/home/lluc/Copy/MUN/phd/results/NL-BELT/Batches21-25combined/R")

mols<-  data.frame(read.csv("mols.csv"))
fames<-    read.csv("fames.csv")
samples<-      read.csv("samples.csv")
resp245<-read.csv("respiration.csv")
resp<-read.csv("allresp.csv")
d13C<-    read.csv("d13C.csv")
amp44<-read.csv("amp44.csv")
nb_C = c(15,15,15,16,16,16,17,18,18,18,19,18,20,NA,22)


#minimum amplitude (44) to accept d13C value
min_amp44 <- 200
#d13C of methanol used for saponification
d13C_meoh <- -50

cond.rel<-samples$trusted_rel==T
cond.abs<-samples$trusted_abs==T
cond.ref <- samples$Region!="Ref"
cond1 <-cond.abs & cond.ref
cond2 <-cond.rel& cond.ref
cond.irms <- samples$trusted_irms==T
cond3 <- cond.irms&cond.ref

sum<-rowSums(mols)

rel<-100*mols/rowSums(mols)
rel.colmeans<-rel
for (i in 1:ncol(rel.colmeans)) {
  rel.colmeans[,i]<-rel[,i]/mean(rel[,i])
}


cy17prec<-(rel$cy17.0/rel$X2.4.16.1w7)
cy19prec<-((rel$cy.19.0+rel$cy.19.0.1)/rel$X18.1w9c)
cyprec<-((rel$cy17.0+rel$cy.19.0+rel$cy.19.0.1)/(rel$X2.4.16.1w7+rel$X18.1w9c))

F2B <-rowSums(rel[,fames$microbial.group=="fungi"])/
    rowSums(rel[,fames$microbial.group=="Gpos"|fames$microbial.group=="Gneg"|fames$microbial.group=="bac"] )
     
F2B.18.2<- rel$X18.2w6.9/ rowSums(rel[,fames$microbial.group=="Gpos"|fames$microbial.group=="Gneg"|fames$microbial.group=="bac"] )

euk.conc <-rowSums(mols[,fames$microbial.group=="fungi"|fames$microbial.group=="euk"])
euk.rel <-rowSums(rel[,fames$microbial.group=="fungi"|fames$microbial.group=="euk"])

allbac.conc <-rowSums(mols[,fames$microbial.group=="Gpos"|fames$microbial.group=="Gneg"|fames$microbial.group=="bac"|fames$microbial.group=="act"] )
allbac.rel <-rowSums(rel[,fames$microbial.group=="Gpos"|fames$microbial.group=="Gneg"|fames$microbial.group=="bac"|fames$microbial.group=="act"] )

Euk2bac <- euk.rel/allbac.rel

fungi3.rel<-rowSums(rel[,fames$microbial.group=="fungi"])
fungi2.rel<-rel$X18.2w6.9 + rel$X18.3w3.6.9
fungi1.rel<- rel$X18.2w6.9

fungi3.conc<-rowSums(mols[,fames$microbial.group=="fungi"])
fungi2.conc<-   mols$X18.2w6.9 + mols$X18.3w3.6.9
fungi1.conc<- mols$X18.2w6.9

monouns.conc<-mols$X15.1+mols$X2.1.16.1a+mols$X2.3.16.1b+mols$X2.4.16.1w7+mols$X16.1c+mols$X18.1b+mols$X18.1c
monouns.rel<-rel$X15.1+rel$X2.1.16.1a+rel$X2.3.16.1b+rel$X2.4.16.1w7+rel$X16.1c+rel$X18.1b+rel$X18.1c


pufa.conc<-rowSums(mols[,fames$saturation=="PUFA"])
pufa.rel<-rowSums(rel[,fames$saturation=="PUFA"])

longeven.rel <-rel$X20.0+rel$X22.0+rel$X24.0
longpufa.rel <- rel$X20.2+rel$X20.3a+rel$X20.3b+rel$X20.4w6.9.12.15+rel$X20.5w3.6.9.12.15+rel$X22.6w3.6.9.12.15.18
ia.rel <- rel$i.15.0+rel$ai.15.0
act.rel<-rowSums(rel[,fames$microbial.group=="act"])

longeven.conc <-mols$X20.0+mols$X22.0+mols$X24.0
longpufa.conc <- mols$X20.2+mols$X20.3a+mols$X20.3b+mols$X20.4w6.9.12.15+mols$X20.5w3.6.9.12.15
ia.conc <- mols$i.15.0+mols$ai.15.0
act.conc<-rowSums(mols[,fames$microbial.group=="act"])

#remove d13C values with amp44 < min
d13C [amp44 < min_amp44] <- NA
tmp <-   data.frame(matrix(ncol=ncol(d13C)-2, nrow=length(unique(d13C[,2]))))
colnames(tmp) <- colnames(d13C)[3:ncol(d13C)]

rownames(tmp) <-    names(tapply(d13C[,3], d13C$X.1, mean))

#calculate mean d13C per sample over multiple injections
for (i in 3:ncol(d13C)) {
  tmp[,i-2] <-tapply(d13C[,i], d13C$X.1, mean)
}



#correct d13C for methanol
meoh_cor_d13C <- tmp
for (i in 1:ncol(tmp))
meoh_cor_d13C[,i] <-  (tmp[,i]*(nb_C[i]+1)-d13C_meoh)/nb_C[i]

#copy d13C for internal standard w/o MeOH correction
meoh_cor_d13C$X20.0EE <- tmp$X20.0EE

#organize d13C values so that they are in the same order as sample info.
d13C_org <- data.frame(matrix(nrow=nrow(samples), ncol=ncol(meoh_cor_d13C)))

colnames(d13C_org) <- colnames(meoh_cor_d13C)


for (i in 1:nrow(samples)) 
  { if (is.na(samples$sample_nr_irms[i])==F)
  d13C_org[i,T] <- meoh_cor_d13C[which(samples$sample_nr_irms[i]==rownames(meoh_cor_d13C)),T]
}

#calculate mass averaged means

# (d13C_org$i15.0*rel$i.15.0+
#   d13C_org$a15.0*rel$i.15.0+
#   d13C_org$X16.0a*rel$X16.0+
#   d13C_org$X16.1w7*rel$X2.4.16.1w7+
#   d13C_org$cy17.0*rel$cy17.0+
#   d13C_org$X18.1w9*rel$X18.1w9c+
#   d13C_org$X18.1w7*rel$X18.1b+
#   d13C_org$X18.2w6*rel$X18.2w6.9+
#   d13C_org$cy19.0*(rel$cy.19.0+rel$cy.19.0.1)+ 
#   d13C_org$X18.3w3*rel$X18.3w3.6.9+
#   d13C_org$X20.0*rel$X20.0+
#   d13C_org$X22.0*rel$X22.0)/
#   (rel$i.15.0+rel$i.15.0+rel$X16.0+rel$X2.4.16.1w7+rel$cy17.0+rel$X18.1w9c+rel$X18.1b+rel$X18.2w6.9+rel$X18.3w3.6.9+rel$X20.0+rel$X22.0)


for (i in 1:ncol(d13C_org)){
  if (i == 1) { plfas_reorg <- data.frame(d13C=d13C_org[cond3,i]-samples$d13C[cond3], Region=samples$Region[cond3], Horizon=samples$Horizon[cond3],Site=samples$Site[cond3],Sample_ID=samples$sample_nr_irms[cond3], plfa=colnames(d13C_org)[i])}
  else { plfas_reorg <- rbind (plfas_reorg, data.frame(d13C=d13C_org[cond3,i]-samples$d13C[cond3], Region=samples$Region[cond3], Horizon=samples$Horizon[cond3], Site=samples$Site[cond3],Sample_ID=samples$sample_nr_irms[cond3],plfa=colnames(d13C_org)[i]))}
}
plfas_reorg <- plfas_reorg[plfas_reorg$plfa!="X20.0EE",T]
plfas_reorg <- plfas_reorg[plfas_reorg$plfa!="X16.0",T]
plfas_reorg <- plfas_reorg[plfas_reorg$plfa!="X15.0",T]

c <- is.na(plfas_reorg$d13C)==F

anova(lm(d13C ~ Horizon * plfa, plfas_reorg))
summary(lm(d13C ~ Horizon * plfa, plfas_reorg))
anova(lm(d13C ~ Horizon * plfa + Region, plfas_reorg))
summary(lm(d13C ~ Horizon * plfa + Region, plfas_reorg))
anova(lm(d13C ~ Horizon * plfa + Region*Horizon + Region*plfa, plfas_reorg))
summary(lm(d13C ~ Horizon * plfa + Region*Horizon + Region*plfa, plfas_reorg))
anova(lm(d13C ~ Horizon * plfa * Region, plfas_reorg))
summary(lm(d13C ~ Horizon * plfa *Region, plfas_reorg))
anova(lm(d13C ~ Horizon * plfa + Region*plfa, plfas_reorg))
summary(lm(d13C ~ Horizon * plfa + Region*plfa, plfas_reorg))

l1 <- lm(d13C ~ Horizon * plfa, plfas_reorg)
l2 <- lm(d13C ~ Horizon * plfa + Region, plfas_reorg)
l3 <- lm(d13C ~ Horizon * plfa + Region*Horizon + Region*plfa, plfas_reorg)
l4 <-lm(d13C ~ Horizon * plfa * Region, plfas_reorg)
l5 <- lm(d13C ~ Horizon * plfa + Region*plfa, plfas_reorg)

anova(l1,l2)
anova(l2,l3)
anova(l2,l4)
anova(l3,l5)

summary(cond1)
summary(cond2)
summary(cond3)

hor.plot(x1, hor=plfas_reorg$Horizon[c], horlev=c("L","F","H","B"), fac=plfas_reorg$Region[c], xlim=c(-1,1), er.type="ci", xlab="residual (permil, 95% CI)", ylab="soil horizon", pch=21:22, pt.bg=2:3)
abline(v=0, col="grey")
abline(v=1, col="grey")
abline(v=-1, col="grey")

lm1<-lm(plfas_reorg$d13C[c] ~ plfas_reorg$Horizon[c] * plfas_reorg$plfa[c])

pdf("residuals.pdf")
plot(lm1$residuals~plfas_reorg$Region[c], xlim=c(0.5,2.5), xlab="region", ylab="residual (d13C ~ horizon * plfa)")
densityplot(~x1, groups=plfas_reorg$Region[c])
densityplot(~x1 | plfas_reorg$plfa, groups=plfas_reorg$Region[c], col=c("blue","green"))
densityplot(~x1 | plfas_reorg$Horizon, groups=plfas_reorg$Region[c], col=c("blue","green"))
dev.off()


x1 <- lm1$residuals
anova(lm(x1~plfas_reorg$Region[c]))
summary(lm(x1~plfas_reorg$Region[c]))


anova(lm(x1~plfas_reorg$Region[c]*plfas_reorg$Horizon[c]))

pdf("hist_residuals_byregion.pdf", height=8, width=5)
par(mfrow=c(3,1))

x2 <- lm1$residuals[plfas_reorg$Region[c]=="GC"]
x3 <- lm1$residuals[plfas_reorg$Region[c]=="ER"]

t.test(x2,x3)
hist(x1, xlim=c(-3,3), prob=T)
curve(dnorm, col=2, mean=mean(x1), sd=sd(x1), add=T)
curve(
  dnorm(x, mean=mean(x1)+1, sd=sd(x1))
      , col=2, , add=T)
?curve
?dnorm

hist(x2, xlim=c(-3,3), prob=T)
curve(dnorm, col=2, mean=mean(x2), sd=sd(x2), add=T)

hist(x3, xlim=c(-3,3), prob=T)
curve(dnorm, col=2, mean=mean(x3), sd=sd(x3), add=T)

densityplot(~x1, groups=plfas_reorg$Region[c])
densityplot(~x1 | plfas_reorg$plfa, groups=plfas_reorg$Region[c], col=c("blue","green"))
densityplot(~x1 | plfas_reorg$Horizon, groups=plfas_reorg$Region[c], col=c("blue","green"))
?densityplot
histogram(~x1,
          panel = function(x1,...){
            panel.histogram(x1,...)
            panel.mathdensity(dmath = dnorm, col = "red",
                              args = list(mean=mean(x1),sd=sd(x1)))
          },
          type="density", tck=0.01)

histogram(~x2,
          panel = function(x2,...){
            panel.histogram(x2,...)
            panel.mathdensity(dmath = dnorm, col = "red",
                              args = list(mean=mean(x2),sd=sd(x2)))
          },
          type="density", tck=0.01)

histogram(~x3,
          panel = function(x3,...){
            panel.histogram(x3,...)
            panel.mathdensity(dmath = dnorm, col = "red",
                              args = list(mean=mean(x3),sd=sd(x3)))
          },
          type="density", tck=0.01)


dev.off()



library(lattice)

hist(x, prob=T)
curve(dnorm, col=2, mean=mean(x), sd=sd(x), add=T)




x <- lm1$residuals[plfas_reorg$Region[c]=="GC"]
histogram(~x,
          panel = function(x,...){
            panel.histogram(x,...)
            panel.mathdensity(dmath = dnorm, col = "red",
                              args = list(mean=mean(x),sd=sd(x)))
          },
          type="density")

x <- lm1$residuals[plfas_reorg$Region[c]=="GC"]
histogram(~x,
          panel = function(x,...){
            panel.histogram(x,...)
            panel.mathdensity(dmath = dnorm, col = "red",
                              args = list(mean=mean(x),sd=sd(x)))
          },
          type="density")




?density

SnowsPenultimateNormalityTest(lm1$residuals)
qqnorm(lm1$residuals)
qqnorm(lm1$residuals[plfas_reorg$Region[c]=="GC"])
qqnorm(lm1$residuals[plfas_reorg$Region[c]=="ER"])
?hist


summary(lm(plfas_reorg$d13C[c] ~ plfas_reorg$Horizon[c] * plfas_reorg$plfa[c]+plfas_reorg$Region[c]))


plf

summary(d13C[c])
################thats it... now some test plots of IRMS data...
dev.off()

pdf("comp_d13c.pdf", height=7, width=5)
par(mfrow=c(4,1), mar=c(2.1,5.1,0.1,2.1) )

hor.plot(d13C_org$X16.0a[cond.irms]-samples$d13C[cond.irms], hor=samples$Horizon[cond.irms], c("L","F","H","B"), fac=samples$Region[cond.irms], pt.bg=2:3, legpl="topright", xlab="d13C 16:0", ylab="Soil horizon", oneway=T, er.type.nested="sd", xlim=c(-8, 7), add=F, lty=2)

hor.plot(d13C_org$X18.1w9[cond.irms]-samples$d13C[cond.irms], hor=samples$Horizon[cond.irms], c("L","F","H","B"), fac=samples$Region[cond.irms], pt.bg=2:3, legpl="topright", xlab="d13C 16:0", ylab="Soil horizon", oneway=T, er.type.nested="sd", xlim=c(-10, 10), add=T)

legend("bottomright", lty=1:2, c("16:0", expression(paste("18:1", omega, "9" ))))

hor.plot(d13C_org$i15.0[cond.irms]-samples$d13C[cond.irms], hor=samples$Horizon[cond.irms], c("L","F","H","B"), fac=samples$Region[cond.irms], pt.bg=2:3, legpl="topright", xlab=expression(paste(Delta, ""^"13", "C" [PLFA - SOM])), ylab="Soil horizon", oneway=T, er.type.nested="sd", xlim=c(-8, 7))

hor.plot(d13C_org$cy17.0[cond.irms]-samples$d13C[cond.irms], hor=samples$Horizon[cond.irms], c("L","F","H","B"), fac=samples$Region[cond.irms], pt.bg=2:3, legpl="topright", xlab="d13C 16:0", ylab="Soil horizon", oneway=T, er.type.nested="sd", xlim=c(-10, 10), add=T, lty=2)

hor.plot(d13C_org$X18.2w6[cond.irms]-samples$d13C[cond.irms], hor=samples$Horizon[cond.irms], c("L","F","H","B"), fac=samples$Region[cond.irms], pt.bg=2:3, legpl="topright", xlab="d13C 16:0", ylab="Soil horizon", oneway=T, er.type.nested="sd", xlim=c(-10, 10), add=T, lty=4)

legend("bottomright", lty=c(1:2,4), c("i15:0","cy17:0", expression(paste("18:2", omega, "6,9" ))))

hor.plot(d13C_org$a15.0[cond.irms]-samples$d13C[cond.irms], hor=samples$Horizon[cond.irms], c("L","F","H","B"), fac=samples$Region[cond.irms], pt.bg=2:3, legpl="topright", xlab=expression(paste(Delta, ""^"13", "C" [PLFA - SOM])), ylab="Soil horizon", oneway=T, er.type.nested="sd", xlim=c(-8, 7))

hor.plot(d13C_org$cy19.0[cond.irms]-samples$d13C[cond.irms], hor=samples$Horizon[cond.irms], c("L","F","H","B"), fac=samples$Region[cond.irms], pt.bg=2:3, legpl="topright", xlab="d13C 16:0", ylab="Soil horizon", oneway=T, er.type.nested="sd", xlim=c(-10, 10), add=T, lty=2)

legend("bottomright", lty=c(1:2,4), c("ai15:0", "cy19:0"))

hor.plot(d13C_org$X18.1w7[cond.irms]-samples$d13C[cond.irms], hor=samples$Horizon[cond.irms], c("L","F","H","B"), fac=samples$Region[cond.irms], pt.bg=2:3, legpl="topright", xlab="d13C 16:0", ylab="Soil horizon", oneway=T, er.type.nested="sd", xlim=c(-8, 7), add=F, lty=1)
hor.plot(d13C_org$X16.1w7[cond.irms]-samples$d13C[cond.irms], hor=samples$Horizon[cond.irms], c("L","F","H","B"), fac=samples$Region[cond.irms], pt.bg=2:3, legpl="topright", xlab="d13C 16:0", ylab="Soil horizon", oneway=T, er.type.nested="sd", xlim=c(-8, 7), add=T, lty=2)
legend("bottomright", lty=c(1:2,4), c(expression(paste("18:1", omega, "7" )),expression(paste("16:1", omega, "7" ))))

dev.off()

hor.plot(d13C_org$X16.1w7[cond.irms]-samples$d13C[cond.irms], hor=samples$Horizon[cond.irms], c("L","F","H","B"), fac=samples$Region[cond.irms], pt.bg=2:3, legpl="topright", xlab="d13C 16:0", ylab="Soil horizon", oneway=T, er.type.nested="sd", xlim=c(-10, 10), add=T)

hor.plot(d13C_org$X16.0[cond.irms]-samples$d13C[cond.irms], hor=samples$Horizon[cond.irms], c("L","F","H","B"), fac=samples$Region[cond.irms], pt.bg=2:3, legpl="topright", xlab="d13C 16:0", ylab="Soil horizon", oneway=T, er.type.nested="sd", xlim=c(-10, 10))

hor.plot(d13C_org$X16.0a[cond.irms]-samples$d13C[cond.irms], hor=samples$Horizon[cond.irms], c("L","F","H","B"), fac=samples$Region[cond.irms], pt.bg=2:3, legpl="topright", xlab="d13C 16:0", ylab="Soil horizon", oneway=T, er.type.nested="sd", xlim=c(-10, 10), add=F, lty=2)

hor.plot(d13C_org$i15.0[cond.irms]-samples$d13C[cond.irms], hor=samples$Horizon[cond.irms], c("L","F","H","B"), fac=samples$Region[cond.irms], pt.bg=2:3, legpl="topright", xlab="d13C 16:0", ylab="Soil horizon", oneway=T, er.type.nested="sd", xlim=c(-10, 10))

hor.plot(d13C_org$a15.0[cond.irms]-samples$d13C[cond.irms], hor=samples$Horizon[cond.irms], c("L","F","H","B"), fac=samples$Region[cond.irms], pt.bg=2:3, legpl="topright", xlab="d13C 16:0", ylab="Soil horizon", oneway=T, er.type.nested="sd", xlim=c(-10, 10), add=F)

hor.plot(d13C_org$X16.1w7[cond.irms]-samples$d13C[cond.irms], hor=samples$Horizon[cond.irms], c("L","F","H","B"), fac=samples$Region[cond.irms], pt.bg=2:3, legpl="topright", xlab="d13C 16:0", ylab="Soil horizon", oneway=T, er.type.nested="sd", xlim=c(-10, 10), add=F)

hor.plot(d13C_org$cy17.0[cond.irms]-samples$d13C[cond.irms], hor=samples$Horizon[cond.irms], c("L","F","H","B"), fac=samples$Region[cond.irms], pt.bg=2:3, legpl="topright", xlab="d13C 16:0", ylab="Soil horizon", oneway=T, er.type.nested="sd", xlim=c(-10, 10), add=F)


hor.plot(d13C_org$X18.1w7[cond.irms]-samples$d13C[cond.irms], hor=samples$Horizon[cond.irms], c("L","F","H","B"), fac=samples$Region[cond.irms], pt.bg=2:3, legpl="topright", xlab="d13C 16:0", ylab="Soil horizon", oneway=T, er.type.nested="sd", xlim=c(-8, 7), add=F)



hor.plot(d13C_org$cy19.0[cond.irms]-samples$d13C[cond.irms], hor=samples$Horizon[cond.irms], c("L","F","H","B"), fac=samples$Region[cond.irms], pt.bg=2:3, legpl="topright", xlab="d13C 16:0", ylab="Soil horizon", oneway=T, er.type.nested="sd", xlim=c(-8, 7), add=F)

hor.plot(d13C_org$X18.3w3[cond.irms]-samples$d13C[cond.irms], hor=samples$Horizon[cond.irms], c("L","F","H","B"), fac=samples$Region[cond.irms], pt.bg=2:3, legpl="topright", xlab="d13C 16:0", ylab="Soil horizon", oneway=T, er.type.nested="sd", xlim=c(-8, 7), add=F)

hor.plot(d13C_org$X20.0[cond.irms]-samples$d13C[cond.irms], hor=samples$Horizon[cond.irms], c("L","F","H","B"), fac=samples$Region[cond.irms], pt.bg=2:3, legpl="topright", xlab="d13C 16:0", ylab="Soil horizon", oneway=T, er.type.nested="sd", xlim=c(-8, 7), add=F)

hor.plot(d13C_org$X22.0[cond.irms]-samples$d13C[cond.irms], hor=samples$Horizon[cond.irms], c("L","F","H","B"), fac=samples$Region[cond.irms], pt.bg=2:3, legpl="topright", xlab="d13C 16:0", ylab="Soil horizon", oneway=T, er.type.nested="sd", xlim=c(-20, 10), add=F)


hor.plot(fungi1.rel[cond.irms]/allbac.rel[cond.irms], hor=samples$Horizon[cond.irms], c("L","F","H","B"), fac=samples$Region[cond.irms], pt.bg=2:3, legpl="topright", xlab="d13C 16:0", ylab="Soil horizon", oneway=F, er.type.nested="sd", xlim=c(0, 1))

hor.plot(fungi1.rel[cond2]/allbac.rel[cond2], hor=samples$Horizon[cond2], c("L","F","H","B"), fac=samples$Region[cond2], pt.bg=2:3, legpl="topright", xlab="d13C 16:0", ylab="Soil horizon", oneway=F, er.type.nested="sd", xlim=c(0, 1))

cond.irms

samples$sample.code[cond.irms]
nrow(d13C_org)

hor.plot(cy19prec, hor=samples$Horizon, c("L","F","H","B"), fac=samples$Region)

pdf("d13c16_0_F2B.pdf", height=5, width=10)
par(mfrow=c(1,2), mar=c(5.1,5.1,4.1,2.1))

x <-fungi2.rel[cond3]/(fungi2.rel[cond3]+allbac.rel[cond3])
y <-(d13C_org$X16.0a[cond3]-samples$d13C[cond3])
corplot(x,y, bg=as.numeric(samples$Horizon[cond3]), pch=as.numeric(samples$Region[cond3])+20, xlab="Fungi : (Fungi + Bacteria)", ylab=expression(paste(Delta, ""^"13", "C" [PLFA(16:0) - SOM])), xlim=c(0,1), ylim=c(-12,2), col=c(1,1))

    lm1 <- lm(y~x)
abline(h=lm1$coefficients[1], lty=3)
abline(h=lm1$coefficients[1]+lm1$coefficients[2], lty=3)

abline(v=1, lty=3)
abline(v=0, lty=3)

y <-(d13C_org$X18.1w9[cond3]-samples$d13C[cond3])

corplot(x, y, bg=as.numeric(samples$Horizon[cond3]), pch=as.numeric(samples$Region[cond3])+20, xlab="Fungi : (Fungi + Bacteria)", ylab=expression(paste(Delta, ""^"13", "C" [PLFA(18:1(w9)) - SOM])), xlim=c(0,1), ylim=c(-8,6), col=c(1,1))
lm1 <- lm(y~x)
abline(h=lm1$coefficients[1], lty=3)
abline(h=lm1$coefficients[1]+lm1$coefficients[2], lty=3)

abline(v=1, lty=3)
abline(v=0, lty=3)

legend("topright", c("L","F","H", "B"), pch=21, pt.bg=c(5,3,4,2), bg="white")

y <-(d13C_org$X18.1w7[cond3]-samples$d13C[cond3])
corplot(x, y, bg=as.numeric(samples$Horizon[cond3]), pch=as.numeric(samples$Region[cond3])+20, xlab="Fungi : (Fungi + Bacteria)", ylab=expression(paste(Delta, ""^"13", "C" [PLFA(18:1(w9)) - SOM])), col=c(1,1))

y <-(d13C_org$cy19.0[cond3]-samples$d13C[cond3])
corplot(x, y, bg=as.numeric(samples$Horizon[cond3]), pch=as.numeric(samples$Region[cond3])+20, xlab="Fungi : (Fungi + Bacteria)", ylab=expression(paste(Delta, ""^"13", "C" [PLFA(18:1(w9)) - SOM])), col=c(1,1))
legend("bottomright", c("L","F","H", "B"), pch=21, pt.bg=c(5,3,4,2), bg="white")
dev.off()

co<-cor.test(d13C_org$X16.0a[cond3]-samples$d13C[cond3], fungi2.rel[cond3]/(fungi2.rel[cond3]+allbac.rel[cond3]))
summary(lm(d13C_org$X16.0a[cond3]-samples$d13C[cond3]~ x))
summary(lm(d13C_org$X16.0[cond3]-samples$d13C[cond3]~ x))
summary(lm(d13C_org$X18.1w9[cond3]-samples$d13C[cond3]~ x))

co<-cor.test(d13C_org$X16.0a[cond3]-samples$d13C[cond3], fungi2.rel[cond3]/(fungi2.rel[cond3]+allbac.rel[cond3]))
lm(d13C_org$X16.0a[cond3]-samples$d13C[cond3]~ x)


x<-(fungi2.rel[cond3]/(fungi2.rel[cond3]+allbac.rel[cond3]))
co
$estimate^2

text(0.3, -4, "asdfas")
paste(co$estimate, siglev(co$p.value))


corplot(d13C_org$X16.0a[cond3], allbac.rel[cond3]/(allbac.rel[cond3]+fungi2.rel[cond3]), bg=as.numeric(samples$Horizon[cond3]), pch=as.numeric(samples$Region[cond3])+20, ylab="Fungi : (Fungi + Bacteria)", xlab=expression(paste(Delta, ""^"13", "C" [PLFA(16:0) - SOM])))

#####################
pdf("pca.pdf")
ord<-rda(rel.colmeans[cond2,T], scale=F)
ord.plot(ord,site.sep2=samples$Horizon[cond2], site.sep1=samples$Site[cond2], pch=21:24, pt.bg=c("green","green","green", "blue","blue","blue"), spe.label.type="text", spe.labels=fames$FAME, cex.leg=.5, spe.mult=1.5, site.cex=1)

par(mfrow=(c(2,2)))
dev.off()

svg("pca2.svg")
ord<-rda(rel[cond2,T], scale=F)
ord.plot(ord,site.sep2=samples$Horizon[cond2], site.sep1=samples$Site[cond2], pch=21:24, pt.bg=c("green","green","green", "blue","blue","blue"), spe.label.type="text", spe.labels=fames$FAME, cex.leg=.5, spe.mult=1.5, site.cex=1)
dev.off()

svg("pca3.svg")
ord<-rda(rel[cond2,T], scale=T)
ord.plot(ord,site.sep2=samples$Horizon[cond2], site.sep1=samples$Site[cond2], pch=21:24, pt.bg=c("green","green","green", "blue","blue","blue"), spe.label.type="text", spe.labels=fames$FAME, cex.leg=.5, spe.mult=1.5, site.cex=1)

horplot(scores(ord, choices=1, display="sites"), )
dev.off()



pdf("groups.pdf", width=12, height=12)
names(list)

par(mfrow=c(3,4))

hor.plot(ia.rel[cond.rel], hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="topright", xlab="I+A 15:0 (mol%,SD)", ylab="Soil horizon", oneway=T, er.type.nested="sd")

hor.plot(ia.conc[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="topright", xlab="I+A 15:0 (mol g-1 d.w.,SD)", ylab="Soil horizon", oneway=T, er.type.nested="sd")

hor.plot(ia.conc[cond.abs]/samples$c_corr[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="topright", xlab="I+A 15:0 (mol g-1 TOC,SD)", ylab="Soil horizon", oneway=T, er.type.nested="sd")

hor.plot(0.000001*ia.conc[cond.abs]*samples$weight[cond.abs]/samples$area[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="topright", xlab="i15:0 + a15:0 (mmol m-2, SD)", ylab="Soil horizon", er.type.nested="sd", oneway=T)

hor.plot( monouns.rel[cond.rel], hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="topright", xlab="Monounsaturated w/o 18:1w9c/t (mol%,SD)", ylab="Soil horizon", oneway=T, er.type.nested="sd")

hor.plot(monouns.conc[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="topright", xlab="Monounsaturated w/o 18:1w9c/t (mol g-1 d.w.,SD)", ylab="Soil horizon", oneway=T, er.type.nested="sd")

hor.plot(monouns.conc[cond.abs]/samples$c_corr[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="topright", xlab="Monounsaturated w/o 18:1w9c/t (mol g-1 TOC,SD)", ylab="Soil horizon", oneway=T, er.type.nested="sd")

hor.plot(0.000001*monouns.conc[cond.abs]*samples$weight[cond.abs]/samples$area[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="topright", xlab="Monounsaturated w/o 18:1w9c/t (mmol m-2, SD)", ylab="Soil horizon", er.type.nested="sd", oneway=T)

hor.plot(fungi2.rel[cond.rel], hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="bottomright", xlab="18:2w6 + 18:3w3 (mol%,SD)", ylab="Soil horizon", oneway=T, er.type.nested="sd")

hor.plot(fungi2.conc[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="bottomright", xlab="18:2w6 + 18:3w3 (mol g-1 d.w.,SD)", ylab="Soil horizon", oneway=T, er.type.nested="sd")

hor.plot(fungi3.conc[cond.abs]/samples$c_corr[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="bottomright", xlab="18:2w6 + 18:3w3 (mol g-1 TOC,SD)", ylab="Soil horizon", oneway=T, er.type.nested="sd")

hor.plot(0.000001*fungi2.conc[cond.abs]*samples$weight[cond.abs]/samples$area[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="bottomright", xlab="18:2w6 + 18:3w3 (mmol m-2, SD)", ylab="Soil horizon", er.type.nested="sd", oneway=T)

#page 2

hor.plot(longeven.rel[cond.rel], hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="topright", xlab="20:0+22:0+24:0 (mol%,SD)", ylab="Soil horizon", oneway=T, er.type.nested="sd")

hor.plot(longeven.conc[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="topright", xlab="20:0+22:0+24:0 (mol g-1 d.w.,SD)", ylab="Soil horizon", oneway=T, er.type.nested="sd")

hor.plot(longeven.conc[cond.abs]/samples$c_corr[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="topright", xlab="20:0+22:0+24:0 (mol g-1 TOC,SD)", ylab="Soil horizon", oneway=T, er.type.nested="sd")

hor.plot(0.000001*longeven.conc[cond.abs]*samples$weight[cond.abs]/samples$area[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="topright", xlab="20:0+22:0+24:0 (mmol m-2, SD)", ylab="Soil horizon", er.type.nested="sd", oneway=T)

hor.plot( longpufa.rel[cond.rel], hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="bottomright", xlab="PUFA 20+ (mol%,SD)", ylab="Soil horizon", oneway=T, er.type.nested="sd")

hor.plot(longpufa.conc[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="bottomright", xlab="PUFA 20+ (mol g-1 d.w.,SD)", ylab="Soil horizon", oneway=T, er.type.nested="sd")

hor.plot(longpufa.conc[cond.abs]/samples$c_corr[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="bottomright", xlab="PUFA 20+ (mol g-1 TOC,SD)", ylab="Soil horizon", oneway=T, er.type.nested="sd")

hor.plot(0.000001*longpufa.conc[cond.abs]*samples$weight[cond.abs]/samples$area[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="topright", xlab="PUFA 20+ (mmol m-2, SD)", ylab="Soil horizon", er.type.nested="sd", oneway=T)

hor.plot(pufa.rel[cond.rel], hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="bottomright", xlab="all PUFA (mol%,SD)", ylab="Soil horizon", oneway=T, er.type.nested="sd")

hor.plot(pufa.conc[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="bottomright", xlab="all PUFA (mol g-1 d.w.,SD)", ylab="Soil horizon", oneway=T, er.type.nested="sd")

hor.plot(pufa.conc[cond.abs]/samples$c_corr[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="bottomright", xlab="all PUFA (mol g-1 TOC,SD)", ylab="Soil horizon", oneway=T, er.type.nested="sd")

hor.plot(0.000001*pufa.conc[cond.abs]*samples$weight[cond.abs]/samples$area[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="bottomright", xlab="all PUFA (mmol m-2, SD)", ylab="Soil horizon", er.type.nested="sd", oneway=T)


dev.off()

pdf("groups2.pdf", width=12, height=12)

par(mfrow=c(3,4))
for (i in c("act","euk","fungi","general","Gneg","Gpos" ))
{
hor.plot(rowSums(rel[cond.rel, fames$microbial.group==i]), hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="bottomright", xlab=paste(i, "(mol%,SD)"), ylab="Soil horizon", oneway=T, er.type.nested="sd")

hor.plot(rowSums(mols[cond.abs, fames$microbial.group==i]), hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="bottomright", xlab=paste(i, "(mol g-1 d.w.,SD)"), ylab="Soil horizon", oneway=T, er.type.nested="sd")

hor.plot(rowSums(mols[cond.abs, fames$microbial.group==i])/samples$c_corr[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="bottomright", xlab=paste(i, "(mol g-1 TOC,SD)"), ylab="Soil horizon", oneway=T, er.type.nested="sd")

hor.plot(0.000001*rowSums(mols[cond.abs, fames$microbial.group==i])*samples$weight[cond.abs]/samples$area[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="bottomright", xlab=paste(i, "(mmol m-2, SD)"), ylab="Soil horizon", er.type.nested="sd", oneway=T)
}
dev.off()

pdf("sumofplfas.pdf", width=7, height=3)
#pdf("sumvsresp.pdf", width=7, height=7)

par(mfrow=c(1,3))

hor.plot(sum[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="bottomright", xlab="Sum of PLFAs (mol g-1 d.w.,SD)", ylab="Soil horizon", oneway=T, er.type.nested="sd")

hor.plot(sum[cond.abs]/samples$c_corr[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="bottomright", xlab="Sum of PLFAs (mol g-1 TOC,SD)", ylab="Soil horizon", oneway=T, er.type.nested="sd", xlim=c(0,100))



attach(resp245)
polygon(
  c((ActrespT2_10_mean[Region=="ER"]+ActrespT2_10_se[Region=="ER"]*sqrt(3))+60,
     ((ActrespT2_10_mean[Region=="ER"]-ActrespT2_10_se[Region=="ER"]*sqrt(3))+60)[4:1]),
  c(4:1,1:4),
  col=rgb(1, 0, 0,0.25), border=1)

polygon(  
  c((ActrespT2_10_mean[Region=="GC"]+ActrespT2_10_se[Region=="GC"]*sqrt(3))+60,
    ((ActrespT2_10_mean[Region=="GC"]-ActrespT2_10_se[Region=="GC"]*sqrt(3))+60)[4:1]),
  c(4:1,1:4),
  col=rgb(0, 1, 0,0.25), border=1)
lines(acc_resp_10_mean[Region=="ER"]/3+60,4:1, col="red", lwd=2)
lines(acc_resp_10_mean[Region=="GC"]/3+60,4:1, col="green", lwd=2)

lines(ActrespT2_10_mean[Region=="ER"]+60,4:1, col="red", lwd=2)
lines(ActrespT2_10_mean[Region=="GC"]+60,4:1, col="green", lwd=2)

resp245

abline(v=60)

hor.plot(0.000001*sum[cond.abs]*samples$weight[cond.abs]/samples$area[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="topright", xlab="Sum of PLFAs (mmol m-2, SD)", ylab="Soil horizon", er.type.nested="sd", oneway=T)
dev.off()

pdf("ratios.pdf", width=7, height=3)

par(mfrow=c(1,3))

hor.plot(F2B[cond.rel], hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="bottomright", xlab="(18:2w6+18:3w3):(Gpos+Gneg+Act)", ylab="Soil horizon", oneway=T, er.type.nested="sd")

hor.plot(Euk2bac[cond.rel], hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="bottomright", xlab="Eukaryots:Bacteria", ylab="Soil horizon", oneway=T, er.type.nested="sd")

hor.plot(cy17prec[cond.rel], hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="topright", xlab="cy17:0/16:1w7", ylab="Soil horizon", oneway=T, er.type.nested="sd")

dev.off()




hor.plot(longeven.rel[cond.rel], hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="bottomright", xlab="even SAFA > 20 (mol%,SD)", ylab="Soil horizon", oneway=T, er.type.nested="sd")

hor.plot(longeven.conc[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="bottomright", xlab="even SAFA >20 (mol g-1 d.w.,SD)", ylab="Soil horizon", oneway=T, er.type.nested="sd")

hor.plot(longeven.conc[cond.abs]/samples$c_corr[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="bottomright", xlab="even SAFA >20 (mol g-1 TOC,SD)", ylab="Soil horizon", oneway=T, er.type.nested="sd")


hor.plot(longpufa[cond.rel], hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="bottomright", xlab="PUFA >20 (mol%,SD)", ylab="Soil horizon", oneway=T, er.type.nested="sd")

hor.plot(longpufa.conc[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="bottomright", xlab="PUFA >20 (mol g-1 d.w.,SD)", ylab="Soil horizon", oneway=T, er.type.nested="sd")

hor.plot(longpufa.conc[cond.abs]/samples$c_corr[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="bottomright", xlab="PUFA >20 (mol g-1 TOC,SD)", ylab="Soil horizon", oneway=T, er.type.nested="sd")



###########
#PLFA per area
##########

pdf("perarea.pdf")
par(mfrow=c(2,2))

Cperm2<-
  hor.plot(0.00001*samples$percent_C[cond.abs]*samples$weight[cond.abs]/samples$area[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="topright", xlab="kg TOC m-2", ylab="Soil horizon", er.type.nested="sd", oneway=T)

tmp<-Cperm2[[1]]
100*tmp[,1:3]/rowSums(tmp[,1:3])

plfaperm2<-
  hor.plot(0.000001*sum[cond.abs]*samples$weight[cond.abs]/samples$area[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="topright", xlab="mmol PLFA m-2", ylab="Soil horizon", er.type.nested="sd", oneway=T)

tmp<-plfaperm2[[1]]
rowSums(tmp[,1:3])
100*tmp[,1:3]/rowSums(tmp[,1:3])
tmp2<-plfaperm2[[2]]
sqrt(rowSums(tmp2[,1:3]^2))

fperm2<-hor.plot(0.000001*fungi2.conc[cond.abs]*samples$weight[cond.abs]/samples$area[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="topright", xlab="18:2w6 + 18:2w3 (mmol m-2, SD)", ylab="Soil horizon", er.type.nested="sd", oneway=T)

tmp<-fperm2[[1]]
rowSums(tmp[,1:3])
100*tmp[,1:3]/rowSums(tmp[,1:3])
tmp2<-plfaperm2[[2]]
sqrt(rowSums(tmp2[,1:3]^2))

hor.plot(0.000001*ia.conc[cond.abs]*samples$weight[cond.abs]/samples$area[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=rep(1,6), nested=samples$Region[cond.abs], pt.bg=2:3, legpl="topright", xlab="i15:0 + a15:0 (mmol m-2, SD)", ylab="Soil horizon", er.type.nested="sd", oneway=T)
dev.off()
#####
#Cy:precursor ratios
#####

hor.plot(cy17prec[cond.rel], hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="topright", xlab="cy17:0/16:1w7 (mol:mol,SD)", ylab="Soil horizon", oneway=T, er.type.nested="sd")

anova(lm(cy17prec[cond2]~samples$Horizon[cond2]))
# Analysis of Variance Table
# 
# Response: cy17prec[cond2]
# Df Sum Sq Mean Sq F value   Pr(>F)    
# samples$Horizon[cond2]  3 1.0272 0.34242  39.466 2.74e-12 ***
#   Residuals              42 0.3644 0.00868                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

lm0<-aov(cy17prec[cond2]~samples$Horizon[cond2])
HSD.test(lm0, "samples$Horizon[cond2]", group=TRUE)
Groups, Treatments and means
# a    B 	 0.6908 
# b 	 L 	 0.4143 
# b 	 H 	 0.3507 
# b 	 F 	 0.3345 

lm0<-aov(cy17prec[cond2]~samples$Horizon[cond2])
HSD.test(lm0, "samples$Horizon[cond2]", group=TRUE)
tmp<-hsd[[5]]
tmp$trt<-factor(tmp$trt, ordered=T, levels=c("L","F","H","B"))
tmp[order(tmp$trt),T]
                
anova(lm(cy17prec[cond2]~samples$Horizon[cond2]*samples$Region[cond2]))
# Analysis of Variance Table
# 
# Response: cy17prec[cond2]
# Df  Sum Sq Mean Sq F value    Pr(>F)    
# samples$Horizon[cond2]                        3 1.02725 0.34242 57.7853 3.164e-14 ***
#   samples$Region[cond2]                         1 0.01946 0.01946  3.2843 0.0778534 .  
# samples$Horizon[cond2]:samples$Region[cond2]  3 0.11976 0.03992  6.7370 0.0009344 ***
#   Residuals                                    38 0.22518 0.00593                      

hor.plot(cy19prec[cond.rel], hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="topright", xlab="cy19:0/18:1w9 (mol:mol,SD)", ylab="Soil horizon", er.type.nested="sd", oneway=T)

anova(lm(F2B[cond2]~samples$Horizon[cond2]))

lm0<-aov(F2B.18.2[cond2]~samples$Horizon[cond2])
HSD.test(lm0, "samples$Horizon[cond2]", group=TRUE)


anova(lm(cy19prec[cond2]~samples$Horizon[cond2]*samples$Region[cond2]))
# Analysis of Variance Table
# 
# Response: cy19prec[cond2]
# Df Sum Sq Mean Sq  F value    Pr(>F)    
# samples$Horizon[cond2]                        3 38.835 12.9449 223.7732 < 2.2e-16 ***
#   samples$Region[cond2]                         1  1.100  1.1004  19.0221  9.53e-05 ***
#   samples$Horizon[cond2]:samples$Region[cond2]  3  1.363  0.4544   7.8553 0.0003363 ***
#   Residuals                                    38  2.198  0.0578           


hor.plot(cyprec[cond.rel], hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="bottomleft", xlab="(cy17:0+cy19:0)/(16:1w7+18:1w9) (mol:mol,SD)", ylab="Soil horizon")

anova(lm(F2B[cond2]~samples$Horizon[cond2]))

lm0<-aov(F2B.18.2[cond2]~samples$Horizon[cond2])
HSD.test(lm0, "samples$Horizon[cond2]", group=TRUE)


anova(lm(cyprec[cond2]~samples$Horizon[cond2]*samples$Region[cond2]))
# Analysis of Variance Table
# 
# Response: cyprec[cond2]
# Df  Sum Sq Mean Sq  F value    Pr(>F)    
# samples$Horizon[cond2]                        3 13.2438  4.4146 222.2535 < 2.2e-16 ***
#   samples$Region[cond2]                         1  0.5598  0.5598  28.1807 5.046e-06 ***
#   samples$Horizon[cond2]:samples$Region[cond2]  3  0.3772  0.1257   6.3293  0.001374 ** 
#   Residuals                                    38  0.7548  0.0199                       

#####
#Fungi:Bacteria ratios
#####

hor.plot(F2B[cond.rel], hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="bottomright", xlab="fungi:bacteria (mol:mol,SD)", ylab="Soil horizon")

anova(lm(F2B[cond2]~samples$Horizon[cond2]))
# Analysis of Variance Table
# 
# Response: F2B[cond2]
# Df Sum Sq Mean Sq F value    Pr(>F)    
# samples$Horizon[cond2]  3 5.2514 1.75047  167.35 < 2.2e-16 ***
#   Residuals              42 0.4393 0.01046                      

lm0<-aov(F2B[cond2]~samples$Horizon[cond2])
HSD.test(lm0, "samples$Horizon[cond2]", group=TRUE)
Groups, Treatments and means
# a    L 	 0.9146 
# b 	 F 	 0.2723 
# b 	 H 	 0.1546 
# c 	 B 	 0.02967 

anova(lm(F2B[cond2]~samples$Horizon[cond2]*samples$Region[cond2]))
# Analysis of Variance Table
# 
# Response: F2B[cond2]
# Df Sum Sq Mean Sq F value  Pr(>F)    
# samples$Horizon[cond2]                        3 5.2514 1.75047 188.049 < 2e-16 ***
#   samples$Region[cond2]                         1 0.0206 0.02063   2.216 0.14484    
# samples$Horizon[cond2]:samples$Region[cond2]  3 0.0650 0.02165   2.326 0.09011 .  



hor.plot(F2B.18.2[cond.rel], hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="bottomright", xlab="fungi (18:2w6,9 only):bacteria (mol:mol,SD)", ylab="Soil horizon", oneway=T)

anova(lm(F2B.18.2[cond2]~samples$Horizon[cond2]))
# Analysis of Variance Table
# 
# Response: F2B.18.2[cond2]
# Df Sum Sq Mean Sq F value    Pr(>F)    
# samples$Horizon[cond2]  3 5.2514 1.75047  167.35 < 2.2e-16 ***
#   Residuals              42 0.4393 0.01046                      
lm0<-aov(F2B.18.2[cond2]~samples$Horizon[cond2])
HSD.test(lm0, "samples$Horizon[cond2]", group=TRUE)
# Groups, Treatments and means
# a    L 	 0.9146 
# b 	 F 	 0.2723 
# b 	 H 	 0.1546 
# c 	 B 	 0.02967 


anova(lm(F2B.18.2[cond2]~samples$Horizon[cond2]*samples$Region[cond2]))
# Analysis of Variance Table
# 
# Response: F2B.18.2[cond2]
#   Df Sum Sq Mean Sq F value  Pr(>F)    
#   samples$Horizon[cond2]                        3 5.2514 1.75047 188.049 < 2e-16 ***
#   samples$Region[cond2]                         1 0.0206 0.02063   2.216 0.14484    
#   samples$Horizon[cond2]:samples$Region[cond2]  3 0.0650 0.02165   2.326 0.09011 .  
#   Residuals                                    38 0.3537 0.00931                   


hor.plot(Euk2bac[cond.rel], hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="bottomright", xlab="fungi (18:2w6,9 only):bacteria (mol:mol,SD)", ylab="Soil horizon", oneway=T)

anova(lm(Euk2bac[cond2]~samples$Horizon[cond2]))
# Analysis of Variance Table
# 
# Response: Euk2bac[cond2]
# Df  Sum Sq Mean Sq F value    Pr(>F)    
# samples$Horizon[cond2]  3 15.4147  5.1382  156.51 < 2.2e-16 ***
#   Residuals              42  1.3789  0.0328                      

lm0<-aov(Euk2bac[cond2]~samples$Horizon[cond2])
HSD.test(lm0, "samples$Horizon[cond2]", group=TRUE)
# Groups, Treatments and means
# a    L 	 1.683 
# b 	 F 	 0.6202 
# b 	 H 	 0.4516 
# c 	 B 	 0.1347 

anova(lm(Euk2bac[cond2]~samples$Horizon[cond2]*samples$Region[cond2]))
# Analysis of Variance Table
# 
# Response: Euk2bac[cond2]
# Df  Sum Sq Mean Sq  F value  Pr(>F)    
#   samples$Horizon[cond2]                        3 15.4147  5.1382 188.6274 < 2e-16 ***
#   samples$Region[cond2]                         1  0.0547  0.0547   2.0090 0.16451    
#   samples$Horizon[cond2]:samples$Region[cond2]  3  0.2890  0.0963   3.5369 0.02358 *  
#   Residuals                                    38  1.0351  0.0272                     



########
#2-way anova for all FAME, rel abundance
########




aov<-anova(lm(rel[cond2,1]~samples$Horizon[cond2]*samples$Region[cond2]))
aov$'Pr(>F)'





lm0<-aov(cy17prec[cond2]~samples$Horizon[cond2])
HSD.test(lm0, "samples$Horizon[cond2]", group=TRUE)

anova(lm(cy17prec[cond2]~samples$Region[cond2]))
aov1<-a
lm1<-aov(cy17prec[cond2]~samples$Horizon[cond2]*samples$Region[cond2])
tuk1<-TukeyHSD(lm1)

str(tuk1)
tuk1[[3]]


plot(tuk1)
HSD.test(lm0, "samples$Horizon[cond2]", group=TRUE)
HSD.test(lm1, "samples$Horizon[cond2]", group=TRUE)
?HSD.test
library("multcomp")
library("multcompView")

multcompLetters(extract_p(tuk1))


tuk <-   glht(lm0, linfct = mcp('samples$Horizon[cond2]'= "Tukey"))
tuk.cld<-cld(tuk)
plot(tuk.cld)

par(mfrow=c(2,2))
cond3<-cond.rel&samples$Horizon=="L"
ord<-rda(rel[cond3,T], scale=T)
ord.plot(ord,site.sep1=samples$Site[cond3], site.sep2=samples$Horizon[cond3], pch=c(21), pt.bg=c(rep(1,3),rep(grey(.6),3)), spe.label.type="text", spe.labels=fames$FAME, cex.leg=.5, spe.mult=2)
title("L horizon")

cond3<-cond.rel&samples$Horizon=="F"
ord<-rda(rel[cond3,T], scale=T)
ord.plot(ord,site.sep1=samples$Site[cond3], site.sep2=samples$Horizon[cond3], pch=c(21), pt.bg=c(rep(1,3),rep(grey(.6),3)), spe.label.type="text", spe.labels=fames$FAME, cex.leg=.5, spe.mult=2)
title("F horizon")

cond3<-cond.rel&samples$Horizon=="H"
ord<-rda(rel[cond3,T], scale=T)
ord.plot(ord,site.sep1=samples$Site[cond3], site.sep2=samples$Horizon[cond3], pch=c(21), pt.bg=c(rep(1,3),rep(grey(.6),3)), spe.label.type="text", spe.labels=fames$FAME, cex.leg=.5, spe.mult=2)
title("H horizon")


cond3<-cond.rel&samples$Horizon=="B"
ord<-rda(rel[cond3,T], scale=T)
ord.plot(ord,site.sep1=samples$Site[cond3], site.sep2=samples$Horizon[cond3], pch=c(21), pt.bg=c(rep(1,3),rep(grey(.6),3)), spe.label.type="text", spe.labels=fames$FAME, cex.leg=.5, spe.mult=2)
title("B horizon")




dev.off()


ord<-rda(rel[cond.rel,T]~samples$Region[cond.rel]*samples$Horizon[cond.rel])

plot(ord)

ord.plot(ord,site.sep1=samples$Horizon[cond.rel], site.sep2=samples$Site[cond.rel], pch=c(21,21,21,22,22,22,24), pt.bg=1:6, spe.label.type="text", spe.labels=fames$FAME, cex.leg=.5)



warnings()

samples$Horizon

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

hor.plot(  rowSums(mols[cond.rel,fames$microbial.group=="fungi"])/
             (rowSums(mols[cond.rel,fames$microbial.group=="Gneg"])+
                rowSums(mols[cond.rel,fames$microbial.group=="Gpos"])), 
           hor=samples$Horizon[cond.rel],  c("L","F","H","B"), fac=samples$Site[cond.rel],
           lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="bottomright", 
           xlab="Fungi:Bacteria (mol:mol, SE)", ylab="Soil horizon", er.type="se", lwd=2)

hor.plot(  rowSums(mols[cond.rel,fames$microbial.group=="Gpos"])/
             rowSums(mols[cond.rel,fames$microbial.group=="Gneg"]),
           hor=samples$Horizon[cond.rel],  c("L","F","H","B"), fac=samples$Site[cond.rel],
           lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="bottomright", 
           xlab="G+/G- (mol:mol, SE)", ylab="Soil horizon", er.type="se", lwd=2)

hor.plot(rowSums(mols[cond.rel,fames$microbial.group=="act"])
           ,
           hor=samples$Horizon[cond.rel],  c("L","F","H","B"), fac=samples$Site[cond.rel],
           lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="bottomright", 
           xlab="actinomycetes (mol:mol, SE)", ylab="Soil horizon", er.type="se", lwd=2)

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

unsaturation16.0<-(rel$X2.1.16.1a+rel$X2.3.16.1b+rel$X2.4.16.1w7+rel$X16.1c)/rel$X16.0

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

hor.plot(((rel$cy.19.0a+rel$cy.19.0b+rel$cy17.0a)/(rel$X2.4.16.1w7+rel$X18.1w9c))[cond.rel],
         hor=samples$Horizon[cond.rel],  c("L","F","H","B"), fac=samples$Site[cond.rel],
         lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="topright", 
         xlab="cyclopropyl:precursor (mol:mol, SE)", ylab="Soil horizon", er.type="se", lwd=2)

hor.plot(unsaturation16.0[cond.rel],
         hor=samples$Horizon[cond.rel],  c("L","F","H","B"), fac=samples$Site[cond.rel],
         lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="bottomleft", 
         xlab="sum(16:1)/16:0 (mol:mol, SE)", ylab="Soil horizon", er.type="se", lwd=2)

hor.plot(
  (rel$X18.1w9c/rel$X18.1b)[cond.rel],
         hor=samples$Horizon[cond.rel],  c("L","F","H","B"), fac=samples$Site[cond.rel],
         lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="bottomright", 
         xlab="18:1w9c/18:1w7, mol:mol, SE", ylab="Soil horizon", er.type="se", 
  lwd=2)

hor.plot(
  (rel$X18.3w3.6.9/rel$X18.2w6.9)[cond.rel],
  hor=samples$Horizon[cond.rel],  c("L","F","H","B"), fac=samples$Site[cond.rel],
  lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="bottomright", 
  xlab="18:3w3c/18:2w6, mol:mol, SE", ylab="Soil horizon", er.type="se", 
  lwd=2)

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

  hor.plot(rel$X22.0[cond.rel], hor=samples$Horizon[cond.rel],  c("L","F","H","B"), fac=samples$Site[cond.rel], lty=rep(1,6), nested=samples$Region[cond.rel], pt.bg=2:3, legpl="bottomright", 
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




summary(sum)

hor.plot(sum[cond2], hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=rep(1), nested=samples$Region)

hor.plot(sum[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=1, nested=samples$Region[cond.abs], er.type.nested="sd", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="Sum of PLFA (nmol g-1 d.w.)", ylab="Soil horizon")

hor.plot(mols$X18.2w6.9[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=1, nested=samples$Region[cond.abs], er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="18:2w6,9 (nmol g-1 d.w.)", ylab="Soil horizon")

hor.plot((mols$X18.2w6.9+mols$X18.3w3.6.9)[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=1, nested=samples$Region[cond.abs], er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="18:2w6,9 (nmol g-1 d.w.)", ylab="Soil horizon")

hor.plot((mols$X18.1w9c+mols$X18.2w6.9+mols$X18.3w3.6.9)[cond.abs], hor=samples$Horizon[cond.abs], c("L","F","H","B"), fac=samples$Site[cond.abs], lty=1, nested=samples$Region[cond.abs], er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="18:2w6,9 (nmol g-1 d.w.)", ylab="Soil horizon")

hor.plot(
  rowSums(mols[cond.rel, fames$microbial.group=="fungi"])/
    (rowSums(mols[cond.rel, fames$microbial.group=="Gpos"]+rowSums(mols[cond.rel, fames$microbial.group=="Gneg"]))), 
  
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="Fungi:Bacteria (SE)", ylab="Soil horizon")

hor.plot(
  mols$i.15.0[cond.rel]/mols$ai.15.0[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="Fungi:Bacteria (SE)", ylab="Soil horizon")

hor.plot(
  mols$cy17.0[cond.rel]/mols$X2.4.16.1w7[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="Fungi:Bacteria (SE)", ylab="Soil horizon")

hor.plot(
  mols$cy.19.0.1[cond.rel]/mols$X18.1w9c[cond.rel], 
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
  rel$X20.0[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="20:0 (mol%, SE)", ylab="Soil horizon",legsize=0.8)

hor.plot(
  rel$X22.0[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="22:0 (mol%, SE)", ylab="Soil horizon",legsize=0.8)

hor.plot(
  rel$X23.0_[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="24:0 (mol%, SE)", ylab="Soil horizon",legsize=0.8)

hor.plot(
  rel$X24.0[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="24:0 (mol%, SE)", ylab="Soil horizon",legsize=0.8)

hor.plot(
  rel$X25.0[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="24:0 (mol%, SE)", ylab="Soil horizon",legsize=0.8)

hor.plot(
  rel$X10.Me.16.0[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="20:3 (mol%, SE)", ylab="Soil horizon",legsize=0.8)

hor.plot(
  rel$X10.Me17.0[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="20:3 (mol%, SE)", ylab="Soil horizon",legsize=0.8)

hor.plot(
  rel$X2.1.16.1a[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="20:3 (mol%, SE)", ylab="Soil horizon",legsize=0.8)

hor.plot(
  rel$X2.3.16.1b[cond.rel], 
  hor=samples$Horizon[cond.rel], c("L","F","H","B"), fac=samples$Site[cond.rel], lty=1, nested=samples$Region[cond.rel] , er.type.nested="se", pch=c(21,22), pt.bg=c(2,3), legpl="bottomright", lwd=2, xlab="20:3 (mol%, SE)", ylab="Soil horizon",legsize=0.8)

hor.plot(
  rel$X16.1c[cond.rel], 
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
  
attach(resp)
resp_act[is.na(resp_act)]<-0
  days<-as.numeric(harvest)
cond<-temperature=="10"
  timeseries(resp_act[cond], xfac=days[cond], sepfac=region[cond], legend="topright")