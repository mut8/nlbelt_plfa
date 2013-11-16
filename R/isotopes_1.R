setwd("/home/lluc/Dropbox/MUN/phd/results/NL-BELT/Batches21-25combined/R")

source('~/R_functions/functions.R')
source('~/R_functions/plfaplotfunction.R')

pch<-c(22,23)
colscale2<-c("red", "blue")

iso<-read.csv('isotopes_first.csv')

var<-iso$SOM
hor.plot(var,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$region, pch=pch, pt.bg=colscale2, xlim=c(-35,10), xlab='d13C', ylab='horizon', legpl='topright')

pdf('isos.pdf')
par(mfrow=c(3,3))

# var<-iso$i15.0
# hor.plot(var,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$region, pch=pch, pt.bg=colscale2, xlim=c(-35,-20), xlab='d13C', ylab='horizon', legpl='topright', main="SOM")
# 
var<-iso$i15.0-iso$SOM
hor.plot(var,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$region, pch=pch, pt.bg=colscale2, xlim=c(-5,5), xlab='d13C', ylab='horizon', legpl='topright', nested=F, main="iu15:0")


var<-iso$a15.0-iso$SOM
hor.plot(var,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$region, pch=pch, pt.bg=colscale2, xlim=c(-5,5), xlab='d13C', ylab='horizon', legpl='topright', nested=F, main="a15:0")

var<-iso$X16.0-iso$SOM
hor.plot(var,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$region, pch=pch, pt.bg=colscale2, xlab='d13C', ylab='horizon', legpl='topright', nested=F, xlim=c(-9,1), main="16:0")

var<-iso$X18.1w7-iso$SOM
hor.plot(var,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$region, pch=pch, pt.bg=colscale2, xlab='d13C', ylab='horizon', legpl='topright', nested=F, xlim=c(-8,2), main="18:1w7")

var<-iso$X18.1w9-iso$SOM
hor.plot(var,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$region, pch=pch, pt.bg=colscale2, xlab='d13C', ylab='horizon', legpl='topright', nested=F, xlim=c(-6,4), main="18:1w9")

var<-iso$X18.2w6-iso$SOM
hor.plot(var,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$region, pch=pch, pt.bg=colscale2, xlab='d13C', ylab='horizon', legpl='topright', nested=F, xlim=c(-9,1), main="18:2w6")

var<-iso$X18.3w6-iso$SOM
hor.plot(var,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$region, pch=pch, pt.bg=colscale2, xlab='d13C', ylab='horizon', legpl='topright', nested=F, xlim=c(-8,2), main="18:3w6")

var<-iso$X18.3w3-iso$SOM
hor.plot(var,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$region, pch=pch, pt.bg=colscale2, xlab='d13C', ylab='horizon', legpl='topright', nested=F, xlim=c(-12,-2), main="18:3")

# var<-iso$X20.4w6 -iso$SOM
# hor.plot(var,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$region, pch=pch, pt.bg=colscale2, xlab='d13C', ylab='horizon', legpl='topright', nested=F, xlim=c(-12,-2), main="20:4w6")

dev.off()

pdf('isos2.pdf')
par(mfrow=c(3,3))

# var<-iso$i15.0
# hor.plot(var,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$region, pch=pch, pt.bg=colscale2, xlim=c(-35,-20), xlab='d13C', ylab='horizon', legpl='topright', main="SOM")
# 
colscale3<-c(1,1,1,2,2,2)
pch2=c(1,1,1,16,16,16)

var<-iso$i15.0-iso$SOM
hor.plot(var,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$site, pch=pch2, col=colscale3, pt.bg=colscale3, xlim=c(-5,5), xlab='d13C', ylab='horizon', legpl='topright', nested=F, main="i15:0", lty=1, lwd=2)


var<-iso$a15.0-iso$SOM
hor.plot(var,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$site, pch=pch2, col=colscale3, pt.bg=colscale3, xlim=c(-5,5), xlab='d13C', ylab='horizon', legpl='topright', nested=F, main="a15:0", lty=1, lwd=2)

var<-iso$X16.0-iso$SOM
hor.plot(var,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$site, pch=pch2, col=colscale3, pt.bg=colscale3, xlim=c(-9,1), xlab='d13C', ylab='horizon', legpl='topright', nested=F, main="16:0", lty=1, lwd=2)
abline(v=0, col="grey")

var<-iso$X18.1w7-iso$SOM
hor.plot(var,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$site, pch=pch2, col=colscale3, pt.bg=colscale3, xlim=c(-8,2), xlab='d13C', ylab='horizon', legpl='topright', nested=F, main="18:1w7", lty=1, lwd=2)
abline(v=0, col="grey")

var<-iso$X18.1w9-iso$SOM
hor.plot(var,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$site, pch=pch2, col=colscale3, pt.bg=colscale3, xlim=c(-6,4), xlab='d13C', ylab='horizon', legpl='topright', nested=F, main="18:1w7", lty=1, lwd=2)
abline(v=0, col="grey")

var<-iso$X18.2w6-iso$SOM
hor.plot(var,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$site, pch=pch2, col=colscale3, pt.bg=colscale3, xlim=c(-9,1), xlab='d13C', ylab='horizon', legpl='topright', nested=F, main="18:2w6", lty=1, lwd=2)
abline(v=0, col="grey")

var<-iso$X18.3w6-iso$SOM
hor.plot(var,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$region, pch=pch, pt.bg=colscale2, xlab='d13C', ylab='horizon', legpl='topright', nested=F, xlim=c(-8,2), main="18:3w6")

var<-iso$X18.3w3-iso$SOM
hor.plot(var,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$region, pch=pch, pt.bg=colscale2, xlab='d13C', ylab='horizon', legpl='topright', nested=F, xlim=c(-12,-2), main="18:3")


pch=c(21,22)



var<-iso$SOM
hor.plot(var,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$region, pch=pch, pt.bg=colscale2, xlim=c(-40,-20), xlab='d13C', ylab='horizon', legpl='topright', sig=F, lty=3)

var<-iso$X18.2w6
hor.plot(var,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$region, pch=pch, pt.bg=colscale2, xlim=c(-9,1), xlab='d13C', ylab='horizon', legpl='topright', nested=F, main="18:2w6", lty=1, lwd=1, addlines=T, sig=F)


par(mfrow=c(3,3))
var<-iso$SOM
xlab=expression(paste(delta, ""^{13}, C))

var2<-iso$i15.0
main="i15:0"
hor.plot(var,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$region, pch=c(1,1), pt.bg=colscale2, col=colscale2, xlim=c(-40,-20), xlab=xlab, ylab='horizon', legpl='topright', sig=F, lty=3)
hor.plot(var2,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$region, pch=pch, pt.bg=colscale2, xlim=c(-9,1), xlab='d13C', ylab='horizon', legpl='topright', nested=F, main="18:2w6", lty=1, lwd=2, addlines=T, sig=F)
legend("topleft", lty=c(3,1), c("SOM", main))

var2<-iso$a15.0
main2="a15:0"
hor.plot(var,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$region, pch=c(1,1), pt.bg=colscale2, col=colscale2, xlim=c(-40,-20), xlab=xlab, ylab='horizon', legpl='topright', sig=F, lty=3)
hor.plot(var2,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$region, pch=pch, pt.bg=colscale2, xlim=c(-9,1), xlab='d13C', ylab='horizon', legpl='topright', nested=F, main="18:2w6", lty=1, lwd=2, addlines=T, sig=F)
legend("topleft", lty=c(3,1), c("SOM", main))

var2<-iso$X15.0
main="15:0"
hor.plot(var,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$region, pch=c(1,1), pt.bg=colscale2, col=colscale2, xlim=c(-40,-20), xlab=xlab, ylab='horizon', legpl='topright', sig=F, lty=3)
hor.plot(var2,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$region, pch=pch, pt.bg=colscale2, xlim=c(-9,1), xlab='d13C', ylab='horizon', legpl='topright', nested=F, main="18:2w6", lty=1, lwd=2, addlines=T, sig=F)
legend("topleft", lty=c(3,1), c("SOM", main))

var2<-iso$X18.1w7
main="18:1w7"
hor.plot(var,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$region, pch=c(1,1), pt.bg=colscale2, col=colscale2, xlim=c(-40,-20), xlab=xlab, ylab='horizon', legpl='topright', sig=F, lty=3)
hor.plot(var2,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$region, pch=pch, pt.bg=colscale2, xlim=c(-9,1), xlab='d13C', ylab='horizon', legpl='topright', nested=F, main="18:2w6", lty=1, lwd=2, addlines=T, sig=F)
legend("topleft", lty=c(3,1), c("SOM", main))


var2<-iso$X16.0
main="16:0"
hor.plot(var,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$region, pch=c(1,1), pt.bg=colscale2, col=colscale2, xlim=c(-40,-20), xlab=xlab, ylab='horizon', legpl='topright', sig=F, lty=3)
hor.plot(var2,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$region, pch=pch, pt.bg=colscale2, xlim=c(-9,1), xlab='d13C', ylab='horizon', legpl='topright', nested=F, main="18:2w6", lty=1, lwd=2, addlines=T, sig=F)
legend("bottomright", lty=c(3,1), c("SOM", main))

var2<-iso$X18.1w9c
main="18:1w9"
hor.plot(var,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$region, pch=c(1,1), pt.bg=colscale2, col=colscale2, xlim=c(-40,-20), xlab=xlab, ylab='horizon', legpl='topright', sig=F, lty=3)
hor.plot(var2,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$region, pch=pch, pt.bg=colscale2, xlim=c(-9,1), xlab='d13C', ylab='horizon', legpl='topright', nested=F, main="18:2w6", lty=1, lwd=2, addlines=T, sig=F)
legend("bottomright", lty=c(3,1), c("SOM", main))

dev.off()


svg("isotope_example2.svg")
var1 <- iso$i15.0-iso$SOM
var2 <- iso$X16.0-iso$SOM
var3 <- iso$X18.1w9c-iso$SOM
var4 <- iso$X18.2w6-iso$SOM

pch<-c(21,22)
colscale2<-1
lty<-c(1,3)

xlab=expression(paste(delta, ""^{13}, C))
main=""
hor.plot(var1,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$region, pch=pch, pt.bg=grey(.7), col=colscale2, xlim=c(-10,8), xlab=xlab, ylab='horizon', sig=F, lty=lty, lwd=2)

hor.plot(var2,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$region, pch=pch, col=colscale2, pt.bg=grey(.5), xlim=c(-9,1), nested=F, main="18:2w6", lty=lty, lwd=2, addlines=T, sig=F)

hor.plot(var3,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$region, pch=pch, col=colscale2, pt.bg=grey(.3), xlim=c(-9,1),  nested=F, main="18:2w6", lty=lty, lwd=2, addlines=T, sig=F)

hor.plot(var4,hor=iso$horizon, horlev=c("L", "F", "H", "B"), fac=iso$region, pch=pch, col=colscale2, pt.bg=grey(0), xlim=c(-9,1),  nested=F, main="18:2w6", lty=lty, lwd=2, addlines=T, sig=F)
legend("bottomright", col=c(grey(.7), grey(.5), grey(.3), grey(0)), c("i15:0","16:0","18:1w9c", "18:2w6"), lty=1, pch=16)
legend("topright", col=1, lty=c(1,3), pt.bg=1, lwd=2, pch=pch, c("cold", "warm"))
dev.off()

svg("isotope_example1.svg", height=3.5, width=4.5)
var1 <- iso$i15.0[iso$region=="GC"]-iso$SOM[iso$region=="GC"]
var2 <- iso$X16.0[iso$region=="GC"]-iso$SOM[iso$region=="GC"]
var3 <- iso$X18.1w9c[iso$region=="GC"]-iso$SOM[iso$region=="GC"]
var4 <- iso$X18.2w6[iso$region=="GC"]-iso$SOM[iso$region=="GC"]
pch<-c(16,17)
colscale2<-1

xlab=expression(paste(Delta,delta, ""^{13}, "C"[PLFA-SOM]))
main=""
hor.plot(var1,hor=iso$horizon[iso$region=="GC"], horlev=c("L", "F", "H", "B"), fac=iso$region[iso$region=="GC"], pch=pch, pt.bg=1, col=colscale2, xlim=c(-8,8), xlab=xlab, ylab='horizon', sig=F, lty=1, lwd=2)

hor.plot(var2,hor=iso$horizon[iso$region=="GC"], horlev=c("L", "F", "H", "B"), fac=iso$region[iso$region=="GC"], pch=pch, col=colscale2, pt.bg=2, xlim=c(-9,1), nested=F, main="18:2w6", lty=1, lwd=2, addlines=T, sig=F)

hor.plot(var3,hor=iso$horizon[iso$region=="GC"], horlev=c("L", "F", "H", "B"), fac=iso$region[iso$region=="GC"], pch=pch, col=colscale2, pt.bg=3, xlim=c(-9,1),  nested=F, main="18:2w6", lty=1, lwd=2, addlines=T, sig=F)

hor.plot(var4,hor=iso$horizon[iso$region=="GC"], horlev=c("L", "F", "H", "B"), fac=iso$region[iso$region=="GC"], pch=pch, col=colscale2, pt.bg=4, xlim=c(-9,1),  nested=F, main="18:2w6", lty=1, lwd=2, addlines=T, sig=F)

legend("bottomright", col=1:4, c("i15:0","16:0","18:1w9c", "18:2w6"), lty=1, pch=16, cex=.7)
dev.off()