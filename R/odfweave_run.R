install.packages("odfWeave")
install.packages("XML")
library("XML")
library("odfWeave")
getImageDefs()

setwd("/home/lluc/Dropbox/MUN/phd/results/NL-BELT/Batches21-25combined/R/")
infile<-"NLBELT_plots.odt"
outfile<-"NLBELT_plots_out.odt"

imageDefs <- getImageDefs()
imageDefs$dispWidth <- 7
imageDefs$dispHeight<- 7
imageDefs$type="eps"
imageDefs$device="postscript"
imageDefs$plotWidth<-7
imageDefs$plotHeight<-7 
imageDefs$args = list(
  horizontal = FALSE, 
  onefile = FALSE, 
  paper = "special")
setImageDefs(imageDefs)
dir<-getwd()
odfWeave(infile, outfile, workDir = paste(getwd(), "/tmp", sep=""))
warnings()
?odfWeave

pwd()
getwd()
sqrt(1:10)

?RweaveOdf
warnings()
ord.plot(ord=rda(rsim), site.sep1=samples$type, site.sep2=samples$days, col=1, pt.bg=colscale, pch=21:24, spe.labels=peaks$origin, spe.label.type="text")

dev.off()

peaks
setwd("~/ligpaper2") 
warnings()