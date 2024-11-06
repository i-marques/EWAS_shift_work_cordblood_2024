################################################################################
#
#   METAL output file processing - ShiftWork MA - Irene Marques - 24.02.2022
#
################################################################################

setwd("Z:/1_My_Projects_/2.ShiftWork EWAS/3.MA/24.02.2022/M3/Postp25.02correct")
## name of METAL output file
infile<-"Z:/1_My_Projects_/2.ShiftWork EWAS/3.MA/24.02.2022/M3/MA.MaternalShiftWork.M3.IreneAug2024.1.txt"

## name of postprocess files
filepref<-"MA_M3_MatShiftWork"

##################### READ FILES ---------------------

# Metal file
df<-read.table(infile, header = T, strings = F, dec=".", sep="\t")

# read crossreactive probes: combination of Naeem and Chen 
list<-read.csv("Z:/1_My_Projects_/2.ShiftWork EWAS/3.MA/crossreactiveprobes.csv")

# read Illumina annotation
annotation <-"Z:/1_My_Projects_/2.ShiftWork EWAS/3.MA/450k_annotationfile_v1_2.csv"
annotation<-read.csv(annotation, header = T, strings = F)
annotation<-annotation[,c(1:2,11:13,17:19,22:33)]
colnames(annotation)[2] <- "MarkerName"

# FLAG PROBES WITH SNPS DO THIS AFTER FINAL MA
flagged<-"Z:/1_My_Projects_/2.ShiftWork EWAS/3.MA/flaggedprobes.csv"
flagged<-read.csv(flagged, header = T, strings = F)
colnames(flagged)[1] <- "MarkerName"

#################### QC -------------------------------

#removing probes present in only one cohort (HetDf = 0)  
df<-df[which(df$HetDf > 0),] 

# remove crossreactive probes
df = subset(df, !(df$MarkerName %in% list$MarkerName))

# add annotation columns
df<-merge(df, annotation, by.x = "MarkerName", by.y = "MarkerName", all.x = T, all.y = F)
dim(df)
#[1] 440191     30

# remove probes in Chr X and Y
df<-df[which(df$CHR!="X"&df$CHR!="Y"),]
dim(df)
#[1] 429959     30

# add "flag" for probes with SNPs
df<-merge(df, flagged, by.x = "MarkerName", by.y = "MarkerName", all.x = T, all.y = F)

# Add Benjamini-Hochberg p-value correction to data
df$P.value.FDR <- p.adjust(df$P.value, method = "BH")
sum(df$P.value.FDR < 0.05)
#[1] 3

# save QCed MA file
outfile <- paste(filepref, ".csv", sep = "")
write.csv(df, file = outfile, row.names = F, quote = F)

########## Calculate lambda

lambda <- qchisq(median(df$P.value, na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)

outfile<- paste(filepref,".lambda.txt", sep = "")
write.table(lambda,file=outfile, row.names=F, quote=F)

#################### Plots -------------------------------
#read QCed MA file if needed
#df <- read.csv("MA_M3_MatShiftWork.csv")

### Manhattan plot ### as figure a) ---------------------
### function to create plot
#Chr <- c(1:22,"X","Y")
Chr <- c(1:22)
GWplot<-function(data,P,Chr,title){
        print(table(data$CHR))
        par(mar=c(5,5,4,2));
        phy.max<-tapply(data$MAPINFO, data$CHR,max,na.rm=T)
        cumlen=0
        for(i in Chr){

        data[data$CHR==i,"loc"]<-data[data$CHR==i,"MAPINFO"]+cumlen
        cumlen<-cumlen+phy.max[i]
        }
        phy.med<-tapply(data$loc, data$CHR,median,na.rm=T)[Chr]
        print(phy.med)
        data$mlgP<--log(data[,P], base=10)
        plot(data[,"loc"],data[,"mlgP"],type="n", xaxt="n", xlab="Chromosome",
            ylab="-log10(P)" , cex.axis=1.5, cex.lab=1.5, xlim=c(0,max(data$loc,na.rm=T)))
        # next line sets figure title aligned to the left (convenient if using a) B)...
        mtext(title, side = 3, line = 2, adj = 0, cex = 1.5)
        axis(side=1, at=phy.med[c(1:19)], labels=Chr[c(1:19)],
         tick=T,cex.axis=1.3,las=3)
        ### adjust 20:22 20:25 to  if using X and Y chr
        axis(side=1, at=phy.med[c(20:22)], labels=Chr[c(20:22)],
         tick=T,cex.axis=1.3,las=3,)
        for(i in Chr){
          if(which(Chr==i) %in% seq(2,30,2)) col="blue" else col="red"
          points(data[data$CHR==i,"loc"],data[data$CHR==i,"mlgP"],col=col,pch=20,cex=0.5)
          #abline(7,0,col="grey", lty=2) ### Bonferroni Correction
          abline(6.45,0,col="grey", lty=3, lwd=3) ### FDR correction --- needs to be adjusted
        }
}
outfile<- paste(filepref, "ManhattanPlot.png", sep = "_")
png(outfile, width =16,height =10,units = "in", res = 300)
GWplot(df,"P.value",Chr,title=" a)") #change title if not using several plots in same figure
dev.off()

### Volcano plot ### as figure b) ---------------------
outfile<- paste(filepref, "Volcano.png", sep = "_")
png(outfile,
    units="cm", 
    width=60,
    height=40,
    pointsize=25, 
    res=300) #300 is usually the requested resulution for publication

with(df, plot((Effect), -log10(P.value), 
                  pch=20, 
                  col="black",
                  main="",
                  xlab="Effect size", 
                  ylab="-log10(P)",
                  cex.axis=1.2,
                  cex.lab=1.2,
                  abline(h=-log10(3.5e-07), col="grey",lty=3, lwd=3)
))
#change title if not using several plots in same figure
mtext("b)", side = 3, line = 2, adj = 0, cex = 1.5)  # Title aligned to the left

dev.off()

### QQ plot ### independent figure  ---------------------

observed <- sort(df$P.value)
observed<- as.numeric(observed)
lobs <- -(log10(observed))

expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))

outfile <- paste(filepref, "QQplot.png", sep = "_") 
png(outfile, width =16,height =10,units = "in", res = 300)
par(cex = 1.3, cex.axis = 1.3, cex.lab = 1.15 )
plot(c(0,30), c(0,30), col="red", lwd=4, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,6), ylim=c(0,7), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=23, cex=.4, bg="black")
dev.off()

############################### QC and Plots done






