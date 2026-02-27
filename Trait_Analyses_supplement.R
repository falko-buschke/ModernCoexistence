install.packages(c("bipartate","rgl"))
library(bipartite)
library(rgl)


# Load the data
comm <- read.table("Raw_data/Gorongosa_data.txt", sep="\t", header=T, stringsAsFactors=T)

# Load the data
coexist <- read.csv("Processed_data/MCT_combined_data.csv")

# Create a vector of species names
sp.name <- comm$Species
R.vect <- (comm$Rmax - comm$M)#/365
# The consumption matrix: the proportional dietary contribution of plant species, multiplied by the daily matabolic intake requirement
cons.mat <- comm[,5:dim(comm)[2]]*(0.05*comm$BM^0.77)
# A vector of species' body masses
BM.vect <- comm$BM
#Consumption
cons <- (0.05*comm$BM^0.77)

D.mat.m <- (dfun((cons.mat))[[2]])


sub <- coexist[which(coexist$PA=="Gorongosa"),]
mergedf=merge(data.frame('key'=sp.name,sp.name),data.frame('key'=sub$Species,sub),
	by='key',all=T, sort=F)

status <- mergedf$Role 


png(filename="Figures/Supp/Figure S4.png",width=20,height=12,units="cm",res=300)

par(mai=c(.55,.55,0.08,0.03))
par(mfcol=c(2,3))
plot(R.vect,mergedf$N, pch=16, col=as.factor(status), log="",cex=1.5,
	las=1, xlab="Intrinsic growth rate", ylab="Niche differences", cex.axis=1.1, cex.lab= 1.3, mgp=c(2.6,0.6,0))
#text(R.vect,mergedf$N ,sp.name,pos=3, cex=0.6)
mtext(paste(round(cor(R.vect,mergedf$N),3),if(cor.test(R.vect,mergedf$N)$p.value < 0.05){"*"},sep=""),
	cex=1.2, side = 3, adj = 0.95, line = -1.5)

plot(R.vect,mergedf$F, pch=16, col=as.factor(status), log="",cex=1.5,
	las=1, xlab="Intrinsic growth rate", ylab="Fitness differences", cex.axis=1.1, cex.lab= 1.3, mgp=c(2.6,0.6,0))
#text(R.vect,mergedf$F ,sp.name,pos=3, cex=0.6)
mtext(paste(round(cor(R.vect,mergedf$F),3),if(cor.test(R.vect,mergedf$F)$p.value < 0.05){"*"},sep=""),
	cex=1.2, side = 3, adj = 0.95, line = -1.5)

plot(cons,mergedf$N, pch=16, col=as.factor(status), log="",cex=1.5,
	las=1, xlab="Daily dietary requirements (kg)", ylab="Niche differences", cex.axis=1.1, cex.lab= 1.3, mgp=c(2.6,0.6,0))
#text(cons,mergedf$N ,sp.name,pos=3, cex=0.6)
mtext(paste(round(cor(cons,mergedf$N),3),if(cor.test(cons,mergedf$N)$p.value < 0.05){"*"},sep=""),
	cex=1.2, side = 3, adj = 0.95, line = -1.5)

plot(cons,mergedf$F, pch=16, col=as.factor(status), log="",cex=1.5,
	las=1, xlab="Daily dietary requirements (kg)", ylab="Fitness differences", cex.axis=1.1, cex.lab= 1.3, mgp=c(2.6,0.6,0))
#text(cons,mergedf$F ,sp.name,pos=3, cex=0.6)
mtext(paste(round(cor(cons,mergedf$F),3),if(cor.test(cons,mergedf$F)$p.value < 0.05){"*"},sep=""),
	cex=1.2, side = 3, adj = 0.95, line = -1.5)

plot(D.mat.m,mergedf$N, pch=16, col=as.factor(status), log="",cex=1.5,
	las=1, xlab="Diet specialisation (Bluthgen 'd)", ylab="Niche differences", cex.axis=1.1, cex.lab= 1.3, mgp=c(2.6,0.6,0))
#text(D.mat.m,mergedf$N ,sp.name,pos=3, cex=0.6)
mtext(paste(round(cor(D.mat.m,mergedf$N),3),if(cor.test(D.mat.m,mergedf$N)$p.value < 0.05){"*"},sep=""),
	cex=1.2, side = 3, adj = 0.95, line = -1.5)

plot(D.mat.m,mergedf$F, pch=16, col=as.factor(status), log="",cex=1.5,
	las=1, xlab="Diet specialisation (Bluthgen 'd)", ylab="Fitness differences", cex.axis=1.1, cex.lab= 1.3, mgp=c(2.6,0.6,0))
#text(D.mat.m,mergedf$F ,sp.name,pos=3, cex=0.6)
mtext(paste(round(cor(D.mat.m,mergedf$F),3),if(cor.test(D.mat.m,mergedf$F)$p.value < 0.05){"*"},sep=""),
	cex=1.2, side = 3, adj = 0.95, line = -1.5)

dev.off()


png(filename="Figures/Supp/Figure S1.png",width=20,height=8,units="cm",res=300)

par(mai=c(.3,.65,0.1,0.1))

par(mfrow=c(1,3))
boxplot(R.vect~status,ylab="Intrinsic growth rate", names=c("Coexist","Excluded"),las=1, xlab="", col="lightgrey")
stripchart(R.vect~status,method = "jitter",jitter=0.2, pch = 16, cex=1.5 ,col = rgb(1,0,0,0.75),vertical = TRUE, add = TRUE)
mtext(paste("W = ",round(wilcox.test(R.vect~status)$statistic,3),if(wilcox.test(R.vect~status)$p.value < 0.05){"*"}else("(n.s.)"),sep=""),
	cex=0.9, side = 3, adj = 0.05, line = -1.5)

boxplot(cons~status,ylab="Daily dietary requirements", names=c("Coexist","Excluded"),las=1, xlab="", col="lightgrey")
stripchart(cons~status,method = "jitter",jitter=0.2, pch = 16, cex=1.5 ,col = rgb(1,0,0,0.75),vertical = TRUE, add = TRUE)
mtext(paste("W = ",round(wilcox.test(cons~status)$statistic,3),if(wilcox.test(cons~status)$p.value < 0.05){"*"}else("(n.s.)"),sep=""),
	cex=0.9, side = 3, adj = 0.05, line = -1.5)

boxplot(D.mat.m~status,ylab="Dietary specialisation (Bluthgen 'd)", names=c("Coexist","Excluded"),las=1, xlab="", col="lightgrey")
stripchart(D.mat.m~status,method = "jitter",jitter=0.2, pch = 16, cex=1.5 ,col = rgb(1,0,0,0.75),vertical = TRUE, add = TRUE)
mtext(paste("W = ",round(wilcox.test(D.mat.m~status)$statistic,3),if(wilcox.test(D.mat.m~status)$p.value < 0.05){"*"}else("(n.s)"),sep=""),
	cex=0.9, side = 3, adj = 0.05, line = -1.5)
dev.off() 


############################################################################################

# Load the data
comm <- read.table("Raw_data/Serengeti_data.txt", sep="\t", header=T, stringsAsFactors=T)

## Create a vector of species names
sp.name <- comm$Species
sp.name <- c("Impala", "Hartebeest", "Wildebeest", "Topi", "Plains zebra", "Thomson's gazelle", "Grant's gazelle", "Buffalo")

R.vect <- (comm$Rmax - comm$M)#/365
# The consumption matrix: the proportional dietary contribution of plant species, multiplied by the daily matabolic intake requirement
cons.mat <- comm[,5:dim(comm)[2]]*(0.05*comm$BM^0.77)
# A vector of species' body masses
BM.vect <- comm$BM
#Consumption
cons <- (0.05*comm$BM^0.77)

D.mat.m <- (dfun((cons.mat))[[2]])


sub <- coexist[which(coexist$PA=="Serengeti"),]
mergedf=merge(data.frame('key'=sp.name,sp.name),data.frame('key'=sub$Species,sub),
	by='key',all=T, sort=F)

status <- mergedf$Role 


png(filename="Figures/Supp/Figure S5.png",width=20,height=12,units="cm",res=300)

par(mai=c(.55,.55,0.08,0.03))
par(mfcol=c(2,3))
plot(R.vect,mergedf$N, pch=16, col=as.factor(status), log="",cex=1.5,
	las=1, xlab="Intrinsic growth rate", ylab="Niche differences", cex.axis=1.1, cex.lab= 1.3, mgp=c(2.6,0.6,0))
#text(R.vect,mergedf$N ,sp.name,pos=3, cex=0.6)
mtext(paste(round(cor(R.vect,mergedf$N),3),if(cor.test(R.vect,mergedf$N)$p.value < 0.05){"*"},sep=""),
	cex=1.2, side = 3, adj = 0.95, line = -1.5)

plot(R.vect,mergedf$F, pch=16, col=as.factor(status), log="",cex=1.5,
	las=1, xlab="Intrinsic growth rate", ylab="Fitness differences", cex.axis=1.1, cex.lab= 1.3, mgp=c(2.6,0.6,0))
#text(R.vect,mergedf$F ,sp.name,pos=3, cex=0.6)
mtext(paste(round(cor(R.vect,mergedf$F),3),if(cor.test(R.vect,mergedf$F)$p.value < 0.05){"*"},sep=""),
	cex=1.2, side = 3, adj = 0.95, line = -1.5)

plot(cons,mergedf$N, pch=16, col=as.factor(status), log="",cex=1.5,
	las=1, xlab="Daily dietary requirements (kg)", ylab="Niche differences", cex.axis=1.1, cex.lab= 1.3, mgp=c(2.6,0.6,0))
#text(cons,mergedf$N ,sp.name,pos=3, cex=0.6)
mtext(paste(round(cor(cons,mergedf$N),3),if(cor.test(cons,mergedf$N)$p.value < 0.05){"*"},sep=""),
	cex=1.2, side = 3, adj = 0.95, line = -1.5)

plot(cons,mergedf$F, pch=16, col=as.factor(status), log="",cex=1.5,
	las=1, xlab="Daily dietary requirements (kg)", ylab="Fitness differences", cex.axis=1.1, cex.lab= 1.3, mgp=c(2.6,0.6,0))
#text(cons,mergedf$F ,sp.name,pos=3, cex=0.6)
mtext(paste(round(cor(cons,mergedf$F),3),if(cor.test(cons,mergedf$F)$p.value < 0.05){"*"},sep=""),
	cex=1.2, side = 3, adj = 0.95, line = -1.5)

plot(D.mat.m,mergedf$N, pch=16, col=as.factor(status), log="",cex=1.5,
	las=1, xlab="Diet specialisation (Bluthgen 'd)", ylab="Niche differences", cex.axis=1.1, cex.lab= 1.3, mgp=c(2.6,0.6,0))
#text(D.mat.m,mergedf$N ,sp.name,pos=3, cex=0.6)
mtext(paste(round(cor(D.mat.m,mergedf$N),3),if(cor.test(D.mat.m,mergedf$N)$p.value < 0.05){"*"},sep=""),
	cex=1.2, side = 3, adj = 0.95, line = -1.5)

plot(D.mat.m,mergedf$F, pch=16, col=as.factor(status), log="",cex=1.5,
	las=1, xlab="Diet specialisation (Bluthgen 'd)", ylab="Fitness differences", cex.axis=1.1, cex.lab= 1.3, mgp=c(2.6,0.6,0))
#text(D.mat.m,mergedf$F ,sp.name,pos=3, cex=0.6)
mtext(paste(round(cor(D.mat.m,mergedf$F),3),if(cor.test(D.mat.m,mergedf$F)$p.value < 0.05){"*"},sep=""),
	cex=1.2, side = 3, adj = 0.95, line = -1.5)

dev.off()


png(filename="Figures/Supp/Figure S2.png",width=20,height=8,units="cm",res=300)

par(mai=c(.3,.65,0.1,0.1))

par(mfrow=c(1,3))
boxplot(R.vect~status,ylab="Intrinsic growth rate", names=c("Coexist","Excluded"),las=1, xlab="", col="lightgrey")
stripchart(R.vect~status,method = "jitter",jitter=0.2, pch = 16, cex=1.5 ,col = rgb(1,0,0,0.75),vertical = TRUE, add = TRUE)
mtext(paste("W = ",round(wilcox.test(R.vect~status)$statistic,3),if(wilcox.test(R.vect~status)$p.value < 0.05){"*"}else("(n.s.)"),sep=""),
	cex=0.9, side = 3, adj = 0.95, line = -1.5)

boxplot(cons~status,ylab="Daily dietary requirements", names=c("Coexist","Excluded"),las=1, xlab="", col="lightgrey")
stripchart(cons~status,method = "jitter",jitter=0.2, pch = 16, cex=1.5 ,col = rgb(1,0,0,0.75),vertical = TRUE, add = TRUE)
mtext(paste("W = ",round(wilcox.test(cons~status)$statistic,3),if(wilcox.test(cons~status)$p.value < 0.05){"*"}else("(n.s.)"),sep=""),
	cex=0.9, side = 3, adj = 0.95, line = -1.5)

boxplot(D.mat.m~status,ylab="Dietary specialisation (Bluthgen 'd)", names=c("Coexist","Excluded"),las=1, xlab="", col="lightgrey")
stripchart(D.mat.m~status,method = "jitter",jitter=0.2, pch = 16, cex=1.5 ,col = rgb(1,0,0,0.75),vertical = TRUE, add = TRUE)
mtext(paste("W = ",round(wilcox.test(D.mat.m~status)$statistic,3),if(wilcox.test(D.mat.m~status)$p.value < 0.05){"*"}else("(n.s)"),sep=""),
	cex=0.9, side = 3, adj = 0.95, line = -1.5)
dev.off() 



############################################################################################

# Load the data
comm <- read.table("Raw_data/Laikipia_data.txt", sep="\t", header=T, stringsAsFactors=T)

# Create a vector of species names
sp.name <- comm$Species
sp.name <- c("Impala", "Hartebeest", "Grevy's zebra", "Plains zebra", "Waterbuck", "Dik-dik", 
	"Grant's gazelle", "Oryx", "Klipspringer","Buffalo")

R.vect <- (comm$Rmax - comm$M)#/365
# The consumption matrix: the proportional dietary contribution of plant species, multiplied by the daily matabolic intake requirement
cons.mat <- comm[,5:dim(comm)[2]]*(0.05*comm$BM^0.77)
# A vector of species' body masses
BM.vect <- comm$BM
#Consumption
cons <- (0.05*comm$BM^0.77)

D.mat.m <- (dfun((cons.mat))[[2]])


sub <- coexist[which(coexist$PA=="Laikipia"),]
mergedf=merge(data.frame('key'=sp.name,sp.name),data.frame('key'=sub$Species,sub),
	by='key',all=T, sort=F)

status <- mergedf$Role 


png(filename="Figures/Supp/Figure S6.png",width=20,height=12,units="cm",res=300)

par(mai=c(.55,.55,0.08,0.03))
par(mfcol=c(2,3))
plot(R.vect,mergedf$N, pch=16, col=as.factor(status), log="",cex=1.5,
	las=1, xlab="Intrinsic growth rate", ylab="Niche differences", cex.axis=1.1, cex.lab= 1.3, mgp=c(2.6,0.6,0))
#text(R.vect,mergedf$N ,sp.name,pos=3, cex=0.6)
mtext(paste(round(cor(R.vect,mergedf$N),3),if(cor.test(R.vect,mergedf$N)$p.value < 0.05){"*"},sep=""),
	cex=1.2, side = 3, adj = 0.95, line = -1.5)

plot(R.vect,mergedf$F, pch=16, col=as.factor(status), log="",cex=1.5,
	las=1, xlab="Intrinsic growth rate", ylab="Fitness differences", cex.axis=1.1, cex.lab= 1.3, mgp=c(2.6,0.6,0))
#text(R.vect,mergedf$F ,sp.name,pos=3, cex=0.6)
mtext(paste(round(cor(R.vect,mergedf$F),3),if(cor.test(R.vect,mergedf$F)$p.value < 0.05){"*"},sep=""),
	cex=1.2, side = 3, adj = 0.95, line = -1.5)

plot(cons,mergedf$N, pch=16, col=as.factor(status), log="",cex=1.5,
	las=1, xlab="Daily dietary requirements (kg)", ylab="Niche differences", cex.axis=1.1, cex.lab= 1.3, mgp=c(2.6,0.6,0))
#text(cons,mergedf$N ,sp.name,pos=3, cex=0.6)
mtext(paste(round(cor(cons,mergedf$N),3),if(cor.test(cons,mergedf$N)$p.value < 0.05){"*"},sep=""),
	cex=1.2, side = 3, adj = 0.95, line = -1.5)

plot(cons,mergedf$F, pch=16, col=as.factor(status), log="",cex=1.5,
	las=1, xlab="Daily dietary requirements (kg)", ylab="Fitness differences", cex.axis=1.1, cex.lab= 1.3, mgp=c(2.6,0.6,0))
#text(cons,mergedf$F ,sp.name,pos=3, cex=0.6)
mtext(paste(round(cor(cons,mergedf$F),3),if(cor.test(cons,mergedf$F)$p.value < 0.05){"*"},sep=""),
	cex=1.2, side = 3, adj = 0.95, line = -1.5)

plot(D.mat.m,mergedf$N, pch=16, col=as.factor(status), log="",cex=1.5,
	las=1, xlab="Diet specialisation (Bluthgen 'd)", ylab="Niche differences", cex.axis=1.1, cex.lab= 1.3, mgp=c(2.6,0.6,0))
#text(D.mat.m,mergedf$N ,sp.name,pos=3, cex=0.6)
mtext(paste(round(cor(D.mat.m,mergedf$N),3),if(cor.test(D.mat.m,mergedf$N)$p.value < 0.05){"*"},sep=""),
	cex=1.2, side = 3, adj = 0.95, line = -1.5)

plot(D.mat.m,mergedf$F, pch=16, col=as.factor(status), log="",cex=1.5,
	las=1, xlab="Diet specialisation (Bluthgen 'd)", ylab="Fitness differences", cex.axis=1.1, cex.lab= 1.3, mgp=c(2.6,0.6,0))
#text(D.mat.m,mergedf$F ,sp.name,pos=3, cex=0.6)
mtext(paste(round(cor(D.mat.m,mergedf$F),3),if(cor.test(D.mat.m,mergedf$F)$p.value < 0.05){"*"},sep=""),
	cex=1.2, side = 3, adj = 0.95, line = -1.5)

dev.off()


png(filename="Figures/Supp/Figure S3.png",width=20,height=8,units="cm",res=300)

par(mai=c(.3,.65,0.1,0.1))

par(mfrow=c(1,3))
boxplot(R.vect~status,ylab="Intrinsic growth rate", names=c("Coexist","Excluded"),las=1, xlab="", col="lightgrey")
stripchart(R.vect~status,method = "jitter",jitter=0.2, pch = 16, cex=1.5 ,col = rgb(1,0,0,0.75),vertical = TRUE, add = TRUE)
mtext(paste("W = ",round(wilcox.test(R.vect~status)$statistic,3),if(wilcox.test(R.vect~status)$p.value < 0.05){"*"}else("(n.s.)"),sep=""),
	cex=0.9, side = 3, adj = 0.95, line = -1.5)

boxplot(cons~status,ylab="Daily dietary requirements", names=c("Coexist","Excluded"),las=1, xlab="", col="lightgrey")
stripchart(cons~status,method = "jitter",jitter=0.2, pch = 16, cex=1.5 ,col = rgb(1,0,0,0.75),vertical = TRUE, add = TRUE)
mtext(paste("W = ",round(wilcox.test(cons~status)$statistic,3),if(wilcox.test(cons~status)$p.value < 0.05){"*"}else("(n.s.)"),sep=""),
	cex=0.9, side = 3, adj = 0.05, line = -1.5)

boxplot(D.mat.m~status,ylab="Dietary specialisation (Bluthgen 'd)", names=c("Coexist","Excluded"),las=1, xlab="", col="lightgrey")
stripchart(D.mat.m~status,method = "jitter",jitter=0.2, pch = 16, cex=1.5 ,col = rgb(1,0,0,0.75),vertical = TRUE, add = TRUE)
mtext(paste("W = ",round(wilcox.test(D.mat.m~status)$statistic,3),if(wilcox.test(D.mat.m~status)$p.value < 0.05){"*"}else("(n.s)"),sep=""),
	cex=0.9, side = 3, adj = 0.95, line = -1.5)
dev.off() 

