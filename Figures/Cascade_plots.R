# Load the data for Gorongosa
g.prob <- read.table("Processed_data/Gorongosa_SpProbIter100.txt", sep="\t", header=T, stringsAsFactors=T)
g.rich <- read.table("Processed_data/Gorongosa_SpRichIter100.txt", sep="\t", header=T, stringsAsFactors=T)

# Load the data for Serengeti
s.prob <- read.table("Processed_data/Serengeti_SpProbIter100.txt", sep="\t", header=T, stringsAsFactors=T)
s.rich <- read.table("Processed_data/Serengeti_SpRichIter100.txt", sep="\t", header=T, stringsAsFactors=T)

# Load the data for Laikipia
l.prob <- read.table("Processed_data/Laikipia_SpProbIter100.txt", sep="\t", header=T, stringsAsFactors=T)
l.rich <- read.table("Processed_data/Laikipia_SpRichIter100.txt", sep="\t", header=T, stringsAsFactors=T)


# Install and load the required packages
#install.packages(c("vioplot","RColorBrewer"))
library(vioplot)
library(RColorBrewer)

########################################################################################
# Set up the plot and dimensions
png(filename="Figures/Complete_cascade.png",width=24,height=18,units="cm",res=300)

# Define plot margins
#Set panel outline. The top panel is used for the legend
m <- (matrix(c(1,2,3,4,4,4,5,6,7), nrow = 3, ncol = 3, byrow = TRUE))

# Set heights for plot panels
layout(mat = m,heights = c(0.44,0.12,0.44))
# Set margins
par(mai=c(.7,.7,0.1,0.1))

########################################################################################
# Make the base violin plot of the simulated data (replicates are iterations)
plot(0,0,type="n", las=1, ylab="Species richness", xlab="Number of resources",
	,cex.axis=1.1, cex.lab= 1.3, mgp=c(2.4,0.6,0), ylim=c(0,10), xlim=c(-0.15,14.15), xaxt='n')
# Add horizontal axis
axis(1, at =seq(0,14,by=2), labels=c("0","20","40","60","80","100","120","140"),cex=1.1)

# Add violin plot
X <- vioplot(g.rich, add=T, col=rgb(0.8,0.8,0.8,1), border=NA,rectCol="black", lineCol="black", 
	colMed=rgb(0.7,0,0,1),wex=1.2)
# Add points for the minimum and maximum resource level
points(c(0,14.4),c(0,8),pch=16,col=rgb(0.7,0,0,1))

# Add points to show the median values
#points(c(0,seq(1,14,by=1),14.4),c(0,X$median,8),pch=1,col=rgb(0,0,0,1), cex=1.3)
# Label the panel
mtext("(a) Gorongosa",cex=1.1, side = 3, adj = 0.03, line = -2,font=1)

# Add a prediction line by defining independent (IV) and dependent (DV) variables and modelling an asymptotic relationship
IV <- rep(1:14,100)
DV <- c(t(g.rich))
nlreg <- nls(DV~a*IV^b,start=list(a=8,b=0.5))

# Predict new data from model parameters
xvals <- seq(0,14.4,l=100)
pred.val <- predict(nlreg, list(IV = xvals) , se.fit = TRUE,interval = "confidence",level = 0.95)
# Add lines to plot
lines(xvals, pred.val, lwd=2, col=rgb(0.8,0,0,0.7), lty=1)

########################################################################################

# Make the base violin plot of the simulated data (replicates are iterations)
plot(0,0,type="n", las=1, ylab="Species richness", xlab="Number of resources",
	,cex.axis=1.1, cex.lab= 1.3, mgp=c(2.4,0.6,0), ylim=c(0,8), xlim=c(-0.15,10), xaxt='n')
# Add horizontal axis
axis(1, at =seq(0,10,by=2), labels=c("0","20","40","60","80","100"),cex=1.1)

# Add violin plot
X <- vioplot(s.rich, add=T, col=rgb(0.8,0.8,0.8,1), border=NA,rectCol="black", lineCol="black", 
	colMed=rgb(0.7,0,0,1),wex=1.2)
# Add points for the minimum and maximum resource level
points(c(0,9.1),c(0,4),pch=16,col=rgb(0.7,0,0,1))

# Add points to show the median values
#points(c(0,seq(1,14,by=1),14.4),c(0,X$median,8),pch=1,col=rgb(0,0,0,1), cex=1.3)
# Label the panel
mtext("(b) Serengeti",cex=1.1, side = 3, adj = 0.03, line = -2,font=1)

# Add a prediction line by defining independent (IV) and dependent (DV) variables and modelling an asymptotic relationship
IV <- rep(1:9,100)
DV <- c(t(s.rich))
nlreg <- nls(DV~a*IV^b,start=list(a=8,b=0.5))

# Predict new data from model parameters
xvals <- seq(0,9.1,l=100)
pred.val <- predict(nlreg, list(IV = xvals) , se.fit = TRUE,interval = "confidence",level = 0.95)
# Add lines to plot
lines(xvals, pred.val, lwd=2, col=rgb(0.8,0,0,0.7), lty=1)


########################################################################################

# Make the base violin plot of the simulated data (replicates are iterations)
plot(0,0,type="n", las=1, ylab="Species richness", xlab="Number of resources",
	,cex.axis=1.1, cex.lab= 1.3, mgp=c(2.4,0.6,0), ylim=c(0,10), xlim=c(-0.15,13), xaxt='n')
# Add horizontal axis
axis(1, at =seq(0,12,by=2), labels=c("0","20","40","60","80","100","120"),cex=1.1)

# Add violin plot
X <- vioplot(l.rich, add=T, col=rgb(0.8,0.8,0.8,1), border=NA,rectCol="black", lineCol="black", 
	colMed=rgb(0.7,0,0,1),wex=1.2)
# Add points for the minimum and maximum resource level
points(c(0,12.1),c(0,6),pch=16,col=rgb(0.7,0,0,1))

# Add points to show the median values
#points(c(0,seq(1,14,by=1),14.4),c(0,X$median,8),pch=1,col=rgb(0,0,0,1), cex=1.3)
# Label the panel
mtext("(c) Laikipia",cex=1.1, side = 3, adj = 0.03, line = -2,font=1)

# Add a prediction line by defining independent (IV) and dependent (DV) variables and modelling an asymptotic relationship
IV <- rep(1:12,100)
DV <- c(t(l.rich))
nlreg <- nls(DV~a*IV^b,start=list(a=8,b=0.5))

# Predict new data from model parameters
xvals <- seq(0,12.6,l=100)
pred.val <- predict(nlreg, list(IV = xvals) , se.fit = TRUE,interval = "confidence",level = 0.95)
# Add lines to plot
lines(xvals, pred.val, lwd=2, col=rgb(0.8,0,0,0.7), lty=1)


########################################################################################
########################################################################################

# Create a species legend using a distinguishable colour ramp
cols2 <- c(brewer.pal(n = 10, name = "Paired"),brewer.pal(n = 10, name = "Paired"))
par(mai=c(0,.7,0,0))

# Make a blank plot
plot(0,0,axes=FALSE,xlab="",ylab="", type=n)

# Create a vector of species names for the legend
sp.names <-c("Buffalo", "Grevy's zebra", "Hartebeest", "Oribi", "Oryx",
	"Plain's zebra", "Reedbuck", "Sable", "Topi", "Waterbuck",
	"Wildebeest", "Grant's gazelle", "Impala", "Thomson's gazelle",
	"Bushbuck", "Dik-dik", "Eland", "Klipspringer", "Kudu", "Nyala")

legend("bottom", ncol=5,  lty=rep(c(1,3),each=10), col=cols2, bty='n',legend=sp.names,cex=1.2, lw=2)



#######################################################################################
# Reset margins
par(mai=c(.7,.7,0.1,0.1))

# A vector so that the species correspond to the colors in the legend
g.cols <- c(13,3,11,8,10,4,7,1,20,15,19)
# Are the lines solid or dashed
dash <- ifelse(g.cols>10,3,1)

# Make the plot
plot(0,0,type="n",xlab="Number of resources", ylab="Presence probability", xlim=c(0,145), 
	ylim=c(0,100),las=1 ,cex.axis=1.1, cex.lab= 1.3, mgp=c(2.4,0.6,0))

# The persistence of the full stable community
pres <- c(0,100,0,100,100,100,0,100,100,100,100)
# Add the lines for each species
for (nsp in 1:11){
	lines(c(0,seq(10,140,by=10),144),c(0,g.prob[,nsp],pres[nsp]), col=cols2[g.cols[nsp]], lwd=2, lty=dash[nsp])
}
# Label the panel
mtext("(d)",cex=1.1, side = 3, adj = 0.03, line = -2,font=1)
###########################################################################################
# A vector so that the species correspond to the colors in the legend
s.cols <- c(13,3,11,9,6,14,12,1)
# Are the lines solid or dashed
dash <- ifelse(s.cols>10,3,1)

# Make the plot
plot(0,0,type="n",xlab="Number of resources", ylab="Presence probability", xlim=c(0,100), 
	ylim=c(0,100),las=1, cex.axis=1.1, cex.lab= 1.3, mgp=c(2.4,0.6,0))

# The persistence of the full stable community
pres <- c(100,0,100,100,0,0,100,0)
# Add the lines for each species
for (nsp in 1:8){
	lines(c(0,seq(10,90,by=10),91),c(0,s.prob[,nsp],pres[nsp]), col=cols2[s.cols[nsp]], lwd=2, lty=dash[nsp])
}
# Label the panel
mtext("(e)",cex=1.1, side = 3, adj = 0.03, line = -2,font=1)
###########################################################################################
# A vector so that the species correspond to the colors in the legend
l.cols <- c(13,3,2,6,10,16,12,5,18,1,17,19)
# Are the lines solid or dashed
dash <- ifelse(l.cols>10,31)

# Make the plot
plot(0,0,type="n",xlab="Number of resources", ylab="Presence probability", xlim=c(0,120), 
	ylim=c(0,100),las=1,cex.axis=1.1, cex.lab= 1.3, mgp=c(2.4,0.6,0))

# The persistence of the full stable community
pres <- c(0,100,100,0,100,0,100,0,0,100,100,0)
# Add the lines for each species
for (nsp in 1:12){
	lines(c(0,seq(10,120,by=10),121),c(0,l.prob[,nsp],pres[nsp]), col=cols2[l.cols[nsp]], lwd=2, lty=dash[nsp])
}
# Label the panel
mtext("(f)",cex=1.1, side = 3, adj = 0.03, line = -2,font=1)
###########################################################################################

# Close plot device and save to file
dev.off()


