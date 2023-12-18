setwd("C:/Users/Falko/Documents/Standalone Research/Mammal Isotopes/Most recent manuscript and files")

# Load the data
comm <- read.csv("Processed_data/MCT_combined_data.csv")


# Define the colour ramp
col.ramp <- c(rgb(0,0.7,0,1), rgb(0.6,0,0,1), rgb(0,0,.8,1))

# Annualised or daily rates (change this parameter for more intutitve interpretation)
rate <- "annual" # or daily

# Correction factor for growth rate
yr <- ifelse(rate=="annual", 365,1)
# Interval range for horizontal axis
int.l <- ifelse(rate=="annual", -1.4, -0.004)
int.h <- ifelse(rate=="annual", .55, 0.002)
int <- c(int.l,int.h)

# Set plot dimensions
png(filename="Figures/Growth rates.png",width=14,height=28,units="cm",res=300)

# Set plot margins
par(mfrow=c(3,1))
par(mai=c(0.5,1.1,0.3,0.15))


########################################################################################
# Use a subset of data from Gorongosa
G.sub <- comm[comm$PA=="Gorongosa",]

# Creat a sub=matrix with growth rates
G.gr <- as.matrix(G.sub[,c(5,6,7)])*yr
rownames(G.gr) <- G.sub[,2]


# Create bar chart
barplot(height=t(G.gr[order(G.gr[,1]),]), beside=T, horiz=T, las=1, col= rep(col.ramp, dim(G.gr)[1]),
	xlim=int,xlab="Growth rates",cex.axis=1.1, cex.lab= 1.3, mgp=c(2.4,0.6,0))

# add zero-line and bounding box
abline(v=0); box()
# Add axis ticks and guidelines
axis(2,at =seq(2.5,42.5,by=4), labels=rep("",11))
abline(h=seq(0.5,44.5,by=4), col="lightgrey")

# Add panel label
mtext("(a) Gorongosa",cex=1.25, side = 3, adj = -0.25, line = 0.8,font=1)

#######################################################################################3
########################################################################################
# Use a subset of data from Serengeti
S.sub <- comm[comm$PA=="Serengeti",]

# Creat a sub=matrix with growth rates
S.gr <- as.matrix(S.sub[,c(5,6,7)])*yr
rownames(S.gr) <- S.sub[,2]


# Create bar chart
barplot(height=t(S.gr[order(S.gr[,1]),]), beside=T, horiz=T, las=1, col= rep(col.ramp, dim(S.gr)[1]),
	xlim=int,xlab="Growth rate",cex.axis=1.1, cex.lab= 1.3, mgp=c(2.4,0.6,0))

# add zero-line and bounding box
abline(v=0); box()
# Add axis ticks and guidelines
axis(2,at =seq(2.5,30.5,by=4), labels=rep("",8))
abline(h=seq(0.5,36.5,by=4), col="lightgrey")

# Add panel label
mtext("(b) Serengeti",cex=1.25, side = 3, adj = -0.25, line = 0.8,font=1)

#######################################################################################3
########################################################################################
# Use a subset of data from Laikipia
L.sub <- comm[comm$PA=="Laikipia",]

# Creat a sub=matrix with growth rates
L.gr <- as.matrix(L.sub[,c(5,6,7)])*yr
rownames(L.gr) <- L.sub[,2]


# Create bar chart
barplot(height=t(L.gr[order(L.gr[,1]),]), beside=T, horiz=T, las=1, col= rep(col.ramp, dim(L.gr)[1]),
	xlim=int,xlab="Growth rate",cex.axis=1.1, cex.lab= 1.3, mgp=c(2.4,0.6,0))

# add zero-line and bounding box
abline(v=0); box()
# Add axis ticks and guidelines
axis(2,at =seq(2.5,46.5,by=4), labels=rep("",12))
abline(h=seq(0.5,48.5,by=4), col="lightgrey")

# Add panel label
mtext("(c) Laikipia",cex=1.25, side = 3, adj = -0.25, line = 0.8,font=1)

# Add legend
legend("bottomleft", pch=22, col="black",pt.bg=col.ramp,bg="white",
	legend=c("Intrinsic growth rate","Invasion growth rate","No-niche growth rate" ),pt.cex=1.5, cex=1.25, pt.lw=0.5)

#######################################################################################3

# Close plot device and save plot
dev.off()

