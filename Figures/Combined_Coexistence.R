setwd("C:/Users/Falko/Documents/Standalone Research/Mammal Isotopes/Most recent manuscript and files")

# Load the data
comm <- read.csv("Processed_data/MCT_combined_data.csv")

# Total number of unique species
tot.S <- length(unique(comm$Species))

# Define a colour ramp with shades of green for grazers, shades of red for mixed feeders, and shades of blue for browsers
cols2 <- c(colorRampPalette(c(rgb(0,0.4,0,0.25),rgb(0.3,0.9,0.3,0.25)),interpolate="linear")(11),
			colorRampPalette(c(rgb(0.6,0,0,0.25),rgb(1,0.5,0,0.25)),interpolate="linear")(3),
			colorRampPalette(c(rgb(0,0.4,1,0.25),rgb(0,0,0.5,0.25)),interpolate="linear")(6))


# Save to file in the 'Figures' directory
png(filename="Figures/Combined Coexistence.png",width=22,height=7,units="cm",res=300)

# Set plot panels and margins
par(mfrow=c(1,3))
par(mai=c(0.55,0.55,0.075,0.075))

##########################################################################################################
# Create blank axis plots
plot(0,0, xlim=c(0,1), ylim=c(0,1),type="n",las=1,
	ylab="Fitness differences",xlab="Niche differences",cex.axis=1, cex.lab= 1.15, mgp=c(2.4,0.6,0))
# A grey polygon for the zone of stable coexistence
polygon(c(-5,5,5),c(-5,5,-5),col=rgb(0,0,0,0.075), border=NA)

# Use a subset of data from Gorongosa
G.sub <- comm[comm$PA=="Gorongosa",]
# Plot the species as points (point characters denote whether species can coexist stably)
points (G.sub$N,G.sub$F,pch=G.sub$Point,bg=cols2[G.sub$Code], cex=1.2, col="black",lw=0.5)
# Label the points
text (G.sub$N,G.sub$F,G.sub$Species, pos=c(4,4,2,2,3,1,4,4,2,2,3), cex=0.55, offset=0.3)
# Label the plot panel
mtext("(a) Gorongosa",cex=0.8, side = 3, adj = 0.03, line = -1.5,font=1)

# Create a vector of species names for the legend
sp.names <-c("Buffalo", "Grevy's zebra", "Hartebeest", "Oribi", "Oryx",
	"Plain's zebra", "Reedbuck", "Sable", "Topi", "Waterbuck",
	"Wildebeest", "Grant's gazelle", "Impala", "Thomson's gazelle",
	"Bushbuck", "Dik-dik", "Eland", "Klipspringer", "Kudu", "Nyala")

# Add a legen to the first panel
legend("bottomleft", pch=21, col="black",pt.bg=cols2,legend=sp.names,pt.cex=1, cex=0.65, pt.lw=0.5)

##########################################################################################################
# Create blank axis plots
plot(0,0, xlim=c(0,1), ylim=c(0,1),type="n",las=1,
	ylab="Fitness differences",xlab="Niche differences",cex.axis=1, cex.lab= 1.15, mgp=c(2.4,0.6,0))
# A grey polygon for the zone of stable coexistence
polygon(c(-5,5,5),c(-5,5,-5),col=rgb(0,0,0,0.075), border=NA)

# Use a subset of data from Serengeti
S.sub <- comm[comm$PA=="Serengeti",]
# Plot the species as points (point characters denote whether species can coexist stably)
points (S.sub$N,S.sub$F,pch=S.sub$Point,bg=cols2[S.sub$Code], cex=1.2, , col="black",lw=0.5)
# Label the points
text (S.sub$N,S.sub$F,S.sub$Species, pos=c(4,4,1,2,2,2,2,2), cex=0.55, offset=0.3)
# Label the panel
mtext("(b) Serengeti",cex=0.8, side = 3, adj = 0.03, line = -1.5,font=1)


##########################################################################################################
# Create blank axis plots
plot(0,0, xlim=c(0,1), ylim=c(0,1),type="n",las=1,
	ylab="Fitness differences",xlab="Niche differences",cex.axis=1, cex.lab= 1.15, mgp=c(2.4,0.6,0))
# A grey polygon for the zone of stable coexistence
polygon(c(-5,5,5),c(-5,5,-5),col=rgb(0,0,0,0.075), border=NA)

# Use a subset of data from Laikipia
L.sub <- comm[comm$PA=="Laikipia",]
# Plot the species as points (point characters denote whether species can coexist stably)
points (L.sub$N,L.sub$F,pch=L.sub$Point,bg=cols2[L.sub$Code], cex=1.2, , col="black",lw=0.5)
# Label the points
text (L.sub$N,L.sub$F,L.sub$Species, pos=c(4,4,4,2,4,3,2,3,4,2,3,2),cex=0.55, offset=0.3)
# Label the panel
mtext("(c) Laikipia",cex=0.8, side = 3, adj = 0.03, line = -1.5,font=1)

# Close plot device and save to file
dev.off()
