# Load the data
gor.pw <- read.table("Pairwise/Gorongosa.txt", sep="\t", header=T, stringsAsFactors=T)

comm<- read.csv("Processed_data/MCT_combined_data.csv")

png(filename="Figures/Figure 4.png",width=16,height=24,units="cm",res=300)
#pdf(file = "Manuscript/Submission Files/Resubmission/Figure4.pdf", height = 9.44882, width = 6.29921) 

par(mfrow=c(3,2))

# Use a subset of data from Gorongosa
G.sub <- comm[comm$PA=="Gorongosa",]

# Create a sub-matrix with niche differences
par(mai=c(0.6,1.2,0.1,0))


stripchart(gor.pw$N[gor.pw$N!=0] ~gor.pw$Invader[gor.pw$N!=0],              # Data
           method = "stack", # Random noise
           pch = 16,          # Pch symbols
           col = rgb(0,0,0,0.5),           # Color of the symbol
           vertical = FALSE,   # Vertical mode
           add = FALSE,
           cex=1,
           las=1, xlim=c(0,1),
           xlab="Niche differences",cex.axis=1.1, cex.lab= 1.3, mgp=c(2.4,0.6,0)) 

gor.mean <- xtabs(gor.pw$N~gor.pw$Invader) /(length(unique(gor.pw$Invader))-1)         

stripchart(gor.mean~ names(gor.mean),            # Data
           method = "stack", # Random noise
           pch = 16,          # Pch symbols
           col = rgb(0,0,1,0.85),           # Color of the symbol
           vertical = FALSE,   # Vertical mode
           add = TRUE,
           cex=2,
           las=1) 

stripchart(G.sub$N~G.sub$Species,              # Data
           method = "stack", # Random noise
           pch = "+",          # Pch symbols
           col = rgb(1,0,0,0.85),           # Color of the symbol
           vertical = FALSE,   # Vertical mode
           add = TRUE,
           cex=2.5,
           las=1) 
abline(h=seq(0.5,44.5,by=1), col="lightgrey")
mtext("(a)",cex=1., side = 3, adj = 0.05, line = -1.8,font=1)

###############################################################################
par(mai=c(0.6,0.2,0.1,0.2))

stripchart( gor.pw$F[gor.pw$F!=0]~gor.pw$Invader[gor.pw$F!=0],              # Data
           method = "stack", # Random noise
           pch = 16,          # Pch symbols
           col = rgb(0,0,0,0.5),           # Color of the symbol
           vertical = FALSE,   # Vertical mode
           add = FALSE,
           cex=1,
           las=1,
           xlab="Fitness differences",yaxt="n",
           xlim=c(-2,1),cex.axis=1.1, cex.lab= 1.3, mgp=c(2.4,0.6,0))  

RY <- gor.pw$RY ; RY[which(RY<0)]<-0

F.dif <- -gor.pw$RY/(gor.pw$F-1)	

gor.F <- 1 - xtabs(F.dif[gor.pw$F!=0]~gor.pw$Invader[gor.pw$F!=0] )        

stripchart(-gor.F/(1-gor.F) ~ names(gor.F ),            # Data
           method = "stack", # Random noise
           pch = 16,          # Pch symbols
           col = rgb(0,0,1,0.85),           # Color of the symbol
           vertical = FALSE,   # Vertical mode
           add = TRUE,
           cex=2,
           las=1) 

stripchart(G.sub$F~G.sub$Species,              # Data
           method = "stack", # Random noise
           pch = "+",          # Pch symbols
           col = rgb(1,0,0,0.85),           # Color of the symbol
           vertical = FALSE,   # Vertical mode
           add = TRUE,
           cex=2.5,
           las=1) 
abline(h=seq(0.5,44.5,by=1), col="lightgrey")
mtext("(b)",cex=1., side = 3, adj = 0.05, line = -1.8,font=1)

###############################################################################

#############################################################################################3
##############################################################################################
################################################################################################
# Load the data
ser.pw <- read.table("Pairwise/Serengeti.txt", sep="\t", header=T, stringsAsFactors=T)

# Use a subset of data from Gorongosa
S.sub <- comm[comm$PA=="Serengeti",]

# Create a sub-matrix with niche differences
par(mai=c(0.6,1.2,0.1,0))

stripchart(ser.pw$N[ser.pw$N!=0] ~ser.pw$Invader[ser.pw$N!=0],              # Data
           method = "stack", # Random noise
           pch = 16,          # Pch symbols
           col = rgb(0,0,0,0.5),           # Color of the symbol
           vertical = FALSE,   # Vertical mode
           add = FALSE,
           cex=1,
           las=1, xlim=c(0,1),
           xlab="Niche differences",cex.axis=1.1, cex.lab= 1.3, mgp=c(2.4,0.6,0)) 

ser.mean <- xtabs(ser.pw$N~ser.pw$Invader) /(length(unique(ser.pw$Invader))-1)         

stripchart(ser.mean~ names(ser.mean),            # Data
           method = "stack", # Random noise
           pch = 16,          # Pch symbols
           col = rgb(0,0,1,0.85),           # Color of the symbol
           vertical = FALSE,   # Vertical mode
           add = TRUE,
           cex=2,
           las=1) 

stripchart(S.sub$N~S.sub$Species,              # Data
           method = "stack", # Random noise
           pch = "+",          # Pch symbols
           col = rgb(1,0,0,0.85),           # Color of the symbol
           vertical = FALSE,   # Vertical mode
           add = TRUE,
           cex=2.5,
           las=1) 
abline(h=seq(0.5,44.5,by=1), col="lightgrey")

mtext("(c)",cex=1., side = 3, adj = 0.05, line = -1.8,font=1)

###############################################################################

par(mai=c(0.6,0.2,0.1,0.2))

stripchart( ser.pw$F[ser.pw$F!=0]~ser.pw$Invader[ser.pw$F!=0],              # Data
           method = "stack", # Random noise
           pch = 16,          # Pch symbols
           col = rgb(0,0,0,0.5),           # Color of the symbol
           vertical = FALSE,   # Vertical mode
           add = FALSE,
           cex=1,
           las=1,
           xlab="Fitness differences",yaxt="n",
           xlim=c(-2,1),cex.axis=1.1, cex.lab= 1.3, mgp=c(2.4,0.6,0))  

F.dif <- -ser.pw$RY/(ser.pw$F-1)	

ser.F <- 1 - xtabs(F.dif[ser.pw$F!=0]~ser.pw$Invader[ser.pw$F!=0] )        


stripchart(-ser.F/(1-ser.F) ~ names(ser.F ),            # Data
           method = "stack", # Random noise
           pch = 16,          # Pch symbols
           col = rgb(0,0,1,0.85),           # Color of the symbol
           vertical = FALSE,   # Vertical mode
           add = TRUE,
           cex=2,
           las=1) 

stripchart(S.sub$F~S.sub$Species,              # Data
           method = "stack", # Random noise
           pch = "+",          # Pch symbols
           col = rgb(1,0,0,0.85),           # Color of the symbol
           vertical = FALSE,   # Vertical mode
           add = TRUE,
           cex=2.5,
           las=1) 
abline(h=seq(0.5,44.5,by=1), col="lightgrey")

mtext("(d)",cex=1., side = 3, adj = 0.05, line = -1.8,font=1)

###############################################################################

#############################################################################################3
##############################################################################################
################################################################################################
# Load the data
lai.pw <- read.table("Pairwise/Laikipia.txt", sep="\t", header=T, stringsAsFactors=T)

# Use a subset of data from Gorongosa
L.sub <- comm[comm$PA=="Laikipia",]

# Create a sub-matrix with niche differences
par(mai=c(0.6,1.2,0.1,0))

stripchart(lai.pw$N[lai.pw$N!=0] ~lai.pw$Invader[lai.pw$N!=0],              # Data
           method = "stack", # Random noise
           pch = 16,          # Pch symbols
           col = rgb(0,0,0,0.5),           # Color of the symbol
           vertical = FALSE,   # Vertical mode
           add = FALSE,
           cex=1,
           las=1, xlim=c(0,1),
           xlab="Niche differences",cex.axis=1.1, cex.lab= 1.3, mgp=c(2.4,0.6,0)) 

lai.mean <- xtabs(lai.pw$N~lai.pw$Invader) /(length(unique(lai.pw$Invader))-1)         

stripchart(lai.mean~ names(lai.mean),            # Data
           method = "stack", # Random noise
           pch = 16,          # Pch symbols
           col = rgb(0,0,1,0.85),           # Color of the symbol
           vertical = FALSE,   # Vertical mode
           add = TRUE,
           cex=2,
           las=1) 

stripchart(L.sub$N~L.sub$Species,              # Data
           method = "stack", # Random noise
           pch = "+",          # Pch symbols
           col = rgb(1,0,0,0.85),           # Color of the symbol
           vertical = FALSE,   # Vertical mode
           add = TRUE,
           cex=2.5,
           las=1) 
abline(h=seq(0.5,44.5,by=1), col="lightgrey")

mtext("(e)",cex=1., side = 3, adj = 0.05, line = -1.8,font=1)

###############################################################################

par(mai=c(0.6,0.2,0.1,0.2))

stripchart( lai.pw$F[lai.pw$F!=0]~lai.pw$Invader[lai.pw$F!=0],              # Data
           method = "stack", # Random noise
           pch = 16,          # Pch symbols
           col = rgb(0,0,0,0.5),           # Color of the symbol
           vertical = FALSE,   # Vertical mode
           add = FALSE,
           cex=1,
           las=1,
           xlab="Fitness differences",yaxt="n",
           xlim=c(-2,1),cex.axis=1.1, cex.lab= 1.3, mgp=c(2.4,0.6,0))  

F.dif <- -lai.pw$RY/(lai.pw$F-1)	

lai.F <- 1 - xtabs(F.dif[lai.pw$F!=0]~lai.pw$Invader[lai.pw$F!=0] )        

stripchart(-lai.F/(1-lai.F) ~ names(lai.F ),            # Data
           method = "stack", # Random noise
           pch = 16,          # Pch symbols
           col = rgb(0,0,1,0.85),           # Color of the symbol
           vertical = FALSE,   # Vertical mode
           add = TRUE,
           cex=2,
           las=1) 

stripchart(L.sub$F~L.sub$Species,              # Data
           method = "stack", # Random noise
           pch = "+",          # Pch symbols
           col = rgb(1,0,0,0.85),           # Color of the symbol
           vertical = FALSE,   # Vertical mode
           add = TRUE,
           cex=2.5,
           las=1) 
abline(h=seq(0.5,44.5,by=1), col="lightgrey")

mtext("(f)",cex=1., side = 3, adj = 0.05, line = -1.8,font=1)

legend("bottomleft", pch=c(16,3),col=c(rgb(0,0,1,0.85),rgb(1,0,0,0.85)),c("Pairwise approx.", "Multispecies"))

###############################################################################
dev.off()


##################################################################
################################################################3
##################################################################3

pt.sz <- 1.2
pt.col.incl <- rgb(0.2,0.8,0.5,0.8)
pt.col.excl <- rgb(0.8,0.2,0.5,1)

png(filename="Figures/Supp/Figure S7.png",width=22,height=8,units="cm",res=300)

par(mfrow=c(1,3))
par(mai=c(0.5,0.5,0.075,0.075))

plot(0,0, xlim=c(0,1), ylim=c(-1,1),type="n",las=1,
	ylab="Fitness differences",xlab="Niche differences",cex.axis=1.1, cex.lab= 1.3, mgp=c(2.4,0.6,0))
polygon(c(-5,5,5),c(-5,5,-5),col=rgb(0,0,0,0.075), border=NA)

excl <- which(gor.pw$N < gor.pw$F)
points(gor.pw$N[gor.pw$N!=0], gor.pw$F[gor.pw$N!=0], col=pt.col.incl, pch=16, cex=pt.sz)
points(gor.pw$N[excl], gor.pw$F[excl], bg=pt.col.excl, pch=21, cex=pt.sz)

labs <- paste (gor.pw$Invader)#,"\n(",gor.pw$Resident,")", sep = "")
text (gor.pw$N[excl],gor.pw$F[excl],labs[excl], cex=0.75, pos=c(2,3,1,3), offset=0.3)
mtext("(a) Gorogosa",cex=1., side = 3, adj = 0.05, line = -1.8,font=1)


##################################################################################
plot(0,0, xlim=c(0,1), ylim=c(-1,1),type="n",las=1,
	ylab="Fitness differences",xlab="Niche differences",cex.axis=1.1, cex.lab= 1.3, mgp=c(2.4,0.6,0))
polygon(c(-5,5,5),c(-5,5,-5),col=rgb(0,0,0,0.075), border=NA)

excl <- which(ser.pw$N < ser.pw$F)
points(ser.pw$N[ser.pw$N!=0], ser.pw$F[ser.pw$N!=0], col=pt.col.incl, pch=16, cex=pt.sz)
points(ser.pw$N[excl], ser.pw$F[excl], bg=pt.col.excl, pch=21, cex=pt.sz)

labs <- paste (ser.pw$Invader)#,"\n(",gor.pw$Resident,")")
text (ser.pw$N[excl],ser.pw$F[excl],labs[excl], cex=0.75, pos=c(1,4,3,3), offset=0.2)

mtext("(b) Serengeti",cex=1., side = 3, adj = 0.05, line = -1.8,font=1)

##################################################################################

plot(0,0, xlim=c(0,1), ylim=c(-1,1),type="n",las=1,
	ylab="Fitness differences",xlab="Niche differences",cex.axis=1.1, cex.lab= 1.3, mgp=c(2.4,0.6,0))
polygon(c(-5,5,5),c(-5,5,-5),col=rgb(0,0,0,0.075), border=NA)

excl <- which(lai.pw$N < lai.pw$F)

points(lai.pw$N[lai.pw$N!=0], lai.pw$F[lai.pw$N!=0], col=pt.col.incl, pch=16, cex=pt.sz)
points(lai.pw$N[excl], lai.pw$F[excl], bg=pt.col.excl, pch=21, cex=pt.sz)


labs <- paste (lai.pw$Invader)#,"\n(",gor.pw$Resident,")")
text (lai.pw$N[excl],lai.pw$F[excl],labs[excl], cex=0.75,, pos=c(2,3), offset=0.2)

mtext("(c) Laikipia",cex=1., side = 3, adj = 0.05, line = -1.8,font=1)

dev.off()
##################################################################################
