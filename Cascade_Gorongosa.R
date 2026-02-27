install.packages(c("bipartitie", "vioplot", "RColorBrewer"))

# Load the data
comm <- read.table("Raw_data/Gorongosa_data.txt", sep="\t", header=T, stringsAsFactors=T)

#Vector of species names
sp.name <- comm$Species

# Total number of plant resources
vegS <- dim(comm)[2]-4

#########################################
# The resource levels to simulate
resource.level <- seq(10,140,by=10)

# The number of interations
no.iter <- 100

# A blank matrix to hold the outputs
rich.mat <- matrix(NA,nrow=no.iter, ncol=length(resource.level))
rownames(rich.mat) <- 1:no.iter
colnames(rich.mat) <- resource.level

####################################################
D.mat.p <- matrix(NA,nrow=no.iter, ncol=length(resource.level))
rownames(D.mat.p) <- 1:no.iter
colnames(D.mat.p) <- resource.level
D.mat.m <- D.mat.p

library(bipartite)
####################################################

# A blank matrix to hold the species presences
sp.prob <- matrix(NA,nrow=length(resource.level),ncol=dim(comm)[1])

# All possible combiantions of species
assemblages <- expand.grid(0:1, 0:1, 0:1, 0:1, 0:1, 0:1, 0:1, 0:1, 0:1, 0:1, 0:1)

# Run a loop for each resource level
for(lev in 1:length(resource.level)) {
	sp.mat <- matrix(0,nrow=no.iter,ncol=dim(comm)[1])

	# Run a loop for each iteration
	for (iter in 1:no.iter) {
		# A vector of intergers species IDs
		sp.id <- 1:dim(comm)[1]
		# The richness of each assemblage
		richness <- rowSums(assemblages)
		# Vector to record whether the assemblage is stable or not
		equilibrium <- rep(NA,dim(assemblages)[1])

		# Sample plant resources to the resource level, without replacement
		res <- sample(1:vegS,resource.level[lev], replace=FALSE)

		# The possible plant species
		veg <- (comm[,5:dim(comm)[2]])
		# The sub-sample of plant resources
		veg.ext <- veg[,res]
		# Standardise proportional consumption, so that each species consumtion sums to one
		veg.stand <- as.matrix(veg.ext/rowSums(veg.ext))
		# If a species has no suitable food, set values to zero (so standardisation does not divide by zero)
		veg.stand[is.nan(veg.stand)] <- 0

		# A vector of maximum daily growth rates (maximum reproduction minus natural mortality)
		R.vect <- (comm$Rmax - comm$M)/365
		# Multiply the standardised consumption proportions by the daily intake
		cons.mat <- veg.stand*(0.05*comm$BM^0.77)

		D.mat.p[iter,lev] <- mean(dfun(t(cons.mat))[[1]])
		D.mat.m[iter,lev] <- mean(dfun((cons.mat))[[1]])

		# Vector of body mass
		BM.vect <- comm$BM
		# The parameter for conversion of 1kg of plant biomass into growth rate
		bi <- 0.1*(1/BM.vect)

		# The species interaction matrix (sum of the product of consumption rates)
		U.mat <- matrix(NA,nrow=dim(cons.mat)[1],ncol=dim(cons.mat)[1])
			for (i in 1:dim(cons.mat)[1]){
				for (j in 1:dim(cons.mat)[1]){
					U.mat[i,j] <- sum(cons.mat[i,]*cons.mat[j,])
				}
			}		
		# Multiple interaction matrix by conversion rate, b_i
		U.mat.b <- t(U.mat * bi)

		# Run a loop for each combination of species
		for (k in 2:dim(assemblages)[1]){
			id <- which(assemblages[k,]==1)
			# If there is only on species, assume it can exist
			if(richness[k]<=1) {equilibrium[k] <- 1} else{
			# If any species in the assembalge have no food resosurces, then the community is not stable
			if(any(colSums(U.mat.b[id,id])==0)) {equilibrium[k] <- 0} else{
			# If the determinant of the consumption matrix is zero, the community is not stable
			if(det(U.mat.b[id,id])==0) {equilibrium[k] <- 0} else{
			# Calcualte equlibrium densities
			N.star <- as.vector((R.vect[id]) %*% solve(U.mat.b[id,id], tol = 1e-170))
			# Condition 1: all species have positive equilibrium densities
			if(all(N.star>0)){
			# Condition 2: None of the other species, not in the commuity, are able to invade
			invaders <- which(assemblages[k,]==0)
			if (length(invaders)>0){
			# Blank vector for ivasion growth rates
			igr <- rep(NA,length(invaders))
			# Loop to calcualte invasion growth rates
				for (n in 1:length(invaders)){
					igr[n] <- R.vect[invaders][n] - sum((U.mat.b)[id,invaders[n]]*N.star)
				}
			equilibrium[k] <- ifelse(any(igr>0),0,1)	
			}
		} else {equilibrium[k] <- 0}
		
		}
	}
}
}
# Print the progress of the loop
print(paste("Iteration =",iter,"; Level =",lev))

# Record the richness for the combination of resource level and simulation iteration
rich.mat[iter,lev] <- max(richness[which(equilibrium==1)])

# Record the assemblage in the stable community for the combination of resource level and simulation iteration
sp.mat[iter,] <- apply(as.matrix(assemblages[which(richness==max(richness[which(equilibrium==1)]) & equilibrium==1),]),2,mean)

}
# Persistence probability of each species
sp.prob[lev,] <- colSums(sp.mat)
}

# Write outputs to file
write.table(rich.mat,file= "Processed_data/Gorongosa_SpRichIter100.txt",quote=T,sep="\t",row.names=F,col.names=T)

colnames(sp.prob) <- sp.name
write.table(sp.prob,file= "Processed_data/Gorongosa_SpProbIter100.txt",quote=T,sep="\t",row.names=F,col.names=T)

library(vioplot)
library(RColorBrewer)

# Set up the plot and dimensions
#png(filename="Gorongosa_cascade.png",width=28,height=12,units="cm",res=300)

# Define plot margins
par(mfrow=c(1,2))
par(mai=c(.8,.8,0.1,0.1))

# Make the base violin plot of the simulated data (replicates are iterations)
X <- vioplot(rich.mat, ylim=c(0,11),las=1, main="", ylab="Maximum species richness", col=rgb(0.8,0.8,0.8,1),
	xlab="Number of resources",cex.axis=1, cex.lab= 1,border=NA, xaxt="n",wex=1.2, rectCol="black",
	lineCol="black", colMed="red", cexMed=1.5, xlim=c(-.15,14.5))
axis(1, at =seq(0,14,by=2), labels=c("0","20","40","60","80","100","120","140"))
points(c(0,seq(1,14,by=1),14.4),c(0,X$median,8),pch=16,col=rgb(0.7,0,0,1), cex=1.3)
mtext("(a)",cex=1.3, side = 3, adj = 0.03, line = -1.25,font=1)

IV <- c(0,seq(1,14,by=1),14.4)
DV <- c(0,X$median,8)
nlreg <- nls(DV~a*IV^b,start=list(a=8,b=0.5))
xvals <- seq(0,14.4,l=100)
lines(xvals, coef(nlreg)[1]*xvals^coef(nlreg)[2], lwd=2, col=rgb(0.8,0,0,0.7), lty=1)


IV <- c(0,seq(1,14,by=1),14.4)
DV <- c(0,X$q1,8)
nlreg <- nls(DV~a*IV^b,start=list(a=8,b=0.5))
xvals <- seq(0,14.4,l=100)
lines(xvals, coef(nlreg)[1]*xvals^coef(nlreg)[2], lwd=1.5, col=rgb(0.8,0,0,0.7), lty=2)

IV <- c(0,seq(1,14,by=1),14.4)
DV <- c(0,X$q3,8)
nlreg <- nls(DV~a*IV^b,start=list(a=8,b=0.5))
xvals <- seq(0,14.4,l=100)
lines(xvals, coef(nlreg)[1]*xvals^coef(nlreg)[2], lwd=1.5, col=rgb(0.8,0,0,0.7), lty=2)


plot(0,0,type="n",xlab="Number of resources", ylab="Presence probability", xlim=c(0,140), 
	ylim=c(0,1.15),las=1, cex.axis=1,cex.lab=1)

#display.brewer.all()
pres <- c(0,1,0,1,1,1,0,1,1,1,1)
cols2 <- brewer.pal(n = 11, name = "Paired")
for (nsp in 1:11){
	lines(c(0,resource.level,vegS),c(0,sp.prob[,nsp]/no.iter,pres[nsp]), col=cols2[nsp], lwd=2)
	points(c(0,resource.level,vegS),c(0,sp.prob[,nsp]/no.iter,pres[nsp]), col=cols2[nsp], pch=16,cex=0.75)
}
legend("topright",cex=0.7,ncol=4,lwd=2,legend=sp.name,col=cols2,lty=1, bty="n")
mtext("(b)",cex=1.3, side = 3, adj = 0.03, line = -1.25,font=1)

#dev.off()


res.mat <- matrix(rep(resource.level,no.iter), nrow=no.iter,byrow=T)

cols2 <- brewer.pal(n = 12, name = "Paired")
cols2 <- c(cols2,cols2)

png(filename="Figures/Supp/Figure S14.png",width=34,height=20,units="cm",res=300)
par(mai=c(.7,.7,0.25,0.1))

par(mfrow=c(3,5))
for (j in 1:14){
  
  trend2 <- glm(as.vector(rich.mat[which(rich.mat[,j]>1),j])~as.vector(D.mat.p[which(rich.mat[,j]>1),j]),family="poisson")
  
 plot(D.mat.p[,j],rich.mat[,j],  pch=21, bg=cols2[res.mat[,j]/10], cex=1.5,las=1, ylab="Herbivore richness", 
       xlab="Mean Bluthgen d' (Plants)", ylim=c(0,10), xlim=c(0,0.5), main=paste("Number of plant resources =",mean(res.mat[,j])))
  
  lines(seq(0,0.6,l=100), (exp(coef(trend2)[1] + (coef(trend2)[2]*seq(0,0.6,l=100)))),lwd=2,lty=ifelse(coef(summary(trend2))[2,4]<0.05,1,2), col=cols2[res.mat[,j]/10])
  
}
dev.off()


png(filename="Figures/Supp/Figure S11.png",width=34,height=20,units="cm",res=300)
par(mai=c(.7,.7,0.25,0.1))

par(mfrow=c(3,5))
for (j in 1:14){ 
  
  trend2 <- glm(as.vector(rich.mat[which(rich.mat[,j]>1),j])~as.vector(D.mat.m[which(rich.mat[,j]>1),j]),family="poisson")
  
  plot(D.mat.m[,j],rich.mat[,j],  pch=21, bg=cols2[res.mat[,j]/10], cex=1.5,las=1, ylab="Herbivore richness", 
       xlab="Mean Bluthgen d'(Herbivores)", ylim=c(0,10), xlim=c(0,0.5), main=paste("Number of plant resources =",mean(res.mat[,j])))
  
  lines(seq(0,0.6,l=100), (exp(coef(trend2)[1] + (coef(trend2)[2]*seq(0,0.6,l=100)))),lwd=2,lty=ifelse(coef(summary(trend2))[2,4]<0.05,1,2), col=cols2[res.mat[,j]/10])
  
}
dev.off()


png(filename="Figures/Supp/Figure S8.png",width=28,height=14,units="cm",res=300)
par(mfrow=c(1,2))
par(mai=c(.8,.8,0.1,0.1))
vioplot(D.mat.p~res.mat, las=1, ylab="Mean Bluthgen 'd (Plants)", xlab="Number of plant resources")
mtext("(a) Plant specialisation",cex=1.3, side = 3, adj = 0.95, line = -1.25,font=1)
vioplot(D.mat.m~res.mat, las=1, ylab="Mean Bluthgen 'd (Herbivores)", xlab="Number of plant resources")

mtext("(b) Herbivore specialisation",cex=1.3, side = 3, adj = 0.95, line = -1.25,font=1)
dev.off()

