setwd("C:/Users/Falko/Documents/Standalone Research/Mammal Isotopes/Most recent manuscript and files")


# Load the data
comm <- read.table("Raw_data/Laikipia_data.txt", sep="\t", header=T, stringsAsFactors=T)

# Create a vector of species names
sp.name <- comm$Species


#########################################
# Create a grid with every possible combination of 12 species 
assemblages <- expand.grid(0:1, 0:1, 0:1, 0:1, 0:1, 0:1, 0:1, 0:1, 0:1, 0:1, 0:1, 0:1)

# A numeric vector with integer species IDs
sp.id <- 1:dim(comm)[1]
# The species richness of every combination of species
richness <- rowSums(assemblages)
# A blank vector to track whether all species in an assemblage have positive equilibrium densities
equilibrium <- rep(NA,dim(assemblages)[1])


# The maximum possible (daily) growth rate: Max birth rate minus natural mortality rate
R.vect <- (comm$Rmax - comm$M)/365
# The consumption matrix: the proportional dietary contribution of plant species, multiplied by the daily matabolic intake requirement
cons.mat <- comm[,5:dim(comm)[2]]*(0.05*comm$BM^0.77)
# A vector of species' body masses
BM.vect <- comm$BM

# Parameter b_i; the efficiency at which a kg of plant consumption is converted into growth rate (assuming 0.1 metabolic scaling across trophic levels)
bi <- 0.1*(1/BM.vect)

# Just a time to process how long the simulation takes
ptm <- proc.time()

# The matrix of the sum of the products of consumption between species i and j.
U.mat <- matrix(NA,nrow=dim(cons.mat)[1],ncol=dim(cons.mat)[1])
	for (i in 1:dim(cons.mat)[1]){
		for (j in 1:dim(cons.mat)[1]){
			U.mat[i,j] <- sum(cons.mat[i,]*cons.mat[j,])
		}
	}

# The consumption matrix, multiplied by parameter b_i; written to a file in Intermediary results for validation purposes
U.mat.b <- (U.mat * bi)
colnames(U.mat.b) <- sp.name
rownames(U.mat.b) <- sp.name
write.table(U.mat.b,file= "Intermediate_data/U_matrix_Laikipia.txt",quote=T,sep="\t",row.names=T,col.names=NA)


# A loop to explore whether combinations of species can coexist stably

for (k in 2:dim(assemblages)[1]){
	id <- which(assemblages[k,]==1)

# Calculate the equilibrium densities
	N.star <- solve(U.mat.b[id,id],R.vect[id])
# First condition: if all species have positive equilibrium densities
	if(all(N.star>0)){

# Second condition: community is not stable if the excluded species can invade (positive invasion growth rate)
	invaders <- which(assemblages[k,]==0)
	if (length(invaders)>0){
		igr <- rep(NA,length(invaders))
		# Calculate invasion growth rate for all excluded species
		for (n in 1:length(invaders)){
			igr[n] <- R.vect[invaders[n]] - bi[invaders[n]]*sum((t(t(cons.mat[id,])*unlist(cons.mat[invaders[n],]))*N.star))

		}
		# If assemblage has positive equilbrium density and none of the remaining species can invase, the assemblage is stable
		equilibrium[k] <- ifelse(any(igr>0),0,1)	
		}
	} 
	else{equilibrium[k] <- 0} # Otherwise, assemblage is unstable
# A progress tracker
print(k)
}
# Report how long the algorithm ran
proc.time() - ptm

# Determin the species richness of the stable community with the highest richness (assuming there are more than 1 stable assemblages)
rich.max <- max(richness[which(equilibrium==1)])

##########################################################################################
##########################################################################################
# Calculating the niche- and fitness differences

# Identify the stable commuity with the highest species richeness
stabcom <- which(richness==rich.max & equilibrium==1)[1]

residents <- which(assemblages[stabcom,]==1) # A vector of residents in the stable community
invaders <- which(assemblages[stabcom,]==0) # A vector of potential (unsuccessful) invaders

# Equilbrium density for stable community
N.star <- solve(U.mat.b[residents,residents],R.vect[residents])

# Conversion factors: (1) based on consumption rates, and (2) based on rations of total metabolic requirments
	c.fact <- matrix(NA,ncol=dim(cons.mat)[1],nrow=dim(cons.mat)[1])
	c.fact2 <- matrix(NA,ncol=dim(cons.mat)[1],nrow=dim(cons.mat)[1])

	for (i in 1:dim(cons.mat)[1]){
		for (j in 1:dim(cons.mat)[1]){
			c.fact[j,i] <- sqrt(sum(cons.mat[j,]*cons.mat[j,])/sum(cons.mat[i,]*cons.mat[i,]))
			c.fact2[j,i] <- (0.05*BM.vect[j]^0.77)/(0.05*BM.vect[i]^0.77)
		}
	}

# Save a plot of the two conversion factors. This is only an intermediate validation step to ensure the C-factors make sense
png(filename="Intermediate_data/C-factors_Laikipia.png",width=16,height=16,units="cm",res=300)
	plot(c.fact2,c.fact, ylab="Calculated c-factors",log="xy" , las=1,xlab="Total dietary consumption ratios"); abline(a=0,b=1, col="red")
dev.off()
# Save C-factors to intermediate file for validation checks
colnames(c.fact) <- sp.name
rownames(c.fact) <- sp.name
write.table(round(c.fact,3),file= "Intermediate_data/c_conversion_Laikipia.txt",quote=T,sep="\t",row.names=T,col.names=NA)

# Blank vectors to save ivasion growth rates (igr) and no-niche growth rate (nngr)
	igr <- rep(NA,length(invaders))
	nngr <- rep(NA,length(invaders))

# For each potential invader, calcualte invasion and no-niche growth rates
	for (n in 1:length(invaders)){
		igr[n] <- R.vect[invaders[n]] - 
			bi[invaders[n]]*
			sum((t(t(cons.mat[residents,])*unlist(cons.mat[invaders[n],]))*N.star))

		nngr[n] <- R.vect[invaders[n]] - 
			bi[invaders[n]]*
			sum(colSums(matrix(rep(unlist(cons.mat[invaders[n],])^2,length(residents)),nrow=length(residents), byrow=F)*
			N.star*c.fact[residents,invaders[n]]))
	}

# Calculate the niche- and fitness difference for the potential invaders.

N.inv <- (igr-nngr)/(R.vect[invaders]-nngr)
F.inv <- (-nngr/R.vect[invaders])/(1-(nngr/R.vect[invaders]))

#######################################################
# To calcualte the niche- and fitness differences of the resident community, we:
# 1. Incrementally remove each species
# 2. Re-calcuate the equilbrium densities for the remaining assemblage
#	2.1. If any species in the sub-cummunity has a negative equilibrium density (due to priority effects), set its density to zero 
# 3. Calcculate N and F the same way as before

# Blank vectors for N and F
N.res <- rep(NA,length(residents))
F.res <- rep(NA,length(residents))

for (m in 1:length(residents)){
# Incrementally remove one species as an invader
	residents.sub <- residents[-m]
	invader.sub <- residents[m]
# Recalcuate equilibrium densities for sub-community
	N.star.sub <- solve(U.mat.b[residents.sub,residents.sub],R.vect[residents.sub])
	# If any species have negative denisties, replace with a zero (only happens in rare instances where there are priorty effects)
	N.star.sub[which(N.star.sub<0)] <- 0

	# Calculate invasion grwoth rate
	igr <- R.vect[invader.sub] - 
		bi[invader.sub]*
		sum((t(t(cons.mat[residents.sub,])*unlist(cons.mat[invader.sub,]))*N.star.sub))

	#Calcualte non-niche growth rate
	nngr <- R.vect[invader.sub] - 
		bi[invader.sub]*
		sum((matrix(rep(unlist(cons.mat[invader.sub,])^2,length(residents.sub)),nrow=length(residents.sub), byrow=F)*
			N.star.sub*c.fact[residents.sub,invader.sub]))

# Calcualte niche- and fitness differneces	
N.res[m] <- (igr-nngr)/(R.vect[invader.sub]-nngr)
F.res[m] <- (-nngr/R.vect[invader.sub])/(1-(nngr/R.vect[invader.sub]))

}




# Make a plot of the results 
cols2 <- colorRampPalette(c(rgb(0,0,0.5,0.25),rgb(0,0.8,0.8,0.25),rgb(1,0.5,0,0.25)),interpolate="linear")(dim(comm)[1])

png(filename="Processed_data/Laikipia_coexistence.png",width=16,height=16,units="cm",res=300)
par(mai=c(0.8,0.8,0.075,0.075))

plot(N.inv,F.inv, xlim=c(0,1), ylim=c(0,1),type="n",las=1,
	ylab="Fitness differences",xlab="Niche differences",cex.axis=1.1, cex.lab= 1.3, mgp=c(2.0,0.6,0))
polygon(c(-5,5,5),c(-5,5,-5),col=rgb(0,0,0,0.075), border=NA)

points (N.inv,F.inv,pch=15,col=cols2[invaders], cex=1.3)
points (N.inv,F.inv,pch=0, cex=1.3)

text (N.inv,F.inv,sp.name[invaders], pos=c(2,2,4,2,2,1,1,3), cex=0.6)

points (N.res,F.res,pch=16,col=cols2[residents], cex=1.3)
points (N.res,F.res,pch=1, cex=1.3)
text (N.res,F.res,sp.name[residents], pos=c(4,4,4,2,4,3), cex=0.6)

legend("topleft", pch=21, col="black",pt.bg=cols2,legend=sp.name,pt.cex=1.3, bty="n",bg='n')
text(0.4,0.05,"Laikipia", cex=1.5)
dev.off()

# Correlation test for niche- and fitness differences
cor.test(c(N.inv,N.res),c(F.inv,F.res))
