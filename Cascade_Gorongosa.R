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
		res <- sample(1:vegS,resource.level[lev], replace=F)

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
					igr[n] <- R.vect[invaders][n] - sum(t(U.mat.b)[id,invaders[n]]*N.star)
				}
			equilibrium[k] <- ifelse(any(igr>0),0,1)	
			}
		} else {equilibrium[k] <- 0}
		# Print the progress of the loop
		print(k)
		}
	}
}
}

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

