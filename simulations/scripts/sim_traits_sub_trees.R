require(phytools)
tree = read.tree("ubd.tre")

for (i in 1:100) {
	#samp = sample(tree$tip.label,5)
	#outtr = paste(paste(toString(i),"-",sep=""),paste(paste(samp,collapse=","),".tre",sep=""),sep="")
	outphy = paste(toString(i),".phy",sep="")
	traits=data.frame(fastBM(tree,sig2=1,nsim=50))
	write.table(traits,sep="\t",quot=F,col.names=F,file=outphy)
}

