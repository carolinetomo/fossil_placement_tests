require(phytools)
tree = read.tree("ubd.tre")

for (i in 1:100) {
	outphy = paste(toString(i),".phy",sep="")
	traits=data.frame(fastBM(tree,sig2=1,nsim=50))
	r = multi2di(di2multi(tree,tol=10.),rand=T)
	for (i in 1:length(r$edge.length)) {r$edge.length[i] = rexp(1)}
	traits1=data.frame(fastBM(r,sig2=1,nsim=10))
	join =merge(traits,traits1,by = 0)
	write.table(join,sep="\t",quot=F,col.names=F,row.names=F,file=outphy)
}

