import os,sys,dendropy
from dendropy.calculate import treecompare

if len(sys.argv) < 3:
    print "usage: "+sys.argv[0]+ " <true tree directory> <inferred tree directory>"
    sys.exit(0)

ttpath = sys.argv[1]
itdir = sys.argv[2]+"/"
rfout = open("ALL.euclidean.unwt.rfdist","w")
rfout.write("trait_type\tunweighted_rf\tweighted_rf\teuclidean_dist\n")

for j in os.listdir(itdir):
    if ".phy.MAP" in j:  
        spls = j.strip().split(".")
        nm = spls[0]
        tree = dendropy.Tree()
        tns = dendropy.TaxonNamespace()
        tt = tree.get_from_path(ttpath,"newick",taxon_namespace=tns)
        it = tree.get_from_path(itdir+j,"newick",taxon_namespace=tns)
        tt.encode_bipartitions()
        it.encode_bipartitions()
        rfout.write(nm+"\t"+str(treecompare.symmetric_difference(tt,it))+"\t"+str(treecompare.weighted_robinson_foulds_distance(tt,it))+"\t"+str(treecompare.euclidean_distance(tt,it))+"\n")


