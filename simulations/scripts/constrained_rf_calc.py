import os,sys,dendropy
from dendropy.calculate import treecompare

if len(sys.argv) < 3:
    print "usage: "+sys.argv[0]+ " <true tree directory> <inferred tree directory>"
    sys.exit(0)

tt = sys.argv[1]
itdir = sys.argv[2]+"/"
rfout = open("ALL.euclidean.unwt.rfdist","w")
rfout.write("trait_type\tunweighted_rf\tweighted_rf\teuclidean_dist\n")

for j in os.listdir(itdir):
    if ".tre" in j:  
        spls = j.strip().split(".")
        num = spls[0]
        tree = dendropy.Tree()
        tns = dendropy.TaxonNamespace()
        tt = tree.get_from_path(tt,"newick",taxon_namespace=tns)
        try:
            it = tree.get_from_path(itdir+num+".MCC.tre","newick",taxon_namespace=tns)
        except:
            continue
        tt.encode_bipartitions()
        it.encode_bipartitions()
        rfout.write(num+"\t"+str(treecompare.symmetric_difference(tt,it))+"\t"+str(treecompare.weighted_robinson_foulds_distance(tt,it))+"\t"+str(treecompare.euclidean_distance(tt,it))+"\n")


