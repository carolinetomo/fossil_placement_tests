import os,sys
from mandos import *

def mrca_tips(tt,tip_name,exclude_list):
    clade = []
    for i in tt.leaves():
        #print i.label,tip_name
        if i.label == tip_name:
            for j in i.parent.leaves():
                if j.label != tip_name and j.label not in exclude_list:
                    clade.append(j.label)
    return clade

def tips2node(t,tip_name,clade,exclude_list):
    found = False
    for i in t.iternodes():
        check_clade = []
        for j in i.iternodes():
            if j.label not in exclude_list and j.label != tip_name and j.label != "":
                try:
                    float(j.label)
                except:
                    check_clade.append(j.label)
        if set(clade)==set(check_clade):
            found = True
            return i
    if found == False:
        print "ERROR: couldn't find equivalent node between true and reconstructed trees"
        sys.exit(0) 

def node_dist(tt,it,tip_name,exclude):
    n1 = None
    n2 = None
    clade = mrca_tips(tt,tip_name,exclude)
    truepar = tips2node(it,tip_name,clade,exclude)
    for i in it.iternodes():
        if i.label == tip_name:
            recon = i
    if truepar == recon.parent:
        return 0
    cur2=recon.parent
    cur1=truepar
    pathlen=0
    #print cur1.get_newick_repr()
    #print cur2.get_newick_repr()
    if truepar in [i for i in recon.parent.iternodes()]:
        while 1:
            cur1 = cur1.parent
            if cur1 == recon.parent:
                return pathlen
            pathlen+=1

    elif recon.parent in [i for i in truepar.parent.iternodes()]:
        while 1:
            cur2 = cur2.parent
            if cur2 ==  truepar.parent:
                return pathlen
            pathlen+=1
    else:
        while 1:
            if cur1 == cur2:
                return pathlen
            if cur2 != it:
                cur2 = cur2.parent
                pathlen+=1
            if cur1 != it:
                cur1 = cur1.parent
                pathlen+=1
    if pathlen > 100000:
        print "failed to find MRCA of true and placement nodes"
        sys.exit(0)
    

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "usage: "+sys.argv[0]+ " <inferred tree directory> <true tree>"
        sys.exit(0)

    itdir = sys.argv[1]+"/"
    rfout = open("ALL.euclidean.unwt.rfdist","w")
    rfout.write("trait_type\tunweighted_rf\tweighted_rf\teuclidean_dist\n")
    count = 0
    allDistances = 0.
    dirls = os.listdir(itdir)
    for j in dirls:
        if "MAP.tre" in j:
            spls = j.strip().split(".")
            num = spls[0]
            exclude = "t33,t10,t4,t13,t41".split(",")  #spls[1].split(".")[0].split(",")
            tt = tree_utils2.read_tree(sys.argv[2])
            try:
                itpath = itdir+num+".MAP.tre"
                it = tree_utils2.read_tree(itpath)
            except:
                print "couldn't read tree", itpath
                sys.exit(0)
            dist = []
            #print itpath
            for f in exclude:
                cur = node_dist(tt,it,f,exclude)
                dist.append(cur)
            sum = 0.
            for i in dist:
                sum+=float(i)
            sum=sum/float(len(dist))
            allDistances+=sum
            #print dist
    print "mean distance: ",allDistances/len(dirls)
