import sys,os

if len(sys.argv) < 2:
    print "usage: python " + sys.argv[0]+ " <input directory>"
    sys.exit(0)

indir = sys.argv[1]+"/"

for i in os.listdir(indir):
    spls = i.strip().split(".")
    if spls[-1] == "t":
        treesfl = indir+i
        outtre = indir+".".join(spls[0:-1])+".MAP.tre"
        #cmd = "sumtrees.py --summary-target=mcct --burnin 200 -F newick "+treesfl+" -o "+outtre
        cmd = "sumMAP.py "+treesfl + " 200"
        os.system(cmd)
