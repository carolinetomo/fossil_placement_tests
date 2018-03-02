import os,sys

if len(sys.argv) != 4:
    print "usage: "+sys.argv[0]+ " traitdir ntaxa ntraits"
    sys.exit(0)

indir = sys.argv[1]+"/"
outdir = "headers/"

for fl in os.listdir(indir):
    cur = open(indir+fl,"r")
    outfl = open(outdir+fl,"w")
    outfl.write(sys.argv[2]+"\t"+sys.argv[3]+"\n")
    for line in cur:
        outfl.write(line)
    outfl.close()
