import sys
import os

if len(sys.argv) != 1:
    print "usage: "+sys.argv[0]+" <fully sampled tree> <fully sampled alignment>"
    sys.exit(0)

alndir = "traits/"
curtree = "36.tre"

count = 0
for fl in os.listdir(alndir):
    spls = fl.strip().split("-")
    pre = fl.strip().split(".")
    rep = spls[0]
    fos = "t33,t10,t4,t13,t41"#spls[1].strip().split(".")[0]#.replace("_",",")
    curaln = alndir+fl
    #curtree = treedir+pre[0]+".tre"
    cmd = "maru -t "+curtree+" -m "+curaln+" -fos "+fos
    cmd += " -gen 1000000 -samp 1000 -pr 250000 -bl 1 -T 1 -W 20 -w flat -o "
    cmd += rep #pre[0] 
    if count % 10 == 0:
        print count
    count += 1
    print(cmd)
    os.system(cmd)

