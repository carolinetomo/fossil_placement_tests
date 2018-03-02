#!/usr/bin/env python2
import sys
import mandos

def replace_string_index(string,index,newval):
    stringlist = []
    for i in string:
        stringlist.append(i)
    stringlist[index] = newval
    return "".join(stringlist)

def decomp_tree(tree,posdic):
    biparts = {} 
    for node in tree.iternodes(order=1):
        if len(node.children) ==0 :
            continue
        if node != tree:
            bs = make_empty_bitstring(tree)            
            l = [i.label for i in node.leaves()] 
            #r = [i.label for i in tree.iternodes() if i.label not in l]
            for tip in l:
                bs = replace_string_index(bs,posdic[tip],"1")
            biparts[bs] =node
    return biparts

def make_empty_bitstring(tree):
    bs = ""
    for node in tree.iternodes():
        if node.istip:
            bs+="0"
    return bs

def get_positions(tree):
    posdic = {}
    count = 0
    for node in tree.iternodes():
        if node.istip:
            posdic[node.label] = count
            count+=1
    return posdic

def flip_bits(bitstring):
    newbs = ""
    for i in bitstring:
        if i == "0":
            newbs+="1"
        elif i == "1":
            newbs+="0"
    return newbs

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print "usage: "+sys.argv[0]+" <posterior trees> <burnin>"
        sys.exit(0)

    burnin = int(sys.argv[2])
    all_trees = []
    tcounts = []
    nstrings= []
    int_lens = []
    tip_lens = {}
    trfl = open(sys.argv[1],"r")
    count = 0
    posdic = {}
    for line in trfl:
        #if count > 2000:
        #    continue       
        if count == 0:
            count +=1
            curtree = mandos.tree_reader2.read_tree_string(line.strip())
            posdic = get_positions(curtree)
            all_trees.append(decomp_tree(curtree,posdic))
            tcounts.append(1)
            nstrings.append(line.strip())
            continue
        if count < burnin:
            count += 1
            if count %1000 == 0:
                print "burning in (ignoring) tree",count
            continue
        curtree = mandos.tree_reader2.read_tree_string(line.strip())
        dtree = decomp_tree(curtree,posdic)
        #print dtree
        #print all_trees
        found = False
        for ind,seen in enumerate(all_trees):
            good = True
            for bp in dtree.keys():
                try:
                    seen[bp]
                except:
                    try:
                        seen[flip_bits(bp)]
                    except:
                        good = False
                        break
                if good == False:
                    break
            if good == True:
                found = True
                tcounts[ind] += 1
        if found == False:
            all_trees.append(dtree)
            nstrings.append(line.strip())
            tcounts.append(1)
        #print len(all_trees)
        #print tcounts

        if count%1000 == 0:
            print "processing tree",count
        count += 1

    high = 0
    highind = 0
    trc = 0
    for ind,val in enumerate(tcounts):
        trc += int(val)
        if val > high:
            high = val
            highind = ind
    
    print "mapping support values to MAP topology"
    bestbp = all_trees[highind]
    nodecounts = {} 

    for ind,t in enumerate(all_trees):
        for mapbp in bestbp:
            found = False           
            good = True
            try:
                t[mapbp]
            except:
                try:
                    t[flip_bits(mapbp)]
                except:
                    good = False
            if good == True:
                found = True
                try:
                    nodecounts[mapbp] += tcounts[ind]
                except:
                    nodecounts[mapbp] = tcounts[ind]
    
    #print nodecounts
    besttree = mandos.tree_reader2.read_tree_string(nstrings[highind])
    dtree = decomp_tree(besttree,posdic)
    for i in dtree.keys():
        pp = nodecounts[i]/float((count-1)-burnin)
        if pp > 1.0: #fix when pp goes > 1.0 due to rounding error
            pp = 1.0
        dtree[i].label = str(pp)

    flnm = sys.argv[1].split(".")[0]+".MAP.tre"
    outtr=open(flnm,"w")
    outtr.write(besttree.get_newick_repr(True)+";\n")
    outtr.close()
    print "wrote results to "+flnm
    #print tcounts

