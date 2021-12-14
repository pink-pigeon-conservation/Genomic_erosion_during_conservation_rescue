# -*- coding: utf-8 -*-
"""
Script to get summary statistics for tree file
"""
import msprime, pyslim, tskit
import numpy as np
import argparse
import allel
import re

print(msprime.__version__)
print(pyslim.__version__)
print(tskit.__version__)
print(allel.__version__)

#### Get arguments
parser = argparse.ArgumentParser()
# tree name
parser.add_argument('--tree',nargs='+',type=str)
# tree name
parser.add_argument('--genTime',nargs='+',type=str)
# seed
parser.add_argument('--seed',nargs='+',type=int)
# generation to print out in slim output
parser.add_argument('--gen',nargs='+',type=str)
# neutral mID to print out in slim output
parser.add_argument('--mID',nargs='+',type=str)
# mutation rate
parser.add_argument('--U',nargs='+',type=str)
# proportion of neutral mutations
parser.add_argument('--neutP',nargs='+',type=str)
# recombination rates file
parser.add_argument('--recRates',nargs='+',type=str)
# recombination positions file
parser.add_argument('--recPos',nargs='+',type=str)
# Ne for recapitation
parser.add_argument('--NeRatio',nargs='+',type=str)
# Ne for recapitation
parser.add_argument('--Outprefix',nargs='+',type=str)

# Ne for recapitation
#parser.add_argument('--sFactor',nargs='+',type=str)


args = parser.parse_args()
TREE=str(args.tree[0])
SEED=int(args.seed[0])
gen=int(args.gen[0])
U=float(args.U[0])
neutP=float(args.neutP[0])
recRates=str(args.recRates[0])
recPos=str(args.recPos[0])
NeRatio=float(args.NeRatio[0])
mID=str(args.mID[0])
genTime=float(args.genTime[0])
Outprefix=str(args.Outprefix[0])
#sFactor=int(args.sFactor[0])


orig_ts = pyslim.load(TREE)
print(f"The orig_ts sequence has {orig_ts.num_trees} trees on a genome of length {orig_ts.sequence_length},"
      f" {orig_ts.num_individuals} individuals, {orig_ts.num_samples} 'sample' genomes,"
      f" and {orig_ts.num_mutations} mutations.")


#U = (2*G)*mu*P
# mu=Uprop/(2*G*P)
#G=geneRegion
#P=proportion of neutral mutations
NoIndv=len(orig_ts.individuals_alive_at(max(orig_ts.individual_times)))
G=orig_ts.sequence_length
mu=U/(2*G)

GENS=orig_ts.slim_generation

orig_max_roots = max(t.num_roots for t in orig_ts.trees())

if orig_max_roots>1:
    #### Recapitate Tree
    # Recombination Map

    with open(recRates) as f:
        RATES=f.readlines()
    RATES = [float(x.strip()) for x in RATES]
    RATES.append(0.0)
    RATES = [x for x in RATES]


    with open(recPos) as f:
        POSITIONS=f.readlines()
    POSITIONS = [int(x.strip()) for x in POSITIONS]
    POSITIONS.insert(0,0)
    POSITIONS[-1]=int(orig_ts.sequence_length)

    recombination_map=msprime.RecombinationMap(POSITIONS, RATES, num_loci=None)

    print('Recapitating '+TREE)
    rts = orig_ts.recapitate(recombination_map = recombination_map, Ne=int(NoIndv/NeRatio))

else:
    print('NO recapitation!! '+TREE)
    orig_max_roots = max(t.num_roots for t in orig_ts.trees())
    print(f"Without recapitation, the max number of roots was {orig_max_roots}")
    rts=orig_ts


# how many roots the tree has before and after recapitation?
recap_max_roots = max(t.num_roots for t in rts.trees())
print(f"Before recapitation, the max number of roots was {orig_max_roots}, "
      f"and after recapitation, it was {recap_max_roots}.")


print('Mutating '+TREE)
muNeutral=(mu*neutP)
mrts = pyslim.SlimTreeSequence(msprime.mutate(rts, rate=(muNeutral)/genTime, random_seed=SEED, keep=True))
print("# mutations: "+str(mrts.num_mutations))
#mrts=orig_ts

t=0
subsample = mrts.individuals_alive_at(t)
subsample_nodes=[mrts.individual(x).nodes for x in subsample] #node numbers for sampled individuals
subsample_nodes=[a for b in subsample_nodes for a in b] #flatten the list
subsample_nodes=np.sort(np.array(subsample_nodes))
o=mrts.simplify(subsample_nodes)
genotypes=np.array(o.genotype_matrix())
mut_count=np.sum(genotypes, axis = 1)

pi_diver=float(o.diversity())
Nc=o.num_individuals
Ne=int(pi_diver/(4*mu))

pi_diver = "{:.2e}".format(pi_diver)


FILEname=Outprefix+"_pyslimParams"+"_gensRun"+str(GENS)+"_genTime"+str(genTime)+"_genOut"+str(gen)+"_mID"+str(mID)+"_Glen"+str(G)+"_U"+str(U)+"_rho"+str(RATES[1])+"_neutP"+str(neutP)+"_pyslimResults"+"_pi"+str(pi_diver)+"_Nc"+str(Nc)+"_Ne" +str(Ne)+"_mutated_fullOutput.txt"
print("output to: "+str(FILEname))
print("### output header")
myfile = open(FILEname, 'w')
FIRST="#OUT: "+str(gen)+" A "+FILEname
myfile.write("%s\n" % FIRST)
myfile.write("%s\n" % "Version: 4")
print("### output populations")
myfile.write("%s\n" % "Populations:")
# if more pops or non-hermafrodite change this!!
POPinfo="p1 "+str(o.num_individuals)+" H"
myfile.write("%s\n" % POPinfo)
print("### output Mutations")
myfile.write("%s\n" % "Mutations:")

count=-1
for m in o.variants():
    count=count+1
    ID=str(m.site.id)
    ID2=str(count)
    POS=str(round(m.site.position))
    POP="p1"
    ORG=str(gen-1)
    COUNT=str(mut_count[count])
    if(len(m.site.mutations[0].metadata)>0):
        mType=str(1)
        s=str(0)
        h=str(0.5)
    else:
        mType=str(1)
        s=str(0)
        h=str(0.5)
    MUT=ID+" "+ID2+" "+mID+" "+POS+" "+s+" "+h+" "+POP+" "+ORG+" "+COUNT
    myfile.write("%s\n" % MUT)

print("### output Individuals")

myfile.write("%s\n" % "Individuals:")

indv=[]
count=-1
for i in mrts.individuals_alive_at(t):
    count=count+1
    INDV=mrts.individual(i)
    ID=INDV.id
    NODES=INDV.nodes
    if len(NODES) ==2:
        AGE=INDV.metadata['age']
        INDV="p1:i"+str(ID)+" H "+"p1:"+str(NODES[0])+" "+"p1:"+str(NODES[1])+" "+str(AGE)
        myfile.write("%s\n" % INDV)

print("### output Genomes")
myfile.write("%s\n" % "Genomes:")
mut=[]
count=-1
for SNP in genotypes:
    count=count+1
    #print(SNP)
    SNP=['NA' if x == 0 else count for x in SNP]
    mut.append(SNP)

mut=np.transpose(mut)
count=-1
for line in mut:
    count=count+1
    line=line[line != 'NA']
    line=line.tolist()
    line=' '.join(line)
    line="p1"+":"+str(count)+" "+"A"+" "+line
    myfile.write("%s\n" % line)

myfile.close()