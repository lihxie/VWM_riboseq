#!/usr/bin/python
#ReorderByList.py
#This script re-orders the featureCounts table by a given list
#This is part of the PipeRiboMNase.sh pipeline, to speed up this matching step
#Version: Yu Sun, 2021-11-13

import sys

def Counter():
    fi=open(sys.argv[1],'r')
    flist=open(sys.argv[2],'r')
    fo=open(sys.argv[3],'w')
    
    ExpDict={}
    for line in fi:
        CurrLine=line.strip().split()
        ExpDict[CurrLine[0]]=CurrLine

    for gene in flist:
        CurrGene=gene.strip()
        FullLine=ExpDict[CurrGene]
        fo.write("\t".join(FullLine)+"\n")

    fi.close()
    flist.close()
    fo.close()

if len(sys.argv) != 4:
    print("This script re-orders the featureCounts table first column by a given list")
    print("This is part of the PipeRiboMNase.sh pipeline, to speed up this matching step")
    print("Usage: [ReorderByList.py] [FeatureCountsTable|at least 2 columns] [List] [Output]")
else:
    Counter()
