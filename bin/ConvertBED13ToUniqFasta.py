#!/usr/bin/python
#ConvertBED13ToUniqFasta.py
#This script converts bed13 to readid-sequence FASTA format, removing duplicates by checking read id.
#This is part of the PipeRiboMNase pipeline.
#Version: Yu Sun, 2021-11-09

import sys

def Counter():
    fi=open(sys.argv[1],'r')
    fo=open(sys.argv[2],'w')
    
    ReadID_Set=set()

    for line in fi:
        CurrLine=line.strip().split()
        ReadID=CurrLine[3]
        Sequence=CurrLine[12]
        if ReadID not in ReadID_Set:
            ReadID_Set.add(ReadID)
            fo.write(">"+ReadID+"\n")
            fo.write(Sequence+"\n")

    fi.close()
    fo.close()

if len(sys.argv) != 3:
    print("This script converts bed13 to readid-sequence FASTA format, removing duplicates by checking read id.")
    print("This is part of the PipeRiboMNase pipeline.")
    print("Usage: [ConvertBED13ToUniqFasta.py] [Data.bed13] [Data.bed13.fa]")
else:
    Counter()
