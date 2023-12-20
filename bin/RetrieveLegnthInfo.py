#!/usr/bin/python
#RetrieveLegnthInfo.py
#This script takes a table as a length matching information as input, and add length info to the last column
#You can designate which column as the key to match
#This is part of the PipeRiboMNase pipeline
#Version: Yu Sun, 2021-11-10

import sys

def Counter():
    fi=open(sys.argv[1],'r')
    KeyColumn=int(sys.argv[2])
    Matching=open(sys.argv[3],'r')
    fo=open(sys.argv[4],'w')
    
    MatchingDict={}
    for line in Matching:
        Currline=line.strip().split()
        MatchingDict[Currline[0]]=Currline[1]

    for line in fi:
        CurrData=line.strip().split()
#        if CurrData[KeyColumn-1] in MatchingDict.keys():    #This if statement will be very slow!
        try:
            RetrieveLen=MatchingDict[CurrData[KeyColumn-1]]  #Try to catch this error
            CurrData.append(RetrieveLen)
            fo.write("\t".join(CurrData)+"\n")
        except KeyError:
            pass

    fi.close()
    fo.close()

if len(sys.argv) != 5:
    print("This script takes a table as a length matching information as input, and add length info to the last column")
    print("You can designate which column as the key to match")
    print("This is part of the PipeRiboMNase pipeline")
    print("Usage: [RetrieveLegnthInfo.py] [Data] [KeyColumn] [LengthTable] [Output]")
else:
    Counter()
