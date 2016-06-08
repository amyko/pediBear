import numpy as np
import matplotlib.pyplot as plt
import pdb
import os.path



def makeIDFile(inPath, outPath, pedPath):

    #open files
    infile = open(inPath)
    outfile = open(outPath, 'wb')
    outfile.write("index\tname\tsex\tsampled\n")
    i=0
    
    #set
    myset = set()

    for line in infile:
        
        fields = line.split()

        sex = fields[4]
        if(fields[4]=='1'):
            sex = 1
        elif(fields[4]=='2'):
            sex = 0
        else: pdb.set_trace()
        


        outfile.write("%d\t%s\t%d\t%d\n" %(i, fields[1], sex, 1))    
        
        i+=1
        myset.add(fields[1])


        
    infile.close()
    

    
    #all other individuals
    infile = open(pedPath)
    infile.readline()
    for line in infile:
        
        fields = line.split()
        
        for j in range(0,3):
        
            if(fields[j] in myset or fields[j]=='0'): continue

        
            sex = 0
            if(fields[3]=='M'): 
                sex = 1
            elif(fields[3]=='F'):
                sex = 0
            else: pdb.set_trace()
            
            outfile.write("%d\t%s\t%d\t%d\n" %(i, fields[j], sex, 0))    
            
            i+=1
            myset.add(fields[j])

        
    infile.close()
    outfile.close()


if __name__ == "__main__":

    #file names
    inPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/jays/102jays.ped"
    outPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/jays/test"
    pedPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/jays/pedigree.txt"

    makeIDFile(inPath, outPath, pedPath)