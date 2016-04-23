import numpy as np
import matplotlib.pyplot as plt
import pdb
import os.path



def makeMapFile(inPath, outPath):

    #open files
    outfile = open(outPath,'wb')
    marker = 1
    
    for chrom in range(1,23):
    
        #open files
        infile = open(inPath+str(chrom))
        #skip header
        infile.readline()
        
        for line in infile:
            
            #pos
            pos = int(line.split()[0])
            pos = pos*1.3/1e6
            
            #write
            outfile.write('%d\tm%d\t%f\n' %(chrom, marker, pos))
            marker+=1
            
            
        infile.close()
    
    outfile.close()
    
    
def makeFreqFile(inPath, outPath):

    #open files
    outfile = open(outPath,'wb')
    marker = 1
    
    for chrom in range(1,23):
    
        # open file
        infile = open(inPath+str(chrom))
        
        #skip header
        infile.readline() 
          
        for line in infile:
            
            #freq
            fields = line.split()
            f1 = float(fields[6])
            f2 = float(fields[9])
            
            #write
            outfile.write('M\tm%d\n' %marker)
            outfile.write('%f\t%f\t%f\t%f\n' %(f1,0,0,f2))
            
            marker+=1
            
            
        infile.close()

    outfile.close()


def makeDatFile(inPath, outPath):

        #open files
    outfile = open(outPath,'wb')
    marker = 1
    
    for chrom in range(1,23):
    
        # open file
        infile = open(inPath+str(chrom))
        
        #skip header
        infile.readline() 
          
        for line in infile:
            
            #write
            outfile.write('M\tm%d\n' %marker)
            
            marker+=1
            
            
        infile.close()

    outfile.close()
    
    return marker
    


def makePedFile(inPath, outPath, numSnp):
    
    #open file
    outfile = open(outPath, "wb")
    
    # missing individuals
    for indiv in range(1,6):
    
        #write info
        outfile.write("%d\t%d\t%d\t%d\t%d\t" %(1,indiv,parents[indiv][0],parents[indiv][1],indiv%2))
    
        #write geno
        for i in range(0,numSnp):
            outfile.write("%d/%d\t" %(0,0))
        outfile.write("\n")
        
    
    # typed individuals
    indivIndex = 0
    for indiv in range(6,9):
        
        #write info
        outfile.write("%d\t%d\t%d\t%d\t%d\t" %(1,indiv,parents[indiv][0],parents[indiv][1],indiv%2))
        
        #write geno
        for chrom in range(1,23):
            
            infile = open(inPath+str(chrom))
            infile.readline()
            
            for line in infile:
                
                geno = line.split()[indivIndex]
                
                if geno[0]=='A':
                    a1 = 1
                else:
                    a1 = 4
                if geno[1]=='A':
                    a2 = 1
                else:
                    a2 = 4
                    
                outfile.write("%d/%d\t" %(a1,a2))
                
            infile.close()
            
        outfile.write("\n")
        indivIndex+=1
            
    outfile.close()
                    
                
        


if __name__ == "__main__":

    #file names
    infoPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/unrelated/msprime.info.pruned."
    genoPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/simulations/sim.test.geno."
    outPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/merlin/test2"

        
    makeMapFile(infoPath, outPath+str(".map"))
    makeFreqFile(infoPath, outPath+str(".freq"))
    numSnp = makeDatFile(infoPath, outPath+str(".dat"))
    
    #parents = {1:[0,0], 2:[0,0], 3:[1,2], 4:[1,2], 5:[1,2], 6:[0,3], 7:[0,4], 8:[0,5], 9:[0,6], 10:[0,7], 11:[0,8]}
    parents = {1:[0,0], 2:[0,1], 3:[0,1], 4:[0,1], 5:[0,2], 6:[0,3], 7:[0,4], 8:[0,5]}
    makePedFile(genoPath, outPath+str(".2.ped"), numSnp)

