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
            outfile.write('F\t%f\t%f\t%f\t%f\n' %(f1,0,0,f2))
            
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
    


def makePedFile(inPath, outPath, numSnp, nUnrelated, nSampled):
    
    #open file
    outfile = open(outPath, "wb")
    
    # missing individuals
    for indiv in range(1,nUnrelated+1):
    
        #write info
        sex = 0
        #sex = indiv%2
        outfile.write("%d\t%d\t%d\t%d\t%d\t" %(1,indiv,parents[indiv][0],parents[indiv][1],sex))
    
        #write geno
        for i in range(1,numSnp):
            outfile.write("%d/%d\t" %(0,0))
        outfile.write("\n")
        
    
    # typed individuals
    indivIndex = 0
    for indiv in range(nUnrelated+1,nUnrelated+nSampled+1):
        
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



def getLkhd(inPath):                    
                
    #open file
    infile = open(inPath)
    
    lkhd = 0
    
    for line in infile:
        
        fields = line.split()
        
        if(len(fields)==0): continue
        
        if(fields[0]=='lnLikelihood'):
            lkhd += float(fields[-1])


    return lkhd

if __name__ == "__main__":

    #file names
    infoPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/unrelated/msprime.info.pruned."
    genoPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/simulations/sim.test.geno."
    outPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/merlin/test2"

        
    #makeMapFile(infoPath, outPath+str(".map"))
    #makeFreqFile(infoPath, outPath+str(".freq"))
    #numSnp = makeDatFile(infoPath, outPath+str(".dat"))
    
    #parents = {1:[0,0], 2:[0,0], 3:[0,0], 4:[1,2], 5:[0,0], 6:[1,2], 7:[0,0], 8:[1,2], 9:[0,0], 10:[3,4], 11:[0,0], 12:[5,6], 13:[0,0], 14:[7,8], 15:[9,10], 16:[11,12], 17:[13,14]}
    parents = {1:[0,0], 2:[0,0], 3:[0,0], 4:[0,0], 5:[0,0], 6:[1,2], 7:[0,0], 8:[1,3], 9:[0,0], 10:[1,4], 11:[0,0], 12:[9,10], 13:[5,6], 14:[7,8], 15:[11,12]}
    #parents = {1:[0,0], 2:[0,0], 3:[0,0], 4:[1,2], 5:[0,0], 6:[1,2], 7:[0,0], 8:[3,4], 9:[0,0], 10:[5,6], 11:[7,8], 12:[9,10]}
    #makePedFile(genoPath, outPath+str(".ped"), numSnp, 12, 3)


    inPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/merlin/ped2.lkhd"
    print(getLkhd(inPath))
