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
            dist = pos*1.3/1e6
            
            #write
            outfile.write('%d\tm%d\t%f\t%d\n' %(chrom, marker, dist, pos))
            marker+=1
            
            
        infile.close()
    
    outfile.close()
    
    


def makePedFile(inPath, outPath, numIndiv):
    
    #open file
    outfile = open(outPath, "wb")
    

    
    # typed individuals  
    for indiv in range(0, numIndiv):
        
        
        #write info
        outfile.write("%d\t%d\t%d\t%d\t%d\t%d\t" %(21+indiv,21+indiv,0,0,1,-9))
        
        #write geno
        for chrom in range(1,23):
            
            geno = np.loadtxt(inPath+str(chrom), dtype=str, skiprows=1, delimiter='\t', usecols=[indiv+3])

            #pdb.set_trace()

            for g in geno:
                    
                outfile.write("%s\t%s\t" %(g[0],g[1]))
                
            
        outfile.write("\n")
            
    outfile.close()





if __name__ == "__main__":

    #file names
    unrelPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/unrelated/msprime.geno.pruned."
    genoPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/simulations/sim.test.geno.error."
    outPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/simulations/test12"

        
    #makeMapFile(infoPath, outPath+str(".map"))
    #makeFreqFile(infoPath, outPath+str(".frq"))
    
    makePedFile(unrelPath, outPath+".all.ped", 100)
    