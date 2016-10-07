import numpy as np
import pdb
import os

def makeHap(nSites):
    
    toReturn = np.empty([2, nSites],str)
    
    for i in range(0, nSites):
        
        if(np.random.random() < .5):
            nuc1 = 'A'
        else:
            nuc1 = 'T'
        if(np.random.random() < .5):
            nuc2 = 'A'
        else:
            nuc2 = 'T'
            
        toReturn[0,i] = nuc1
        toReturn[1,i] = nuc2
    
    return toReturn

def makeChild(mom, dad):
    
    nSites = mom.shape[1]
    lam = .0013
    
    toReturn = np.empty([2, nSites], str)

    # choose current haplotype at first snp
    currHapMom = int(np.random.random() < .5)
    currHapDad =  int(np.random.random() < .5)
    toReturn[0,0] = mom[currHapMom, 0]
    toReturn[1,0] = dad[currHapDad, 0]
    
    # remaining sites
    n = 0
    for i in range(1,nSites):
        
        #recombine for mom hap
        if(np.random.poisson(lam) % 2 == 1):
            currHapMom = 1 - currHapMom
            n+=1
        
        if(np.random.poisson(lam) % 2 == 1):
            currHapDad = 1 - currHapDad
            n+=1               
    
        toReturn[0,i] = mom[currHapMom, i]
        toReturn[1,i] = dad[currHapDad, i]
    
    print n

    return toReturn


def writePedLine(outfile, name, pa, ma, sex, hap):
    
    outfile.write("fam1\t%s\t%s\t%s\t%d\t-1\t" %(name, pa, ma, sex))
    
    for i in range(0, hap.shape[1]):
        
        outfile.write("%s %s\t" %(hap[0,i], hap[1,i]))
    
    outfile.write("\n")
    
    
def makeMapFile(fileName, nSites):
    
    outfile = open(fileName,'wb')
    
    for i in range(0,nSites):
        outfile.write("1\t%d\t%f\n" %(i, i*100000*1.3e-6))
        
        
def makeFreqFile(fileName, nSites):
    
    outfile = open(fileName,'wb')
    
    for i in range(0,nSites):
        outfile.write("M\t%d\n" %i)
        outfile.write("F\t.5\t.5\t0\t0\n")
        
    outfile.close()
    
def makeDatFile(fileName, nSites):
    
    outfile = open(fileName,'wb')
    
    for i in range(0,nSites):
        outfile.write("M\t%d\n" %i)
        
    outfile.close()


def makeInfoFile(fileName, nSites):
    
    outfile = open(fileName,'wb')

    #header
    outfile.write("CHROM\tPOS\tLDPOS\tU\tu\tV\tv\tpA\tpC\tpG\tpT\tA\tB\tC\tD\n")

    for i in range(0,nSites):
        
        outfile.write("1\t%d\t-1\tX\tX\tA\tT\t.5\t.5\t.5\t.5\t-1\t-1\t-1\t-1\n" %(i*100000))
        
    outfile.close()
        

if __name__=='__main__':
    
    dir = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/simulations/"
    
    #make files
    makeInfoFile(dir + "nucFam.info", 3000)
    #makeMapFile(dir + "nucFam.merlin.map", 3000)
    #makeDatFile(dir + "nucFam.dat", 3000)
    #makeFreqFile(dir + "nucFam.freq", 3000)
    pdb.set_trace()
    
    
    #outfile
    outfile = open(dir + "nucFam.ped",'wb')
    
    #make parents
    mom = makeHap(3000)
    dad = makeHap(3000)
    writePedLine(outfile, 'dad', '0', '0', 1, dad)
    writePedLine(outfile, 'mom', '0', '0', 2, mom)
    
    
    #make children
    for i in range(1,5):
        child = makeChild(mom, dad)
        writePedLine(outfile, 'c%d' %i, 'dad', 'mom', 1, child)
