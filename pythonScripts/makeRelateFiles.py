import numpy as np
import matplotlib.pyplot as plt
import pdb
import os.path



def makeGenoFile(pairPath, genoPath, outPath):

    #open files
    outfile = open(outPath,'wb')
    
    
    #write pair info
    for i in range(0,2):
        genos = np.loadtxt(pairPath, dtype=str, delimiter='\t', skiprows=1, usecols=[i])
        
        for g in genos:
            
            if(g=='AA'): outfile.write("%d\t" %1)
            elif(g=='AT' or g=='TA'): outfile.write("%d\t" %2)
            else: outfile.write("%d\t" %3)
        
    
        outfile.write("\n")
    
    #concatenate 50 individuals for allele freq estimation
    for i in range(50,100):
        
        genos = np.loadtxt(genoPath, dtype=str, delimiter='\t', skiprows=1, usecols=[i])
        #print len(genos)
        
        for g in genos:
            
            if(g=='AA'): outfile.write("%d\t" %1)
            elif(g=='AT' or g=='TA'): outfile.write("%d\t" %2)
            else: outfile.write("%d\t" %3)
        
    
        outfile.write("\n")
    

    outfile.close()



def makeChrFile(outPath, numSnps):


    outfile = open(outPath,'wb')
    
    for i in range(0,numSnps):
        outfile.write("%d\n" %1)


def makeIndivFile(outPath, numIndiv):

    outfile = open(outPath,'wb')
    
    for i in range(0,numIndiv):
        if(i<2): outfile.write("%d\t" %1)
        
        else: outfile.write("%d\t" %0)
    
    outfile.write('\n')
    outfile.close()
    
    
def pedToDat(pedPath, altPath, outPath):
    
    outfile = open(outPath, 'wb')
    
    #selected individuals
    infile = open(pedPath)
    for line in infile:
        geno = line.split()[6:]
        
        for i in range(0, len(geno)/2):
            a1 = geno[2*i]
            a2 = geno[2*i+1]
            
            g = 2
            
            if(a1=='A' and a2=='A'): g = 1
            elif(a1=='T' and a2=='T'): g = 3
            
            outfile.write('%d\t' %g)
        
        outfile.write("\n")

    infile.close()

    # alt reference pop
    infile = open(altPath)
    for line in infile:
        geno = line.split()[6:]
        
        for i in range(0, len(geno)/2):
            a1 = geno[2*i]
            a2 = geno[2*i+1]
            
            g = 2
            if(a1=='A' and a2=='A'): g = 1
            elif(a1=='T' and a2=='T'): g = 3
            
            outfile.write('%d\t' %g)
        
        outfile.write("\n")
            
                
    infile.close()
    outfile.close()      
        
        

    

if __name__ == "__main__":
        
    testName = "test"
    
        
    #file names
    indivPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/ibd/"+testName+".indiv"
    chrPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/ibd/"+testName+".chr"
    pairPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/ibd/"+testName+".geno."
    genoPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/unrelated/msprime.geno.pruned.1"
    outPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/ibd/"+testName+".geno.concat."

    #makeChrFile(chrPath, 3541)
    #makeIndivFile(indivPath, 52)

    pedPath = "/Users/kokocakes/Google Drive/Research/pediBear/data/simulations/genotypes/test12.pruned.10k."
    altPath = "/Users/kokocakes/Google Drive/Research/pediBear/data/simulations/genotypes/test12.all.pruned.10k.ped"
    outPath = "/Users/kokocakes/Google Drive/Research/pediBear/data/simulations/relateFiles/geno."

    for t in range(0,100):
        print(t)
        pedToDat(pedPath+str(t)+".ped", altPath, outPath+str(t))

            
    
        
    pdb.set_trace()
    
    fixedTPR = [ 0.59877667,  0.59409326,  0.58564486,  0.57774069,  0.5718019 , 0.56074547,  0.55058438,  0.53579102,  0.52217112,  0.49701309]   
    fixedFPR = [ 0.02156979,  0.02140035,  0.01144539,  0.0112495 ,  0.01097371, 0.01083524,  0.01069637,  0.01043933,  0.0102061 ,  0.0100275 ]
    estTPR = [ 0.60419083,  0.59980281,  0.59423987,  0.58732693,  0.5786517 ,0.57043254,  0.56017463,  0.5461609 ,  0.53211058,  0.50726091]
    estFPR = [ 0.03188398,  0.03157868,  0.03141643,  0.03122449,  0.02129386, 0.02098951,  0.0207332 ,  0.02055863,  0.02035207,  0.0101476 ]
    

    fig = plt.figure(facecolor='white')
    plt.plot(fixedFPR, fixedTPR, color='blue', label='relationship known')
    plt.plot(estFPR, estTPR, color='red', label='relationship unknown')
    plt.scatter(fixedFPR, fixedTPR, color='blue')
    plt.scatter(estFPR, estTPR, color='red')
    
    plt.title("ROC")
    plt.ylabel("true positive rate")
    plt.xlabel("false positive rate")
    plt.legend(loc='lower right')
    
    
    plt.show()
    
    pdb.set_trace()  
