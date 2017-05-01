import numpy as np
import matplotlib.pyplot as plt
import pdb
import os.path


#make geno file from tped file
def makeGenoFile(tpedPath, refPath, outPath, outPath2, t):


    os.system("/Users/kokocakes/Google\ Drive/Research/pediBear/plink_mac/plink --tfile %s --recode --out %stemp" %(tpedPath, outPath))


    
    #open files
    infile = open(outPath2+"temp.ped")
    outfile = open(outPath2+"geno", 'wb')
    
    
    #write sampled genotypes
    for line in infile:
        fields = line.split()
 
        toWrite = ""
        
        for i in range(3, len(fields)/2):
            
            g1 = fields[2*i]
            g2 = fields[2*i+1]
            
            if((g1=='A' and g2=='T') or (g1=='T' and g2=='A')):
                toWrite += "2 "
            elif(g1=='A' and g2=='A'):
                toWrite += "1 "
            elif(g1=='T' and g2=='T'):
                toWrite += "3 "
                
        
        outfile.write(toWrite+"\n")
            
    infile.close()
   
    
    
    #concatenate 100 individuals for allele freq estimation

    #open files
    infile = open(altPath)

    #write sampled genotypes
    for line in infile:
        fields = line.split()
        
        toWrite = ""
        
        for i in range(3, len(fields)/2):
            g1 = fields[2*i]
            g2 = fields[2*i+1]
            
            if((g1=='A' and g2=='T') or (g1=='T' and g2=='A')):
                toWrite += "2 "
            elif(g1=='A' and g2=='A'):
                toWrite += "1 "
            elif(g1=='T' and g2=='T'):
                toWrite += "3 "
                
        
        outfile.write(toWrite+"\n")
            
    infile.close()
    outfile.close()



def makePosFile(inPath, outPath):


    outfile = open(outPath,'wb')
    chr = np.loadtxt(inPath,usecols=[2],dtype=float)
    
    for x in chr:
        outfile.write(str(x*100)+"\n")
        
    outfile.close()

    
def makeChrFile(inPath, outPath):


    outfile = open(outPath,'wb')
    chr = np.loadtxt(inPath,usecols=[0],dtype=int)
    
    for x in chr:
        outfile.write(str(x)+"\n")
        
    outfile.close()


def makeIndivFile(outPath, numIndiv, sampleNum):

    outfile = open(outPath,'wb')
    
    for i in range(0,numIndiv):
        if(i<sampleNum): outfile.write("%d\t" %1)
        
        else: outfile.write("%d\t" %0)
    
    outfile.write('\n')
    outfile.close()
    
    
def makeIBDFile(inPath, outPath, numIndiv):
    
    outfile = open(outPath)
    infile = open(inPath)
    
    oufile.write("FID1\tIID1\tFID2\tIID2\tIBD0\tIBD1\tIBD2\tPI_HAT\n")
    
    for i in range(1, numIndiv+1):
        for j in range(i+1, numIndiv+1):
            
            k1 = float(infile.readline().split()[0])
            k2 = float(infile.readline().split()[0])
            k0 = 1 - ki - k2
            pi_hat = .5*k1 + k2
            
            outfile.write("%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\n" %(i,i,j,j,k0, k1, k2, pi_hat))
            
            
    outfile.close()
    infile.close()
    


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
        
    testName = "sim4"
    
    myDir = os.path.expanduser('~') + "/Google\ Drive/Research/pediBear/data/simulations/simPed4/"
    myDir2 = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/simulations/simPed4/"
    hmmPath = os.path.expanduser('~') + "/Google\ Drive/Research/tools/relate-master/src/relateHMM"

    altPath = "/Users/kokocakes/Google Drive/Research/pediBear/data/unrelated/msprime.unrel.try2.pruned.fixed.ped"
    outPath = myDir + "relateFiles/"
    outPath2 = myDir2 + "relateFiles/"
    
    #chr
    #makePosFile(myDir2+"sim4.0.tped", outPath2+"pos")
    #makeChrFile(myDir2+"sim4.0.tped", outPath2+"chr")
    makeIndivFile(outPath2+"indiv", 118,18)

    for t in range(0,1):
        print(t)
        tpedPath = myDir + testName+".%d" %t
        
        # make input file
        makeGenoFile(tpedPath, altPath, outPath, outPath2, t)
        
        pdb.set_trace()
        
        # run relate; call from data director
        os.system(hmmPath + "-g geno -p pos -c chr -o options -d indiv -post postout -k kout")
            
        pdb.set_trace()
        # make ibd file
        
        
        
        # remove outfiles
        os.system('rm postout')
        os.sytem('rm kout')
    
    
    '''
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
    '''
