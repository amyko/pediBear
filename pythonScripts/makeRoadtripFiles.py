import numpy as np
import matplotlib.pyplot as plt
import pdb
import os.path
import random



def makeGenoFile(inPath, outPath):

    #open files
    infile = open(inPath)
    outfile = open(outPath,'wb')
    
    #skip header
    infile.readline()
    
    #write pair info
    for line in infile:
        
        fields = line.split()
        
        for g in fields:
            
            if(g=='AA'): outfile.write("%d\t" %2)
            elif(g=='AT' or g=='TA'): outfile.write("%d\t" %1)
            else: outfile.write("%d\t" %0)
        
    
        outfile.write("\n")
    

    
    infile.close()
    outfile.close()


def makePed():

    outfile = open(os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/assoc/noPedigree.ped", 'wb')

    fam = 1
    for i in range(1,101):              
                
        outfile.write("%d\t%d\t%d\t%f\n" %(fam, 1, 1,0))
        fam+=1

#pedigree
def makePed2():
    
    outfile = open(os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/assoc/pedigree.ped", 'wb')

    for fam in range(1,21):              


        for i in range(1,6):
        
            kin = 0
        
            for j in range(i,6):
        
                if(i!=j):
                    kin = 0.015625
        
                outfile.write("%d\t%d\t%d\t%f\n" %(fam, i,j, kin))
                
                
                
def makePheno():
    
    infile = open(os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/assoc/msprime.info.pruned.1")
    
    infile.readline()
    
    lineNum = 0
    for line in infile:
        lineNum+=1
        fields = line.split()
        freq = float(fields[6])
        if (freq < .6 and freq > .4): 
            print freq
            print lineNum
            break
        
        
    infile.close()
        
        
def makePheno2(inPath, outPath, targetLine, r1, r2, e, t):
    
    outfile = open(outPath+'noPedigree.'+str(t), 'wb')
    outfile2 = open(outPath+'pedigree.'+str(t), 'wb')
    infile = open(inPath)
    
    #find geno line
    lineNum = 0
    geno = ''
    for line in infile:
        if(lineNum==targetLine):
            geno = line.split()
            break
        lineNum+=1
        
    infile.close()
    

    #write file
    fam1 = 1
    fam2 = 1
    indiv = 1
    n = 0
    for g in geno:

        environ = np.random.normal(fam2*e, .01)
        #environ = 0

        prob = r1 + environ
        affected = 1
        
        if(g=='2'):
            prob = r2 + environ
            
        if(random.random() < prob):
            affected = 2
            n+=1
        
        #no pedigree
        outfile.write("%d\t%d\t%d\n" %(fam1,1,affected))
        fam1+=1
        
        #pedigree
        outfile2.write("%d\t%d\t%d\n" %(fam2,indiv,affected))
        
        
        #increment
        if(indiv%5==0):
            indiv = 1
            fam2+=1
        else:
            indiv+=1
        
    
    outfile.close()
    outfile2.close()
    print n
    
    
    
def makePheno3(inPath, outPath, targetLine, r1, r2, e, t):
    
    outfile = open(outPath+'noPedigree.'+str(t), 'wb')
    outfile2 = open(outPath+'pedigree.'+str(t), 'wb')
    infile = open(inPath)
    
    #find geno line
    lineNum = 0
    geno = ''
    for line in infile:
        if(lineNum==targetLine):
            geno = line.split()
            break
        lineNum+=1
        
    infile.close()
    

    #write file
    fam1 = 1
    fam2 = 1
    indiv = 1
    n = 0
    
    for g in geno:

        #compute environmental factor
        mean = (np.random.normal(fam2*e,.01) + np.random.normal(fam2*e,.01))/2
        mean = (mean + np.random.normal(mean,.01))/2
        environ = (mean + np.random.normal(mean,.01))/2

        prob = r1 + environ
        affected = 1
        
        if(g=='2'):
            prob = r2 + environ
            
        if(random.random() < prob):
            affected = 2
            n+=1
        
        #no pedigree
        outfile.write("%d\t%d\t%d\n" %(fam1,1,affected))
        fam1+=1
        
        #pedigree
        outfile2.write("%d\t%d\t%d\n" %(fam2,indiv,affected))
        
        
        #increment
        if(indiv%5==0):
            indiv = 1
            fam2+=1
        else:
            indiv+=1
        
    
    outfile.close()
    outfile2.close()
    print n
    
    

if __name__ == "__main__":
        
    
    #file names
    inPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/assoc/geno."
    genoPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/assoc/genofile."
    phenoPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/assoc/pheno."

    #param
    targetLine = 30
    r1 = .4
    r2 = .05
    e = .03

    #makePheno()
    #pdb.set_trace()
    makePed2()
    for t in range(0,20):
        #print(t)
        #makeGenoFile(inPath+str(t), genoPath+str(t))
        makePheno2(genoPath+str(t), phenoPath, targetLine, r1, r2, e, t)

