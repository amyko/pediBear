import numpy as np
import pdb
import os

def writeTruePath(outPath, fam1, fam2):

    #open outfile
    outfile = open(outPath,'wb')

    #write path
    for i in range(0,len(fam1)):
        for j in range(i+1, len(fam1)):
            outfile.write('%d\t%d\t%d\t%d\t%d\n' %(fam1[i],fam1[j],3,3,2))
            
    for i in range(0,len(fam2)):
        for j in range(i+1, len(fam2)):
            outfile.write('%d\t%d\t%d\t%d\t%d\n' %(fam2[i],fam2[j],3,3,2))


    for i in range(0,len(fam1)):
        for j in range(0, len(fam2)):
            outfile.write('%d\t%d\t%d\t%d\t%d\n' %(fam1[i],fam2[j],4,4,2))



def concatUnrelated(truePath, outPath):

    outfile = open(outPath,'wb')
    infile = open(truePath)
    
    #copy family relationships
    for line in infile:
        outfile.write(line)
        
    for i in range(0,10):
        for j in range(10,20):
            outfile.write("%d\t%d\t%d\t%d\t%d\n" %(i,j,0,0,0))
            
    for i in range(10,20):
        for j in range(i+1,20):
            outfile.write("%d\t%d\t%d\t%d\t%d\n" %(i,j,0,0,0))
            
    infile.close()
    outfile.close()



if __name__ == "__main__":

    #file names
    testName = "test11"
    truePath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/simulations/results/" + testName + ".true"
    outPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/simulations/results/" + testName + ".true2"
    
    concatUnrelated(truePath, outPath)
    
    #data
    #fam1 = [0,1,2,3,4]
    #fam2 = [5,6,7,8,9]

    #write
    #writeTruePath(outPath, fam1, fam2)