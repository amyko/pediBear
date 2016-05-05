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


if __name__ == "__main__":

    #file names
    testName = "test8"
    outPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/simulations/results/" + testName + ".true"
    
    #data
    fam1 = [0,1,2,3,4]
    fam2 = [5,6,7,8,9]

    #write
    writeTruePath(outPath, fam1, fam2)