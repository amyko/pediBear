import numpy as np
import matplotlib.pyplot as plt
import pdb
import os.path


class myPath(object):
    
    def __init__(self, up, down, visit):
        self.up = up
        self.down = down
        self.visit = visit
        
    def __hash__(self):
        return hash((self.up, self.down, self.visit))

    def __eq__(self, other):
        return (self.up, self.down, self.visit) == (other.up, other.down,other.visit)

    def __ne__(self, other):
        # Not strictly necessary, but to avoid having both x==y and x!=y
        # True at the same time
        return not(self == other)


def getResult(nIndiv, truePath, saPath, otherPath, otherPath2, pathToOmega):
    
    nPairs = nIndiv*(nIndiv-1)/2
    
    #get true path
    trueDict = getTruePath(truePath)

    saErrors = getMeanErrorSortByMeisosis(saPath, trueDict, pathToOmega, nIndiv)
    otherErrors = getMeanErrorSortByMeisosis(otherPath, trueDict, pathToOmega, nIndiv)
    otherErrors2 = getMeanErrorSortByMeisosis(otherPath2, trueDict, pathToOmega, nIndiv)
    
    return saErrors, otherErrors, otherErrors2


def getMeanErrorSortByMeisosis(inPath, trueDict, pathToOmega, nIndiv):

    #dictionary to store accuracy rates
    accDict = dict()

    #read data
    infile = open(inPath)
    mylist = []


    
    for line in infile:

        fields = line.split()
        if(fields[0]=='>'): 
            iter = float(fields[1])
            continue
    
        if(iter in mylist): 
            continue
    
        i = int(fields[0])
        j = int(fields[1])

        
        acc = 1 - float(fields[2])
        
        #acc = float(fields[2])


        k1k2 = pathToOmega[trueDict[(i,j)]]
        key = .25*k1k2[0] + .5*k1k2[1]
        
 
    
        if key in accDict:
            accDict[key].append(acc)
        else:
            accDict[key] = [acc]
            
 
            
    

    infile.close()


    keys = -np.sort(-np.array(accDict.keys()))
    keys = keys[0:-1]
    print keys

    
    data = []
    data.append(accDict[0.0])
    
    
    
    for k in keys:
        data.append(accDict[k])

    means = [np.mean(x) for x in data]
    

    return means


def getTruePath(truePath):
    
    infile = open(truePath)
    
    trueDict = dict()
    
    for line in infile:

        fields = line.split()

        i = int(fields[0])
        j = int(fields[1])
        up = int(fields[2])
        down = int(fields[3])
        numVisit = int(fields[4])
    

        trueDict[(i,j)] = myPath(up, down, numVisit)
    
    
    infile.close()
    
    return trueDict


def path2k1k2(pathToOmega):

    path2k1k2 = dict()
    
    infile = open(pathToOmega)
    
    for line in infile:
        
        fields = line.split()

        
        key = myPath(int(fields[0]), int(fields[1]), int(fields[2]))
        key2 = myPath(int(fields[1]), int(fields[0]), int(fields[2]))
        k1 = float(fields[5].split('/')[0]) / int(fields[5].split('/')[1])
        k2 = 0
        if(key.up==1 and key.down==1 and key.visit==2):
            k2 = .25
        
        
        path2k1k2[key] = [k1, k2]
        path2k1k2[key2] = [k1,k2]
        
        
    infile.close()
    
    return path2k1k2


if __name__ == "__main__":

    #path2omega
    omegaPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/simulations/pathToOmega.txt" 
    pathToOmega = path2k1k2(omegaPath)
    
    #dir
    resultDir = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/simulations/results/"

    
    #test A
    nIndiv = 20
    saPath = resultDir + "sim6.sa.revised.try2.mapAcc"
    #otherPath2 = resultDir + "test2.mapAcc"
    otherPath = resultDir + "sim6.padre.plink.thresh1.mapAcc"
    otherPath2 = resultDir + "sim6.padre.relate.thresh1.mapAcc"
    truePath = resultDir + "sim6.true"
    saA, otherA, other2A = getResult(nIndiv, truePath, saPath, otherPath, otherPath2, pathToOmega)

    #pdb.set_trace()

    #PLOT
    xdata = [i for i in range(1,len(saA)+1)]

    tickMarks = ['0', '1/4', '1/8', '1/16', '1/32', '1/64', '1/128', '1/256']
    #tickMarks = ['0', '1/4', '1/8', '1/16', '1/32', '1/64']

    print saA
    print otherA
    print other2A
    
    #plot
    fig = plt.figure(facecolor='white')
    
    plt.scatter(xdata, saA, color='blue', label='CLAPPER')
    plt.scatter(xdata, otherA, color='black', marker='>', label='PPP')
    plt.scatter(xdata, other2A, color='red', marker='*', label='RPP')
   
    plt.xticks(xdata,tickMarks)
    plt.legend(loc='upper left')
    plt.xlim([0,len(xdata)+1])
    plt.ylim([0,1])
    plt.xlabel("True kinship coefficient ($\phi$)")
    plt.ylabel("Average error rate ($ \\bar e $)")
    #plt.ylabel("Kinship coefficient distance ($ d $)")


    
    plt.show()   