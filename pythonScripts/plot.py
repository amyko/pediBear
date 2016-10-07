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



def getMeanError(inPath, nIndiv):

    #dictionary to store accuracy rates
    accDict = dict()

    #read data
    infile = open(inPath)
    
    for line in infile:
        
        #pdb.set_trace()
    
        fields = line.split()
        if(fields[0]=='>'): continue
    
        i = int(fields[0])
        j = int(fields[1])
        acc = 1 - float(fields[2])
    
        if (i,j) in accDict:
            accDict[(i,j)].append(acc)
        else:
            accDict[(i,j)] = [acc]
    
    
    infile.close()
    
    #make box plot
    data = []
    keys = []
    for i in range(0,nIndiv):
        for j in range(i+1,nIndiv):
            key = (i,j)
            data.append(accDict[key])
            keys.append(key)

    means = [np.mean(x) for x in data]

    return data, means, keys



def getMeanErrorSortByMeisosis(inPath, trueDict, pathToOmega, nIndiv):

    #dictionary to store accuracy rates
    accDict = dict()

    #read data
    infile = open(inPath)
    mylist = [34,89]
    #mylist = []
    
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
        acc = np.sqrt(float(fields[2]))


        k1k2 = pathToOmega[trueDict[(i,j)]]
        key = .25*k1k2[0] + .5*k1k2[1]

    
        if key in accDict:
            accDict[key].append(acc)
        else:
            accDict[key] = [acc]
            
 
            
    

    infile.close()


    keys = -np.sort(-np.array(accDict.keys()))
    keys = keys[0:-1]
    
    
    data = []
    data.append(accDict[0.0])

    
    for k in keys:
        data.append(accDict[k])


    means = [np.mean(x) for x in data]
   # pdb.set_trace()

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

    nIndiv = 20
    nPairs = nIndiv*(nIndiv-1)/2


    #file names
    resultDir = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/simulations/results/" 
    testName = "test12"
    truePath = resultDir + testName + ".true"
    omegaPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/simulations/pathToOmega.txt" 
    
    #accuracy
    priorPath = resultDir + testName + ".pruned.2k.prior.4gen.mapAcc"
    noPriorPath = resultDir + testName + ".pruned.2k.noPrior.4gen.mapAcc"
    pairwisePath = resultDir + testName + ".pruned.2k.prior2.4gen.mapAcc"
    #plinkPath = resultDir + testName + ".pruned.10k.primus.mapAcc"

    #kinship distance
    priorPath = resultDir + testName + ".pruned.2k.prior.4gen.kinshipDist"
    noPriorPath = resultDir + testName + ".pruned.2k.noPrior.4gen.kinshipDist"
    pairwisePath = resultDir + testName + ".pruned.2k.prior2.4gen.kinshipDist"
    
    #get true path and path2Omega
    trueDict = getTruePath(truePath)
    pathToOmega = path2k1k2(omegaPath)
    
    
    priorMeans = getMeanErrorSortByMeisosis(priorPath, trueDict, pathToOmega, nIndiv)
    noPriorMeans = getMeanErrorSortByMeisosis(noPriorPath, trueDict, pathToOmega, nIndiv)
    pairwiseMeans = getMeanErrorSortByMeisosis(pairwisePath, trueDict, pathToOmega, nIndiv)
    xdata = [i for i in range(0,len(priorMeans))]

    print(priorMeans)
    print(pairwiseMeans)

    tickMarks = ['0', '1/4', '1/8', '1/16', '1/32', '1/64', '1/128', '1/256']

    
    #plot
    fig = plt.figure(facecolor='white')
    ax = fig.add_subplot(111)

    #noPriorMeans[5] = .46
    plt.scatter(xdata, priorMeans, color='blue', label='Poi(n)')
    plt.scatter(xdata, pairwiseMeans, color='red', marker='^', label='Poi(n/2)') 
    plt.scatter(xdata, noPriorMeans, color='green', marker='>', label='no prior')
    
    plt.legend(loc='upper left')
    #plt.ylim([-.1,1.1])
    plt.xlim([-1,len(xdata)+1])
    plt.xlabel("True kinship coefficient")
    plt.ylabel("Distance")
    ax.set_xticks(xdata)
    ax.set_xticklabels(tickMarks)
    plt.setp(tickMarks)

    plt.show()    