import numpy as np
import matplotlib.pyplot as plt
import pdb
import os.path

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



def getMeanErrorSortByMeisosis(inPath, truePath, nIndiv):
    
    #get true path
    trueDict = getTruePath(truePath)

    
    #dictionary to store accuracy rates
    accDict = dict()

    #read data
    infile = open(inPath)
    
    for line in infile:

        fields = line.split()
        if(fields[0]=='>'): continue
    
        i = int(fields[0])
        j = int(fields[1])
        acc = 1 - float(fields[2])
    
    
        key = trueDict[(i,j)][0] + trueDict[(i,j)][1]
    
    
        if key in accDict:
            accDict[key].append(acc)
        else:
            accDict[key] = [acc]
    
    
    infile.close()


    #make box plot
    data = []
    keys = np.sort(accDict.keys())
    temp = keys[0]
    keys = keys[1:]
    keys = np.append(keys, [temp])
    
    
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
    

        trueDict[(i,j)] = [up, down, numVisit]
    
    
    infile.close()
    
    return trueDict


if __name__ == "__main__":

    #file names
    testName = "test11"
    mcmcPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/simulations/results/" + testName + ".mcmc.map.acc"
    pairwisePath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/simulations/results/" + testName + ".pairwise.map.acc"
    truePath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/simulations/results/" + testName + ".true"
    nIndiv = 20
    nPairs = nIndiv*(nIndiv-1)/2

    
    #get means
   # xdata = [i for i in range(1,nPairs+1)]
    #mcmcData, mcmcMeans, keys = getMeanError(mcmcPath, nIndiv)
    #pairData, pairwiseMeans, keys = getMeanError(pairwisePath, nIndiv)
    
    mcmcMeans = getMeanErrorSortByMeisosis(mcmcPath, truePath, nIndiv)
    pairwiseMeans = getMeanErrorSortByMeisosis(pairwisePath, truePath, nIndiv)
    xdata = [i for i in range(0,len(mcmcMeans))]
    
    
    #plot
    plt.figure()
   #plt.boxplot(mcmcData)
    plt.scatter(xdata, mcmcMeans, color='blue', label='mcmc mean error')
    plt.scatter(xdata, pairwiseMeans, color='red', label='pairwise mean error')
    plt.legend()
    #plt.ylim([-.1,1.1])
    plt.xlim([-1,len(xdata)+1])
    plt.xlabel("relationship category (in increasing distance)")
    plt.ylabel("error rate")
    plt.title("error rate for 100 simulations; " + testName)
    plt.show()    