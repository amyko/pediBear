import numpy as np
import matplotlib.pyplot as plt
import pdb
import os.path

def getMeanErrorSortByMeisosis(inPath, truePath, omegaPath, nIndiv):
    
    #get omega
    omegaDict = getPathToOmega(omegaPath)

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
        
        if(i>=nIndiv or j>=nIndiv): continue
        
        acc = 1 - float(fields[2])

        myPath = (trueDict[(i,j)][0], trueDict[(i,j)][1], trueDict[(i,j)][2])
        key = omegaDict[myPath]
        
        if(key==.03125):
            print (i,j)
    
        if key in accDict:
            accDict[key].append(acc)
        else:
            accDict[key] = [acc]
    

    infile.close()
    
    
    keys = [-x for x in accDict.keys()]
    
    
    
    data = []
    keys = -np.sort(keys)

    for k in keys:
        data.append(accDict[k])


    means = [np.mean(x) for x in data]

    return keys, means



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


def getPathToOmega(myPath):
    
    infile = open(myPath)
    
    myDict = dict()
    
    for line in infile:

        fields = line.split()

        up = int(fields[0])
        down = int(fields[1])
        numVisit = int(fields[2])
 

        ibd2 = float(fields[6].split('/')[0]) / float(fields[6].split('/')[1])
        ibd1 = float(fields[5].split('/')[0]) / float(fields[5].split('/')[1])

        myDict[(up,down,numVisit)] = ibd2 + .5*ibd1
        myDict[(down,up,numVisit)] = ibd2 + .5*ibd1
    
    
    infile.close()
    
    return myDict
    


if __name__ == "__main__":

    #file names
    resultDir = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/jays/" 
    testName = "75jays"
    mcmcPath = resultDir + testName + ".sa.map.acc.50indiv.0"
    pairwisePath = resultDir + testName + ".pruned.50_1.pairwise.map.acc"
    truePath = resultDir + testName + ".true"
    omegaPath = resultDir + "pathToOmega.txt"
    nIndiv = 30
    nPairs = nIndiv*(nIndiv-1)/2

    
    xdata, pairwiseMeans = getMeanErrorSortByMeisosis(pairwisePath, truePath, omegaPath, nIndiv)
    print xdata
    print pairwiseMeans
    xdata, mcmcMeans = getMeanErrorSortByMeisosis(mcmcPath, truePath, omegaPath, nIndiv)    
    print xdata
    print mcmcMeans


    tickMarks = ['%.4f' %i for i in xdata]
    xdata = range(1,len(xdata)+1)


    
    #plot
    fig = plt.figure(facecolor='white')
    ax = fig.add_subplot(111)
    #fig, ax = plt.figure(facecolor='white')
   #plt.boxplot(mcmcData)
    plt.scatter(xdata, mcmcMeans, color='blue', label='simulated annealing')
    plt.scatter(xdata, pairwiseMeans, color='red', label='pairwise', marker='^')

    plt.legend(loc='upper left')
    #plt.ylim([-.1,1.1])
    plt.xlim([-1,len(xdata)+1])
    plt.xlabel("true kinship coefficient x 2")
    plt.ylabel("error rate")
    ax.set_xticks(xdata)
    ax.set_xticklabels(tickMarks)
    plt.setp(tickMarks)

    #plt.title("error rate for 100 simulations; " + testName)
    plt.show() 