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

        if(trueDict[(i,j)][2]==1): 
            key = 4

        else:
            key = trueDict[(i,j)][0] + trueDict[(i,j)][1]

    
    
        if key in accDict:
            accDict[key].append(acc)
        else:
            accDict[key] = [acc]
    
    print accDict.keys()
    infile.close()


    #make box plot
    data = []
    keys = np.sort(accDict.keys())
    #temp = keys[0]
    #keys = keys[1:]
    
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
    resultDir = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/simulations/results/" 
    testName = "test12"
    mcmcPath = resultDir + testName + ".sa.map.acc.prior.10k"
    pairwisePath = resultDir + testName + ".pruned.10k.pairwise.map"
    testName = "test12"
    plinkPath = resultDir + testName + ".primus.plink.map.acc"
    relatePath = resultDir + testName + ".primus.relate.map.acc"
    truePath = resultDir + testName + ".true"
    nIndiv = 20
    nPairs = nIndiv*(nIndiv-1)/2

    
    
    mcmcMeans = getMeanErrorSortByMeisosis(mcmcPath, truePath, nIndiv)
    pairwiseMeans = getMeanErrorSortByMeisosis(pairwisePath, truePath, nIndiv)
    #plinkMeans = getMeanErrorSortByMeisosis(plinkPath, truePath, nIndiv)
    #relateMeans = getMeanErrorSortByMeisosis(relatePath, truePath, nIndiv)
    xdata = [i for i in range(0,len(mcmcMeans))]

    #tickMarks = [1.0/(4*2**i) for i in range(0,7)]
    tickMarks = ['0', '1/4', '1/8', '1/16', '1/32', '1/64', '1/128', '1/256']

    #pdb.set_trace()
    print mcmcMeans
   # print plinkMeans
    #print relateMeans
    print pairwiseMeans

    
    #plot
    fig = plt.figure(facecolor='white')
    ax = fig.add_subplot(111)
    #fig, ax = plt.figure(facecolor='white')
   #plt.boxplot(mcmcData)
    plt.scatter(xdata, mcmcMeans, color='blue', label='simulated annealing')
    plt.scatter(xdata, pairwiseMeans, color='red', label='pairwise', marker='^')
    plt.show()
    pdb.set_trace()
    
    plt.scatter(xdata, plinkMeans, color='magenta', label='PLINK + PRIMUS', marker='>')
    plt.scatter(xdata, relateMeans, color='green', label='RELATE + PRIMUS', marker='p')
    plt.legend(loc='upper left')
    #plt.ylim([-.1,1.1])
    plt.xlim([-1,len(xdata)+1])
    plt.xlabel("true kinship coefficient")
    plt.ylabel("error rate")
    ax.set_xticks(xdata)
    ax.set_xticklabels(tickMarks)
    plt.setp(tickMarks)

    #plt.title("error rate for 100 simulations; " + testName)
    plt.show()    