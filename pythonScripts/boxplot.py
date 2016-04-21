import numpy as np
import matplotlib.pyplot as plt
import pdb
import os.path




def getMeanError(inPath):

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
    data = [accDict[key] for key in accDict.keys()]
    means = [np.mean(x) for x in data]

    return data, means


if __name__ == "__main__":

    #file names
    testName = "test2"
    mcmcPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/simulations/results/" + testName + ".kinship.acc"
    pairwisePath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/simulations/results/" + testName + ".kinship.pairwise.acc"
    nIndiv = 6
    nPairs = nIndiv*(nIndiv-1)/2

    
    #get means
    xdata = [i for i in range(1,nPairs+1)]
    mcmcData, mcmcMeans = getMeanError(mcmcPath)
    pairData, pairwiseMeans = getMeanError(pairwisePath)

    
    #plot
    plt.figure()
    plt.boxplot(mcmcData)
    plt.scatter(xdata, mcmcMeans, color='blue', label='mcmc mean error')
    plt.scatter(xdata, pairwiseMeans, color='red', label='pariwise mean error')
    plt.legend()
    #plt.ylim([-.1,1.1])
    plt.xlabel("pair")
    plt.ylabel("error rate")
    plt.title("error rate for 100 simulations; " + testName)
    plt.show()    
    