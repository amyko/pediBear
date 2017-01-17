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


def getResult(nIndiv, truePath, saPath, otherPath, pathToOmega):
    
    nPairs = nIndiv*(nIndiv-1)/2
    
    #get true path
    trueDict = getTruePath(truePath)

    saErrors = getMeanErrorSortByMeisosis(saPath, trueDict, pathToOmega, nIndiv)
    otherErrors = getMeanErrorSortByMeisosis(otherPath, trueDict, pathToOmega, nIndiv)
    
    
    return saErrors, otherErrors


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
    testName = "sim2"
    saPath = resultDir + testName + ".n1.mapAcc"
    otherPath = resultDir + testName + ".5.pairwise.mapAcc"
    otherPath = resultDir + "sim2.mcmc.3chains.mapAcc"
    truePath = resultDir + "sim2.true"
    saB, otherB = getResult(nIndiv, truePath, saPath, otherPath, pathToOmega)


    #test B
    nIndiv = 20
    saPath = resultDir + "sim1.n1.mapAcc"
    otherPath = resultDir + "sim1.5.pairwise.mapAcc"
    otherPath = resultDir + "sim1.mcmc.3chains.mapAcc"
    truePath = resultDir + "test12.true"
    saA, otherA = getResult(nIndiv, truePath, saPath, otherPath, pathToOmega)
    
    
    
    #test C
    nIndiv = 18
    saPath = resultDir + "sim4.n1.mapAcc"
    otherPath = resultDir + "sim4.5.pairwise.mapAcc"
    truePath = resultDir + "sim4.true"
    saC, otherC = getResult(nIndiv, truePath, saPath, otherPath, pathToOmega)
    
    
 
    #PLOT
    xdata = [i for i in range(1,len(saA)+1)]

    tickMarks = ['0', '1/4', '1/8', '1/16', '1/32', '1/64', '1/128', '1/256']

    
    #plot
    fig = plt.figure(facecolor='white')
    ax = fig.add_subplot(111)
    ax3 = fig.add_subplot(313)
    ax1 = fig.add_subplot(311, sharex=ax3)
    ax2 = fig.add_subplot(312, sharex=ax3)

    # Turn off axis lines and ticks of the big subplot
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')


    ax1.scatter(xdata, saA, color='blue', label='SA')
    ax1.scatter(xdata, otherA, color='red', marker='>', label='Pairwise')
    ax2.scatter(xdata, saB, color='blue')
    ax2.scatter(xdata, otherB, color='red', marker='>')
    ax3.scatter(xdata, saC, color='blue')
    ax3.scatter(xdata, otherC, color='red', marker='>')
    
    #ax1.boxplot(saA)
    #ax2.boxplot(saB)
    #ax3.boxplot(saC)
    #ax1.set_ylim([-.01,.04])
    #ax2.set_ylim([-.01,.04])
    #ax3.set_ylim([-.01,.14])
    
    
    ax1.text(.98, .95, "A", transform=ax1.transAxes,fontsize=16, fontweight='bold', va='top', ha="right")
    ax2.text(.98, .95, "B", transform=ax2.transAxes,fontsize=16, fontweight='bold', va='top', ha='right')
    ax3.text(.98, .95, "C", transform=ax3.transAxes,fontsize=16, fontweight='bold', va='top', ha='right')
    
    
    
    ax1.legend(loc='upper left')
    plt.xlim([0,len(xdata)+1])
    ax.set_xlabel("True kinship coefficient ($\phi$)", fontsize=12)
    ax.set_ylabel("Average error rate ($ \\bar e $)", fontsize=12)
    #fig.text(0.04, 0.5, 'Kinship coefficient distance ($d$)', va='center', rotation='vertical')
    ax3.set_xticks(xdata)
    ax3.set_xticklabels(tickMarks)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    
    plt.show()   