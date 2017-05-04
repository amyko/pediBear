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


def getResult(nIndiv, truePath, saPath, pathToOmega):
    
    nPairs = nIndiv*(nIndiv-1)/2
    
    #get true path
    trueDict = getTruePath(truePath)

    saErrors = getMeanErrorSortByMeisosis(saPath, trueDict, pathToOmega, nIndiv)
 
    return saErrors


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

        
        #acc = 1 - float(fields[2])
        
        acc = float(fields[2])


        k1k2 = pathToOmega[trueDict[(i,j)]]
        key = .25*k1k2[0] + .5*k1k2[1]
        
        if(key==0): continue
        
        acc = acc / key
        
  
    
        if key in accDict:
            accDict[key].append(acc)
        else:
            accDict[key] = [acc]
            
 
            
    

    infile.close()


    keys = -np.sort(-np.array(accDict.keys()))
    #keys = keys[0:-1]
    print keys

    
    data = []
    #data.append(accDict[0.0])
    
    
    
    for k in keys:
        data.append(accDict[k])

    #means = [np.mean(x) for x in data]
    

    return data


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
    saPath = resultDir + "sim6.sa.revised.try2.kinshipDist"
    plinkPath = resultDir + "sim6.padre.plink.thresh2.kinshipDist"
    relatePath = resultDir + "sim6.padre.relate.thresh2.kinshipDist"
    truePath = resultDir + "sim6.true"
    sa = getResult(nIndiv, truePath, saPath, pathToOmega)
    plink = getResult(nIndiv, truePath, plinkPath, pathToOmega)
    relate = getResult(nIndiv, truePath, relatePath, pathToOmega)

    #pdb.set_trace()

    #PLOT
    xdata = [i for i in range(1,len(sa)+1)]


    tickMarks = ['1/4', '1/8', '1/16', '1/32', '1/64', '1/128', '1/256']
    #tickMarks = [ '1/4', '1/8', '1/16', '1/32', '1/64']

    
    #plot
    fig = plt.figure(facecolor='white', figsize=[20,5])
    ax = fig.add_subplot(111)
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132, sharey=ax1)
    ax3 = fig.add_subplot(133, sharey=ax1)
    
    # Turn off axis lines and ticks of the big subplot
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

    bp1 = ax1.boxplot(sa, patch_artist=True)
    bp2 = ax2.boxplot(plink, patch_artist=True)
    bp3 = ax3.boxplot(relate, patch_artist=True)
    
    
    #color
    plt.setp(bp1['medians'], color='magenta')
    plt.setp(bp1['whiskers'], color='blue')
    plt.setp(bp1['fliers'], color='blue')
    plt.setp(bp1['boxes'], color='blue', alpha = .5)
    plt.setp(bp2['medians'], color='magenta')
    plt.setp(bp2['whiskers'], color='black')
    plt.setp(bp2['fliers'], color='black')
    plt.setp(bp2['boxes'], color='black', alpha = .5)
    plt.setp(bp3['medians'], color='magenta')
    plt.setp(bp3['whiskers'], color='red')
    plt.setp(bp3['fliers'], color='red')
    plt.setp(bp3['boxes'], color='red', alpha = .5)
    
    ax1.text(.05, .98, "CLAPPER", transform=ax1.transAxes,fontsize=12, fontweight='bold', va='top', ha="left")
    ax2.text(.05, .98, "PPP", transform=ax2.transAxes,fontsize=12, fontweight='bold', va='top', ha='left')
    ax3.text(.05, .98, "RPP", transform=ax3.transAxes,fontsize=12, fontweight='bold', va='top', ha='left')
    
    

    
    ax3.set_xticks(xdata)
    ax1.set_xticklabels(tickMarks)
    ax2.set_xticklabels(tickMarks)
    ax3.set_xticklabels(tickMarks)
    plt.xlim([0,len(xdata)+1])
    plt.ylim([-0.1,8])
    ax.set_xlabel("True kinship coefficient ($\phi$)")
    ax.set_ylabel("Kinship coefficient distance ($ d $)")
    

    plt.show()  