import numpy as np
import matplotlib.pyplot as plt
import pdb
import os.path

def getLkhds(inPath, sampleRate, plotInterval):

    infile = open(inPath)
    
    n=0
    lkhds = []
    
    for line in infile:
        
        fields = line.split()
        
        if(fields[0]=='>'):
            
            #increment
            n+=1
            
            if(n%plotInterval==0):
            
                lkhd = float(fields[1])
                
                lkhds.append(lkhd)
            



    return lkhds


def getPosteriorProp(inPath, plotInterval):

    infile = open(inPath)
    
    n=0
    ratio = []
    
    for line in infile:
        
        fields = line.split()
        
        if(n%plotInterval==0):
            
            r = float(fields[1])
                
            ratio.append(r)
            
        n+=1



    return ratio


def getDistFromTruth(inPath, plotInterval):

    infile = open(inPath)
    
    n=0
    dist = []
    
    for line in infile:
        
        fields = line.split()
        
        if(n%plotInterval==0):
            
            r = float(fields[0])
                
            dist.append(r)
            
        n+=1



    return dist



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
    samplePath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/simulations/results/mcmc.sample.dist."
    convPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/simulations/results/mcmc.sample.conv"
    truePath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/simulations/results/" + testName + ".true"
    nIndiv = 20
    nPairs = nIndiv*(nIndiv-1)/2

    
    # get likelihood
    xdata = [i for i in range(1,20000000, 100000)]
    lkhd1 = getDistFromTruth(samplePath+str(0), 10)
    lkhd2 = getDistFromTruth(samplePath+str(1), 10)
    lkhd3 = getDistFromTruth(samplePath+str(2), 10)
    lkhd4 = getDistFromTruth(samplePath+str(3), 10)
    
    #conv = getPosteriorProp(convPath, 100)
    #expected = [4.98 for i in range(len(xdata))]
    #lkhd1 = getLkhds(samplePath+str(0), 50, 10)
    #lkhd2 = getLkhds(samplePath+str(1), 50, 10)
    #lkhd3 = getLkhds(samplePath+str(2), 50, 10)
    #lkhd4 = getLkhds(samplePath+str(3), 50, 10)
    
    #pdb.set_trace()
    
    
    #plot
    plt.figure()
    plt.plot(xdata, lkhd1, xdata, lkhd2, xdata, lkhd3, xdata, lkhd4)
    #plt.plot(xdata, conv, xdata, expected)
    plt.ylim([-.0001,.001])
    #plt.xlim([-1,len(xdata)+1])
    plt.xlabel("# iterations")
    plt.ylabel("distance from truth")
    plt.title("convergence")
    plt.show()   
