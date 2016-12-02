import numpy as np
import matplotlib.pyplot as plt
import pdb
import os.path

def getLkhds(inPath, plotInterval):

    infile = open(inPath)
    
    n=0
    lkhds = []
    xdata = []
    
    for line in infile:
        
        fields = line.split()
        
        if(fields[0]=='>'): continue
        
        #if(int(fields[0]) < 1e7): continue
            
            #increment
        n+=1
            
        if(n%plotInterval==0):
            
            lkhd = float(fields[1])
                
            lkhds.append(lkhd)
            xdata.append(int(fields[0]))
            



    return xdata,lkhds


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
    testName = "100tasiilaq.admixed0.05.pruned0.05.n2.sampleDepth2."
    samplePath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/simulations/results/mcmc.sample.dist."
    convPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/inuits/"+testName
    truePath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/simulations/results/" + testName + ".true"
    nIndiv = 20
    nPairs = nIndiv*(nIndiv-1)/2

    
    # get likelihood
    xdata, lkhd1 = getLkhds(convPath+str(0)+".lkhd", 2)
    xdata, lkhd2 = getLkhds(convPath+str(1)+".lkhd", 2)
    xdata, lkhd3 = getLkhds(convPath+str(2)+".lkhd", 2)
    xdata, lkhd4 = getLkhds(convPath+str(3)+".lkhd", 2)
    xdata, lkhd5 = getLkhds(convPath+str(4)+".lkhd", 2)
    #pdb.set_trace()

    
    #pdb.set_trace()
    
    
    #plot
    plt.figure(facecolor='white')
    plt.plot(xdata, lkhd1, xdata, lkhd2, xdata, lkhd3, xdata, lkhd4, xdata, lkhd5)
    plt.plot(xdata, lkhd1)
    #plt.plot(xdata, conv, xdata, expected)
    plt.ylim([min(lkhd1)-100, max(lkhd1)+100])
    #plt.xlim([-1,len(xdata)+1])
    plt.xlabel("Iteration")
    plt.ylabel("Composite likelihood")
    #plt.title("convergence")
    plt.show()   
