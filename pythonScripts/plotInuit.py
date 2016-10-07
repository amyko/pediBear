import numpy as np
import pdb
import os
import matplotlib.pyplot as plt

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


def plotConv():
    
    dir = '/Users/kokocakes/Google Drive/Research/pediBear/data/inuits/mcmc.2k.newPriorN2.sampleDepth2'
    
    plt.figure()
    
    for i in range(0,5):
            
        myPath = dir + "."+str(i)+".lkhd"
        
        y = np.loadtxt(myPath, dtype=float, usecols=[1], skiprows=1)
        x = range(0,len(y))
        
        plt.plot(x,y)
        
    plt.show()
    
    pdb.set_trace()
    
    

if __name__=='__main__':
    
    
    plotConv()
    
    #get path to kinship
    pathToOmega = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/simulations/pathToOmega.txt" 
    myDict = path2k1k2(pathToOmega)
    
    #plot k1 vs k2
    k1 = []
    k2 = []
    
    infile = open('/Users/kokocakes/Google Drive/Research/pediBear/data/inuits/mcmc.sample.0')
    infile.readline()
    for line in infile:
        print line
        fields = line.split()
        if(fields[0]=='>'): break
        key = key = myPath(int(fields[2]), int(fields[3]), int(fields[4]))
        print (key.up, key.down, key.visit)
        k1k2 = myDict[key]
        
        k1.append(k1k2[0])
        k2.append(k1k2[1])
        
    
    
    infile.close()
    
    
    #plot
    plt.hist(k1, 100)
    plt.show()
    pdb.set_trace()
    plt.scatter(k1,k2)
    plt.xlabel('k1')
    plt.ylabel('k2')
    plt.show()
    
    
        
    #pdb.set_trace()
        
    
    
    