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
    testName = "test10"
    inPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/simulations/results/test.conv."

    
    #get data
    x = np.loadtxt(inPath+str(0),usecols=[0], dtype=int)
    y0 = np.loadtxt(inPath+str(0),usecols=[1], dtype=float)
    y1 = np.loadtxt(inPath+str(1),usecols=[1], dtype=float)
    
    #pdb.set_trace()
    
    #plot
    plt.figure()
    plt.plot(x[0::50],y0[0::50])
    plt.plot(x[0::50],y1[0::50])
    #plt.plot(x,y1)
    plt.show()