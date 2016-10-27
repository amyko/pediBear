import numpy as np
import matplotlib.pyplot as plt
import pdb

def meanError(inPath):
    
    infile = open(inPath)
    
    numCorrect = 0
    nTotal = 0
    
    for line in infile:
        
        fields = line.split()
        
        if(fields[0])=='>': continue
        
        numCorrect += float(fields[2])
        nTotal+=1
        
    return numCorrect/float(nTotal)


def meanDist(inPath):
    
    infile = open(inPath)
    
    numCorrect = 0
    nTotal = 0
    
    for line in infile:
        
        fields = line.split()
        
        if(fields[0])=='>': continue
        
        numCorrect += np.sqrt(float(fields[2]))
        nTotal+=1
        
    return numCorrect/float(nTotal)



if __name__=='__main__':
    
    
    dir = '/Users/kokocakes/Google Drive/Research/pediBear/data/simulations/results/'
    lengths = [10,20,30,40]
    r_sqrs = ['.013', '.025', '.050', '.075', '.100']
    
    data = np.zeros([4,5])
    
    for i in range(0,4):
        len = lengths[i]
        for j in range(0,5):
            
            r_sqr = r_sqrs[j]
            
            inPath = dir+ 'sim3.%dmorgan.0%s.pairwise.kinshipDist' %(len,r_sqr)
            
            data[i,j] = meanDist(inPath)

       
    #normalize data
    data_norm = (data-data.min()) / (data.max() - data.min())

    pdb.set_trace()
    
    #plot
    fig, ax = plt.subplots()
    fig.patch.set_facecolor('white')
    heatmap = ax.pcolor(data, cmap=plt.cm.Blues, alpha = .8)
    plt.colorbar(heatmap)
    
    #ax.set_frame_on(False)
    ax.set_yticks(np.arange(data.shape[0]) + 0.5, minor=False)
    ax.set_xticks(np.arange(data.shape[1]) + 0.5, minor=False)
    xlabels = ['10','20','30','40']
    ylabels = ['.013', '.025', '.050', '.075', '.100']
    ax.set_xticklabels(ylabels, minor=False)
    ax.set_yticklabels(xlabels, minor=False)
    ax.grid(False)
    plt.xlabel('LD threshold ($r^2$)')
    plt.ylabel('Recombination length (Morgan)')
    
    
    ax = plt.gca()
    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    
    plt.show()
    #fig.savefig('prune.png', transparent=True)
    