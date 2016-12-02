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


    print (inPath, numCorrect)
        
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


def getResult(dir, testName):
    
    lengths = [10,20,30,40]
    r_sqrs = ['.025', '.050', '.075', '.100']
    
    data = np.zeros([4,4])
    
    for i in range(0,4):
        len = lengths[i]
        for j in range(0,4):
            
            r_sqr = r_sqrs[j]
            
            inPath = dir+ '%s.%dmorgan.0%s.pairwise.mapAcc' %(testName,len,r_sqr)
            
            data[i,j] = meanError(inPath)
            
            
    return data


if __name__=='__main__':
    
    
    dir = '/Users/kokocakes/Google Drive/Research/pediBear/data/simulations/results/'

    thirdCousin = getResult(dir, 'cousins4')
    secondCousin = getResult(dir, 'cousins3')
    unrelated = getResult(dir, 'unrel')
    
    
    secondCousin = np.array([[ 0.3 ,  0.38,  0.36,  0.48],
       [ 0.34,  0.46,  0.52,  0.4 ],
       [ 0.36,  0.54,  0.4 ,  0.32],
       [ 0.52,  0.48,  0.46,  0.38]])
    
    thirdCousin = np.array([[ 0.12,  0.2 ,  0.4 ,  0.18],
       [ 0.1 ,  0.18,  0.2 ,  0.14],
       [ 0.24,  0.3 ,  0.32,  0.16],
       [ 0.28,  0.32,  0.32,  0.24]])

       
    #normalize data
    #data_norm = (data-data.min()) / (data.max() - data.min())

    #pdb.set_trace()

    #plot
    fig = plt.figure(facecolor='white', figsize=[18,5])
    ax = fig.add_subplot(111)
    ax3 = fig.add_subplot(133)
    ax1 = fig.add_subplot(131, sharex=ax3, sharey=ax3)
    ax2 = fig.add_subplot(132, sharex=ax3, sharey=ax3)
    
    # Turn off axis lines and ticks of the big subplot
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

    heatmap = ax1.pcolor(unrelated, cmap=plt.cm.Blues, alpha = .8, vmin=0, vmax=1)
    heatmap = ax2.pcolor(thirdCousin, cmap=plt.cm.Blues, alpha = .8, vmin=0, vmax=1)
    heatmap = ax3.pcolor(secondCousin, cmap=plt.cm.Blues, alpha = .8, vmin=0, vmax=1)
    cbar_ax = fig.add_axes([0.92, 0.15, 0.05, 0.7])
    fig.colorbar(heatmap, cax=cbar_ax)
    
    #ax.set_frame_on(False)
    ax1.set_yticks(np.arange(unrelated.shape[0]) + 0.5, minor=False)
    ax1.set_xticks(np.arange(unrelated.shape[1]) + 0.5, minor=False)
    xlabels = ['10','20','30','40']
    ylabels = ['.025', '.050', '.075', '.100']
    ax1.set_xticklabels(ylabels, minor=False)
    ax1.set_yticklabels(xlabels, minor=False)
    ax1.grid(False)
    ax.set_xlabel('LD threshold ($r^2$)')
    ax.set_ylabel('Recombination length (Morgan)')
    ax1.set_title("Unrelated")
    ax2.set_title("Third cousins")
    ax3.set_title("Second cousins")
    
    
    ax = plt.gca()
    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    
    plt.show()
    #fig.savefig('prune.png', transparent=True)
    