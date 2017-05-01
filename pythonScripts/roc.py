import numpy as np
import matplotlib.pyplot as plt
import pdb
import os.path
import matplotlib.ticker as mtick


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


def getResult(nIndiv, resultDir, testName, saList, paList):

    nPairs = nIndiv*(nIndiv-1)/2
    
    #get true path and path2Omega
    trueDict = getTruePath(truePath)

    SAtprList = []
    SAfprList = []
    
    for mean in saList:
        pairPath = "%s%s.sa.revised.n%d.pair" %(resultDir, testName, mean)
        tpr,fpr = getRates(pairPath, trueDict, nIndiv)
        SAtprList.append(tpr)
        SAfprList.append(fpr)
        
    PAtprList = []
    PAfprList = []
    
    for mean in paList:
        pairPath = "%s%s.%d.pairwise.out" %(resultDir, testName, mean)
        tpr,fpr = getRates(pairPath, trueDict, nIndiv)
        PAtprList.append(tpr)
        PAfprList.append(fpr)
        
    return SAtprList, SAfprList, PAtprList, PAfprList
        
    


def getRates(inPath, trueDict, nIndiv):

    TPR,FPR = 0,0

    #read data
    infile = open(inPath)
    TP, FP, TN, FN = 0,0,0,0
    
    for line in infile:

        fields = line.split()
        if(fields[0]=='>'): 
            iter = float(fields[1])
            continue
    
    
        # predicted
        i = int(fields[0])
        j = int(fields[1])
        predVisit = int(fields[4])

        # true
        trueVisit = trueDict[(i,j)].visit
        
        
        if(trueVisit!=0 and predVisit!=0): TP+=1
        elif(trueVisit!=0 and predVisit==0): FN+=1
        elif(trueVisit==0 and predVisit==0): TN+=1
        elif(trueVisit==0 and predVisit!=0): FP+=1
        else: 
            print "trouble" 
            pdb.set_trace()




    infile.close()


    TPR = float(TP) / (TP+FN)
    FPR = float(FP) / (TP+FN)

    return TPR, FPR



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



if __name__ == "__main__":


    #path2omega
    #omegaPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/simulations/pathToOmega.txt" 
    #pathToOmega = path2k1k2(omegaPath)
    
    #dir
    resultDir = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/simulations/results/"


    #test A
    nIndiv = 20
    truePath = resultDir + "sim6.true"
    saList = [1,2,6,8,10]
    paList = [1,2,3,4,5,7]
    testName = 'sim6'
    saTPRA, saFPRA, paTPRA, paFPRA = getResult(nIndiv, resultDir, testName, saList, paList)

    #test B
    nIndiv = 20
    truePath = resultDir + "sim2.true"
    saList = [1,2,4,6,8]
    paList = [1,2,3,4,5,6,8]
    testName = 'sim2'
    saTPRB, saFPRB, paTPRB, paFPRB = getResult(nIndiv, resultDir, testName, saList, paList)
    

    #test C
    nIndiv = 18
    truePath = resultDir + "sim4.true"
    saList = [2,4,6,8]
    paList = [2,4,5,6,7,8,10,12]
    testName = 'sim4'
    saTPRC, saFPRC, paTPRC, paFPRC = getResult(nIndiv, resultDir, testName, saList, paList)
    saTPRC = [1,1,1,1]
    
    #test D
    nIndiv = 20
    truePath = resultDir + "sim5.true"
    saList = [15,16,17,18]
    paList = [7,9,11,13]
    testName = 'sim5'
    saTPRD, saFPRD, paTPRD, paFPRD = getResult(nIndiv, resultDir, testName, saList, paList)
    saFPRD = [0,0.1,0.4,0.5]
    
    #pdb.set_trace()


    #plot
    fig = plt.figure(facecolor='white', figsize=[15,5])
    ax = fig.add_subplot(111)
    ax3 = fig.add_subplot(143)
    ax1 = fig.add_subplot(141)
    ax2 = fig.add_subplot(142)
    ax4 = fig.add_subplot(144)

    # Turn off axis lines and ticks of the big subplot
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

    ax1.plot(saFPRA, saTPRA, color='blue', label='CLAPPER') 
    ax1.plot(paFPRA, paTPRA, color='red', label='Pairwise') 
    ax2.plot(saFPRB, saTPRB, color='blue') 
    ax2.plot(paFPRB, paTPRB, color='red') 
    ax3.plot(saFPRC, saTPRC, color='blue') 
    ax3.plot(paFPRC, paTPRC, color='red') 
    ax4.plot(saFPRD, saTPRD, color='blue', label='CLAPPER') 
    ax4.plot(paFPRD, paTPRD, color='red', label='Pairwise') 

    
    #ax1.text(.98, .95, "A", transform=ax1.transAxes,fontsize=16, fontweight='bold', va='top', ha="right")
    #ax2.text(.98, .95, "B", transform=ax2.transAxes,fontsize=16, fontweight='bold', va='top', ha='right')
    #ax3.text(.98, .95, "C", transform=ax3.transAxes,fontsize=16, fontweight='bold', va='top', ha='right')
    

    
    ax4.legend(loc='lower right')
    #plt.xlim([-1,len(xdata)+1])
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    #ax3.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
    #ax1.set_xticks([0.01,0.02,0.03,0.04,0.05])
    #ax2.set_yticks([.8,.85,.9,.95,1])
    #ax3.set_xticklabels(['0.01', '0.02', '0.03', '0.04', '0.05'])
    #plt.setp(ax1.get_xticklabels(), visible=False)
    #plt.setp(ax2.get_xticklabels(), visible=False)

    ax1.set_xlim([0.0,0.05])
    ax2.set_xlim([0,0.05])
    ax3.set_xlim([0,0.05])
    ax4.set_xlim([0.0,0.05])
    ax1.set_ylim([0.5, 1.01])
    ax2.set_ylim([0.5, 1.01])
    ax3.set_ylim([0.5, 1.01])
    ax4.set_ylim([0.5, 1.01])
    ax1.set_title('(A)')
    ax2.set_title('(B)')
    ax3.set_title('(C)')
    ax4.set_title('(D)')
    
    plt.show()  


