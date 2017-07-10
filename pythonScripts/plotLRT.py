import numpy as np
import matplotlib.pyplot as plt
import pdb


def lkhd(fileName, up, down, nVisit):
    
    infile = open(fileName)
    
    toReturn = np.array([])
    
    process = False
    for line in infile:
        
        fields = line.split()

        
        if(fields[0]=='>'):
            
            if(fields[1]==up and fields[2]==down and fields[3]==nVisit):
                process = True
                continue
            else:
                if(process==True):
                    break

        if(process==False): continue
        
        if(process==True):
            toReturn = np.append(toReturn, float(fields[2]))
            
            if(len(toReturn)==3242): print (fields[0], fields[1])
            
        
    infile.close()
            
    return toReturn
    


if __name__=="__main__":
    
    
    fileName = "/Users/kokocakes/Google Drive/Research/pediBear/data/simulations/genotypes/test12.all.pruned.10k.dep.pairwise"

    lkhd_000 = lkhd(fileName, '0','0','0')
    lkhd_551 = lkhd(fileName, '2','2','2')
    diff = lkhd_000 - lkhd_551
    print min(diff)

    
    fileName = '/Users/kokocakes/Google Drive/Research/pediBear/data/simulations/'
    lkhd_000 = np.loadtxt(fileName+"unrel.unlinked.merlin.lkhd",dtype=float,usecols=[2])
    lkhd_551 = np.loadtxt(fileName+'thirdCousin.unlinked.merlin.lkhd',dtype=float, usecols=[2])
    diff= lkhd_000 - lkhd_551
    
    lkhd_000 = np.loadtxt(fileName+"unrel.merlin.lkhd",dtype=float,usecols=[2])
    lkhd_551 = np.loadtxt(fileName+'thirdCousin.merlin.lkhd',dtype=float, usecols=[2])
    diff2 = lkhd_000 - lkhd_551
    
    fig = plt.figure(facecolor='white')
    plt.ylabel("Frequency")
    plt.xlabel("L(unrelated) - L(third cousins)", style='italic')
    weights = np.ones_like(diff)/len(diff)
    weights2 = np.ones_like(diff2)/len(diff2)
    bins = np.linspace(-10, 5, 100)
    plt.hist(diff,bins,color='blue', label='Unlinked markers', weights=weights)
    plt.hist(diff2,bins,color='green', label='Linked markers', weights=weights2)
    plt.legend(loc='upper left')
    plt.show()
    
    pdb.set_trace()
    
    
    #plot likelihood surface comparison
    x = [1,2,3,4,5]
    tickMarks = ['.07','.20','.40','.67','1']
    full = [-18716.360,-17410.471,-14504.051,-11584.392,-8545.180]
    pair = [-18716.362569, -18507.8639172,-18118.430455,-17453.2292502,-16946.5014504]
    cond = [-18716.362569,-17673.869308999998,-15726.701997999999,-12400.695975000002,-9867.056975999993]
    
    
    fig = plt.figure(facecolor='white')
    ax = fig.add_subplot(111)
    
    plt.plot(x,full,color='red', label='Full likelihood')
    plt.plot(x,cond,color='blue', label='Composite likelihood A')
    plt.plot(x,pair,color='green', label='Composite likelihood B')
    
    #plt.plot(x,full,color='red')
    #plt.plot(x,cond,color='blue')
    #plt.plot(x,pair,color='green')
    
    plt.legend(loc='upper left')
    plt.xlabel("Proportion of correct pairwise relationships")
    plt.ylabel("Log likelihood")
    ax.set_xticks(x)
    ax.set_xticklabels(tickMarks)
    plt.setp(tickMarks)
    
    plt.show()
    
    pdb.set_trace()
    
    


    
    