import numpy as np
import pdb
import random



# sample individuals from first 3 generations
def sample(N, n, d):

    #init samples
    samples = set()
    numSampled = 0
    totalUnits = 0.0
    for x in range(1,d+2):
        totalUnits += x
    
    
    #sample
    while(numSampled < n):
        
        # choose depth
        p = np.random.random()
        depth = -1
        cdf = 0
        for x in range(1,d+2):
            
            cdf += x/totalUnits
            
            if(p < cdf):
                depth = d-x+1
                break
        
        
        #choose individual
        indiv = np.random.randint(0,N+1)
        id = "%d_%d" %(depth, indiv)
        if(id in samples): continue
        
        numSampled+=1
        
        print id
        

    return samples
    
    
def polygamousPedigree(samples):
    
    
    

if __name__ == "__main__":

    N = 200
    n = 20
    d = 2

    sample(N, n, d)