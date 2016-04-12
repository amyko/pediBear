import numpy as np
import matplotlib.pyplot as plt
import pdb
import os.path

testName = "test3"
inPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/simulations/results/" + testName + ".kinship.acc"
nIndiv = 6
nPairs = nIndiv*(nIndiv-1)/2


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
    acc = float(fields[2])

    if (i,j) in accDict:
        accDict[(i,j)].append(acc)
    else:
        accDict[(i,j)] = [acc]


infile.close()


#make box plot
data = [accDict[key] for key in accDict.keys()]
means = [np.mean(x) for x in data]
xdata = [i for i in range(1,nPairs+1)]

#pdb.set_trace()

plt.figure()
plt.boxplot(data)
plt.scatter(xdata, means)
plt.ylim([-.1,1.1])
plt.xlabel("pair")
plt.ylabel("accuracy")
plt.title("pairwise accuracy for 100 simulations; " + testName)
plt.show()
