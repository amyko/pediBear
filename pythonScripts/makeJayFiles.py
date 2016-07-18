import numpy as np
import matplotlib.pyplot as plt
import pdb
import os.path

badInd = [31,32,39,40,41,98,6,7,8,9,10,11,12,21,22,26,27,34,48,49,50,52,53,87,88,90,91] #inbred or looped

def makeIDFile(inPath, outPath, pedPath):

    indexPath =  os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/jays/102jays.index2name"

    #build dictionary
    names = np.loadtxt(indexPath, dtype=str, delimiter='\t', usecols=[1], skiprows=1)
    indices = np.loadtxt(indexPath, dtype=int, delimiter='\t', usecols=[0], skiprows=1)
    name2index = dict(zip(names,indices))


    #open files
    infile = open(inPath)
    outfile = open(outPath, 'wb')
    outfile.write("index\tname\tsex\tsampled\n")
    i=0
    
    #set
    myset = set()

    #typed individuals
    for line in infile:
        
        fields = line.split()
        
        if(name2index[fields[1]] in badInd): cotinue

        sex = fields[4]
        if(fields[4]=='1'):
            sex = 1
        elif(fields[4]=='2'):
            sex = 0
        else: pdb.set_trace()
        


        outfile.write("%d\t%s\t%d\t%d\n" %(i, fields[1], sex, 1))    
        
        i+=1
        myset.add(fields[1])


        
    infile.close()
    

    
    #all other individuals
    infile = open(pedPath)
    infile.readline()
    for line in infile:
        
        fields = line.split()
        
        for j in range(0,3):
        
            if(fields[j] in myset or fields[j]=='0' or name2index[fields[j]] in badInd): continue

        
            sex = 0
            if(fields[3]=='M'): 
                sex = 1
            elif(fields[3]=='F'):
                sex = 0
            else: pdb.set_trace()
            
            outfile.write("%d\t%s\t%d\t%d\n" %(i, fields[j], sex, 0))    
            
            i+=1
            myset.add(fields[j])

        
    infile.close()
    outfile.close()



def removeBadIndiv(inPath, outPath, indexPath):
    
    #name2Index
    ids = np.loadtxt(indexPath, dtype=int, skiprows=1, usecols=[0], delimiter='\t')
    names = np.loadtxt(indexPath, dtype=str, skiprows=1, usecols=[1], delimiter='\t')

    
    name2index = dict(zip(names,ids))

    
    #clean data
    outfile = open(outPath, 'wb')
    infile = open(inPath)
    
    for line in infile:
        fields = line.split()
        #pdb.set_trace()
        if(name2index.get(fields[1]) in badInd): continue


        for x in fields:
            outfile.write("%s " %x)
        outfile.write("\n")
        
        #outfile.write(line)
        
    
    #pdb.set_trace()
    infile.close()
    outfile.close()


def makeGenoFiles(inPath, outPath, posPath):

    #open files
    chroms = np.loadtxt(posPath, dtype=int, delimiter='\t', usecols=[0])
    pos = np.loadtxt(posPath, dtype=int, delimiter='\t', usecols=[3])
    infile = open(inPath)
    outfile = open(outPath+str(1), 'wb')
    outfile.write("POS\tA1\tA2\t")
    for i in range(0,75):
        outfile.write("%d\t" %i)
    outfile.write('\n')

    
    currChrom = 0
    col = 6
    chromNum = 1
    
    #pdb.set_trace()
    
    for i in range(0,len(chroms)):
        
        chrom = chroms[i]
        
        #open new file
        if(chrom!=currChrom and chrom!=1):
            
            print chrom
            outfile.close()
            chromNum+=1
            outfile = open(outPath+str(chromNum), 'wb')
            
            #write header
            outfile.write("POS\tA1\tA2\t")
            for j in range(0,75):
                outfile.write("%d\t" %j)
            outfile.write('\n')
                
            currChrom = chrom
        
        #get SNPS for this position
        snps = np.loadtxt(inPath, dtype=str, delimiter=' ', usecols=[col, col+1])
        col+=2
        
        missingData = False
        for snp in snps:
            if(snp[0]=='0' or snp[1]=='0'): 
                missingData = True
                break
        if(missingData==True): continue
        
        
        #get alleles
        a1 = snps[0][0]
        a2 = ''
        for s in snps:
            for j in range(0,2):
                if(s[j]!=a1): 
                    a2 = s[j]
                    break
                
        if(a2==''): 
            print("monomorphic")
            continue
        
        #write pos, a1, a2
        outfile.write("%d\t%s\t%s\t" %(pos[i], a1, a2))
        
        #write genotype
        for snp in snps:
            if((snp[0]!=a1 and snp[0]!=a2) or (snp[1]!=a1 and snp[1]!=a2)): 
                print("not biallelic")
                pdb.set_trace()
            outfile.write("%s%s\t" %(snp[0],snp[1]))
        outfile.write("\n")
        
    infile.close()
        
        
def makeMapFile(genoPath, outPath):

    outfile = open(outPath, 'wb')
    
    snpID = 1

    for chrom in range(1,35):
        
        infile = open(genoPath+str(chrom))
        infile.readline()
        
        for line in infile:
            
            fields = line.split()
            
            outfile.write("%d\tsnp%d\t%d\t%s\n" %(chrom, snpID, 0, fields[0]))
            snpID+=1
        infile.close()
        
    outfile.close()
        


def replacePedID(mydir):
    
    infile = open(mydir+"75jays.pruned.50_1.ped")
    outfile = open(mydir + "75jays.pruned.50_1.founders.ped", 'wb')
    
    names = np.loadtxt(mydir+"75jays.index2name", dtype=str, delimiter='\t', skiprows=1, usecols=[1])
    ids = np.loadtxt(mydir+"75jays.index2name", dtype=int, delimiter='\t', skiprows=1, usecols=[0])
    
    name2id = dict(zip(names, ids))
    
    for line in infile:
        
        fields = line.split()
        
        #change name
        fields[1] = str(name2id[fields[1]] + 1)
        
        #change parents
        fields[2] = '0'
        fields[3] = '0'
        
        for x in fields:
            outfile.write("%s " %x)
        outfile.write("\n")
        
    infile.close()
    outfile.close()
    
def makeSexFile(mydir):
    
    outfile = open(mydir+"75jays.founders.sex", 'wb')
   
    sexes = np.loadtxt(mydir+"75jays.index2name", dtype=int, delimiter='\t', skiprows=1, usecols=[2])
    
    i=1
    while(i<76):
        sex = sexes[i-1]
        outfile.write("1\t%d\t%d\n" %(i,sex))
        i+=1
    outfile.close() 

if __name__ == "__main__":

    mydir = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/jays/"
    
    #file names
    inPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/jays/75jays.pruned.1000_15.ped"
    outPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/jays/75jays.pruned.1000_15.geno."
    posPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/jays/102jays.pruned.1000_15.map"
    indexPath = os.path.expanduser('~') + "/Google Drive/Research/pediBear/data/jays/102jays.index2name"

    #removeBadIndiv(inPath, outPath, indexPath)

    #makeGenoFiles(inPath, outPath, posPath)
    
    #replacePedID(mydir)
    makeSexFile(mydir)
    
    #makeMapFile(outPath, posPath)
