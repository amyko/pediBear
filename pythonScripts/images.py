import matplotlib.pyplot as plt
from PIL import Image
import matplotlib.image as mpimg
import numpy as np
import pdb

def loadImage(imagePath):

    img = Image.open(imagePath)
    #img = img.resize((w,h), Image.ANTIALIAS)
    w,h = img.size
    img = img.crop((50, 300, w-50, h-300))
    return img

def plotSimPed():
    
    dir = '/Users/kokocakes/Google Drive/Research/pediBear/manuscript/'
    
    simA = loadImage(dir+"simPed1.tiff")
    simB = loadImage(dir+"simPed2.tiff")
    simC = loadImage(dir+"simPed4.tiff")
    
    fig = plt.figure(facecolor='white', figsize=[20,5])
    ax1 = fig.add_subplot(1,3,1)
    ax2 = fig.add_subplot(1,3,2)
    ax3 = fig.add_subplot(1,3,3)
    #f, (ax1,ax2,ax3) = plt.subplots(1,3)
    ax1.imshow(simA)
    ax2.imshow(simB)
    ax3.imshow(simC)
    
    ax1.axis('off')
    ax2.axis('off')
    ax3.axis('off')
    ax1.set_title('(A)', fontsize=26)
    ax2.set_title('(B)', fontsize=26)
    ax3.set_title('(C)', fontsize=26)
    plt.show()
    
def plotGreenland():
    
    dir = '/Users/kokocakes/Google Drive/Research/pediBear/manuscript/'

    img = loadImage(dir+"greenland.n2.tiff")
    fig = plt.figure(facecolor='white', figsize=[15,10])
    plt.imshow(img)
    plt.axis('off')
    plt.show()

if __name__=='__main__':
    
    plotSimPed()
    
    pdb.set_trace()
    
    dir = '/Users/kokocakes/Google Drive/Research/pediBear/manuscript/'

    img = loadImage(dir+"greenland.1844snps.n2.tiff")
    fig = plt.figure(facecolor='white', figsize=[12,8])
    plt.imshow(img)
    plt.axis('off')
    plt.show()
    