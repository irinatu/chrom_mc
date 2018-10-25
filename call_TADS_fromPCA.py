import matplotlib
matplotlib.use('Agg')
from matplotlib.mlab import PCA
import numpy as np
import optparse, scipy.ndimage
from sys import argv
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

def clip(arr, stddevs=10):
    arr = np.ma.masked_invalid(arr)
    mean = np.mean(arr)
    stddev = np.var(arr) ** 0.5
    np.clip(arr, 0, mean + stddevs * stddev, out=arr)
    arr = np.ma.filled(arr, 0)
    return arr

def clip_and_blur(arr, stddevs=5, blur=1): # function from Krzysiek script
    arr = np.ma.masked_invalid(arr)
    mean = np.mean(arr)
    stddev = np.var(arr) ** 0.5
    np.clip(arr, 0, mean + stddevs * stddev, out=arr)
    arr = np.ma.filled(arr, 0)
    scipy.ndimage.gaussian_filter(arr, blur, output=arr)
    np.clip(arr, mean * 0.01, mean + stddevs * stddev, out=arr)
    return arr
  
def dist_normalization(ma):
    l = ma.shape[0]
    row = l
    col = 0
    #print ma
    for di in range(l):
        mean = ma.diagonal(di).sum()/ma.diagonal(di).shape[0]
        if mean == 0.0:
            #print "Zero", mean, di
            continue
        #print mean, di
        #print ma.diagonal(di)
        np.fill_diagonal(ma[0:row, col:l], ma.diagonal(di)/float(mean))
        row -= 1
        col += 1
        if row == 1: break
    ma = np.triu(ma).T + np.triu(ma)
    np.fill_diagonal(ma, ma.diagonal()/2.0)
    return ma
    
def principle_component(mac):
    zmiennosc = np.std(mac, axis=0) # find columns with the same values
    if len(zmiennosc) == mac.shape[0]: 
        mac1 = np.delete(mac, np.where(zmiennosc==0.0)[0].tolist(), axis =1) # delete columns with the same values
        #print np.std(mac1, axis=0), mac1.shape
    else:
        pass 
        #print len(zmiennosc), arr_nor.shape[0]
    
    results = PCA(mac1)
    return results
    
def TADSdef(pca):
    zero_crossings = np.where(np.diff(np.sign(pca)))[0] #indexes of the sign henging
    start = 1
    nr = 1
    for i in zero_crossings:
        if i+1 != start:
            print "%i\t2L\t%i\t%i" %(nr, start, i+1)
            nr +=1
            start=i+2
    if i != len(pca)-1:
        print "%i\t2L\t%i\t%i" %(nr, start, len(pca_res)-1)    
    
if __name__=="__main__":
    
    optparser = optparse.OptionParser(usage = "%prog [<options>]")
    optparser.add_option('-m', type = "string", default = "", dest="Matrix", help = "Numpy matrix in npy format (real) - will be cliped and blured")
    optparser.add_option('-c', type = "string", default = "", dest="Matrix2", help = "Second Numpy matrix in npy format (simulation)")
    (opts,args) = optparser.parse_args()
    if len(argv) ==1:
        print optparser.format_help()
        exit(1)
    
    arr = np.load(opts.Matrix)
    arr = clip_and_blur(arr)
    arr_nor = dist_normalization(arr)
    pca_res = principle_component(arr_nor).Y[:,0]
    y_ki = savgol_filter(pca_res, 51,3)
    
    arr2 = np.load(opts.Matrix2)
    arr_nor2 = dist_normalization(arr2)
    pca_res2 = principle_component(arr_nor2).Y[:,0]
    y_ki2 = savgol_filter(pca_res2, 51,3)

    plt.plot(pca_res, linewidth=2, label=opts.Matrix.split(".")[0])
    plt.plot(pca_res2, linewidth=2, label=opts.Matrix2.split(".")[0])
    plt.plot(y_ki, label=opts.Matrix.split(".")[0]+"_smooth" )
    plt.plot(y_ki2, label=opts.Matrix2.split(".")[0]+"_smooth" )
    plt.legend()
    plt.show()
    
    print "FIRST"
    TADSdef(pca_res)
    print "SECOND"
    TADSdef(pca_res2)
    
      
    #print len(pca_res)
    #zero_crossings = np.where(np.diff(np.sign(pca_res)))[0] #indexes of the sign henging
    #start = 0
    #nr = 1
    #for i in zero_crossings:
    #    print "%i\t2L\t%i\t%i" %(nr, start, i)
    #    nr +=1
    #    start=i+1
    #if i != len(pca_res)-1:
    #    print "%i\t2L\t%i\t%i" %(nr, start, len(pca_res)-1)
    
