from matplotlib.mlab import PCA
import numpy as np
import optparse
from sys import argv

def clip(arr, stddevs=10):
  arr = np.ma.masked_invalid(arr)
  mean = np.mean(arr)
  stddev = np.var(arr) ** 0.5
  np.clip(arr, 0, mean + stddevs * stddev, out=arr)
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
    
    
if __name__=="__main__":
    
    optparser = optparse.OptionParser(usage = "%prog [<options>]")
    optparser.add_option('-m', type = "string", default = "", dest="Matrix", help = "Numpy matrix in npy format")
    (opts,args) = optparser.parse_args()
    if len(argv) ==1:
        print optparser.format_help()
        exit(1)
    
    arr = clip(np.load(opts.Matrix))
    arr_nor = dist_normalization(arr)
    pca_res = principle_component(arr_nor).Y[:,0]
    #print len(pca_res)
    zero_crossings = np.where(np.diff(np.sign(pca_res)))[0] #indexes of the sign henging
    start = 0
    nr = 1
    for i in zero_crossings:
        print "%i\t2L\t%i\t%i" %(nr, start, i)
        nr +=1
        start=i+1
    if i != len(pca_res)-1:
        print "%i\t2L\t%i\t%i" %(nr, start, len(pca_res)-1)
    
