import numpy as np
our = np.load('HiC.npy')
four = np.load('4_sum_normed.npy')
our.shape
four.shape
our = our[12:257,12:257]
four = four[12:257,12:257]
four = np.nan_to_num(four)
import scipy.stats
scipy.stats.pearson(four.reashape(60025), our.reshape(60025))
scipy.stats.pearsonr(four.reashape(60025), our.reshape(60025))
scipy.stats.pearsonr(four.reshape(60025), our.reshape(60025))
four[:100,:100].shape
scipy.stats.pearsonr(four[:100,:100].reshape(10000), our[:100,:100].reshape(10000))
four.shape
four
lour = np.log10(our)
lour
lfour = np.log10(four)
our
np.log10([1e-15, -3.])
np.log10([1e-15, 0.])
lfour
np.nan_to_num(lfour)
np.nan_to_num(our)
lfour
our = np.nan_to_num(our)
lfour = np.nan_to_num(lfour)
lour = np.nan_to_num(lour)
scipy.stats.pearsonr(lfour.reshape(60025), lour.reshape(60025))
scipy.stats.pearsonr(lfour[:100,:100].reshape(10000), lour[:100,:100].reshape(10000))
np.amax(lfour)
np.amin(lfour)
np.amin(lour)
np.amax(lour)
np.log10(four)
a = np.log10(four)
np.amin(a)
np.where(lfour == -1.7976931348623157e+308)
np.where(lfour == -1.7976931348623157e+308).shape
len(np.where(lfour == -1.7976931348623157e+308))
len(np.where(lfour == -1.7976931348623157e+308)[0])
four
te = np.load('4_sum_normed.npy')
te[12:257]
te[13:257]
te[14]
te[12]
:
te[12:257,12:257]
np.isnan(te[12:257,12:257])
np.where(np.isnan(te[12:257,12:257])=='True')
np.where(np.isnan(te[12:257,12:257])=='True').shape
len(np.where(np.isnan(te[12:257,12:257])=='True'))
len(np.where(np.isnan(te[12:257,12:257])=='True')[0])
np.isinf(te[12:257,12:257])
len(np.where(np.isinf(te[12:257,12:257])=='True')[0])
np.where(np.isnan(te[12:257,12:257])=='True')[0]
np.where(np.isinf(te[12:257,12:257])=='True')[0]
np.where(np.isinf(te[12:257,12:257]))[0]
np.where(np.isnan(te[12:257,12:257]))[0]
log(0)
import math
log10(0)
scipy.stats.pearsonr(four.reshape(60025), our.reshape(60025))
plt
import matplotlib.pyplot as plt
plt.plot(four.reshape(60025), our.reshape(60025))
plt.show()
four.reshape(60025)
plt.plot([1,2,3,4],[1,2,3,4])
plt.show()
plt.plot([1,2,3,4],[1,2,3,4],'ro')
plt.show()
plt.plot(four.reshape(60025), our.reshape(60025), 'ro')
plt.show()
plt.show()
plt.plot(four.reshape(60025), our.reshape(60025), 'ro')
plt.show()
plt.show()
plt.plot(four.reshape(60025), our.reshape(60025), 'ro')
plt.xlim(xmax=5000)
plt.plot(four.reshape(60025), our.reshape(60025), 'ro')
plt.xlim(xmax=5000)
plt.show()
plt.plot(four.reshape(60025), our.reshape(60025), 'ro')
plt.xlim(xmax=1000)
plt.show()
plt.xlim(xmax=1000)
plt.plot(four.reshape(60025), our.reshape(60025), 'ro')
plt.xlim(xmax=400)
plt.show()
plt.plot(four.reshape(60025), our.reshape(60025), 'ro')
plt.xlim(xmax=100)
plt.show()
np.amin(lour)
np.amax(lour)
np.amax(lfour)
np.amin(lfour)
a = 0.9e-200
a
lfour[lfour<0.]=0.1e-100
np.amin(lfour)
plt.plot(four.reshape(60025), our.reshape(60025))
scipy.stats.pearsonr(lfour.reshape(60025), lour.reshape(60025))
np.amin(four)
np.log10??
np.log10(0)
plt.plot(lfour.reshape(60025), lour.reshape(60025))
plt.show()
plt.plot(lfour.reshape(60025), lour.reshape(60025),'ro')
plt.show()
plt.plot(lfour.reshape(60025), lour.reshape(60025),'ro')
scipy.stats.pearsonr(four.reshape(60025), our.reshape(60025))
%save corr_HiCmaps.py
%history -f corr_HiCmaps.py
