import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
#from pylab import *

inp = raw_input("Put pendolino output: ")
inp_list = inp.split()
out = inp_list[0].split('.')[0]

for file in inp_list:
    col = np.genfromtxt(file, delimiter=' ', skip_header  = 3, skip_footer = 3, usecols = (1, 3, 5, 8)) # 1- iteracja, 3 - zaakcept krok, 5 - energy, 8 - gyration
    try:
        columns
        columns = np.vstack((columns, col))
    except NameError:
        columns = col

   
#columns = np.vstack(columns, col)

fig = plt.figure()
#x_sm = np.linspace(columns[:50000,1].min(), columns[:50000, 1].max(), 50)
#y_sm = interp1d(columns[:50000,1], columns[:50000, 2], x_sm)
#plt.plot(columns[:50000,1], columns[:50000, 2], "-", x_sm, y_sm, "--", linewidth =0.5)
#plt.scatter(columns[:50000,1], columns[:50000, 2], linewidth = 0.5)
plt.plot( columns[:, 2], "-", linewidth =0.5)
plt.axis([0,len(columns[:, 2]), columns[:, 2].min(), columns[:, 2].max()])
plt.title("Energy")
plt.show()
#plt.scatter(columns[:,1], columns[:, 2], linewidth = 0.5)
#plt.show()
fig.savefig(out+"_energy.png")
fig = plt.figure()
plt.plot( columns[:, 3], "r--", linewidth =0.5)
plt.axis([0,len(columns[:, 3]), 0, columns[:, 3].max()])
dlugosc = len(columns[:, 3])/5
start_mean = len(columns[:, 3]) - dlugosc
print np.mean(columns[start_mean:, 3])
plt.title("Gyration")
plt.show()
fig.savefig(out+"_gyrat.png")


