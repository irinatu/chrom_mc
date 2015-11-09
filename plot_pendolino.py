import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
#from pylab import *

inp = raw_input("Put pendolino output: ")
out = inp.split('.')[0]

columns = np.genfromtxt(inp, delimiter=' ', skip_header  = 3, skip_footer = 3, usecols = (1, 3, 5, 8)) # 1- iteracja, 3 - zaakcept krok, 5 - energy, 8 - gyration

fig = plt.figure()
#x_sm = np.linspace(columns[:50000,1].min(), columns[:50000, 1].max(), 50)
#y_sm = interp1d(columns[:50000,1], columns[:50000, 2], x_sm)
#plt.plot(columns[:50000,1], columns[:50000, 2], "-", x_sm, y_sm, "--", linewidth =0.5)
#plt.plot(columns[:50000,1], columns[:50000, 2], "-", linewidth =0.5)
plt.scatter(columns[:50000,1], columns[:50000, 2], linewidth = 0.5)
plt.title("Energy")
plt.show()
fig.savefig(out+"_energy.png")

plt.plot(columns[:,1], columns[:, 3], "-", linewidth =0.5)
plt.title("Gyration")
plt.show()
fig.savefig(out+"_gyrat.png")


