import numpy as np
import matplotlib.pyplot as plt



columns = np.genfromtxt('test.out', delimiter=' ', skip_header  = 3, skip_footer = 1, usecols = (1, 3, 5, 8)) # 1- iteracja, 3 - zaakcept krok, 5 - energy, 8 - gyration

fig = plt.figure()
plt.plot(columns[:,1], columns[:, 2], "-", linewidth =2)
plt.show()