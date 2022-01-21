import numpy as np
import matplotlib.pyplot as plt


file = r'/Users/zuzannaj/Documents/Spring2022/scripts/tests/test.py2022-01-20 13:26:15.899376boxMC.boxMC'
data = np.loadtxt(file)
data = data[1:]
vol = data[::2]
p = data[1::2]

plt.plot(vol,p)
plt.show()