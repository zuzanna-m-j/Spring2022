import matplotlib.pyplot as plt
import numpy as np
import math
file = r'/scripts/tests/test.py2022-01-20 09:50:04.384755.sdf'
def extrapolate(s, dx, xmax, degree=5):
  s = s[:int(math.ceil(xmax/dx))]
  x = np.arange(0,xmax,dx)
  x += dx/2
  x = x[:len(s)]
  p = np.polyfit(x, s, degree)
  return np.polyval(p, x),x

h = np.loadtxt(file,skiprows=1)
for i in range(1):
    d = h[0,1:200]
    plt.plot(d,label=i)

plt.legend()
plt.show()