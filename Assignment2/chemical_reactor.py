import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import scipy.optimize as sci

k = 0.01
depth = [0, 0.5, 1, 1.5, 2, 2.5, 3]
temp = [70, 70, 55, 22, 13, 10, 10]

cs = CubicSpline(depth, temp, bc_type='clamped')
x = np.linspace(0,3)

# Task 1
# plt.plot(depth, temp, 'o', label='data')
# plt.plot(x, cs(x), label='Interpolated curve')

# Task 2
tdepth = sci.fsolve(cs, 1, args=2, xtol=10e-6)
plt.axvline(x=tdepth, linestyle='dashdot', label='Thermocline depth', lw=0.5)

# Task 3
plt.plot(x,-cs(x,1),label='Flux')

# Show plot
plt.legend(loc='lower left', ncol=2)
x_label = plt.xlabel('\nDepth (centimeter)')
y_label = plt.ylabel('\nTemperature (celsius)')
plt.show()