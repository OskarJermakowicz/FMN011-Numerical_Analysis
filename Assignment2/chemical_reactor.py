import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import scipy.optimize as sci

depth = [0, 0.5, 1, 1.5, 2, 2.5, 3]
temp = [70, 70, 55, 22, 13, 10, 10]

cs = CubicSpline(depth, temp, bc_type='clamped')
x = np.linspace(0,3)

# Task 1
plt.plot(depth, temp, 'o', label='data')
plt.plot(x, cs(x), label='Interpolated curve')

# Task 2
tdepth = sci.fsolve(cs, 1, args=2, xtol=10e-6)
plt.axvline(x=tdepth, linestyle='dashdot', label='Thermocline depth', lw=0.5)

# Task 3
plt.plot(x,-cs(x,1),label='Flux')

# Task 4
plt.plot(x, cs(x, 1), 'b', label='First derivative')
plt.plot(x, cs(x, 2), 'r', label='Second derivative')

# Task 5
t = sci.fsolve(cs, 1.07874, args=2, xtol=10e-6)
print("Approximation of depth at 50C:", t[0])
print("Approximation of temperature at 1.7m:", cs(1.7))

plt.axvline(x=1.7, linestyle='dashdot', label='Approx. of 1.7m', lw=0.5, color='g')
plt.axhline(y=cs(1.7), linestyle='dashdot', lw=0.5, color='g')
plt.axhline(y=50, linestyle='dashdot', label='Approx. of 50C', lw=0.5, color='b')
plt.axvline(x=1.07874, linestyle='dashdot', lw=0.5, color='b')

# Show plot
plt.legend(loc='lower left', ncol=2)
x_label = plt.xlabel('\nDepth (meter)')
y_label = plt.ylabel('\nTemperature (celsius)')
plt.show()