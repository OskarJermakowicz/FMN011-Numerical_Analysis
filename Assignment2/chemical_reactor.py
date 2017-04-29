import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

def task1():
    k = 0.01
    depth = [0, 0.5, 1, 1.5, 2, 2.5, 3]
    temp = [70, 70, 55, 22, 13, 10, 10]

    cs = CubicSpline(depth, temp, bc_type='clamped')
    x = np.linspace(0,3)
    plt.plot(depth, temp, 'o', label='data')
    plt.plot(x, cs(x), 'r', label='Interpolated curve')
    plt.legend(loc='lower left', ncol=2)
    plt.show()

print("-- Task 1:")
task1()