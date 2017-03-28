import numpy as np
import math

def geval(L, b, d, i):
    P_i = (1/(2*b))*(b**2 + L[2*i-2]**2 - L[2*i-1]**2)
    h_i = math.sqrt(L[2*i-2]**2 - P_i**2)

    if i == 1:
        X_i = (math.sqrt(3)/6)*(2*b+d-3*P_i)
        Y_i = (1/2)*(d+P_i)
    elif i == 2:
        X_i = (-1)*(math.sqrt(3)/6)*(b+2*d)
        Y_i = (1/2)*(b-2*P_i)
    elif i == 3:
        X_i = (-1)*(math.sqrt(3)/6)*(b-d-3*P_i)
        Y_i = (-1)*(1/2)*(b+d-P_i)

    return [h_i, X_i, Y_i]



def task1():
    #L = [3,3,3,3,3,3]
    L = [8,8,8,8,8,8]
    b = 15
    d = 1
    i = [1,2,3]
    for index in i:
        print(geval(L, b, d, index))

print("Running...")
task1()