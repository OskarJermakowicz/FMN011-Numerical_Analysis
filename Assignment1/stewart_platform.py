import numpy as np
import math

# Returns a list of three elements, for each i, which includes h, x and y for each.
def geval(L, b, d):
    res = []

    for i in [1,2,3]:
        P_i = (1/(2*b))*(b**2 + L[2*i-2]**2 - L[2*i-1]**2)
        h_i = math.sqrt(L[2*i-2]**2 - P_i**2)
        res.append([h_i, 0, 0])

    # XP1 and YP1
    res[0][1] = (math.sqrt(3) / 6) * (2 * b + d - 3 * P_i)
    res[0][2] = (1 / 2) * (d + P_i)

    # XP2 and YP2
    res[1][1] = (-1) * (math.sqrt(3) / 6) * (b + 2 * d)
    res[1][2] = (1 / 2) * (b - 2 * P_i)

    # XP3 and YP3
    res[2][1] = (-1) * (math.sqrt(3) / 6) * (b - d - 3 * P_i)
    res[2][2] = (-1) * (1 / 2) * (b + d - P_i)

    return res

# Determines the left hand side of equations for calculating the X-coordinates of the vertices in the upper frame. Ideally the left hand side should equal zero.
def stf(a, gev, xt):
    equations = [0, 0, 0]

    equations[0] = a**2 + 2*xt[0]*xt[1] - 2*xt[0]*(gev[0][1] + math.sqrt(3)*(gev[0][2]-gev[1][2])) - 2*gev[1][1]*xt[1] - ((math.sqrt(3)*gev[0][1]-gev[0][2])**2 + (gev[0][0]**2 + gev[1][0]**2) - 4*gev[0][1]**2 - gev[1][1]**2) + 2*math.sqrt((gev[0][0]**2 - 4*(xt[0]-gev[0][1])**2)*(gev[1][0]**2 - (xt[1] - gev[1][1])**2))
    equations[1] = a**2 - 4*xt[0]*xt[2] - 2*xt[0]*(gev[0][1] - 3*gev[2][1] + math.sqrt(3)*(gev[0][2]-gev[2][2])) - 2*xt[2]*((-3)*gev[0][1] + gev[2][1] + math.sqrt(3)*(gev[0][2] - gev[2][2])) - ((math.sqrt(3)*(gev[0][1] + gev[2][1]) - gev[0][2] + gev[2][2])**2 + (gev[0][0]**2 + gev[2][0]**2) - 4*gev[0][1]**2 - 4*gev[2][1]**2) + 2*math.sqrt((gev[0][0]**2 - 4*(xt[0] - gev[0][1])**2)*(gev[2][0]**2 - 4*(xt[2] - gev[2][1])**2))
    equations[2] = a**2 + 2*xt[1]*xt[2] - 2*xt[2]*(gev[2][1] + math.sqrt(3)*(gev[1][2] - gev[2][2])) - 2*gev[1][1]*xt[1] - ((math.sqrt(3)*gev[2][1] - gev[1][2] + gev[2][2])**2 + (gev[1][0]**2 + gev[2][0]**2) - gev[1][1]**2 - 4*gev[2][1]**2) + 2*math.sqrt((gev[1][0]**2 - (xt[1] - gev[1][1])**2)*(gev[2][0] - 4*(xt[2] - gev[2][1])**2))

    return equations


def task1():
    #L = [3, 3, 3, 3, 3, 3]
    L = [8, 8, 8, 8, 8, 8]
    b = 15
    d = 1

    for i in geval(L, b, d):
        print(i)

def task2():
    L = [11.5, 11.5, 11.5, 11.5, 11.5, 11.5]
    b = 15
    d = 1
    a = 10
    xt = [2, 2, 2]

    print(stf(a, geval(L, b, d), xt))

print("-- Task 1:")
task1()

print("\n-- Task 2:")
task2()