import math as m
import numpy as np
import matplotlib.pyplot as plt

M = .00     # maximum camber [% chord]
P = .00     # maximum camber position [% chord]
t = .12     # thickness [% chord]
a = 15.     # half-angle (only for diamond airfoil)
n = 200     # number of panels

# IMPORTANT
# - points are in counter-clockwise order
# - there are no repeated points
# - lines formed by points do not intersect
# - custom points are in (x y) format in mesh/body.txt
# e.g. mesh/body.txt:
# x1 y1
# x2 y2
# x3 y3
# ...

def airfoil(M, P, t, n):
    x = np.linspace(0, 1, n//2)
    y = np.zeros(n//2)
    for i in range(len(x)):
        if x[i] < P:
            y[i] = (M/P**2)*(2*P*x[i] - x[i]**2)
        else:
            y[i] = (M/(1 - P)**2)*(1 - 2*P + 2*P*x[i] - x[i]**2)
    a0, a1, a2, a3, a4 = .2969, -.1260, -.3516, .2843, -.1036
    yt = t/.2*(a0*x**.5 + a1*x + a2*x**2 + a3*x**3 + a4*x**4)
    yu, yl = y + yt, y - yt
    X = np.concatenate((np.flip(x[:]), x[1:-1]))
    Y = np.concatenate((np.flip(yu[:]), yl[1:-1]))
    XY = np.column_stack((X, Y))
    np.savetxt('mesh/body.txt', XY, fmt='%8.5f')
    custom()

def diamond(a, n):
    le = (0, 0)
    te = (1, 0)
    top = (0.5, 0.5*np.tan(np.radians(a)))
    btm = (0.5, -0.5*np.tan(np.radians(a)))
    def interpolate(start, end, n):
        return np.linspace(start, end, n, endpoint=False)
    te_to_top = interpolate(te, top, n//4)
    top_to_le = interpolate(top, le, n//4)
    le_to_btm = interpolate(le, btm, n//4)
    btm_to_te = interpolate(btm, te, n//4)
    XY = np.vstack([te_to_top, top_to_le, le_to_btm, btm_to_te])
    np.savetxt('mesh/body.txt', XY, fmt='%8.5f')
    custom()

def custom():
    XY = np.loadtxt('mesh/body.txt')
    X, Y = XY[:, 0], XY[:, 1]
    xmin, xmax = np.min(X), np.max(X)
    ymin, ymax = np.min(Y), np.max(Y)
    X = (X - xmin)/(xmax - xmin)
    Y = (Y - ymin)/(xmax - xmin)
    Y = Y - .5*(np.max(Y) - np.min(Y))
    XY = np.column_stack((X, Y))
    np.savetxt('mesh/body.txt', XY, fmt='%8.5f')
    plt.title('Body Preview', fontweight='bold')
    plt.plot(X, Y)
    plt.axis('equal')
    plt.show()

# airfoil(M, P, t, n)
diamond(a, n)
