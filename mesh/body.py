import os
import math as m
import numpy as np
import matplotlib.pyplot as plt

M = .00             # maximum camber [% chord]
P = .00             # maximum camber position [% chord]
t = .12             # thickness [% chord]
a = m.radians(5.)   # angle of attack [deg]
ha = 15.            # half-angle (only for diamond airfoil)
n = 200             # number of panels

def airfoil(M, P, t, n):
    x = (1 - np.cos(np.linspace(0, m.pi, int(n/2))))/2  # cosine spacing
    # x = np.linspace(0, 1, n//2) # linear spacing
    y = np.zeros(n//2)
    for i in range(len(x)):
        if x[i] < P:
            y[i] = (M/P**2)*(2*P*x[i] - x[i]**2)
        else:
            y[i] = (M/(1 - P)**2)*(1 - 2*P + 2*P*x[i] - x[i]**2)
    a0, a1, a2, a3, a4 = .2969, -.1260, -.3516, .2843, -.1036
    yt = t/.2*(a0*x**.5 + a1*x + a2*x**2 + a3*x**3 + a4*x**4)
    yu, yl = y + yt, y - yt
    xa = np.concatenate((x[1:-1], np.flip(x[:])))
    ya = np.concatenate((yu[1:-1], np.flip(yl[:])))
    X = xa*m.cos(a) + ya*m.sin(a)
    Y = -xa*m.sin(a) + ya*m.cos(a)
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
    le_to_top = interpolate(le, top, n//4)
    top_to_te = interpolate(top, te, n//4)
    te_to_btm = interpolate(te, btm, n//4)
    btm_to_le = interpolate(btm, le, n//4)
    XY = np.vstack([le_to_top, top_to_te, te_to_btm, btm_to_le])
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

if not os.path.exists('mesh/body.txt'):
    # airfoil(M, P, t, n)
    diamond(ha, n)
else:
    custom()
