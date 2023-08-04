import math as m
import numpy as np

# airfoil parameters

M = .00  # maximum camber [% chord]
P = .00  # maximum camber position [% chord]
t = .12  # thickness [% chord]
C = 1.  # chord length [m]

n = 100  # number of panels
a = m.radians(0.)  # angle of attack [rad]
xle = 24.5  # leading edge x-coordinate [m]
yle = 25.  # leading edge y-coordinate [m]

L = 50.  # domain length [m]
H = 50.  # domain height [m]

sft = 25.  # size field thickness
sfin = .004  # size field internal element size
sfout = .4  # size field external element size
xmin = 24.  # size field minimum x-coordinate
xmax = 26.  # size field maximum x-coordinate
ymin = 24.  # size field minimum y-coordinate
ymax = 26.  # size field maximum y-coordinate

# program start

# create airfoil

# x = (1 - np.cos(np.linspace(0, m.pi, int(n/2))))/2  # cosine spacing
x = np.linspace(0, 1, int(n/2))
y = np.zeros(int(n/2))

for i in range(len(x)):
    if x[i] < P:
        y[i] = (M/P**2)*(2*P*x[i] - x[i]**2)
    else:
        y[i] = (M/(1 - P)**2)*(1 - 2*P + 2*P*x[i] - x[i]**2)

a0 = .2969
a1 = -.1260
a2 = -.3516
a3 = .2843
a4 = -.1036
yt = t/.2*(a0*x**.5 + a1*x + a2*x**2 + a3*x**3 + a4*x**4)
yu = y + yt
yl = y - yt

xa = C*np.concatenate((np.flip(x[:]), x[1:-1]))
ya = C*np.concatenate((np.flip(yu[:]), yl[1:-1]))
X = xa*m.cos(a) + ya*m.sin(a) + xle
Y = -xa*m.sin(a) + ya*m.cos(a) + yle

# create .geo file

geo = open('airfoil.geo', 'w')

# domain points

geo.write('\n// domain points\n\n')
geo.write('Point(%d) = {%.10f, %.10f, 0, 1.0};\n' % (1, 0., 0.))
geo.write('Point(%d) = {%.10f, %.10f, 0, 1.0};\n' % (2, L, 0.))
geo.write('Point(%d) = {%.10f, %.10f, 0, 1.0};\n' % (3, L, H))
geo.write('Point(%d) = {%.10f, %.10f, 0, 1.0};\n' % (4, 0., H))

# airfoil points

geo.write('\n// airfoil points\n\n')
for i in range(len(X)):
    geo.write('Point(%d) = {%.10f, %.10f, 0, 1.0};\n' % (i + 5, X[i], Y[i]))

# domain lines

geo.write('\n// domain lines\n\n')
geo.write('Line(1) = {1, 2};\n')  # bottom
geo.write('Line(2) = {2, 3};\n')  # right
geo.write('Line(3) = {3, 4};\n')  # top
geo.write('Line(4) = {4, 1};\n')  # left

# airfoil lines

geo.write('\n// airfoil lines\n\n')
for i in range(len(X) - 1):
    geo.write('Line(%d) = {%d, %d};\n' % (i + 5, i + 5, i + 6))
geo.write('Line(%d) = {%d, %d};\n' % (len(X) + 4, len(X) + 4, 5))

# curve loops and plane surface

geo.write('\n// curve loops and plane surface\n\n')
geo.write('Curve Loop(1) = {1, 2, 3, 4};\n')
airfoil = ', '.join(map(str, list(range(5, len(X) + 5))))
geo.write('Curve Loop(2) = {' + airfoil + '};\n')
geo.write('Plane Surface(1) = {1, 2};\n')

# physical groups

geo.write('\n// physical groups\n\n')
geo.write('Physical Curve("FREESTREAM", %d) = {1, 2, 3, 4};\n' % (len(X) + 5))
geo.write('Physical Curve("SLIPWALL", %d) = {' % (
    len(X) + 6) + airfoil + '};\n')

# size field

geo.write('\n// size field\n\n')
geo.write('Field[1] = Box;\n')
geo.write('Field[1].Thickness = %.10f;\n' % sft)
geo.write('Field[1].VIn = %.10f;\n' % sfin)
geo.write('Field[1].VOut = %.10f;\n' % sfout)
geo.write('Field[1].XMin = %.10f;\n' % xmin)
geo.write('Field[1].XMax = %.10f;\n' % xmax)
geo.write('Field[1].YMin = %.10f;\n' % ymin)
geo.write('Field[1].YMax = %.10f;\n' % ymax)
geo.write('Background Field = 1;\n')
