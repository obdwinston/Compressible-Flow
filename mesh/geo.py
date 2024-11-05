import numpy as np

R = 25.         # domain radius
esf = 4.        # element size factor
sfr = 1.        # size field radius
sft = 25.       # size field thickness
vin = .004      # size field internal element size
vout = .4       # size field external element size
xc = 25.        # size field centre x-coordinate
yc = 25.        # size field centre x-coordinate

XY = np.loadtxt('mesh/body.txt')
X, Y = XY[:, 0] + R - .5, XY[:, 1] + R

geo = open('mesh/mesh.geo', 'w')

geo.write('\n// domain points\n\n') # counter-clockwise
geo.write('Point(%d) = {%.10f, %.10f, 0, 1.0};\n' % (1, R, R))
geo.write('Point(%d) = {%.10f, %.10f, 0, 1.0};\n' % (2, 2*R, R))
geo.write('Point(%d) = {%.10f, %.10f, 0, 1.0};\n' % (3, R, 2*R))
geo.write('Point(%d) = {%.10f, %.10f, 0, 1.0};\n' % (4, 0., R))
geo.write('Point(%d) = {%.10f, %.10f, 0, 1.0};\n' % (5, R, 0.))

geo.write('\n// body points\n\n') # clockwise
for i in range(len(X)):
    geo.write('Point(%d) = {%.10f, %.10f, 0, 1.0};\n' % (i + 6, X[i], Y[i]))

geo.write('\n// domain lines\n\n')
geo.write('Circle(1) = {2, 1, 3};\n')
geo.write('Circle(2) = {3, 1, 4};\n')
geo.write('Circle(3) = {4, 1, 5};\n')
geo.write('Circle(4) = {5, 1, 2};\n')

geo.write('\n// body lines\n\n')
for i in range(len(X) - 1):
    geo.write('Line(%d) = {%d, %d};\n' % (i + 5, i + 6, i + 7))
geo.write('Line(%d) = {%d, %d};\n' % (len(X) + 4, len(X) + 5, 6))

geo.write('\n// curve loops\n\n')
geo.write('Curve Loop(1) = {1, 2, 3, 4};\n')
body = ', '.join(map(str, list(range(5, len(X) + 5))))
geo.write('Curve Loop(2) = {' + body + '};\n')

geo.write('\n// plane surface\n\n') # domain boundary before body boundary
geo.write('Plane Surface(1) = {1, 2};\n')

geo.write('\n// physical groups\n\n')
geo.write('Physical Curve("FREESTREAM", %d) = {1, 2, 3, 4};\n' % (len(X) + 5))
geo.write('Physical Curve("SLIPWALL", %d) = {' % (len(X) + 6) + body + '};\n')

geo.write('\n// size field\n\n')
geo.write('Mesh.MeshSizeFactor = %.10f;\n' % esf)
geo.write('Field[1] = Ball;\n')
geo.write('Field[1].Radius = %.10f;\n' % sfr)
geo.write('Field[1].Thickness = %.10f;\n' % sft)
geo.write('Field[1].VIn = %.10f;\n' % vin)
geo.write('Field[1].VOut = %.10f;\n' % vout)
geo.write('Field[1].XCenter = %.10f;\n' % xc)
geo.write('Field[1].YCenter = %.10f;\n' % yc)
geo.write('Background Field = 1;\n')
