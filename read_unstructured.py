import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import moviepy.editor as mpy

t = 5.
dt = 1e-4
nint = 100

folder = 'data'
file = 'density'
movie = 'movie.mp4'

x1 = -1. # left domain limit
x2 = 1. # right domain limit
y1 = -1. # bottom domain limit
y2 = 1. # top domain limit
res = 1000 # interpolation resolution
# minval = .5 # minimum display value
# maxval = 2.5 # maximum display value

plt.figure(figsize=(2*6.4, 2*4.8))
N = int(t/dt) + 1
for n in range(1, N):
    if n % nint == 0:
        print(n, N)

        read = open(os.path.join(sys.path[0], folder, file + '{:010d}.txt'.format(n))).readlines()
        
        values = []
        points = []
        for line in read:
            ccx = float(line.split()[0])
            ccy = float(line.split()[1])
            
            values.append(line.split()[2])
            points.append([ccx, ccy])

        plt.title('Density at t = %.3f' % (n*dt))
        plt.xlabel('x')
        plt.ylabel('y')
        grid_x, grid_y = np.mgrid[x1:x2:complex(0, res), y1:y2:complex(0, res)]
        grid_z = griddata(points, values, (grid_x, grid_y), method='linear')
        plt.imshow(grid_z.T, extent=(x1,x2,y1,y2), origin='lower', cmap='jet')
        # plt.imshow(grid_z.T, extent=(x1,x2,y1,y2), origin='lower', vmin=minval, vmax=maxval, cmap='jet')
        plt.colorbar()
        plt.axis('equal')
        plt.savefig(os.path.join(sys.path[0], folder, file + '{:010d}.png'.format(n)))
        plt.clf()

image = []
for n in range(1, N):
    if n % nint == 0:
        image.append(os.path.join(sys.path[0], folder, file + '{:010d}.png'.format(n)))
clip = mpy.ImageSequenceClip(image, fps=int(len(image)/5))
clip.write_videofile(movie)
