import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import moviepy.editor as mpy

nx = 200
ny = 200
nint = 10
x = np.linspace(-1., 1., nx)
y = np.linspace(-1., 1., ny)
dx = x[1] - x[0]
dy = y[1] - y[0]
dt = .1*min(dx, dy)
t = 10.
g = 1.4

folder = 'data'
file = 'density'
movie = 'movie.mp4'

plt.figure(figsize=(2*6.4, 2*4.8))
N = int(t/dt) + 1
for n in range(1, N):
    if n % nint == 0:
        print(n, N)

        read = open(os.path.join(sys.path[0], folder, file + '{:010d}.txt'.format(n)), 'r').readlines()
        
        u = np.zeros((nx, ny))
        for line in read:
            i = int(line.split()[0])
            j = int(line.split()[1])
            u[i - 1, j - 1] = float(line.split()[2])

        plt.title('Density at t = %.3f' % (n*dt))
        plt.xlabel('x')
        plt.ylabel('y')
        xx, yy = np.meshgrid(x, y)
        cmap = plt.get_cmap('jet')
        levels = np.linspace(np.min(np.min(u)), np.max(np.max(u)), 100)
        plt.contourf(xx, yy, np.transpose(u), cmap=cmap, levels=levels)
        plt.colorbar()
        plt.axis('equal')
        plt.savefig(os.path.join(sys.path[0], folder, file + '{:010d}.png'.format(n)))
        plt.clf()

image = []
for n in range(1, N):
    if n % nint == 0:
        image.append(os.path.join(sys.path[0], folder, file + '{:010d}.png'.format(n)))
clip = mpy.ImageSequenceClip(image, fps=int(len(image)/t))
clip.write_videofile(movie)
