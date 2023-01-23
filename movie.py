import os
import math as m
import numpy as np
import matplotlib.pyplot as plt
import moviepy.editor as mpy

nx = 200
ny = 200
x = np.linspace(-1., 1., nx)
y = np.linspace(-1., 1., ny)
dx = x[1] - x[0]
dy = y[1] - y[0]
dt = .1*min(dx, dy)
t = 10.
g = 1.4

data = 'data'
figure = 'figure'

plt.figure(figsize=(2*6.4, 2*4.8))
Nt = int(t/dt)
for nt in range(Nt):
    print(nt)

    U = np.load(os.path.join(data, 'U{:d}.npy'.format(nt)))

    plt.title('Density at t = %.3f' % (nt*dt))
    plt.xlabel('x')
    plt.ylabel('y')
    xx, yy = np.meshgrid(x, y)
    cmap = plt.get_cmap('jet')
    levels = np.linspace(np.min(np.min(U[:, :, 0])), np.max(np.max(U[:, :, 0])), 100)
    plt.contourf(xx, yy, np.transpose(U[:, :, 0]), cmap=cmap, levels=levels)
    plt.colorbar()
    plt.axis('equal')
    plt.savefig(os.path.join(figure, 'U{:d}.png'.format(nt)))
    plt.clf()

image = []
for nt in range(Nt):
    image.append(os.path.join(figure, 'U{:d}.png'.format(nt)))
clip = mpy.ImageSequenceClip(image, fps=int(Nt/10))
clip.write_videofile('video.mp4')
