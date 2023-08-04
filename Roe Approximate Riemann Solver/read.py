import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import moviepy.editor as mpy

folder = 'data/'  # data folder
variable = 'M'  # data variable

xmin = 24.  # size field minimum x-coordinate
xmax = 26.  # size field maximum x-coordinate
ymin = 24.  # size field minimum y-coordinate
ymax = 26.  # size field maximum y-coordinate

write = 100  # write interval
frames = 37594 # number of frames
duration = 25.  # animation duration
resolution = 100  # figure resolution

# program start

def in_field(px, py):
    if (xmin < px < xmax) and (ymin < py < ymax):
        return True
    else:
        return False

plt.figure(figsize=(2*6.4, 2*4.8))
grid_x, grid_y = np.mgrid[xmin:xmax:complex(xmin, resolution*(xmax - xmin)), ymin:ymax:complex(ymin, resolution*(ymax - ymin))]

# plot M

airfoil_x = []
airfoil_y = []
lines = open(folder + 'Cp.txt').readlines()
for line in lines:
    px = float(line.split()[0])
    py = float(line.split()[1])
    airfoil_x.append(px)
    airfoil_y.append(py)
px = airfoil_x[0]
py = airfoil_y[0]
airfoil_x.append(px)
airfoil_y.append(py)

points = []
indices = []
lines = open(folder + variable + '{:010d}.txt'.format(write)).readlines()
i = 0
for line in lines:
    px = float(line.split()[0])
    py = float(line.split()[1])
    if in_field(px, py):
        points.append([px, py])
        indices.append(i)
    i += 1

for it in range(1, frames + 1):
    if it % write == 0:
        values = []
        lines = open(folder + variable + '{:010d}.txt'.format(it)).readlines()
        for i in range(len(indices)):
            line = lines[indices[i]]
            values.append(float(line.split()[2]))
        grid_z = griddata(points, values, (grid_x, grid_y), method='linear')

        plt.title('Mach Number')
        plt.xlabel('x')
        plt.ylabel('y')
        cmap = plt.get_cmap('jet')
        levels = np.linspace(min(values), max(values), resolution)
        plt.contourf(grid_x, grid_y, grid_z, cmap=cmap, levels=levels)
        plt.fill(airfoil_x, airfoil_y, zorder=10)
        plt.colorbar()
        plt.axis('equal')
        plt.savefig(folder + variable + '{:010d}.png'.format(it))
        plt.clf()

        print(it)

image = []
for it in range(1, frames + 1):
    if it % write == 0:
        image.append(folder + variable + '{:010d}.png'.format(it))
clip = mpy.ImageSequenceClip(image, fps=int(frames/write/duration))
clip.write_videofile(folder + 'animation.mp4')

# plot Cp

points = []
values = []
lines = open(folder + 'Cp.txt').readlines()
for line in lines:
    points.append(float(line.split()[0]))
    values.append(-float(line.split()[2]))

offset = min(points)
for i in range(len(points)):
    points[i] -= offset

plt.title('Pressure Coefficient')
plt.xlabel('c')
plt.ylabel('Cp (-)')
plt.xlim([0., 1.])
plt.plot(points, values)
plt.grid('on')
plt.savefig(folder + 'Cp.png')
plt.show()
