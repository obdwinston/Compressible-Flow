import numpy as np
import matplotlib.pyplot as plt
import moviepy.editor as mpy
from scipy.interpolate import griddata

xb = [24., 26.]     # grid x-bounds
yb = [24., 26.]     # grid y-bounds
fs = [8, 6]        # figure size
res = 200           # grid resolution
fps = 25            # frames per second

# read configuration

with open("config.txt", "r") as file:
    for line in file:
        if line.startswith("M="):
            M = float(line.split()[1])
        elif line.startswith("g="):
            g = float(line.split()[1])
        elif line.startswith("nt="):
            nt = int(line.split()[1])
        elif line.startswith("nw="):
            nw = int(line.split()[1])

# surface plot

nxy = np.loadtxt('data/nxy.txt')
data = np.loadtxt('data/Wn_{:010d}.txt'.format(nt))
cp = -2.*(data[:, 3] - 1./g)/M**2

plt.figure(figsize=(8, 6))
plt.title('Pressure Coefficient', fontweight='bold')
plt.xlabel('$c$')
plt.ylabel('$-C_{p}$')
plt.plot(nxy[:, 0] - 24.5, cp, lw=3, label='Current')
plt.grid('on')
plt.legend()
plt.savefig('data/cp.png')
plt.show()

# field animation

print('Generating animation...')
cxy = np.loadtxt('data/cxy.txt')
x, y = cxy[:, 0], cxy[:, 1]
grid_x, grid_y = np.mgrid[xb[0]:xb[1]:res*1j, yb[0]:yb[1]:res*1j]

within_box = (x > xb[0]) & (x < xb[1]) & (y > yb[0]) & (y < yb[1])
indices = np.where(within_box)[0]

frames = []
for i in range(nw, nt + nw, nw):
    data = np.loadtxt(f'data/Wc_{i:010d}.txt')
    values = np.sqrt(data[indices, 1]**2 + data[indices, 2]**2)
    grid_values = griddata((x[indices], y[indices]), values, (grid_x, grid_y), method='linear')

    vmin, vmax = np.min(values), np.max(values)
    levels = np.linspace(vmin, vmax, res)

    fig, ax = plt.subplots(figsize=fs)
    cf = ax.contourf(grid_x, grid_y, grid_values, levels=levels, cmap='jet')
    cb = fig.colorbar(cf, ax=ax)
    ax.fill(nxy[:, 0], nxy[:, 1])
    ax.set_title('Mach Number', fontweight='bold')
    ax.set_xlim(xb), ax.set_ylim(yb)
    ax.set_aspect('equal')

    fig.savefig(f'data/mach_{i:010d}.png')
    fig.clf()

    frames.append(f'data/mach_{i:010d}.png')
    print('Saving frame %d of %d (%.3f%%) ... Done!' % (i, nt, i/nt*100.))

animation = mpy.ImageSequenceClip(frames, fps=fps)
animation.write_videofile('data/speed.mp4')
print('Done!')
