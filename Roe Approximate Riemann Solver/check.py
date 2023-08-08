import math as m
import numpy as np
import matplotlib.pyplot as plt

mesh = 'airfoil.su2'  # mesh file

# program start

lines = open(mesh, 'r').readlines()

nodes = []  # [[nx, ny], ...]
faces = []  # [[n1, n2], ...]
cells = []  # [[n1, n2, n3], ...]

cell_faces = []  # [[f1, f2, f3], ...]
face_cells = []  # [[c1, c2], ...]
cell_cells = []  # [[c1, c2, c3], ...]
face_faces = []  # [[f1, f2], ...] (local indices)

Sf = []  # face area [Sf1, ...]
nf = []  # face normal [[nfx, nfy], ...]
cf = []  # face centre [[cfx, cfy], ...]
wf = []  # face weighting [wf1, ...]

Vc = []  # cell volume [Vc1, ...]
cc = []  # cell centre [[ccx, ccy], ...]
rc = []  # cell face centre [[[rc1x, rc1y], [rc2x, rc2y], [rc3x, rc3y]], ...]
sc = []  # cell face sign [[sc1, sc2, sc3], ...]

# cells

n_cells = int(lines[1].split()[1])
for i in range(n_cells):
    n1 = int(lines[2 + i].split()[1])
    n2 = int(lines[2 + i].split()[2])
    n3 = int(lines[2 + i].split()[3])
    cells.append([n1, n2, n3])

# nodes

n_nodes = int(lines[2 + n_cells].split()[1])
for i in range(n_nodes):
    nx = float(lines[3 + n_cells + i].split()[0])
    ny = float(lines[3 + n_cells + i].split()[1])
    nodes.append([nx, ny])

# faces

index = 5 + n_cells + n_nodes
n_freestream_faces = int(lines[index].split()[1])
for i in range(n_freestream_faces):
    n1 = int(lines[index + 1 + i].split()[1])
    n2 = int(lines[index + 1 + i].split()[2])
    faces.append([n1, n2])

index += 2 + n_freestream_faces
n_slipwall_faces = int(lines[index].split()[1])
for i in range(n_slipwall_faces):
    n1 = int(lines[index + 1 + i].split()[2])  # swap
    n2 = int(lines[index + 1 + i].split()[1])  # swap
    faces.append([n1, n2])

for i in range(n_cells):
    n1 = cells[i][0]
    n2 = cells[i][1]
    n3 = cells[i][2]

    if [n1, n2] not in faces and [n2, n1] not in faces:
        faces.append([n1, n2])
    if [n2, n3] not in faces and [n3, n2] not in faces:
        faces.append([n2, n3])
    if [n3, n1] not in faces and [n1, n3] not in faces:
        faces.append([n3, n1])
n_faces = len(faces)

start = 0
end = n_freestream_faces
freestream_faces = list(range(start, end))

start += n_freestream_faces
end += n_slipwall_faces
slipwall_faces = list(range(start, end))

start += n_slipwall_faces
end = n_faces
interior_faces = list(range(start, end))
n_interior_faces = len(interior_faces)

boundary_faces = freestream_faces + slipwall_faces
n_boundary_faces = len(boundary_faces)

# cell_faces

for i in range(n_cells):
    n1 = cells[i][0]
    n2 = cells[i][1]
    n3 = cells[i][2]

    cell_faces_list = []
    for j in range(n_faces):
        if [n1, n2] == faces[j] or [n2, n1] == faces[j]:
            cell_faces_list.append(j)
        if [n2, n3] == faces[j] or [n3, n2] == faces[j]:
            cell_faces_list.append(j)
        if [n3, n1] == faces[j] or [n1, n3] == faces[j]:
            cell_faces_list.append(j)
    cell_faces.append(cell_faces_list)

# face_cells

for i in range(n_faces):
    face_cells_list = []
    for j in range(n_cells):
        if i in cell_faces[j]:
            face_cells_list.append(j)
    if len(face_cells_list) == 1:  # boundary face
        face_cells.append([face_cells_list[0], face_cells_list[0]])
    if len(face_cells_list) == 2:  # interior face
        face_cells.append([face_cells_list[0], face_cells_list[1]])

# Sf, nf, cf

for i in range(n_faces):
    n1 = faces[i][0]
    n2 = faces[i][1]
    n1x = nodes[n1][0]
    n1y = nodes[n1][1]
    n2x = nodes[n2][0]
    n2y = nodes[n2][1]

    Sf.append(m.sqrt((n2x - n1x)**2 + (n2y - n1y)**2))
    nf.append([(n2y - n1y)/Sf[i], -(n2x - n1x)/Sf[i]])
    cf.append([.5*(n2x + n1x), .5*(n2y + n1y)])

    c1 = face_cells[i][0]
    c2 = face_cells[i][1]

    face_faces_list = []
    for j in range(3):
        if cell_faces[c1][j] == i:
            face_faces_list.append(j)
    for j in range(3):
        if cell_faces[c2][j] == i:
            face_faces_list.append(j)
    face_faces.append(face_faces_list)

# Vc, cc, rc, sc

for i in range(n_cells):
    n1 = cells[i][0]
    n2 = cells[i][1]
    n3 = cells[i][2]
    n1x = nodes[n1][0]
    n1y = nodes[n1][1]
    n2x = nodes[n2][0]
    n2y = nodes[n2][1]
    n3x = nodes[n3][0]
    n3y = nodes[n3][1]

    Vc.append(.5*abs((n1x*(n2y - n3y) + n2x*(n3y - n1y) + n3x*(n1y - n2y))))
    cc.append([(n1x + n2x + n3x)/3, (n1y + n2y + n3y)/3])

    rc_list = []
    sc_list = []
    cell_cells_list = []
    for j in range(3):
        fj = cell_faces[i][j]
        cfx = cf[fj][0]
        cfy = cf[fj][1]
        ccx = cc[i][0]
        ccy = cc[i][1]

        rc_list.append([cfx - ccx, cfy - ccy])

        if i == face_cells[fj][0]:
            sc_list.append(1)
            cell_cells_list.append(face_cells[fj][1])
        else:
            sc_list.append(-1)
            cell_cells_list.append(face_cells[fj][0])

    rc.append(rc_list)
    sc.append(sc_list)
    cell_cells.append(cell_cells_list)

# wf

for i in range(n_faces):
    c1 = face_cells[i][0]
    c2 = face_cells[i][1]

    d1x = cc[c1][0] - cf[i][0]
    d1y = cc[c1][1] - cf[i][1]
    d2x = cc[c2][0] - cf[i][0]
    d2y = cc[c2][1] - cf[i][1]

    d1 = 1/m.sqrt(d1x*d1x + d1y*d1y)
    d2 = 1/m.sqrt(d2x*d2x + d2y*d2y)

    wf.append(d1/(d1 + d2))

# check mesh

ci = 101  # check cell
fi = 237  # check face

plt.figure(figsize=(8, 8))

for i in range(n_nodes):
    nx = nodes[i][0]
    ny = nodes[i][1]

    plt.text(nx, ny, i, horizontalalignment='center',
            verticalalignment='center')

for i in range(n_faces):
    n1 = faces[i][0]
    n2 = faces[i][1]
    n1x = nodes[n1][0]
    n1y = nodes[n1][1]
    n2x = nodes[n2][0]
    n2y = nodes[n2][1]

    cfx = cf[i][0]
    cfy = cf[i][1]

    plt.text(cfx, cfy, i, horizontalalignment='center',
            verticalalignment='center')

    if i in interior_faces:
        plt.plot([n1x, n2x], [n1y, n2y], c='r', linewidth=2.)
    else:
        plt.plot([n1x, n2x], [n1y, n2y], c='k', linewidth=.5)

    if i == fi:
        c1 = face_cells[fi][0]
        c2 = face_cells[fi][1]
        f1 = face_faces[fi][0]
        f2 = face_faces[fi][1]

        print('Face', fi, c1, c2, f1, f2, cell_faces[c1], cell_faces[c2])

for i in range(n_cells):
    ccx = cc[i][0]
    ccy = cc[i][1]

    plt.text(ccx, ccy, i, horizontalalignment='center',
            verticalalignment='center')

    if i == ci:
        for j in range(3):
            fj = cell_faces[ci][j]
            cj = cell_cells[ci][j]
            c1 = face_cells[fj][0]
            c2 = face_cells[fj][1]

            cfx = cf[fj][0]
            cfy = cf[fj][1]

            plt.arrow(cfx, cfy, sc[ci][j]*nf[fj][0], sc[ci][j]*nf[fj]
                    [1], color='g', head_width=.2)
            plt.arrow(ccx, ccy, rc[ci][j][0], rc[ci]
                    [j][1], color='b', head_width=.2)

            print('Cell', ci, fj, cj, c1, c2, sc[ci][j])

plt.axis('equal')
plt.grid('on')
plt.show()
