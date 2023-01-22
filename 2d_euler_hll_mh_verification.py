import os
import math as m
import numpy as np
import matplotlib.pyplot as plt

nx = 100
ny = 100
x = np.linspace(-1., 1., nx)
y = np.linspace(-1., 1., ny)
dx = x[1] - x[0]
dy = y[1] - y[0]
dt = .1*min(dx, dy)
t = .52
g = 1.4

r = np.zeros((nx, ny))
u = np.zeros((nx, ny))
v = np.zeros((nx, ny))
p = np.zeros((nx, ny))

U = np.zeros((nx, ny, 4))
Un = np.zeros((nx, ny, 4))
Us = np.zeros((nx, ny, 4))
Ue = np.zeros((nx, ny, 4))
Uw = np.zeros((nx, ny, 4))
Un_bar = np.zeros((nx, ny, 4))
Us_bar = np.zeros((nx, ny, 4))
Ue_bar = np.zeros((nx, ny, 4))
Uw_bar = np.zeros((nx, ny, 4))

# initial conditions
Nx = int(nx/2)
Ny = int(ny/2)

r[Nx:, Ny:] = .5313
u[Nx:, Ny:] = 0
v[Nx:, Ny:] = 0
p[Nx:, Ny:] = .4

r[:Nx, Ny:] = 1.
u[:Nx, Ny:] = .7276
v[:Nx, Ny:] = 0
p[:Nx, Ny:] = 1.

r[:Nx, :Ny] = .8
u[:Nx, :Ny] = 0
v[:Nx, :Ny] = 0
p[:Nx, :Ny] = 1.

r[Nx:, :Ny] = 1.
u[Nx:, :Ny] = 0
v[Nx:, :Ny] = .7276
p[Nx:, :Ny] = 1.

U[:, :, 0] = r
U[:, :, 1] = r*u
U[:, :, 2] = r*v
U[:, :, 3] = .5*r*(u*u + v*v) + p/(g - 1)

def get_xflux(i, j, Ue_bar, Uw_bar):
    rL = Ue_bar[i - 1, j, 0]
    uL = Ue_bar[i - 1, j, 1]/rL
    vL = Ue_bar[i - 1, j, 2]/rL
    EL = Ue_bar[i - 1, j, 3]
    pL = (g - 1)*(EL - .5*rL*(uL*uL + vL*vL))

    rR = Uw_bar[i, j, 0]
    uR = Uw_bar[i, j, 1]/rR
    vR = Uw_bar[i, j, 2]/rR
    ER = Uw_bar[i, j, 3]
    pR = (g - 1)*(ER - .5*rR*(uR*uR + vR*vR))

    # compute wave speeds
    aL = m.sqrt(g*pL/rL)
    aR = m.sqrt(g*pR/rR)
    sL = min(uL - aL, uR - aR)
    sR = max(uL + aL, uR + aR)
    
    # compute intercell flux
    UL = np.array([rL, rL*uL, rL*vL, EL])
    UR = np.array([rR, rR*uR, rR*vR, ER])
    FL = np.array([rL*uL, rL*uL*uL + pL, rL*uL*vL, uL*(EL + pL)])
    FR = np.array([rR*uR, rR*uR*uR + pR, rR*uR*vR, uR*(ER + pR)])

    if sL >= 0:
        return FL
    elif sL <= 0 and sR >= 0:
        return (sR*FL - sL*FR + sL*sR*(UR - UL))/(sR - sL)
    elif sR <= 0:
        return FR

def get_yflux(i, j, Un_bar, Us_bar):
    rB = Un_bar[i, j - 1, 0]
    uB = Un_bar[i, j - 1, 1]/rB
    vB = Un_bar[i, j - 1, 2]/rB
    EB = Un_bar[i, j - 1, 3]
    pB = (g - 1)*(EB - .5*rB*(uB*uB + vB*vB))

    rT = Us_bar[i, j, 0]
    uT = Us_bar[i, j, 1]/rT
    vT = Us_bar[i, j, 2]/rT
    ET = Us_bar[i, j, 3]
    pT = (g - 1)*(ET - .5*rT*(uT*uT + vT*vT))

    # compute wave speeds
    aB = m.sqrt(g*pB/rB)
    aT = m.sqrt(g*pT/rT)
    sB = min(vB - aB, vT - aT)
    sT = max(vB + aB, vT + aT)

    # compute intercell flux
    UB = np.array([rB, rB*uB, rB*vB, EB])
    UT = np.array([rT, rT*uT, rT*vT, ET])
    GB = np.array([rB*vB, rB*uB*vB, rB*vB*vB + pB, vB*(EB + pB)])
    GT = np.array([rT*vT, rT*uT*vT, rT*vT*vT + pT, vT*(ET + pT)])

    if sB >= 0:
        return GB
    elif sB <= 0 and sT >= 0:
        return (sT*GB - sB*GT + sB*sT*(UT - UB))/(sT - sB)
    elif sT <= 0:
        return GT

Nt = int(t/dt)
for nt in range(Nt):
    print('%d of %d' % (nt, Nt))
    
    ############
    ## Step 1 ## : Reconstruct data and boundary extrapolated values
    ############

    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            for k in range(4):
                
                # Step 1a : Compute slope vectors

                dn = U[i, j + 1, k] - U[i, j, k]
                ds = U[i, j, k] - U[i, j - 1, k]
                de = U[i + 1, j, k] - U[i, j, k]
                dw = U[i, j, k] - U[i - 1, j, k]
                di = .5*(de + dw)
                dj = .5*(dn + ds)

                # Step 1b : Compute slope limiters

                if dn != 0:
                    rj = ds/dn
                    if rj <= 0:
                        xj = 0
                    elif rj > 0 and rj <= .5:
                        xj = 2*rj
                    elif rj > .5 and rj <= 1.:
                        xj = 1.
                    else:
                        xj = min(rj, 2/(1 + rj), 2)
                else:
                    xj = 0

                if de != 0:
                    ri = dw/de
                    if ri <= 0:
                        xi = 0
                    elif ri > 0 and ri <= .5:
                        xi = 2*ri
                    elif ri > .5 and ri <= 1.:
                        xi = 1.
                    else:
                        xi = min(ri, 2/(1 + ri), 2)
                else:
                    xi = 0
                
                # Step 1c : Compute boundary extrapolated values

                Un[i, j, k] = U[i, j, k] + .5*xj*dj
                Us[i, j, k] = U[i, j, k] - .5*xj*dj
                Ue[i, j, k] = U[i, j, k] + .5*xi*di
                Uw[i, j, k] = U[i, j, k] - .5*xi*di

    # transmissive boundary
    Un[1:-1, 0, :] = Us[1:-1, 1, :]
    Us[1:-1, 0, :] = Us[1:-1, 1, :]
    Un[1:-1, -1, :] = Un[1:-1, -2, :]
    Us[1:-1, -1, :] = Un[1:-1, -2, :]

    Ue[0, 1:-1, :] = Uw[1, 1:-1, :]
    Uw[0, 1:-1, :] = Uw[1, 1:-1, :]
    Ue[-1, 1:-1, :] = Ue[-2, 1:-1, :]
    Uw[-1, 1:-1, :] = Ue[-2, 1:-1, :]
    
    ############
    ## Step 2 ## : Evolve boundary extrapolated values
    ############

    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            rn = Un[i, j, 0]
            un = Un[i, j, 1]/rn
            vn = Un[i, j, 2]/rn
            En = Un[i, j, 3]
            pn = (g - 1)*(En - .5*rn*(un*un + vn*vn))

            rs = Us[i, j, 0]
            us = Us[i, j, 1]/rs
            vs = Us[i, j, 2]/rs
            Es = Us[i, j, 3]
            ps = (g - 1)*(Es - .5*rs*(us*us + vs*vs))

            re = Ue[i, j, 0]
            ue = Ue[i, j, 1]/re
            ve = Ue[i, j, 2]/re
            Ee = Ue[i, j, 3]
            pe = (g - 1)*(Ee - .5*re*(ue*ue + ve*ve))

            rw = Uw[i, j, 0]
            uw = Uw[i, j, 1]/rw
            vw = Uw[i, j, 2]/rw
            Ew = Uw[i, j, 3]
            pw = (g - 1)*(Ew - .5*rw*(uw*uw + vw*vw))

            Gn = np.array([rn*vn, rn*un*vn, rn*vn*vn + pn, vn*(En + pn)])
            Gs = np.array([rs*vs, rs*us*vs, rs*vs*vs + ps, vs*(Es + ps)])
            Fe = np.array([re*ue, re*ue*ue + pe, re*ue*ve, ue*(Ee + pe)])
            Fw = np.array([rw*uw, rw*uw*uw + pw, rw*uw*vw, uw*(Ew + pw)])

            Un_bar[i, j, :] = Un[i, j, :] + .5*(dt/dx)*(Fw - Fe) + .5*(dt/dy)*(Gs - Gn)
            Us_bar[i, j, :] = Us[i, j, :] + .5*(dt/dx)*(Fw - Fe) + .5*(dt/dy)*(Gs - Gn)
            Ue_bar[i, j, :] = Ue[i, j, :] + .5*(dt/dx)*(Fw - Fe) + .5*(dt/dy)*(Gs - Gn)
            Uw_bar[i, j, :] = Uw[i, j, :] + .5*(dt/dx)*(Fw - Fe) + .5*(dt/dy)*(Gs - Gn)

    # transmissive boundary
    Un_bar[1:-1, 0, :] = Us_bar[1:-1, 1, :]
    Us_bar[1:-1, 0, :] = Us_bar[1:-1, 1, :]
    Un_bar[1:-1, -1, :] = Un_bar[1:-1, -2, :]
    Us_bar[1:-1, -1, :] = Un_bar[1:-1, -2, :]

    Ue_bar[0, 1:-1, :] = Uw_bar[1, 1:-1, :]
    Uw_bar[0, 1:-1, :] = Uw_bar[1, 1:-1, :]
    Ue_bar[-1, 1:-1, :] = Ue_bar[-2, 1:-1, :]
    Uw_bar[-1, 1:-1, :] = Ue_bar[-2, 1:-1, :]

    ############
    ## Step 3 ## : Solve Riemann problem at cell interfaces
    ############

    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            rij = r[i, j]
            uij = u[i, j]
            vij = v[i, j]
            pij = p[i, j]
            Eij = .5*rij*(uij*uij + vij*vij) + pij/(g - 1)

            # Step 3a : Solve x-split Riemann problem for x-fluxes

            Fw = get_xflux(i, j, Ue_bar, Uw_bar)
            Fe = get_xflux(i + 1, j, Ue_bar, Uw_bar)

            # Step 3b : Solve y-split Riemann problem for y-fluxes

            Gs = get_yflux(i, j, Un_bar, Us_bar)
            Gn = get_yflux(i, j + 1, Un_bar, Us_bar)
            
            # Step 3c : Update solution for Euler equations

            U[i, j, :] = np.array([rij, rij*uij, rij*vij, Eij]) + (dt/dx)*(Fw - Fe) + (dt/dy)*(Gs - Gn)

    # transmissive boundary
    U[1:-1, 0, :] = U[1:-1, 1, :]
    U[1:-1, -1, :] = U[1:-1, -2, :]
    U[0, 1:-1, :] = U[1, 1:-1, :]
    U[-1, 1:-1, :] = U[-2, 1:-1, :]

    r[:, :] = U[:, :, 0]
    u[:, :] = U[:, :, 1]/r
    v[:, :] = U[:, :, 2]/r
    p[:, :] = (g - 1)*(U[:, :, 3] - .5*r*(u*u + v*v))

plt.figure(figsize=(2*6.4, 2*4.8))
plt.title('Density at t = %.2f' % t)
plt.xlabel('x')
plt.ylabel('y')
xx, yy = np.meshgrid(x, y)
cmap = plt.get_cmap('viridis')
levels = np.linspace(.54, 1.7, 30)
plt.contour(xx, yy, np.transpose(U[:, :, 0]), cmap=cmap, levels=levels)
plt.colorbar()
plt.axis('equal')
plt.show()
