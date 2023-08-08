# Roe Approximate Riemann Solver

## Program Files

| File          | Purpose                                                                     |
| :---:         | :---                                                                        |
| main.f90      | Main program containing functions and subroutines for 2D Euler flow solver. |
| read.py       | Post-processing program to create Mach animation and $C_p$ plot.            |
| airfoil.py    | NACA airfoil generation program to create .geo file for Gmsh.               |
| airfoil.geo   | File containing geometric information to generate mesh in Gmsh.             |
| airfoil.su2   | File containing mesh information to input to main program.                  |
| check.py      | Mesh visualisation program with Matplotlib to verify mesh.                  |

## Program Verification

### NACA Airfoil

https://github.com/obdwinston/Compressible-Flow/assets/104728656/f74f930d-8189-48a7-b372-813bd94c9349

![image](https://github.com/obdwinston/Compressible-Flow/assets/104728656/3ff9166a-3f14-4626-a24d-3da8e4ebc5d0)

### Diamond Airfoil

https://github.com/obdwinston/Compressible-Flow/assets/104728656/67043caa-0926-450b-a79d-1ddc554d24ef

![image](https://github.com/obdwinston/Compressible-Flow/assets/104728656/c6e0221a-a349-4e64-9687-d55e1db7538c)

## Program Theory

![image](https://github.com/obdwinston/Compressible-Flow/assets/104728656/3bd90d02-4c9b-41f5-b896-1a73c90102a4)

![image](https://github.com/obdwinston/Compressible-Flow/assets/104728656/07bc1d6d-7e0d-4e9f-a5e1-6bdcde6603be)

# HLL Approximate Riemann Solver

https://user-images.githubusercontent.com/104728656/215265995-4928992b-e33f-423a-94d8-b56625bb4769.mp4

**Domain Parameters**

| Parameter     | Value                    |
| :---:         | :---:                    |
| $L_x$         | $x \in$ [-1.0, 1.0]      |
| $L_y$         | $y \in$ [-1.0, 1.0]      |
| $L_t$         | $t \in$ [0.0, 10.0]      |
| $D_x$         | $d_x \in$ [-0.60, -0.40] |
| $D_y$         | $d_y \in$ [-0.01, 0.01]  |
| $n_x$         | 200                      |
| $n_y$         | 200                      |

**Initial Conditions**

| Parameter             | Value              |
| :---:                 | :---:              |
| $\tilde{\rho}_\infty$ | 1.0                |
| $\tilde{u}_\infty$    | $M_\infty^x$ = 1.4 |
| $\tilde{v}_\infty$    | $M_\infty^y$ = 0.0 |
| $\tilde{p}_\infty$    | 1.0 / $\gamma$     |

![image](https://user-images.githubusercontent.com/104728656/213981312-e6ba28d0-b8ac-4c0a-b1e3-c3993dc08189.png)

# HLLC Approximate Riemann Solver

https://user-images.githubusercontent.com/104728656/215268445-400de7f4-adea-485c-b2d5-7ada80ffddfb.mp4

**Domain Parameters**

| Parameter     | Value                    |
| :---:         | :---:                    |
| $L_x$         | $x \in$ [-1.0, 1.0]      |
| $L_y$         | $y \in$ [-1.0, 1.0]      |
| $L_t$         | $t \in$ [0.0, 10.0]      |
| $D_x$         | $d_x \in$ [-0.60, -0.40] |
| $D_y$         | $d_y \in$ [-0.01, 0.01]  |
| $n_x$         | 200                      |
| $n_y$         | 200                      |

**Initial Conditions**

| Parameter             | Value              |
| :---:                 | :---:              |
| $\tilde{\rho}_\infty$ | 1.0                |
| $\tilde{u}_\infty$    | $M_\infty^x$ = 1.4 |
| $\tilde{v}_\infty$    | $M_\infty^y$ = 0.0 |
| $\tilde{p}_\infty$    | 1.0 / $\gamma$     |

![image](https://user-images.githubusercontent.com/104728656/218246064-1b87e775-c30a-4d11-9ceb-2375cfed20bd.png)
