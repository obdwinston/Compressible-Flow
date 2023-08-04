# Roe Approximate Riemann Solver

| File          | Purpose                                                                     |
| :---:         | :---                                                                        |
| main.f90      | Main program containing functions and subroutines for 2D Euler flow solver. |
| read.py       | Post-processing program to create Mach animation and $C_p$ plot.            |
| airfoil.py    | NACA airfoil generation program to create .geo file for Gmsh.               |
| airfoil.geo   | File containing geometric information to generate mesh in Gmsh.             |
| airfoil.su2   | File containing mesh information to input to main program.                  |

![image](https://github.com/obdwinston/Compressible-Flow/assets/104728656/67bc44fa-4182-4ab1-8b36-975101d80584)

![image](https://github.com/obdwinston/Compressible-Flow/assets/104728656/4834fa56-b3db-4360-be1a-0568befd4f5f)

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
