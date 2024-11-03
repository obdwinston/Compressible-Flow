## Program Overview

This program implements a Riemann solver for the two-dimensional Euler equations on an unstructured triangular mesh using the finite volume method. The HLLC fluxes are computed by solving the x-split Riemann problem at each face, taking advantage of the rotational invariance of the flux vectors. The solver employs an explicit multi-stage Runge-Kutta temporal discretisation, along with a multi-slope MUSCL gradient reconstruction and van Albada limiter to ensure stability and accuracy, especially near discontinuities.

![](https://github.com/user-attachments/assets/74e6c744-b7af-4530-bac7-a2e9dd163051)

## Program Files

```
root/
│
├── data/                   # saved data folder
│   ├── cp.png              # pressure coefficient plot
│   └── mach.mp4            # mach number animation
│
├── mesh/                   # mesh folder
│   ├── body.txt            # body coordinates (user input)
│   ├── body.py             # generates body.txt
│   ├── geo.py              # generates mesh.geo from body.txt
│   ├── su2.py              # generates mesh.su2 from mesh.geo
│   └── mesh.f90            # generates mesh.txt from mesh.su2
│
├── mods/                   # modules folder
│   ├── mod_mesh.f90        # mesh type and procedures
│   ├── mod_config.f90      # configuration type and procedures
│   ├── mod_solve.f90       # solver procedures
│   ├── mod_flux.f90        # HLLC flux procedures
│   └── mod_utils.f90       # utility procedures
│
├── run.sh                  # script to run program
├── main.f90                # script to run solver
├── read.py                 # script to read and plot saved data
├── config.txt              # configuration for solver (user input)
└── requirements.txt        # dependencies for solver
```

Clone repository:

```bash
git clone https://github.com/obdwinston/Compressible-Flow.git && cd Compressible-Flow
```

Execute program (for macOS users):

```bash
chmod +x run.sh && ./run.sh
```

For Windows users, you need to modify `run.sh` accordingly before executing the program. For custom bodies, coordinates in `mesh/body.txt` should be `x y` space-delimited and in clockwise order, with no repeated points or intersecting lines.

## Solver Verification

### Diamond Airfoil

| Half-Angle | Mach Number | Angle of Attack |
| :--------: | :---------: | :-------------: |
|    15°     |      2      |       0°        |

![](https://github.com/user-attachments/assets/2ba0a703-6139-4788-ba7c-fe2d5112adc3)

https://github.com/user-attachments/assets/08b6ac5a-815d-43f4-acbe-edc91b71cd5b

![](https://github.com/user-attachments/assets/f6672bd8-a343-436c-bd48-3f719abc8828)

### NACA Airfoil

| NACA Designation | Mach Number | Angle of Attack |
| :--------------: | :---------: | :-------------: |
|       0012       |     0.8     |      1.25°      |

https://github.com/user-attachments/assets/2042e650-e90c-4776-9f77-088fa2f2cae7

![](https://github.com/user-attachments/assets/473e00f2-fe8e-4163-a66c-23550a8d2b61)

## Solver Theory

![](https://github.com/user-attachments/assets/fce29d2a-54bb-46a5-9d35-4c071ea22d82)

## References

[1] Toro (2009). _Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction._  
[2] Blazek (2015). _Computational Fluid Dynamics: Principles and Applications._  
[3] Hou et al. (2015). _An Efficient Unstructured MUSCL Scheme for Solving the 2D Shallow Water Equations._  
[4] Curcic (2021). _Modern Fortran: Building Efficient Parallel Applications._  
[5] Anderson (2020). _Modern Compressible Flow with Historical Perspective._  
[6] Pulliam (1986). _Artificial Dissipation Models for the Euler Equations._
