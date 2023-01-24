# HLL Approximate Riemann Solver

## Solver Simulation

In progress.

**Domain Parameters**

| Parameter     | Value                           |
| :---:         | :---:                           |
| $L_x$         | $\tilde{x} \in$ [-1.0, 1.0]     |
| $L_y$         | $\tilde{y} \in$ [-1.0, 1.0]     |
| $L_t$         | $\tilde{t} \in$ [0.0, 10.0]     |
| $D_x$         | $\tilde{d}_x \in$ [-0.05, 0.05] |
| $D_y$         | $\tilde{d}_y \in$ [-0.05, 0.05] |
| $n_x$         | 200                             |
| $n_y$         | 200                             |

**Initial Conditions**

| Parameter             | Value              |
| :---:                 | :---:              |
| $\tilde{\rho}_\infty$ | 1.0                |
| $\tilde{u}_\infty$    | $M_\infty^x$ = 1.0 |
| $\tilde{v}_\infty$    | $M_\infty^y$ = 0.0 |
| $\tilde{p}_\infty$    | 1.0 / $\gamma$     |

## Solver Verification

![image](https://user-images.githubusercontent.com/104728656/213871643-a3c32857-ff93-4d02-932e-e641346562ee.png)

**Domain Parameters**

| Parameter     | Value               |
| :---:         | :---:               |
| $L_x$         | x $\in$ [-1.0, 1.0] |
| $L_y$         | y $\in$ [-1.0, 1.0] |
| $L_t$         | t $\in$ [0.0, 0.52] |
| $n_x$         | 400                 |
| $n_y$         | 400                 |

**Initial Conditions**

| Quadrant      | $\rho_0(x,y)$ | $u_0(x,y)$    | $v_0(x,y)$    | $p_0(x,y)$    |
| :---:         | :---:         | :---:         | :---:         | :---:         |
| x > 0, y > 0  | 0.5313        | 0.0           | 0.0           | 0.4           |
| x < 0, y > 0  | 1.0           | 0.7276        | 0.0           | 1.0           |
| x < 0, y < 0  | 0.8           | 0.0           | 0.0           | 1.0           |
| x > 0, y < 0  | 1.0           | 0.0           | 0.7276        | 1.0           |

## Solver Theory

![image](https://user-images.githubusercontent.com/104728656/213981312-e6ba28d0-b8ac-4c0a-b1e3-c3993dc08189.png)

# HLLC Approximate Riemann Solver

## Solver Simulation

In progress.

**Domain Parameters**

| Parameter     | Value                           |
| :---:         | :---:                           |
| $L_x$         | $\tilde{x} \in$ [-1.0, 1.0]     |
| $L_y$         | $\tilde{y} \in$ [-1.0, 1.0]     |
| $L_t$         | $\tilde{t} \in$ [0.0, 10.0]     |
| $D_x$         | $\tilde{d}_x \in$ [-0.01, 0.01] |
| $D_y$         | $\tilde{d}_y \in$ [-0.01, 0.01] |
| $n_x$         | 200                             |
| $n_y$         | 200                             |

**Initial Conditions**

| Parameter             | Value              |
| :---:                 | :---:              |
| $\tilde{\rho}_\infty$ | 1.0                |
| $\tilde{u}_\infty$    | $M_\infty^x$ = 1.0 |
| $\tilde{v}_\infty$    | $M_\infty^y$ = 0.0 |
| $\tilde{p}_\infty$    | 1.0 / $\gamma$     |

## Solver Theory

![image](https://user-images.githubusercontent.com/104728656/213981397-6b5dfabd-97ae-4c3a-bc0c-508f74820e9b.png)
