# HLL Approximate Riemann Solver

## Solver Simulation

Simulation and troubleshooting in progress.

**Domain Parameters**

| Parameter     | Value                           |
| :---:         | :---:                           |
| $L_x$         | $\tilde{x} \in$ [-1.0, 1.0]     |
| $L_y$         | $\tilde{y} \in$ [-1.0, 1.0]     |
| $L_t$         | $\tilde{t} \in$ [0.0, 5.0]      |
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

![Figure 2](https://user-images.githubusercontent.com/104728656/213871643-a3c32857-ff93-4d02-932e-e641346562ee.png)

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

![Figure 3](https://user-images.githubusercontent.com/104728656/213872397-aaa42bca-244a-4216-bf62-d158abddaddd.png)
