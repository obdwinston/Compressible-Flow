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

## Structured Grid

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

## Unstructured Grid

https://user-images.githubusercontent.com/104728656/219235803-acf27799-3b23-4881-b4fd-9de153abb7c4.mp4

**Domain Parameters**

| Parameter     | Value                        |
| :---:         | :---:                        |
| $L_x$         | $x \in$ [-1.0, 1.0]          |
| $L_y$         | $y \in$ [-1.0, 1.0]          |
| $L_t$         | $t \in$ [0.0, 5.0]           |
| Radius        | $r_c$ = 0.02                 |
| Centre        | ($x_c$, $y_c$) = (-0.5, 0.0) |
| Cells         | $n_{cells}$ ≈ 26000          |

**Initial Conditions**

| Parameter             | Value              |
| :---:                 | :---:              |
| $\tilde{\rho}_\infty$ | 1.0                |
| $\tilde{u}_\infty$    | $M_\infty^x$ = 1.2 |
| $\tilde{v}_\infty$    | $M_\infty^y$ = 0.0 |
| $\tilde{p}_\infty$    | 1.0 / $\gamma$     |

![image](https://user-images.githubusercontent.com/104728656/218762819-8685847a-4d48-48d3-b418-5513fefc0404.png)
