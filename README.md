# HLL Approximate Riemann Solver

**Domain Parameters**

| Parameter     | Value             |
| :---:         | :---:             |
| $L_x$         | x $\in$ [-1, 1]   |
| $L_y$         | y $\in$ [-1, 1]   |
| $L_t$         | t $\in$ [0, 0.52] |
| $n_x$         | 400               |
| $n_y$         | 400               |

**Initial Conditions**

| Quadrant      | $\rho_0(x,y)$ | $u_0(x,y)$    | $v_0(x,y)$    | $p_0(x,y)$    |
| :---:         | :---:         | :---:         | :---:         | :---:         |
| $x>0, y>0$    | 0.5313        | 0.0           | 0.0           | 0.4           |
| $x<0, y>0$    | 1.0           | 0.7276        | 0.0           | 1.0           |
| $x<0, y<0$    | 0.8           | 0.0           | 0.0           | 1.0           |
| $x>0, y<0$    | 1.0           | 0.0           | 0.7276        | 1.0           |

![image](https://user-images.githubusercontent.com/104728656/212928586-139c3e6c-2682-4f25-aff8-03fd80fe8af1.png)
