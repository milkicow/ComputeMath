# 8 lab

## Theory
### Transport equastion


In this lab using 3 difference schemes we solve transport equastion:



<p align="center">
  <img src="images/1.png">
</p>

<p align="center">
  <img src="images/2.png">
</p>

<p align="center">
  <img src="images/3.png">
</p>


### Analytical solution

Compare the computational and analytical solution of the following equation:

$
\begin{equation}
 \begin{cases}
    u_t' + u_x'= tx
    \\
   u(0, x) = \frac{x^3}{12}
   \\
   u(t, 0) = \frac{t^3}{12}  
 \end{cases}
\end{equation}
$

Analytical solution:

<p align="center">
  <img src="images/4.png">
</p>

Where c = 1

<p align="center">
  <img src="images/5.png">
</p>

My computational solution:

<p align="center">
  <img src="images/6.png">
</p>

Prooved by analytical method