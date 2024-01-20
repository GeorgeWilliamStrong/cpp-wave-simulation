# C++ 2D Wave Equation Solver

A minimal 2D wave equation solver for a variable velocity model using the finite difference method.

## Mathematical Formulation

The 2D wave equation is defined as follows:

$$
\frac{\partial^2 u}{\partial t^2} = \gamma\left(\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2}\right)
$$

In equation 1, \(u\) is the displacement, \(t\) is time, \(\gamma=v^{2}\), and \(v\) is the velocity of the medium.

## Theory

A numerical solution to equation 1 was achieved using a finite difference method (FDM) approach. In order to implement the FDM, the space and time domains are discretized as follows:

$$
x_{i} = i \Delta x\ \textup{for}\ i = 0, 1, 2,...,nx-1, \hspace{0.7cm} y_{j} = j \Delta y\ \textup{for}\ j = 0,1,2,...,ny-1, \hspace{0.7cm} t_{l} = l \Delta t\ \textup{for}\ l = 0,1,2,...,nt-1.
$$

Equation 1 was then expressed as

$$
\frac{\partial^2 u}{\partial t^2}=\frac{\partial }{\partial x}\left ( \gamma \frac{\partial u}{\partial x} \right )+\frac{\partial }{\partial y}\left ( \gamma \frac{\partial u }{\partial y} \right ).
$$

The components on the right-hand side of equation 2 were approximated by finite differences in two stages:

1. The outer expression for each component was discretized as follows:

$$
\frac{\partial }{\partial x}\left ( \gamma \frac{\partial u}{\partial x} \right )\approx \frac{1}{\Delta x}\left ( \gamma\frac{\partial u}{\partial x}\bigg\rvert_{x=x_{i}+0.5} - \gamma\frac{\partial u}{\partial x}\bigg\rvert_{x=x_{i}-0.5} \right ).
$$

2. The inner expressions were then approximated by the following:

$$
\begin{align*}
\gamma\frac{\partial u}{\partial x}\bigg\rvert_{x=x_{i}+0.5}&\approx 0.5(\gamma_{i}+\gamma_{i+ 0.5})\frac{u_{i+1}+u_{i}}{\Delta x}, \\
\gamma\frac{\partial u}{\partial x}\bigg\rvert_{x=x_{i}-0.5}&\approx 0.5(\gamma_{i}+\gamma_{i- 0.5})\frac{u_{i}-u_{i-1}}{\Delta x}
\end{align*}
$$

Finally, the term on the left-hand side of equation 2 was approximated by the following finite difference expression:

$$
\frac{\partial^2 u}{\partial t^2}\approx \frac{u_{ij}^{l-1}-2u_{ij}^{l}+u_{ij}^{l+1}}{\Delta t^{2}}.
$$

Substituting the above expressions into equation 2 and rearranging yields the following formula for calculating the displacement, \(u_{ij}\), at the future time-step, \(t=l+1\):

$$
u_{ij}^{l+1}=2u_{ij}^{l}-u_{ij}^{l-1}+\Psi_{ij}
$$

where \(\Psi_{ij}\) was defined as follows:

$$
\begin{align*}
\Psi_{ij}&=\frac{\Delta t^{2}}{\Delta x}\left ( 0.5(\gamma_{ij}+\gamma_{i+1,j})\frac{u_{i+1,j}-u_{ij}}{\Delta x} - 0.5(\gamma_{ij}+\gamma_{i-1,j})\frac{u_{ij}-u_{i-1,j}}{\Delta x} \right )\\
&\quad +\frac{\Delta t^{2}}{\Delta y}\left ( 0.5(\gamma_{ij}+\gamma_{i,j+1})\frac{u_{i,j+1}-u_{ij}}{\Delta y} - 0.5(\gamma_{ij}+\gamma_{i,j-1})\frac{u_{ij}-u_{i,j-1}}{\Delta y} \right ).
\end{align*}
$$

Equation 3 is used to update the inner points of the model for each time increment. Initial conditions are required to define equation 3 at \(t=0\). Firstly, \(u(t=0,x,y)=I(x,y)\) was used to define \(u_{ij}^{l=0}\) and \(I(x,y)\) was chosen to be the Ricker wavelet function:

$$
I(x,y)=\frac{a}{\pi s^{2}}\left ( 1 - 0.5 \left ( \frac{x^{2}+y^{2}}{s^{2}} \right ) \right )e^{-\frac{x^{2}-y^{2}}{2s^{2}}}
$$

where \(a\) is the amplitude and \(s\) is the spread. Secondly, \(\frac{\partial u}{\partial t}(t=0,x,y)=0\) was approximated as follows:

$$
\frac{\partial u}{\partial t}\approx \frac{u_{ij}^{l+1}-u_{ij}^{l-1}}{2\Delta t}=0,
$$

illustrating that \(u_{ij}^{l+1}=u_{ij}^{l-1}\). Substituting \(u_{ij}^{l+1}\) for \(u_{ij}^{l-1}\) in equation 3, rearranging, and utilizing the fact that \(u_{ij}^{l+1}=u_{ij}^{l-1}\), the following expression for the fictitious value of \(u_{ij}^{l-1}\) at \(t=0\) is obtained:

$$
u_{ij}^{l-1}=u_{ij}^{l}+\frac{\Psi_{ij}}{2}
$$

The basic algorithm for the program can then be expressed as follows:

1. Set initial condition for all points: \(u_{ij}^{l=0}=I(x,y)\)
2. Set fictitious prior displacement for inner points: \(u_{ij}^{l-1}=u_{ij}^{l}+\frac{\Psi_{ij}}{2}\)
3. While \(t<nt\):
    a. Update inner points: \(u_{ij}^{l+1}=2u_{ij}^{l}-u_{ij}^{l-1}+\Psi_{ij}\)
    b. Initialise all points for the next time-step: \(u_{ij}^{l-1}=u_{ij}^{l}\), \(u_{ij}^{l}=u_{ij}^{l+1}\), \(t=t+\Delta t\)

## User Instructions

Upon running the program, `main.cpp`, the user will be asked 'manually enter runfile parameters? (y/n): '. If 'y' is entered, the program prompts the user for the required parameters. If 'n' is entered, the program will load the parameters from 'runfile.txt' which are arranged as follows: model dimensions (\(lx\), \(ly\)), velocities (\(v1, v2, v3, v4\)), grid spacing (\(dx\)), time (\(t\)), Courant number (\(c\)), source origin coordinates (\(ox\), \(oy\)), source amplitude (\(a\)), and source spread (\(s\)). Also note that \(dy = dx\), and the program uses the Courant number to set the time increment (\(dt = cdx/v\)). A suitable runfile has been provided - this was used to generate the results in Figure 1. The program writes a file named 'FDWE.txt' containing the 2D array, \(u_{ij}\), at the final specified time, \(t\). A python notebook, `PLOT-FDWE.ipynb`, that can load 'FDWE.txt' and generate plots such as those in Figure 1 has also been provided.

## Results

The program results are demonstrated in Figure 1, using the example runfile at various time intervals. Note that the boundaries are reflecting (\(u=0\)), and the Courant number needs to be chosen carefully to avoid numerical dispersion and instability. As illustrated in Figure 1, this program can be used to model the propagation of seismic P-waves through geological media with variable velocities.

![Velocity model and modelled displacement at increasing times](FDWE2.png)

*Fig. 1: Velocity model (top-left) and modelled displacement, \(u\), at increasing times, \(t\). Reflections from the first velocity contrast can be seen from 4 seconds, and an emergent head-wave can be observed at 10-12 seconds.*

## References

- Langtangen, H. P. (2013). *A Primer on Scientific Programming with Python*. Springer.
