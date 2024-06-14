# Discontinuous Galerkin Method with physics based numerical flux for the numerical modelling of:
A. 1-D elastic wave equation, and
B. Earthquake rupture dynamics with non-linear friction laws

## A. 1-D ealstic wave equation (concise theory)
### Basic Equations

The source-free elastic wave equation in a heterogeneous 1D medium is 

$$\begin{align}
\rho(x)\partial_t v(x,t) -\partial_x \sigma(x,t) & = 0\\
\frac{1}{\mu(x)}\partial_t \sigma(x,t) -\partial_x v(x,t) & = 0  
\end{align}$$

with $\rho(x)$ the density, $\mu(x)$ the shear modulus and $x = [0, L]$. At the boundaries $x = 0, x = L$ we pose the general well-posed linear boundary conditions

$$\begin{equation}
\begin{split}
B_0(v, \sigma, Z_{s}, r_0): =\frac{Z_{s}}{2}\left({1-r_0}\right){v} -\frac{1+r_0}{2} {\sigma} = 0,  \quad \text{at} \quad x = 0, \\
 B_L(v, \sigma, Z_{s}, r_n): =\frac{Z_{s}}{2} \left({1-r_n}\right){v} + \frac{1+r_n}{2}{\sigma} = 0, \quad \text{at} \quad  x = L.
 \end{split}
\end{equation}$$

with the reflection coefficients $r_0$, $r_n$ being real numbers and  $|r_0|, |r_n| \le 1$. 

Note that at $x = 0$,  while  $r_0 = -1$ yields a clamped wall, $r_0 = 0$  yields  an absorbing boundary, and  with $r_0 = 1$  we have a free-surface boundary condition. Similarly, at $x = L$, $r_n = -1$ yields a clamped wall, $r_n = 0$ yields an absorbing boundary, and  $r_n = 1$  gives a free-surface boundary condition.

1. Discretize the spatial domain $x$ into $K$ elements and denote the ${k}^{th}$ element $e^k = [x_{k}, x_{k+1}]$ and the element width 

$$\begin{equation} 
\partial_x{x}_k = x_{k+1}-x_{k}.
\end{equation}$$

Consider two adjacent elements $e^k = [x_{k}, x_{k+1}]$ and  $e^{k+1} = [x_{k+1}, x_{k+2}]$ with an interface at $x_{k+1}$. At the interface we pose the physical conditions for a locked interface

$$\begin{align}
\text{force balance}:  \quad &\sigma^{-} = \sigma^{+} = \sigma, \nonumber \\
\text{no slip}: \quad & [\[ v]\]  = 0,
\end{align}$$

where $$[\[ v]\] = v^{+} - v^{-}$$, and $v^{-}, \sigma^{-}$ and $v^{+}, \sigma^{+}$ are the fields in $e^k = [x_{k}, x_{k+1}]$ and  $e^{k+1} = [x_{k+1}, x_{k+2}]$, respectively. 

2. Within the element derive the weak form of the equation by multiplying both sides by an arbitrary test function and integrating over the element.

3) Next map the $e^k = [x_{k}, x_{k+1}]$ to a reference element $\xi = [-1, 1]$

4) Inside the transformed  element  $\xi \in [-1, 1]$, approximate the solution  and material parameters by a polynomial interpolant,  and write 
\begin{equation}
v^k(\xi, t) = \sum_{j = 1}^{N+1}v_j^k(t) \mathcal{L}_j(\xi), \quad \sigma^k(\xi, t)  = \sum_{j = 1}^{N+1}\sigma_j^k(t) \mathcal{L}_j(\xi),
\end{equation}

\begin{equation}
\rho^k(\xi) = \sum_{j = 1}^{N+1}\rho_j^k \mathcal{L}_j(\xi), \quad \mu^k(\xi) = \sum_{j = 1}^{N+1}\mu_j^k \mathcal{L}_j(\xi),
\end{equation}

where $ \mathcal{L}_j$ is the $j$th interpolating polynomial of degree $N$. If we consider  nodal basis then the interpolating polynomials satisfy $ \mathcal{L}_j(\xi_i) = \delta_{ij}$.

The interpolating nodes $\xi_i$, $i = 1, 2, \dots, N+1$ are the nodes of a Gauss quadrature with

\begin{equation}
 \sum_{i = 1}^{N+1} f(\xi_i)w_i \approx \int_{-1}^{1}f(\xi) d\xi,
\end{equation}

where $w_i$ are quadrature weights.

5) At the element boundaries $\xi = \pm 1$, we generate $\widehat{v}^{k}(\pm 1, t)$ $\widehat{\sigma}^{k}(\pm 1, t)$ by solving a Riemann problem and constraining the solutions against interface and boundary conditions. Then numerical fluctuations $F^k(-1, t)$ and $G^k(1, t)$ are obtained by penalizing hat variables against the incoming characteristics only.

6) Finally, the flux fluctuations are appended to the semi-discrete PDE with special penalty weights and we have 


\begin{equation}
\begin{split}
\frac{d \boldsymbol{v}^k( t)}{ d t} &= \frac{2}{\Delta{x}_k} W^{-1}({\boldsymbol{\rho}}^{k})\left(Q \boldsymbol{\sigma}^k( t) - \boldsymbol{e}_{1}F^k(-1, t)- \boldsymbol{e}_{N+1}G^k(1, t)\right),
\end{split}
\end{equation}

\begin{equation}
\begin{split}
\frac{d \boldsymbol{\sigma}^k( t)}{ d t} &= \frac{2}{\Delta{x}_k} W^{-1}(1/{\boldsymbol{\mu}^{k}})\left(Q \boldsymbol{v}^k( t)  + \boldsymbol{e}_{1}\frac{1}{Z_{s}^{k}(-1)}F^k(-1, t)- \boldsymbol{e}_{N+1}\frac{1}{Z_{s}^{k}(1)}G^k(1, t)\right),
\end{split}
\end{equation}


where 
\begin{align}
\boldsymbol{e}_{1} = [ \mathcal{L}_1(-1), \mathcal{L}_2(-1), \dots,  \mathcal{L}_{N+1}(-1) ]^T, \quad  \boldsymbol{e}_{N+1} = [ \mathcal{L}_1(1), \mathcal{L}_2(1), \dots,  \mathcal{L}_{N+1}(1) ]^T,
\end{align}
and
\begin{align}
G^k(1, t):= \frac{Z_{s}^{k}(1)}{2} \left(v^{k}(1, t)-\widehat{v}^{k}(1, t) \right) + \frac{1}{2}\left(\sigma^{k}(1, t)- \widehat{\sigma}^{k}(1, t)\right), 
\end{align}
\begin{align}
F^{k}(-1, t):= \frac{Z_{s}^{k}(-1)}{2} \left(v^{k}(-1, t)-\widehat{v}^{k}(-1, t) \right) - \frac{1}{2}\left(\sigma^{k}(-1, t)- \widehat{\sigma}^{k}(-1, t)\right).
\end{align}

And the weighted elemental mass matrix $W^N(a)$ and the stiffness matrix $Q^N $ are defined by

\begin{align}
W_{ij}(a) = \sum_{m = 1}^{N+1} w_m \mathcal{L}_i(\xi_m)  {\mathcal{L}_j(\xi_m)} a(\xi_m), \quad Q_{ij} = \sum_{m = 1}^{N+1} w_m \mathcal{L}_i(\xi_m)  {\mathcal{L}_j^{\prime}(\xi_m)}.
\end{align}

7) Time extrapolation can be performed using any stable time stepping scheme like Runge-Kutta or ADER scheme.This notebook implements both Runge-Kutta and ADER schemes for solving the free source version of the elastic wave equation in a homogeneous media. To keep the problem simple, we use as spatial initial condition a Gauss function with half-width $\delta$

\begin{equation}
v(x,t=0)  = e^{-1/\delta^2 (x - x_{o})^2}, \quad \sigma(x,t=0) = 0
\end{equation}

**** Exercises****
1. Lagrange polynomial is used to interpolate the solution and the material parameters. First use polynomial degree 2 and then 6. Compare the simulation results in terms of accuracy of the solution (third and fourth figures give erros). At the end of simulation, time required to complete the simulation is also printed. Also compare the time required to complete both simulations.

2. We use quadrature rules: Gauss-Legendre-Lobatto and Gauss-Legendre. Run simulations once using Lobatto and once using Legendre rules. Compare the difference.

3. Now fix the order of polynomial to be 6, for example. Then use degree of freedom 100 and for another simulation 250. What happpens? Also compare the timre required to complete both simulations.

4. Experiment with the boundary conditions by changing the reflection coefficients $r_0$ and $r_n$.

5. You can also play around with sinusoidal initial solution instead of the Gaussian.

6. Change the time-integrator from RK to ADER. Observe if there are changes in the solution or the CFL number. Vary the polynomial order N.

## B. Earthquake rupture dynamics with non-linear friction laws
### Basic Equations

To begin, consider the domain $\Omega = \Omega_{-}\cup  \Omega_{+}$, with $  \Omega_{-}:= [0, x_0]$,  $\Omega_{+}:= [x_0, L]$, $0<x_0<  L$. 

The sub-domains are governed by the source-free elastic wave equation in a heterogeneous 1D medium is 

\begin{align}
\rho^{\pm}(x)\partial_t v^{\pm}(x,t) -\partial_x \sigma^{\pm}(x,t) & = 0\\
\frac{1}{\mu^{\pm}(x)}\partial_t \sigma^{\pm}(x,t) -\partial_x v^{\pm}(x,t) & = 0  
\end{align}

Here, the field variables and material parameters in the sub-domains  $  \Omega_{\pm}$ with the superscripts $\pm$:  $v^{\pm}$,  $\sigma^{\pm}$, $\rho^{\pm}$,  $\mu^{\pm}$, $Z_s^{\pm}$.

The interface at $x = x_0$ is governed by the general non-linear frictional condition. 

\begin{align}
\text{force balance}:  \quad &\sigma^{-} = \sigma^{+} = \sigma, \nonumber \\
\text{fricition law}: \quad &\sigma =  \sigma_{n}\frac{f\left(\left|[\![ v ]\!]\right|,\psi\right)}{\left|[\![ v ]\!]\right|}[\![ v ]\!]
\end{align}

Here, $\sigma_{n} > 0$ is the effective normal stress,  $f\left(\left|[\![ v ]\!]\right|,\psi\right)$ is a non-linear friction coefficient modeling the frictional response of the interface. In general $f\left(\left|[\![ v ]\!]\right|,\psi\right)\ge 0$ is a monotically increasing positive function, with $f\left(0, \psi\right) = 0 $. Note that with $f\left(\left|[\![ v ]\!]\right|,\psi\right)$ we can describe both rate-and-state, and slip-weakening friction laws.

At the external boundaries $ x = 0, x = L$ we pose the general well-posed linear boundary conditions

\begin{equation}
\begin{split}
B_0(v, \sigma, Z_{s}, r_0): =\frac{Z_{s}}{2}\left({1-r_0}\right){v} -\frac{1+r_0}{2} {\sigma} = 0,  \quad \text{at} \quad x = 0, \\
 B_L(v, \sigma, Z_{s}, r_n): =\frac{Z_{s}}{2} \left({1-r_n}\right){v} + \frac{1+r_n}{2}{\sigma} = 0, \quad \text{at} \quad  x = L.
 \end{split}
\end{equation}

with the reflection coefficients $r_0$, $r_n$ being real numbers and  $|r_0|, |r_n| \le 1$. 

Note that at $x = 0$,  while  $r_0 = -1$ yields a clamped wall, $r_0 = 0$  yields  an absorbing boundary, and  with $r_0 = 1$  we have a free-surface boundary condition. Similarly, at $x = L$, $r_n = -1$ yields a clamped wall, $r_n = 0$ yields an absorbing boundary, and  $r_n = 1$  gives a free-surface boundary condition.

We define the mechanical energy in each subdomain by

\begin{equation}
E^{\pm}(t) = \frac{1}{2}\int_{\Omega_{\pm}}{\left({\rho^{\pm}(x)} |v^{\pm}(x, t)|^2 + \frac{1}{\mu^{\pm}(x)}|\sigma^{\pm}(x, t)|^2\right) dx}.
\end{equation}

The elastic wave equation with the frictional interface condition, satisfies the energy balance

\begin{equation}
\frac{d \left(E^-(t)+E^+(t)\right)}{dt}  = -\sigma [\![ v]\!]  -v^-(0, t)\sigma^-(0, t) + v^+(L, t)\sigma^+(L, t),\le 0.
\end{equation}

With the above friction law the first term in the right hand side of the energy rate, $\sigma [\![ v]\!] = f\left(\left|[\![ v ]\!]\right|,\psi\right){\left|[\![ v ]\!]\right|} \ge 0$, is the work done by friction on the fault which is dissipated as heat. 

#### The discontinuous Galerkin spectral element method (DGSEM)

We will design a provably stable DGSEM obeying the energy balance at the discrete level.

1) Discretize the spatial domain $x$ into $K$ elements and denote the ${k}^{th}$ element $e^k = [x_{k}, x_{k+1}]$ and the element width $\Delta{x}_k = x_{k+1}-x_{k}$. Consider two adjacent elements  $e^k = [x_{k}, x_{k+1}]$ and  $e^{k+1} = [x_{k+1}, x_{k+2}]$ with an interface at $x_{k+1}$. At the internal non-frictional interfaces we pose the physical conditions for a locked interface

\begin{align}
\text{force balance}:  \quad &\sigma^{-} = \sigma^{+} = \sigma, \nonumber \\
\text{no slip}: \quad & [\![ v]\!]  = 0
\end{align}

where $v^{-}, \sigma^{-}$ and $v^{+}, \sigma^{+}$ are the fields in $e^k = [x_{k}, x_{k+1}]$ and  $e^{k+1} = [x_{k+1}, x_{k+2}]$, respectively. 

2) Within the element derive the weak form of the equation by multiplying both sides by an arbitrary test function and integrating over the element.

3) Next map the $e^k = [x_{k}, x_{k+1}]$ to a reference element $\xi = [-1, 1]$

4) Inside the transformed  element  $\xi \in [-1, 1]$, approximate the solution  and material parameters by a polynomial interpolant,  and write 
\begin{equation}
v^k(\xi, t) = \sum_{j = 1}^{N+1}v_j^k(t) \mathcal{L}_j(\xi), \quad \sigma^k(\xi, t)  = \sum_{j = 1}^{N+1}\sigma_j^k(t) \mathcal{L}_j(\xi),
\end{equation}

\begin{equation}
\rho^k(\xi) = \sum_{j = 1}^{N+1}\rho_j^k \mathcal{L}_j(\xi), \quad \mu^k(\xi) = \sum_{j = 1}^{N+1}\mu_j^k \mathcal{L}_j(\xi),
\end{equation}

where $ \mathcal{L}_j$ is the $j$th interpolating polynomial of degree $N$. If we consider  nodal basis then the interpolating polynomials satisfy $ \mathcal{L}_j(\xi_i) = \delta_{ij}$.

The interpolating nodes $\xi_i$, $i = 1, 2, \dots, N+1$ are the nodes of a Gauss quadrature with

\begin{equation}
 \sum_{i = 1}^{N+1} f(\xi_i)w_i \approx \int_{-1}^{1}f(\xi) d\xi,
\end{equation}

where $w_i$ are quadrature weights.

5) At the element boundaries $\xi = \pm 1$, we generate transformed hat-variables $\widehat{v}^{k}(\pm 1, t)$, $\widehat{\sigma}^{k}(\pm 1, t)$ by solving a Riemann problem and constraining the solutions against interface and boundary conditions. These hat-variables encode the solutions of the IBVP at the element boundaries. For non-linear friction laws we must solve a non-linear algebraic problem to compute hat-variables at the element boundaries. Once the hat-variables are computed we communicate them to the element by penalizing hat variables against the incoming characteristics only. This gives the numerical fluctuations $F^k(-1, t)$ and $G^k(1, t)$ at the element boundaries. 

6) Finally, the flux fluctuations are appended to the semi-discrete PDE with special penalty weights and we have 

\begin{equation}
\begin{split}
\frac{d \boldsymbol{v}^k( t)}{ d t} &= \frac{2}{\Delta{x}_k} W^{-1}({\boldsymbol{\rho}}^{k})\left(Q \boldsymbol{\sigma}^k( t) - \boldsymbol{e}_{1}F^k(-1, t)- \boldsymbol{e}_{N+1}G^k(1, t)\right),
\end{split}
\end{equation}

\begin{equation}
\begin{split}
\frac{d \boldsymbol{\sigma}^k( t)}{ d t} &= \frac{2}{\Delta{x}_k} W^{-1}(1/{\boldsymbol{\mu}^{k}})\left(Q \boldsymbol{v}^k( t)  + \boldsymbol{e}_{1}\frac{1}{Z_{s}^{k}(-1)}F^k(-1, t)- \boldsymbol{e}_{N+1}\frac{1}{Z_{s}^{k}(1)}G^k(1, t)\right),
\end{split}
\end{equation}

where 
\begin{align}
\boldsymbol{e}_{1} = [ \mathcal{L}_1(-1), \mathcal{L}_2(-1), \dots,  \mathcal{L}_{N+1}(-1) ]^T, \quad  \boldsymbol{e}_{N+1} = [ \mathcal{L}_1(1), \mathcal{L}_2(1), \dots,  \mathcal{L}_{N+1}(1) ]^T,
\end{align}
and
\begin{align}
G^k(1, t):= \frac{Z_{s}^{k}(1)}{2} \left(v^{k}(1, t)-\widehat{v}^{k}(1, t) \right) + \frac{1}{2}\left(\sigma^{k}(1, t)- \widehat{\sigma}^{k}(1, t)\right), 
\end{align}
\begin{align}
F^{k}(-1, t):= \frac{Z_{s}^{k}(-1)}{2} \left(v^{k}(-1, t)-\widehat{v}^{k}(-1, t) \right) - \frac{1}{2}\left(\sigma^{k}(-1, t)- \widehat{\sigma}^{k}(-1, t)\right).
\end{align}

And the weighted elemental mass matrix $W(a)$ and the stiffness matrix $Q $ are defined by

\begin{align}
W_{ij}(a) = \sum_{m = 1}^{N+1} w_m \mathcal{L}_i(\xi_m)  {\mathcal{L}_j(\xi_m)} a(\xi_m), \quad Q_{ij} = \sum_{m = 1}^{N+1} w_m \mathcal{L}_i(\xi_m)  {\mathcal{L}_j^{\prime}(\xi_m)}.
\end{align}

7) Time extrapolation can be performed using any stable time stepping scheme like Runge-Kutta or ADER scheme.This notebook implements both Runge-Kutta and ADER schemes for solving the elastic wave equation in a hoterogeneous medium with a nonlinear frictional interface, modeled by a linear friction law (LN), slip-weakening friction law (SW) or the rate-and-state friction law (RS).

# References:
This implementation is based on the following paper: Duru et al. (2019), A new discontinuous Galerkin spectral element method for elastic waves with physically motivated numerical fluxes. https://arxiv.org/pdf/1802.06380



