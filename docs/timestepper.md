# Time-stepper details

### Scalar evolution

The advection-diffusion equation for the continous variable is

$$
\partial_t T + \nabla \cdot(\boldsymbol{u} T) = \kappa \nabla^2 T
$$

Rewrite this as linear/nonlinear split:

$$
\partial_t T + \mathcal{L}(T) = \mathcal{N}(\boldsymbol{u}, T)
$$

In Fourier space, this can be expressed

$$
\partial_t \hat{T} + \kappa |\mathbf{k}|^2 \hat{T} = -i\mathbf{k} \widehat{\boldsymbol{u} T}
$$

The linear diffusive terms are treated numerically by a Crank-Nicolson decomposition:

$$
\frac{T^{n+1} - T^{n}}{\Delta t} = -\kappa |\mathbf{k}|^2 \frac{T^n + T^{n+1}}{2}
$$

and the nonlinear terms are treated explicitly using a third-order Runge-Kutta scheme

$$
\frac{T^{l+1} - T^l}{\alpha_l\Delta t} = \beta_l \mathcal{N}^l + \gamma_l \mathcal{N}^{l-1}
$$

Applying the Crank-Nicolson scheme at each RK substep leads to the full equation:

$$
\frac{T^{l+1} - T^l}{\alpha_l\Delta t} + \kappa |\mathbf{k}|^2 \frac{T^l + T^{l+1}}{2} = \beta_l \mathcal{N}^l + \gamma_l \mathcal{N}^{l-1}
$$

Collecting terms to get the next time step in terms of the previous one:

$$
\left[1 + \frac{1}{2}\alpha_l \Delta t\kappa|\mathbf{k}|^2 \right] T^{l+1} = [1 - \frac{1}{2}\alpha_l \Delta t\kappa|\mathbf{k}|^2] T^l + \alpha_l\beta_l \Delta t \mathcal{N}^l + \alpha_l\gamma_l\Delta t \mathcal{n}^l
$$

Mathcing the notation in DIABLO, 

$$
\left[1 + \frac{1}{2}\bar{h}_l \Delta t\kappa|\mathbf{k}|^2 \right] T^{l+1} = [1 - \frac{1}{2}\bar{h}_l \Delta t\kappa|\mathbf{k}|^2] T^l + \bar{h}_l\beta_l \Delta t \mathcal{N}^l + \bar{h}_l\bar{\zeta}_l\Delta t \mathcal{n}^l
$$

$$
[1 + \mathtt{temp1} |\mathbf{k}|^2] \hat{T}^{l+1} = [1 - \mathtt{temp1}|\mathbf{k}|^2] \hat{T}^{l} + \mathtt{temp2}(-i\mathbf{k}\widehat{\mathbf{u}T}) + \mathtt{temp3}(-ik\widehat{\mathbf{u}T})
$$

## Velocity field
The three components of velocity satisfy a similar equation in spectral space, e.g.

$$
\partial_t \hat{u_i} + i \mathbf{k}\cdot \widehat{\boldsymbol{u} u_i} = -ik_i\hat{p} + -\nu |\mathbf{k}|^2 \hat{u_i} + Ri_n \hat{\theta}_n \delta_{ig}
$$

Again, we treat the diffusive (viscous) term semi-implicitly and all other terms explicitly using the Runge-Kutta method.

## Fractional step for pressure evolution
A discretised form of the above equation for the velocity gives an intermediate field $u^*$ which may not be divergence free.
We therefore use a pressure correction $\phi$ such that

$$
\frac{u_i^{l+1} - u_i^*}{\bar{h}^l \Delta t} = -\partial_i \phi
$$

to enforce that the velocity at the next time step is divergence free.
Taking the divergence of the above equation and setting $\nabla \cdot \boldsymbol{u}^{l+1}=0$ gives

$$
\frac{-\partial_i u_i^*}{\bar{h}^l \Delta t} = -\partial_{ii} \phi , \quad \Rightarrow -k_i k_i \hat{\phi} = \frac{ik_i \hat{u}_i^*}{\bar{h}^l \Delta t}
$$

So the two steps are to solve for $\phi$ and then to create the divergence-free field:

$$
\hat{\phi} = -i  \frac{k_i \hat{u}_i^*}{|\mathbf{k}|^2\bar{h}^l\Delta t}, \quad \hat{u}_i^{l+1} = \hat{u}_i^* - ik_i \bar{h}^l \Delta t \hat{\phi}
$$

Finally, we also need to update the pressure field for the next time step using this correction:

$$
\hat{p}^{l+1} = \hat{p}^l + \phi
$$

## Pressure initialisation
At the first time step, if the pressure field is not prescribed, we need to solve for it.
Taking the divergence of the momentum equation and imposing incompressibility tells us that

$$
\nabla^2 p = -\frac{\partial u_i}{\partial x_j} \frac{\partial u_j}{\partial x_i}
$$

We can be efficient with the number of transforms performed to calculate derivatives here by noting that $R_{ij}=R_{ji}$ where $R$ denotes the expression on the right side of the above equation.
It's also useful that

$$
R_{22} = -(\partial_y v)^2 = -(\partial_x u + \partial_z w)^2 = R_{11} + R_{33} - 2\partial_x u \partial_z w
$$

This step is performed in the subroutine `compute_initial_pressure` in the `timestepper` module.