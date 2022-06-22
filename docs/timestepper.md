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