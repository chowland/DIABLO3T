# Validation case
## Diffusion of an oscillatory field

In the absence of advection, the scalars simply satisfy the diffusion equation

$$
\partial_t \theta_n = \kappa_n \nabla^2 \theta_n
$$

We can seek a separable solution of the form

$$
\theta_n(x,y,t) = \sin(kx+ly)\Theta(t)
$$

giving

$$
\nabla^2 \theta_n = -(k^2 + l^2) \theta_n \quad \Rightarrow \Theta'(t) = -\kappa_n(k^2+l^2)\Theta
$$

with the solution

$$
\Theta(t) = \Theta_0 \exp[-\kappa_n(k^2 + l^2) t]
$$

We can check the statistics of this solution by noting the mean square of the scalar must obey

$$
\langle \Theta^2\rangle = \frac{\Theta_0}{2} \exp [-2\kappa_n (k^2 + l^2) t]
$$