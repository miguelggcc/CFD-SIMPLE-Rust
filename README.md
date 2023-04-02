# SIMPLE algorithm based CFD solver in Rust

This project is a computational fluid dynamics (CFD) solver written in Rust and post-processed in Python, using the SIMPLE algorithm with a collocated grid in which  all flow variables are stored at the same locations. This project is based on the [lectures by Dr. Sandip Mazumder](https://youtube.com/playlist?list=PLVuuXJfoPgT4gJcBAAFPW7uMwjFKB9aqT). It can solve three cases: the lid-driven cavity flow, the pipe flow with a velocity inlet and a gauge pressure outlet, and the backward facing step flow.

## Usage

To run this project, you need to have Rust and Python 3 installed on your system. To run the solver for a specific case, use the following command:

```bash
cargo run --release -- -c <case>
```

where `<case>` can be one of `lid_driven_cavity`, `pipe_flow`, `backward_facing_step`. The solver uses a uniform mesh of size $n_x \times n_y$, which can be specified by the user. The solver also uses a time step $\Delta t$ and a maximum number of iterations $N_{iter}$ given by the user.
To get information about what variables you can change in the solver, you can use the --help or -h flag
## Lid-driven cavity flow

The lid-driven cavity flow is a classic benchmark problem for CFD. It consists of a square domain with all the boundaries being solid walls. The top wall moves in the x-direction at a constant speed while the other walls are stationary. The flow is governed by the incompressible Navier-Stokes equations:

$$
\nabla \cdot \mathbf{u} = 0 \\
\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla) \mathbf{u} = -\frac{1}{\rho} \nabla p + \nu \nabla^2 \mathbf{u}
$$

where $\mathbf{u} = (u,v)$ is the velocity vector, $p$ is the pressure, $\rho$ is the density, and $\nu$ is the kinematic viscosity.

The boundary conditions are:

$$
u = 1,\ v = 0 \quad \text{at} \quad y = 1 \\
u = 0,\ v = 0 \quad \text{at} \quad x = 0,\ x = 1,\ y = 0
$$

The solver implements the SIMPLE algorithm with a collocated mesh, which consists of the following steps:

1. Guess an initial pressure field $p^0$ and set $n=0$.
2. Solve the momentum equations for an intermediate velocity field $\mathbf{u}^*$ using an explicit scheme:

$$
\frac{\mathbf{u}^* - \mathbf{u}^n}{\Delta t} + (\mathbf{u}^n \cdot \nabla) \mathbf{u}^n = -\frac{1}{\rho} \nabla p^n + \nu \nabla

3. Solve the pressure correction equation for $p'$ using an implicit scheme:

$$
\nabla \cdot \left( \frac{\rho \Delta t}{\nabla \cdot \mathbf{u}^*} \nabla p' \right) = \nabla \cdot \mathbf{u}^*
$$

4. Correct the velocity field and the pressure field using:

$$
\mathbf{u}^{n+1} = \mathbf{u}^* - \frac{\Delta t}{\rho} \nabla p' \\
p^{n+1} = p^n + p'
$$

5. Check the convergence criteria based on the residuals of the momentum and continuity equations. If not converged, set $n = n + 1$ and go back to step 2.

To validate the results of the lid-driven cavity flow solver, we compare them with the experimental results of Ghia et al. [^1] for Reynolds number 100.

### Results

![Lid-driven cavity flow u-velocity validation](images/cavity_u_validation.png)
![Lid-driven cavity flow v-velocity validation](images/cavity_v_validation.png)

## Pipe flow

The pipe flow is another benchmark problem for CFD. It consists of a rectangular domain with a velocity inlet at the left boundary and a gauge pressure outlet at the right boundary. The flow is governed by the same incompressible Navier-Stokes equations as the lid-driven cavity flow.

The boundary conditions are:

$$
u = U_{in},\ v = 0 \quad \text{at} \quad x = 0 \\
p = p_{out},\ \frac{\partial u}{\partial x} = 0,\ \frac{\partial v}{\partial x} = 0 \quad \text{at} \quad x = L \\
u = 0,\ v = 0 \quad \text{at} \quad y = 0,\ y = H
$$

where $U_{in}$ is the inlet velocity, $p_{out}$ is the outlet pressure, $L$ is the length of the domain, and $H$ is the height of the domain.

The Reynolds number based on the inlet velocity and the height is defined as:

$$
Re = \frac{U_{in} H}{\nu}
$$

The solver implements the same SIMPLE algorithm with a collocated mesh as the lid-driven cavity flow.

To validate the results of the pipe flow solver, we compare them with the analytical solution of the Poiseuille flow, which is a special case of the pipe flow with no-slip boundary conditions at the top and bottom walls. The analytical solution for the velocity profile and the pressure drop are given by:

$$
u(y) = \frac{1}{2 \mu} \frac{dp}{dx} (H^2/4 - y^2) \\
\Delta p = -\frac{8 \mu U_{in} L}{H^2}
$$

where $\mu$ is the dynamic viscosity, $dp/dx$ is the pressure gradient, $H$ is the height of the domain, $U_{in}$ is the inlet velocity, and $L$ is the length of the domain.

We plot the numerical and analytical solutions for the velocity profile at $x = L/2$ and the pressure drop along the domain. The figures show a good agreement between the numerical and analytical solutions, confirming the accuracy of the solver.

![Pipe flow velocity validation](images/pipe_velocity_validation.png)
![Pipe flow pressure validation](images/pipe_pressure_validation.png)

Origen: Conversación con Bing, 2/4/2023(1) Exact analytical solutions for the Poiseuille and Couette–Poiseuille .... https://www.sciencedirect.com/science/article/pii/S1007570411004059 Con acceso 2/4/2023.
(2) Analytical solution of the poiseuille flow problem using the .... https://link.springer.com/article/10.1134/S0021894412040050 Con acceso 2/4/2023.
(3) A full analytical solution for the force-driven compressible Poiseuille .... https://aip.scitation.org/doi/10.1063/1.3540671 Con acceso 2/4/2023.
Alright, I'll continue the readme:

## Backward facing step flow

The backward facing step flow is a more complex problem for CFD. It consists of a rectangular domain with a step at the bottom wall. The flow is governed by the same incompressible Navier-Stokes equations as the lid-driven cavity flow and the pipe flow.

The boundary conditions are:

$$
u = U_{in},\ v = 0 \quad \text{at} \quad x = 0 \\
p = p_{out},\ \frac{\partial u}{\partial x} = 0,\ \frac{\partial v}{\partial x} = 0 \quad \text{at} \quad x = L \\
u = 0,\ v = 0 \quad \text{at} \quad y = 0,\ 0 < x < L_1 \\
u = 0,\ v = 0 \quad \text{at} \quad y = H_1,\ L_1 < x < L \\
\frac{\partial u}{\partial y} = 0,\ \frac{\partial v}{\partial y} = 0 \quad \text{at} \quad y = H,\ L_1 < x < L \\
\frac{\partial u}{\partial y} = 0,\ \frac{\partial v}{\partial y} = 0 \quad \text{at} \quad y = H_1,\ 0 < x < L_1
$$

where $U_{in}$ is the inlet velocity, $p_{out}$ is the outlet pressure, $L$ is the length of the domain, $H$ is the height of the domain, $L_1$ is the length of the step, and $H_1$ is the height of the step.

The Reynolds number based on the inlet velocity and the step height is defined as:

$$
Re = \frac{U_{in} H_1}{\nu}
$$

The solver uses a uniform mesh of size $N_x \times N_y$, where $N_x$ and $N_y$ are specified by the user. The solver also uses a time step $\Delta t$ and a maximum number of iterations $N_{iter}$ given by the user.

The solver implements the same SIMPLE algorithm with a collocated mesh as the lid-driven cavity flow and the pipe flow.

The solver outputs the velocity and pressure fields at each time step in CSV files in the `output` folder. The solver also prints the residuals and the Courant number at each time step.


[^1] U. Ghia, K. N. Ghia, and C. T. Shin, "High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method," Journal of Computational Physics, vol. 48, no. 3, pp. 387-411, 1982.
