# SIMPLE algorithm based CFD solver in Rust

This project is a computational fluid dynamics (CFD) solver written in Rust and post-processed with matplotlib, using the SIMPLE algorithm with a collocated grid to integrate the 2D incompressible steady Navier-Stokes equations. This project is based on the [lectures by Dr. Sandip Mazumder](https://youtube.com/playlist?list=PLVuuXJfoPgT4gJcBAAFPW7uMwjFKB9aqT). It can solve three cases: the lid-driven cavity flow, the pipe flow with a velocity inlet and a gauge pressure outlet, and the backward facing step flow.

## Usage

To run this project, you need to have Rust, Python 3 and matplotlib installed on your system. To run the solver for a specific case, use the following command:

```bash
cargo run --release -- -c <case>
```

where `<case>` can be one of `lid_driven_cavity`, `pipe_flow`, `backward_facing_step`. The solver uses a uniform mesh of size $n_x \times n_y$, which can be specified by the user. To get information about what variables you can change in the solver, you can use the `--help` or `-h` flag.

## Lid-driven cavity flow

The lid-driven cavity flow is a classic benchmark problem for CFD. It consists of a square domain with all the boundaries being solid walls. The top wall moves in the x-direction at a constant speed while the other walls are stationary. The flow is governed by the incompressible Navier-Stokes equations:

$$
\nabla \cdot \mathbf{u} = 0 \\
\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla) \mathbf{u} = -\frac{1}{\rho} \nabla p + \nu \nabla^2 \mathbf{u}
$$

where $\mathbf{u} = (u,v)$ is the velocity vector, $p$ is the pressure, $\rho$ is the density, and $\nu$ is the kinematic viscosity.

The boundary conditions are:

$$
u = 1,\ v = 0 \quad \text{at} \quad y = 1
$$

$$
u = 0,\ v = 0 \quad \text{at} \quad x = 0,\ x = 1,\ y = 0
$$

To validate the results of the lid-driven cavity flow solver, we compare them with the experimental results of Ghia et al. [^1] for Reynolds number 100.

### Results
<p align="middle">
<img align="middle" src="https://user-images.githubusercontent.com/100235899/229377769-c678d206-57b4-490d-be1c-7f23985b3fe1.png" width="500">
<img align="middle" src="https://user-images.githubusercontent.com/100235899/229377781-3fa8f610-eff3-400a-87e2-779cc5ab3fd9.png" width="350">
</p>
<p align="middle">
<img align="middle" src="https://user-images.githubusercontent.com/100235899/229377785-288878e7-f831-412a-93c9-298fc458671d.png" width="350">
<img align="middle" src="https://user-images.githubusercontent.com/100235899/229377966-760cac88-bfd7-4ad6-bcd7-805a0f285b8f.png" width="500">
</p>

## Pipe flow

The pipe flow is another benchmark problem for CFD. It consists of a rectangular domain with a velocity inlet at the left boundary and a gauge pressure outlet at the right boundary. The flow is governed by the same incompressible Navier-Stokes equations as the lid-driven cavity flow.

The boundary conditions are:

$$
u = U_{in},\ v = 0 \quad \text{at} \quad x = 0
$$

$$
p = p_{out} \quad \text{at} \quad x = L
$$

$$
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
u(y) = 1.5 U_{in}(1-y^2/h^2)
$$

$$
\frac{dp}{dx} = -3 \frac{\mu}{h^2}U_{in}
$$

where $\mu$ is the dynamic viscosity, $dp/dx$ is the pressure gradient and $h$ is $H/2$. The figures show a good agreement between the numerical and analytical solutions, confirming the accuracy of the solver.

### Results
<p align="middle">
<img align="middle" src="https://user-images.githubusercontent.com/100235899/229377925-5450b28c-307b-47d7-b3f2-a725b53f1ac1.png" width="600">
<img align="middle" src="https://user-images.githubusercontent.com/100235899/229377933-2c9f70ce-9328-4034-9a8f-3bc24ec8a405.png" width="350">
</p>
<p align="middle">
<img src="https://user-images.githubusercontent.com/100235899/229377937-ba8cd109-ecf9-4162-8769-348b3c6cca2c.png" width="350">
<img src="https://user-images.githubusercontent.com/100235899/229377941-23c1eab2-6f6e-462b-a57f-28bcc33967c4.png" width="350">
</p>


## Backward facing step flow

The backward facing step flow is a more complex problem for CFD. It consists of a rectangular domain with a step at the bottom wall.

The boundary conditions are:

$$
u = U_{in},\ v = 0 \quad \text{at} \quad x = 0 
$$

$$
p = p_{out} \quad \text{at} \quad x = L
$$

$$
u = 0,\ v = 0 \quad \text{at} \quad y = 0,\ 0 < x < L_1
$$

$$
u = 0,\ v = 0 \quad \text{at} \quad y = H_1,\ L_1 < x < L
$$

where $U_{in}$ is the inlet velocity, $p_{out}$ is the outlet pressure, $L$ is the length of the domain, $H$ is the height of the domain, $L_1$ is the length of the step, and $H_1$ is the height of the step.

The Reynolds number based on the inlet velocity and the step height is defined as:

$$
Re = \frac{U_{in} H_1}{\nu}
$$

### Results
<p align="middle">
<img align="middle" src="https://user-images.githubusercontent.com/100235899/229378653-89ebd29d-b123-46a6-a00d-f7837e24ad50.png" width="600">
<img align="middle" src="https://user-images.githubusercontent.com/100235899/229378657-5a16376b-d7bf-41fb-a7a9-a8ea1ff8f2e8.png" width="350">
</p>
<p align="middle">
<img align="middle" src="https://user-images.githubusercontent.com/100235899/229378063-14a6db9d-a107-4109-83c6-11c89d138332.png" width="600">
<img align="middle" src="https://user-images.githubusercontent.com/100235899/229378087-fa05128d-50a4-43b5-ab9c-b22608cf1952.png" width="350">
</p>
<p align="middle">
<img align="middle" src="https://user-images.githubusercontent.com/100235899/229378559-de2ab640-d964-4823-8383-461a96f2b825.png" width="600">
<img align="middle" src="https://user-images.githubusercontent.com/100235899/229378566-af9966fd-57dc-4f04-a281-2cd263b05ec1.png" width="350">
</p>

[^1]: U. Ghia, K. N. Ghia, and C. T. Shin, "High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method," Journal of Computational Physics, vol. 48, no. 3, pp. 387-411, 1982.
