# SIMPLE algorithm based CFD solver in Rust

This project is a computational fluid dynamics (CFD) solver written in Rust and post-processed with matplotlib, using the SIMPLE algorithm with a collocated grid to integrate the 2D incompressible steady Navier-Stokes equations. This project is based on the [lectures by Dr. Sandip Mazumder](https://youtube.com/playlist?list=PLVuuXJfoPgT4gJcBAAFPW7uMwjFKB9aqT). It can solve three cases: lid-driven cavity flow, pipe flow with a velocity inlet/gauge pressure outlet, and backward facing step flow.

## Usage

To run this project, you need to have Rust, Python 3 and matplotlib installed on your system. To run the solver for a specific case, use the following command:

```bash
cargo run --release -- -c <case>
```

where `<case>` can be one of `lid_driven_cavity`, `pipe_flow`, `backward_facing_step`. The solver uses a uniform mesh of size $n_x \times n_y$, which can be specified by the user. To get information about what variables you can change in the solver, you can use the `--help` or `-h` flag.

## Lid-driven cavity flow

The lid-driven cavity flow is a classic benchmark problem for CFD. It consists of a square domain with all the boundaries being solid walls. The top wall moves in the x-direction at a constant speed while the other walls are stationary. The flow is governed by the incompressible Navier-Stokes equations:

$$
\nabla \cdot \mathbf{u} = 0
$$

$$
(\mathbf{u} \cdot \nabla) \mathbf{u} = -\frac{1}{\rho} \nabla p + \nu \nabla^2 \mathbf{u}
$$

where $\mathbf{u} = (u,v)$ is the velocity vector, $p$ is the pressure, $\rho$ is the density, and $\nu$ is the kinematic viscosity.

The boundary conditions are:

$$
u = 1,\ v = 0 \quad \text{at} \quad y = 1
$$

$$
u = 0,\ v = 0 \quad \text{at} \quad x = 0,\ x = 1,\ y = 0
$$

To validate the results of the lid-driven cavity flow solver, they're compared with the experimental data of Ghia et al. [^1] for Reynolds number 100.

### Results
<p align="middle">
<img align="middle" src="https://user-images.githubusercontent.com/100235899/229572745-e5539738-8548-4ee7-b361-dd14d4a2f6f6.svg" width="500">
</p>

<p align="middle">
  <img align="middle" src="https://user-images.githubusercontent.com/100235899/229572848-ae0dfb31-bba6-4213-946a-8d057cdc63ed.png" width="350">
<img align="middle" src="https://user-images.githubusercontent.com/100235899/229575038-2325a136-3db4-4299-8077-b4d000d56beb.png" width="350">
</p>

<p align="middle">
<img align="middle" src="https://user-images.githubusercontent.com/100235899/229581475-6d71f5e5-3aff-48f5-9972-73cfc3be70cd.svg" width="500">
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

To validate the results of the pipe flow solver, they're compared with the analytical solution of the Poiseuille flow, in which the u velocity profile and the pressure drop are given by:

$$
u(y) = 1.5 U_{in}(1-y^2/h^2)
$$

$$
\frac{\partial p}{\partial x} = -3 \frac{\mu}{h^2}U_{in}
$$

where $\mu$ is the viscosity and $h$ is $H/2$. The figures show a good agreement between the numerical and analytical solutions, confirming the accuracy of the solver.

### Results
<p align="middle">
<img align="middle" src="https://user-images.githubusercontent.com/100235899/229579428-475ffa16-94d3-4bdc-b140-e5ba93208390.svg" width="600">
</p>

<p align="middle">
<img align="middle" src="https://user-images.githubusercontent.com/100235899/229579602-a5f8a843-045f-4cac-9579-357d01c08411.png" width="350">
</p>
<p align="middle">
<img src="https://user-images.githubusercontent.com/100235899/229579766-89d5a017-5a71-4f8a-9c67-bd528c9deb30.png" width="300">
  <img src="https://user-images.githubusercontent.com/100235899/229579879-6502113e-1383-43f5-b324-3ad52f101653.png" width="300">
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
u = 0,\ v = 0 \quad \text{at} \quad y = H_1,\ 0 < x < L_1
$$

$$
u = 0,\ v = 0 \quad \text{at} \quad y = 0,\ L_1 < x < L
$$

where $U_{in}$ is the inlet velocity, $p_{out}$ is the outlet pressure, $L$ is the length of the domain, $H$ is the height of the domain, $L_1$ is the length of the step, and $H_1$ is the height of the step.

The Reynolds number based on the inlet velocity and the total height is defined as:

$$
Re = \frac{U_{in} H}{\nu}
$$

### Results
<p align="middle">
<img align="middle" src="https://user-images.githubusercontent.com/100235899/229580299-3055ac31-2814-4780-bf93-ffd870637793.svg" width="600">
<img align="middle" src="https://user-images.githubusercontent.com/100235899/229580400-548a8f22-5c1f-4587-a21f-4e17eb22419d.png" width="350">
</p>

<p align="middle">
<img align="middle" src="https://user-images.githubusercontent.com/100235899/229580687-e8275574-fe76-4b39-8a07-b0ac8f850795.svg" width="600">
<img align="middle" src="https://user-images.githubusercontent.com/100235899/229580775-a3bdc81b-dc04-4a9f-8a57-e0677f95aeb2.png" width="350">
</p>

<p align="middle">

<img align="middle" src="https://user-images.githubusercontent.com/100235899/229581044-676f548b-bd6e-45b0-b599-34760fb8aeb1.svg" width="600">
<img align="middle" src="https://user-images.githubusercontent.com/100235899/229581111-4ebc922b-cb16-4d4d-838d-b40cc115de32.png" width="350">
</p>


[^1]: U. Ghia, K. N. Ghia, and C. T. Shin, "High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method," Journal of Computational Physics, vol. 48, no. 3, pp. 387-411, 1982.
